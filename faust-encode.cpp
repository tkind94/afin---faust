#define __STDC_LIMIT_MACROS
#include <assert.h>
#include <limits.h>
#include <pthread.h>
#include <semaphore.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>

#include <algorithm>
#include <unordered_map>

#include "divsufsort.h"
#include "faust-encode.h"
#include "matching-statistics.h"

Encode::Encode(
    size_t raw_block_size,
    size_t block_size,
    size_t nthreads,
    size_t length_threshold,
    bool dynamic_threshold,
    char *input_path,
    char *output_path)
{
    config.raw_block_size = raw_block_size;
    config.block_size = block_size;
    config.nthreads = nthreads;
    config.length_threshold = length_threshold;
    config.dynamic_threshold = dynamic_threshold;
    input.path = input_path;
    output.path = output_path;

    block.fp = fopen(input.path, "r");

    if (block.fp == NULL) {
        fprintf(stderr, "Error: Can not open %s\n", input.path);
        exit(EXIT_FAILURE);
    }

    char buf[BUFSIZ];
    fgets(buf, BUFSIZ, block.fp);

    read.length = strlen(buf); // Include '\n'

    fseek(block.fp, 0L, SEEK_END);

    input.nbytes = (size_t)ftell(block.fp);

    fseek(block.fp, 0L, SEEK_SET);

    input.nreads = input.nbytes / read.length;

    assert(input.nbytes % read.length == 0); // Assume uniform read length.

    read.processed = 0;

    block.index = NULL;

    block.count = 0;

    block.nreads = config.block_size / read.length;

    block.nbytes = block.nreads * read.length;

    block.raw_nreads = config.raw_block_size / read.length;

    size_t initial_nbytes = block.raw_nreads * read.length;

    input.nblocks = (input.nbytes - initial_nbytes) / block.nbytes;

    if ((input.nbytes - initial_nbytes) % block.nbytes > 0)
        input.nblocks++;

    input.nblocks++; // Make sure to include the initial raw block.

    block.sa = new int32_t[block.nbytes];

    block.text = new uint8_t[block.nbytes + 1];

    encoding = new Encoding[block.nreads];

    pthread_mutex_init(&shared.encoding_mutex, NULL);

    output.stream = new Writer(output_path);
    output.stream->write(input.nreads, 64);
    output.stream->write(block.raw_nreads, 32);
    output.stream->write(block.nreads, 32);
    output.stream->write(read.length - 1, 16);
    output.stream->write(config.length_threshold, 16);
    output.stream->write(dynamic_threshold, 1);

    output.nbits = 161; // 64 + 32 + 32 + 16 + 16 + 1
    output.nomatch = 0;
    output.nfullmatch = 0;
    output.nrcomplement = 0;

    printf("input.path               : %s\n", input.path);
    printf("output.path              : %s\n", output.path);
    printf("input.nbytes             : %zu\n", input.nbytes);
    printf("input.nreads             : %zu\n", input.nreads);
    printf("input.nblocks            : %zu\n", input.nblocks);
    printf("read.length              : %zu\n", read.length - 1); // Excluding \n.
    printf("block.nbytes             : %zu\n", block.nbytes);
    printf("block.nreads             : %zu\n", block.nreads);
    printf("block.raw_nreads         : %zu\n", block.raw_nreads);
    printf("block.initial_nbytes     : %zu\n", initial_nbytes);
    printf("config.raw_block_size    : %zu\n", config.raw_block_size);
    printf("config.block_size        : %zu\n", config.block_size);
    printf("config.length_threshold  : %zu\n", config.length_threshold);
    printf("config.dynamic_threshold : %d\n", config.dynamic_threshold);
    printf("config.nthreads          : %zu\n", config.nthreads);
    printf("\n");
}

Encode::~Encode()
{
    delete output.stream;
    if (block.index)
        delete block.index;
    sem_destroy(&shared.mutex);
    sem_destroy(&shared.nempty);
    sem_destroy(&shared.nstored);
    pthread_mutex_destroy(&shared.encoding_mutex);
    delete[] encoding;
    delete[] block.sa;
    delete[] block.text;
    fclose(block.fp);
}

bool Encode::next_block()
{
    if (block.count > 0) {
        read.processed += block.nreads;

        sem_destroy(&shared.mutex);
        sem_destroy(&shared.nempty);
        sem_destroy(&shared.nstored);
    }

    sem_init(&shared.mutex, 0, 1);

    sem_init(&shared.nempty, 0, MAX_THREADS);

    sem_init(&shared.nstored, 0, 0);

    shared.nread = 0;

    shared.nwrite = 0;

    shared.finished = false;

    block.nreads = config.block_size / read.length;

    char *ptr = (char*)block.text;

    size_t nbytes = 0;
    size_t i = 0;

    while (fgets(ptr, BUFSIZ, block.fp)) {
        size_t length = strlen(ptr);
        assert(length == read.length);

        ptr[length - 1] = '\0'; // Chomp.

        ptr += length;

        nbytes += length;

        i++;

        if (block.count == 0 &&
            i == block.raw_nreads)
            break;

        if (i == block.nreads)
            break;
    }

    if (i == 0)
        return false;

    block.nbytes = nbytes;

    block.nreads = i;

    for (i = 0; i < block.nreads; i++) {
        encoding[i].reference = SIZE_MAX;
        encoding[i].ref_position = 0;
        encoding[i].dst_position = 0;
        encoding[i].length = 0;
        encoding[i].isrcomplement = false;
        encoding[i].nbits = ceil_log2(read.processed + 1) + (2 * (read.length - 1));
    }

    block.count++;

    if (block.count == 1)
        return true; // We not need to build index for the initial raw block.

    divsufsort(block.text, block.sa, (int32_t)block.nbytes);

    if (block.index)
        delete block.index;

    block.index = new MatchingStatistics(block.text, block.sa, block.nbytes);

    return true;
}

static void *produce(void *arg)
{
    assert(arg);

    Encode *ctx = (Encode*)arg;

    FILE *fp = fopen(ctx->input.path, "r");
    assert(fp);

    Progress progress(ctx->read.processed);

    Stopwatch stopwatch;

    while (1) {
        sem_wait(&ctx->shared.nempty);
        sem_wait(&ctx->shared.mutex);

        if (ctx->shared.nread == ctx->read.processed) {
            ctx->shared.finished = true;
            sem_post(&ctx->shared.nstored); // Let consumer threads terminate.
            sem_post(&ctx->shared.nempty);
            sem_post(&ctx->shared.mutex);
            break;
        }

        size_t i = ctx->shared.nread % MAX_THREADS;

        char *ptr = (char*)ctx->shared.reads[i];

        fgets(ptr, MAX_READ_LENGTH, fp);

        size_t length = strlen(ptr);

        ptr[length - 1] = '\0'; // Chomp.

        assert(length == ctx->read.length); // Including \n.

        ctx->shared.nread++;

        sem_post(&ctx->shared.mutex);
        sem_post(&ctx->shared.nstored);

        progress.tick(1);
    }

    progress.done();

    stopwatch.stop();

    fclose(fp);

    return NULL;
}

static void *consume(void *arg)
{
    assert(arg);

    Encode *ctx = (Encode*)arg;

    uint8_t read[MAX_READ_LENGTH] = {0};

    while (1) {
        sem_wait(&ctx->shared.nstored);
        sem_wait(&ctx->shared.mutex);

        if (ctx->shared.finished &&
            ctx->shared.nwrite == ctx->read.processed) {
            sem_post(&ctx->shared.nstored);
            sem_post(&ctx->shared.mutex);
            break;
        }

        size_t i = ctx->shared.nwrite % MAX_THREADS;

        size_t reference = ctx->shared.nwrite;

        size_t n = ctx->read.length - 1;

        memcpy(read, ctx->shared.reads[i], n);

        ctx->shared.nwrite++;

        sem_post(&ctx->shared.mutex);
        sem_post(&ctx->shared.nempty);

        ctx->process_read(read, n, reference, false);

        rcomplement(read, n);

        ctx->process_read(read, n, reference, true);
    }

    return NULL;
}

struct Candidate {
    size_t dst_index;
    size_t dst_position;
    size_t ref_position;
    size_t length;
};

void Encode::process_read(uint8_t *data, size_t n, size_t reference, bool isrcomplement)
{
    assert(data);
    assert(n == read.length - 1);

    struct {
        int pos[MAX_READ_LENGTH];
        int len[MAX_READ_LENGTH];
    } ms = {{0}, {0}};

    // Compute matching statistics.
    int pos = block.index->dollar;
    int len = 0;
    int N = (int)n;
    for (int i = N - 1; i >= 0; i--) {
        block.index->extend_left(pos, len, data[i]);
        if (len && ms.len[i] < len) {
            ms.len[i] = len;
            ms.pos[i] = pos;
        }
    }

    std::unordered_map<size_t, Candidate> candidates;

    size_t computed_threshold = config.length_threshold;
    if (config.dynamic_threshold && (block.count - 1) < (input.nblocks * 0.5))
        computed_threshold = std::max((size_t)16, config.length_threshold / 2);

    for (size_t i = 0; i < n; i++) {
        if (ms.len[i] < (int)computed_threshold)
            continue;

        // Compute [s,e] range in SA from LCP and MS.
        int32_t *lcp = block.index->LCP;

        size_t s = ms.pos[i];
        while(1) {
            if (s == 0) break;
            if (lcp[s] < ms.len[i]) break;
            s--;
        }

        size_t e = ms.pos[i];
        while(1) {
            if (e == block.nbytes - 1) break;
            if (lcp[e + 1] < ms.len[i]) break;
            e++;
        }

        // Update candidates.
        for (size_t j = s; j <= e; j++) {
            size_t k = block.sa[j] / read.length;

	    // Already found a full match.
	    if (encoding[k].length == read.length - 1)
                continue;

	    // Only check against candidates that have an equal or better match length.
	    if (encoding[k].length > (size_t)ms.len[i])
                continue;

            auto it = candidates.find(k);

            if ((it == candidates.end()) ||
                (it != candidates.end() && candidates[k].length < (size_t)ms.len[i]))
                candidates[k] = {
                    block.sa[j] / read.length,
                    block.sa[j] % read.length,
                    i,
                    (size_t)ms.len[i]
                };
        }
    }

    size_t reference_bits = ceil_log2(read.processed + 1);
    size_t position_bits = ceil_log2((read.length - 1) - computed_threshold + 2);

    // Update block encoding data.
    for (auto it = candidates.begin(); it != candidates.end(); it++)
    {
        size_t index = it->first;

        Candidate candidate = it->second;

        size_t read_length = read.length - 1; // Excluding \n.

        // Already found a full match.
        if (encoding[index].length == read_length)
            continue;

        // Only check against candidates that have an equal or better match length.
        if (encoding[index].length > candidate.length)
            continue;

        size_t nbits = reference_bits + 1 + position_bits;

        // For an exact match we store the reference pointer and set pos (read.length+1).
        if (candidate.length < read_length) {
            nbits +=
                position_bits +
                ceil_log2(read_length -
                          std::max(candidate.ref_position,
                                   candidate.dst_position) + 1) +
                (2 * (read_length - candidate.length));
        }

        pthread_mutex_lock(&shared.encoding_mutex);

        if (encoding[index].nbits > nbits) {
            encoding[index].nbits = nbits;
            encoding[index].reference = reference;
            encoding[index].ref_position = candidate.ref_position;
            encoding[index].dst_position = candidate.dst_position;
            encoding[index].length = candidate.length;
            encoding[index].isrcomplement = isrcomplement;

            uint8_t *target = block.text + (index * read.length);
            size_t partition = candidate.dst_position;
            size_t j = 0;

            for (size_t i = 0; i < partition; i++)
                encoding[index].differences[j++] = target[i];

            for (size_t i = partition + candidate.length; i < (read.length - 1); i++)
                encoding[index].differences[j++] = target[i];

            assert(j == (read.length-1) - candidate.length);
        }

        pthread_mutex_unlock(&shared.encoding_mutex);
    }
}

void Encode::write_block()
{
    size_t nomatch = 0;
    size_t nfullmatch = 0;
    size_t nrcomplement = 0;
    size_t nbits = 0;

    size_t computed_threshold = config.length_threshold;
    if (config.dynamic_threshold && (block.count - 1) < (input.nblocks * 0.5))
        computed_threshold = std::max((size_t)16, config.length_threshold / 2);

    size_t reference_bits = ceil_log2(read.processed + 1);
    size_t position_bits = ceil_log2((read.length - 1) - computed_threshold + 2);

    for (size_t i = 0; i < block.nreads; i++) {
        if (encoding[i].length == read.length - 1)
            nfullmatch++;

        if (encoding[i].isrcomplement)
            nrcomplement++;

        nbits += encoding[i].nbits;

        if (encoding[i].reference == SIZE_MAX) {
            nomatch++;

            output.stream->write(read.processed, reference_bits);

            uint8_t *target = block.text + (i * read.length);

            for (size_t j = 0; j < read.length - 1; j++) {
                uint8_t c = nucleotide_to_int[target[j]];
                output.stream->write(c, 2);
            }
        }
        else {
            output.stream->write(encoding[i].reference, reference_bits);
            output.stream->write(encoding[i].isrcomplement, 1);

            if (encoding[i].length == read.length - 1) {
                output.stream->write((read.length - 1) - computed_threshold + 1, position_bits);
                continue;
            }

            size_t length_bits = ceil_log2((read.length - 1) -
                                           std::max(encoding[i].ref_position,
                                                    encoding[i].dst_position) + 1);

            output.stream->write(encoding[i].ref_position, position_bits);
            output.stream->write(encoding[i].dst_position, position_bits);
            output.stream->write(encoding[i].length, length_bits);

            for (size_t j = 0; j < (read.length-1) - encoding[i].length; j++) {
                uint8_t c = nucleotide_to_int[encoding[i].differences[j]];
                output.stream->write(c, 2);
            }
        }
    }

    output.nbits += nbits;
    output.nomatch += nomatch;
    output.nfullmatch += nfullmatch;
    output.nrcomplement += nrcomplement;

    printf("\n");
    vmpeak(stdout);

    printf("\n");
    printf("block.id           : %zu\n", block.count);
    printf("block.nreads       : %zu\n", block.nreads);
    printf("block.nbytes       : %zu\n", block.nbytes);
    printf("block.nomatch      : %zu (%.2f%%)\n", nomatch, (nomatch/(double)block.nreads)*100.0);
    printf("block.fullmatch    : %zu (%.2f%%)\n", nfullmatch, (nfullmatch/(double)block.nreads)*100.0);
    printf("block.nrcomplement : %zu (%.2f%%)\n", nrcomplement, (nrcomplement/(double)block.nreads)*100.0);
    printf("block.nbits        : %zu\n", nbits);
    nbits += 7; // Round up to the nearest byte.
    printf("block.nbytes       : %zu (%.2f MB)\n", nbits / 8, nbits / 8.0 / 1e6);
    printf("block.threshold    : %zu\n", computed_threshold);
    printf("block.compression  : %.2f%%\n", ((nbits/8.0)/block.nbytes)*100.0);
    printf("\n");
    printf("computed_threshold: %zu\n", computed_threshold);
    printf("reference_bits: %zu\n", reference_bits);
    printf("position_bits: %zu\n", position_bits);
    printf("\n");
}

void Encode::write_raw_block()
{
    size_t nbits = (block.nreads * (read.length - 1)) * 2;

    output.nbits += nbits;
    output.nomatch += block.nreads;

    Progress progress(block.nreads);

    for (size_t i = 0; i < block.nreads; i++) {
        uint8_t *target = block.text + (i * read.length);

        for (size_t j = 0; j < read.length - 1; j++) {
            uint8_t c = nucleotide_to_int[target[j]];
            output.stream->write(c, 2);
        }

        progress.tick(1);
    }

    progress.done();

    printf("\n");
    printf("block.id (raw)    : %zu\n", block.count);
    printf("block.nreads      : %zu\n", block.nreads);
    printf("block.nbits       : %zu\n", nbits);
    nbits += 7; // Round up to the nearest byte.
    printf("block.nbytes      : %zu (%.2f MB)\n", nbits / 8, nbits / 8.0 / 1e6);
    printf("block.compression : %.2f%%\n", ((nbits/8.0)/block.nbytes)*100.0);
    printf("\n");
}

void Encode::process_block()
{
    if (block.count == 1) {
        write_raw_block();
        return;
    }

    pthread_t producer;

    pthread_create(&producer, NULL, produce, this);

    pthread_t consumer[MAX_THREADS];

    for (size_t i = 0; i < config.nthreads; i++)
        pthread_create(consumer + i, NULL, consume, this);

    pthread_join(producer, NULL);

    for (size_t i = 0; i < config.nthreads; i++)
        pthread_join(consumer[i], NULL);

    write_block();
}

void Encode::run()
{
    Stopwatch stopwatch;

    while (next_block())
        process_block();

    stopwatch.stop("[runtime]");

    printf("\n");
    vmpeak(stdout);

    printf("\n");
    printf("output.nomatch      : %zu (%.2f%%)\n", output.nomatch, (output.nomatch/(double)input.nreads)*100.0);
    printf("output.fullmatch    : %zu (%.2f%%)\n", output.nfullmatch, (output.nfullmatch/(double)input.nreads)*100.0);
    printf("output.nrcomplement : %zu (%.2f%%)\n", output.nrcomplement, (output.nrcomplement/(double)input.nreads)*100.0);
    printf("output.nbits        : %zu\n", output.nbits);
    printf("output.nbytes       : %zu (%.2f MB)\n", output.nbits / 8, output.nbits / 8.0 / 1e6);
    printf("output.compression  : %.2f%%\n", ((output.nbits/8.0)/input.nbytes)*100.0);
    printf("\n");
}

static size_t s2b(const char *str)
{
    assert(str);

    char *ptr;

    size_t bytes = strtol(str, &ptr, 10);

    switch(*ptr) {
    case 'g':
    case 'G':
        bytes *= 1000; // Fall through.

    case 'm':
    case 'M':
        bytes *= 1000; // Fall through.

    case 'k':
    case 'K':
        bytes *= 1000; // Fall through.

    case 'b':
    case 'B':
    default:
        break;
    }

    return bytes;
}

#define DEFAULT_RAW_BLOCK_SIZE "1g"
#define DEFAULT_BLOCK_SIZE "2g"
#define DEFAULT_LENGTH_THRESHOLD 64
#define DEFAULT_DYNAMIC_THRESHOLD false
#define DEFAULT_THREAD_COUNT 1

static void help(char *path)
{
    assert(path);

    printf("usage: %s -[iorblrth]\n"
           "i|input PATH\n"
           "o|output PATH                  (default: encode.output.TIMESTAMP)\n"
           "r|initial-raw-block-size VALUE (default: %s)\n"
           "b|block-size VALUE             (default: %s)\n"
           "l|length-threshold VALUE       (default: %d)\n"
           "d|dynamic-threshold            (default: %s)\n"
           "t|threads VALUE                (default: %d)\n"
           "h|help\n",
           path,
           DEFAULT_RAW_BLOCK_SIZE,
           DEFAULT_BLOCK_SIZE,
           DEFAULT_LENGTH_THRESHOLD,
           DEFAULT_DYNAMIC_THRESHOLD ? "true" : "fasle",
           DEFAULT_THREAD_COUNT);

    exit(EXIT_SUCCESS);
}

int main(int argc, char **argv)
{
    if (argc == 1) help(argv[0]);

    printf("command:");
    for (int i = 0; i < argc; i++)
        printf(" %s", argv[i]);
    printf("\n\n");

    static const char *optstring = "i:o:r:b:l:dt:h";

    static struct option longopts[] = {
        {"input",                  required_argument, NULL, 'i'},
        {"output",                 required_argument, NULL, 'o'},
        {"initial-raw-block-size", required_argument, NULL, 'r'},
        {"block-size",             required_argument, NULL, 'b'},
        {"length-threshold",       required_argument, NULL, 'l'},
        {"dynamic-threshold",      no_argument,       NULL, 'd'},
        {"threads",                required_argument, NULL, 't'},
        {"help",                   no_argument,       NULL, 'h'},
        { NULL,                    0,                 NULL,  0 }
    };

    char default_output[BUFSIZ];
    time_t t;
    time(&t);
    strftime(default_output, BUFSIZ, "encode.output.%Y%m%d%H%M%S", localtime(&t));

    char *input = NULL;
    char *output = default_output;
    size_t raw_block_size = s2b(DEFAULT_RAW_BLOCK_SIZE);
    size_t block_size = s2b(DEFAULT_BLOCK_SIZE);
    size_t length_threshold = DEFAULT_LENGTH_THRESHOLD;
    size_t dynamic_threshold = DEFAULT_DYNAMIC_THRESHOLD;
    size_t threads = DEFAULT_THREAD_COUNT;

    int opt;
    while ((opt = getopt_long(argc, argv, optstring, longopts, NULL)) != -1) {
        switch (opt) {
        case 'i':
            input = optarg;
            break;

        case 'o':
            output = optarg;
            break;

        case 'r':
            raw_block_size = s2b(optarg);
            break;

        case 'b':
            block_size = s2b(optarg);
            break;

        case 'l':
            length_threshold = s2b(optarg);
            break;

        case 'd':
            dynamic_threshold = true;
            break;

        case 't':
            threads = s2b(optarg);
            break;

        case 'h':
            help(argv[0]);
            break;
        }
    }

    if (input == NULL) {
        fprintf(stderr, "Error: Input file required\n");
        exit(EXIT_FAILURE);
    }

    Encode encode(
        raw_block_size,
        block_size,
        threads,
        length_threshold,
        dynamic_threshold,
        input,
        output
    );

    encode.run();

    return EXIT_SUCCESS;
}
