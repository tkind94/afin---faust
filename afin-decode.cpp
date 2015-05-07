#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <stack>

#include "common.h"

// 10  reference +  is base node (read length * 2 bpb)
// 11  reference +  is internal node (pos,pos,len,non-factor bases)
// 00 !reference + !pop stack once decoded
// 01 !reference +  pop stack once decoded
const size_t reference_bases = 2;
const size_t reference_internal = 3;
const size_t leaf = 0;
const size_t leaf_pop = 1;

int main(int argc, char **argv)
{
    if (argc != 3) {
        fprintf(stderr, "usage: %s input output\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    Reader reader(argv[1]);

    FILE *fp = fopen(argv[2], "w");
    assert(fp);

    size_t nreads = reader.read(64);
    size_t read_length = reader.read(16);
    size_t length_threshold = reader.read(16);

    printf("input            : %s\n", argv[1]);
    printf("output           : %s\n", argv[2]);
    printf("nreads           : %zu\n", nreads);
    printf("read_length      : %zu\n", read_length);
    printf("length_threshold : %zu\n", length_threshold);
    printf("\n");

    std::stack<uint8_t*> references; // :D

    uint8_t read[MAX_READ_LENGTH] = {0};
    uint8_t rcomp[MAX_READ_LENGTH] = {0};
    uint8_t *reference = NULL;
    uint8_t *ptr;

    Stopwatch stopwatch;

    Progress progress(nreads);

    for (size_t i = 0; i < nreads; i++) {
        uint8_t header = reader.read(2);

        if (header == reference_bases) {
            if (reference != NULL)
                delete[] reference;

            while (!references.empty()) {
                reference = references.top();
                delete[] reference;
                references.pop();
            }

            reference = new uint8_t[read_length + 1];
            for (size_t j = 0; j < read_length; j++)
                reference[j] = int_to_nucleotide[reader.read(2)];

            memcpy(rcomp, reference, read_length);

            rcomplement(rcomp, read_length);

            fwrite(reference, 1, read_length, fp);

            fputc(10, fp);

            progress.tick(1);

            continue;
        }

        size_t position_bits = ceil_log2(read_length - length_threshold + 2);

        size_t isrcomplement = reader.read(1);

        ptr = reference;
        if (isrcomplement)
            ptr = rcomp;

        size_t reference_position = reader.read(position_bits);

        if (reference_position == read_length - length_threshold + 1) {
            memcpy(read, ptr, read_length);
        }
        else {
            size_t read_position = reader.read(position_bits);

            size_t length_bits = ceil_log2(read_length -
                                           std::max(reference_position,
                                                    read_position) + 1);

            size_t length = reader.read(length_bits);

            size_t index = 0;

            for (size_t j = 0; j < read_position; j++)
                read[index++] = int_to_nucleotide[reader.read(2)];

            for (size_t j = 0; j < length; j++)
                read[index++] = ptr[reference_position + j];

            for (size_t j = read_position + length; j < read_length; j++)
                read[index++] = int_to_nucleotide[reader.read(2)];
        }

        fwrite(read, 1, read_length, fp);

        fputc(10, fp);

        if (header == reference_internal) {
            references.push(reference);
            reference = new uint8_t[read_length + 1];
            memcpy(reference, read, read_length);
            memcpy(rcomp, reference, read_length);
            rcomplement(rcomp, read_length);
        }
        else if (header == leaf_pop) {
            delete[] reference;
            reference = NULL;
            if (!references.empty()) {
                reference = references.top();
                references.pop();
                memcpy(rcomp, reference, read_length);
                rcomplement(rcomp, read_length);
            }
        }

        progress.tick(1);
    }

    fclose(fp);

    progress.done();

    stopwatch.stop("\n[runtime]");

    if (reference != NULL)
        delete[] reference;

    while (!references.empty()) {
        reference = references.top();
        delete[] reference;
        references.pop();
    }

    printf("\n");
    vmpeak(stdout);

    return EXIT_SUCCESS;
}
