#define __STDC_LIMIT_MACROS
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <time.h>

#include <algorithm>

#include "common.h"

struct Read {
    uint32_t read;
    uint32_t reference;
    uint8_t reference_position;
    uint8_t read_position;
    uint8_t factor_length;
    uint8_t isrcomplement;
    uint8_t bases[MAX_READ_LENGTH / 4];
};

// Note: Read.bases as a static array is not a good idea. If you start
// playing with MAX_READ_LENGTH memory blows up. Ideally you would
// allocated exactly what is needed dynamically.

Reader *input;
Writer *output;
size_t ntrees;
size_t nreads;
size_t raw_nreads;
size_t block_nreads;
size_t read_length;
size_t length_threshold;
size_t dynamic_threshold;
size_t nbits;

Read* reads;
uint32_t *T;

const size_t reference_bases = 2;    // (read length * 2 bpb)
const size_t reference_internal = 3; // pos,pos,len,non-factor bases
const size_t leaf = 0;               // continue once decoded
const size_t leaf_pop = 1;           // pop stack once decoded

void write_reference(uint32_t index) {
    output->write(reference_bases, 2);

    //fprintf(stderr, "bases: ");

    size_t bit = 0;
    for (size_t j = 0; j < read_length; j++) {
        uint8_t base = reads[index].bases[bit / 8] & 0x3;
        reads[index].bases[bit / 8] >>= 2;
        bit += 2;

        //fprintf(stderr, "%c", int_to_nucleotide[base]);

        output->write(base, 2);
    }

    nbits += 2 + (2 * read_length);

    //fprintf(stderr, "\n");
}

void write_internal(size_t type, uint32_t index) {
    output->write(type, 2);

    size_t position_bits = ceil_log2(read_length - length_threshold + 2);

    size_t length_bits = ceil_log2(
        read_length -
        std::max(reads[index].reference_position,
                 reads[index].read_position) + 1);

    output->write(reads[index].isrcomplement, 1);

    output->write(reads[index].reference_position, position_bits);

    nbits += 2 + 1 + position_bits;

    //fprintf(stderr, "internal type: %zu", type);

    if (reads[index].reference_position == read_length - length_threshold + 1) {
        //fprintf(stderr, "\n");
        return;
    }

    output->write(reads[index].read_position, position_bits);

    output->write(reads[index].factor_length, length_bits);

    size_t bit = 0;
    for (size_t j = 0; j < read_length - reads[index].factor_length; j++) {
        uint8_t base = reads[index].bases[bit / 8] & 0x3;
        reads[index].bases[bit / 8] >>= 2;
        bit += 2;

        output->write(base, 2);
    }

    nbits += position_bits + length_bits + (2 * (read_length - reads[index].factor_length));
}

void fn(uint32_t parent, uint32_t child) {
    if (child == UINT32_MAX) return; // Then the child is a leaf.
    while (child < nreads && reads[child].reference == reads[parent].read) {
        if (T[reads[child].read] != UINT32_MAX) // The child is an internal reference.
            write_internal(reference_internal, child);
        else if ((child + 1) < nreads && reads[child + 1].reference == reads[parent].read)
            write_internal(leaf, child);
        else
            write_internal(leaf_pop, child);
        fn(child, T[reads[child].read]);
        child++;
    }
}

int main(int argc, char **argv)
{
    if (argc != 3) {
        fprintf(stderr, "usage: %s input output\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    input = new Reader(argv[1]);

    nreads = input->read(64);
    raw_nreads = input->read(32);
    block_nreads = input->read(32);
    read_length = input->read(16);
    length_threshold = input->read(16);
    dynamic_threshold = input->read(1);

    size_t nblocks = (nreads - raw_nreads) / block_nreads;
    if ((nreads - raw_nreads) % block_nreads > 0)
        nblocks++;

    nblocks++; // Include the initial raw block.

    output = new Writer(argv[2]);

    printf("input             : %s\n", argv[1]);
    printf("output            : %s\n", argv[2]);
    printf("nblocks           : %zu\n", nblocks);
    printf("nreads            : %zu\n", nreads);
    printf("raw_nreads        : %zu\n", block_nreads);
    printf("block_nreads      : %zu\n", block_nreads);
    printf("read_length       : %zu\n", read_length);
    printf("length_threshold  : %zu\n", length_threshold);
    printf("dynamic_threshold : %zu\n", dynamic_threshold);
    printf("\n");

    reads = new Read[nreads];

    T = new uint32_t[nreads];

    printf("\nsizeof(reads): %zu\n", nreads * sizeof(Read));
    printf("sizeof(T): %zu\n\n", nreads * sizeof(*T));

    nbits = 0;

    Stopwatch runtime;

    Stopwatch decoding("[decoding]");

    uint32_t i;

    for (i = 0; i < nreads; i++) {
        reads[i].read = i;
        reads[i].reference = UINT32_MAX;

        T[i] = UINT32_MAX;
    }

    // Read initial block of reads coded in 2 bpb.
    for (i = 0; i < raw_nreads; i++) {
        size_t bit = 0;
        for (size_t j = 0; j < read_length; j++) {
            uint8_t base = input->read(2);
            reads[i].bases[bit / 8] |= (base << (bit % 8));
            bit += 2;
        }
    }

    size_t nprocessed = raw_nreads;

    for (size_t block = 1; block < nblocks; block++) {
        size_t computed_threshold = length_threshold;

        if (dynamic_threshold &&
            block < (nblocks * 0.5))
            computed_threshold = std::max((size_t)16, length_threshold / 2);

        size_t reference_bits = ceil_log2(nprocessed + 1);

        size_t position_bits = ceil_log2(read_length - computed_threshold + 2);

        for (size_t read = 0; read < block_nreads; read++) {
            if (i == nreads) break;

            reads[i].reference = input->read(reference_bits);

            if (reads[i].reference == nprocessed) {
                reads[i].reference = UINT32_MAX;

                size_t bit = 0;
                for (size_t j = 0; j < read_length; j++) {
                    uint8_t base = input->read(2);
                    reads[i].bases[bit / 8] |= (base << (bit % 8));
                    bit += 2;
                }
            }
            else {
                reads[i].isrcomplement = input->read(1);
                reads[i].reference_position = input->read(position_bits);

                if (reads[i].reference_position < read_length - computed_threshold + 1) {
                    reads[i].read_position = input->read(position_bits);

                    size_t length_bits =
                        ceil_log2(read_length -
                                  std::max(reads[i].reference_position,
                                           reads[i].read_position) + 1);

                    reads[i].factor_length = input->read(length_bits);

                    size_t bit = 0;
                    for (size_t j = 0; j < read_length - reads[i].factor_length; j++) {
                        uint8_t base = input->read(2);
                        reads[i].bases[bit / 8] |= (base << (bit % 8));
                        bit += 2;
                    }
                }
                else if (dynamic_threshold && computed_threshold < length_threshold) {
                    // We have to readjust reference position value on
                    // full matches to account for Afin's fixed length
                    // threshold restriction.
                    reads[i].reference_position = read_length - length_threshold + 1;
                }
            }

            i++;

        } // for (read ...

        nprocessed += block_nreads;

    } // for (block ...

    decoding.stop();

    assert(i == nreads);

    Stopwatch sorting("\n[sorting reads]");

    std::sort(reads, reads + nreads, [](const Read& lhs, const Read& rhs) {
        if (lhs.reference == rhs.reference)
            return lhs.read < rhs.read;
        return lhs.reference < rhs.reference;
    });

    sorting.stop();

    Stopwatch construct("\n[constructing trees]");

    // T[ref] -> First position of ref in sorted reads.
    T[reads[0].reference] = 0;
    for (i = 1; i < nreads; i++) {
        if (reads[i].reference == UINT32_MAX) break;
        if (reads[i - 1].reference != reads[i].reference)
            T[reads[i].reference] = i;
    }

    output->write(nreads, 64);
    output->write(read_length, 16);
    output->write(length_threshold, 16);
    nbits += 96;

    // Construct the trees.
    ntrees = 0;
    for (i = 0; i < nreads; i++) {
        if (reads[i].reference != UINT32_MAX) continue;
        write_reference(i);
        fn(i, T[reads[i].read]);
        ntrees++;
    }

    printf("nreads: %zu\n", nreads);
    printf("ntrees: %zu\n", ntrees);
    printf("average tree size: %.2f\n", nreads/(double)ntrees);
    printf("nbits: %zu (%.2f MB)\n", nbits, (nbits+7)/8.0/1e6);

    construct.stop();

    runtime.stop("\n[runtime]");

    delete[] T;
    delete[] reads;
    delete output;
    delete input;

    printf("\n");
    vmpeak(stdout);

    return EXIT_SUCCESS;
}
