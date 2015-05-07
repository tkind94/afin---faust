#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <algorithm>

#include "common.h"

#define add_base(x) do {                        \
        data[bit / 8] |= ((x) << (bit % 8));    \
        bit += 2;                               \
    } while(0)

#define base_at(i) ((data[((i)*2) / 8] >> (((i)*2) % 8)) & 0x3)

int main(int argc, char **argv)
{
    if (argc != 3) {
        fprintf(stderr, "usage: %s input output\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    Reader reader(argv[1]);

    FILE *fp = fopen(argv[2], "w");

    size_t nreads = reader.read(64);
    size_t raw_nreads = reader.read(32);
    size_t block_nreads = reader.read(32);
    size_t read_length = reader.read(16);
    size_t length_threshold = reader.read(16);
    size_t dynamic_threshold = reader.read(1);

    size_t nblocks = (nreads - raw_nreads) / block_nreads;
    if ((nreads - raw_nreads) % block_nreads > 0)
        nblocks++;

    nblocks++; // Include the initial raw block.

    size_t data_nbits = (nreads * read_length) * 2;

    size_t data_nbytes = (data_nbits + 7) / 8;

    printf("input             : %s\n", argv[1]);
    printf("output            : %s\n", argv[2]);
    printf("nreads            : %zu\n", nreads);
    printf("nblocks           : %zu\n", nblocks);
    printf("raw_nreads        : %zu\n", block_nreads);
    printf("block_nreads      : %zu\n", block_nreads);
    printf("read_length       : %zu\n", read_length);
    printf("length_threshold  : %zu\n", length_threshold);
    printf("dynamic_threshold : %zu\n", dynamic_threshold);
    printf("memory usage      : %zu\n", data_nbytes);
    printf("\n");

    uint8_t *data = new uint8_t[data_nbytes];
    size_t bit = 0;
    uint8_t base;

    uint8_t read[MAX_READ_LENGTH];

    Stopwatch runtime;

    Progress progress(nreads);

    // Read initial block of reads coded in 2 bpb.
    for (size_t i = 0; i < raw_nreads; i++) {
        for (size_t j = 0; j < read_length; j++) {
            read[j] = reader.read(2);
            add_base(read[j]);
        }

        for (size_t j = 0; j < read_length; j++)
            fputc(int_to_nucleotide[read[j]], fp);
        fputc(10, fp);

        progress.tick(1);
    }

    assert(bit == raw_nreads * read_length * 2);

    size_t nprocessed = raw_nreads;

    for (size_t block = 1; block < nblocks; block++) {
        size_t computed_threshold = length_threshold;

        if (dynamic_threshold &&
            block < (nblocks * 0.5))
            computed_threshold = std::max((size_t)16, length_threshold / 2);

        size_t reference_bits = ceil_log2(nprocessed + 1);

        size_t position_bits = ceil_log2(read_length - computed_threshold + 2);

        for (size_t i = 0; i < block_nreads; i++) {
            if (nprocessed + i == nreads) break;

            size_t reference = reader.read(reference_bits);

            if (reference == nprocessed) {
                for (size_t j = 0; j < read_length; j++) {
                    read[j] = reader.read(2);
                    add_base(read[j]);
                }
            }
            else {
                size_t isrcomplement = reader.read(1);

                uint8_t rcomplement[MAX_READ_LENGTH] = {0};
                size_t rbit;

                if (isrcomplement) {
                    size_t j = read_length - 1;
                    rbit = 0;
                    while (1) {
                        size_t z = (reference * read_length) + j;
                        base = base_at(z);
                        if (base == 0) base = 3;      // 'A' -> 'T'
                        else if (base == 3) base = 0; // 'T' -> 'A'
                        else if (base == 1) base = 2; // 'C' -> 'G'
                        else if (base == 2) base = 1; // 'G' -> 'C'
                        rcomplement[rbit / 8] |= (base << (rbit % 8));
                        rbit += 2;
                        if (j == 0) break;
                        j--;
                    }
                }

                size_t ref_position = reader.read(position_bits);

                if (ref_position == read_length - computed_threshold + 1) {
                    if (isrcomplement) {
                        for (size_t j = 0; j < read_length; j++) {
                            rbit = j * 2;
                            read[j] = (rcomplement[rbit / 8] >> (rbit % 8)) & 0x3;
                            add_base(read[j]);
                        }
                    }
                    else {
                        for (size_t j = 0; j < read_length; j++) {
                            read[j] = base_at((reference * read_length) + j);
                            add_base(read[j]);
                        }
                    }
                }
                else {
                    size_t dst_position = reader.read(position_bits);

                    size_t length_bits =
                        ceil_log2(read_length -
                                  std::max(ref_position, dst_position) + 1);

                    size_t length = reader.read(length_bits);

                    for (size_t j = 0; j < dst_position; j++) {
                        read[j] = reader.read(2);
                        add_base(read[j]);
                    }

                    if (isrcomplement) {
                        for (size_t j = 0; j < length; j++) {
                            rbit = (ref_position + j) * 2;
                            read[dst_position + j] = (rcomplement[rbit / 8] >> ((rbit % 8))) & 0x3;
                            add_base(read[dst_position + j]);
                        }
                    }
                    else {
                        for (size_t j = 0; j < length; j++) {
                            read[dst_position + j] = base_at((reference * read_length) + ref_position + j);
                            add_base(read[dst_position + j]);
                        }
                    }

                    for (size_t j = dst_position + length; j < read_length; j++) {
                        read[j] = reader.read(2);
                        add_base(read[j]);
                    }
                }
            }

            for (size_t j = 0; j < read_length; j++)
                fputc(int_to_nucleotide[read[j]], fp);
            fputc(10, fp);

            progress.tick(1);

        } // for (read ...

        nprocessed += block_nreads;

    }  // for (block ...

    fclose(fp);

    progress.done();

    assert(bit == read_length * nreads * 2);

    runtime.stop("\n[runtime]");

    delete[] data;

    printf("\n");
    vmpeak(stdout);

    return EXIT_SUCCESS;
}
