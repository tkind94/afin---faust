#ifndef ENCODE_H
#define ENCODE_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <semaphore.h>
#include <pthread.h>

#include "common.h"

#define MAX_THREADS 32

struct MatchingStatistics;

struct Encoding {
    size_t reference;
    size_t ref_position;
    size_t dst_position;
    size_t length;
    size_t nbits;
    bool isrcomplement;
    uint8_t differences[MAX_READ_LENGTH];
};

struct Encode {
    Encoding *encoding;

    struct {
        FILE *fp;
        uint8_t *text;
        int32_t *sa;
        size_t nbytes;
        size_t nreads;
        size_t raw_nreads;
        size_t count;
        MatchingStatistics *index;
    } block;

    struct {
        size_t length;
        size_t processed;
    } read;

    struct {
        char *path;
        size_t nbytes;
        size_t nreads;
        size_t nblocks;
    } input;

    struct {
        char *path;
        size_t nbits;
        size_t nomatch;
        size_t nfullmatch;
        size_t nrcomplement;
        Writer *stream;
    } output;

    struct {
        bool dynamic_threshold;
        size_t length_threshold;
        size_t block_size;
        size_t raw_block_size;
        size_t nthreads;
    } config;

    struct {
        pthread_mutex_t encoding_mutex;
        uint8_t reads[MAX_THREADS][MAX_READ_LENGTH];
        sem_t mutex;
        sem_t nempty;
        sem_t nstored;
        size_t nread;
        size_t nwrite;
        bool finished;
    } shared;

    Encode(size_t, size_t, size_t, size_t, bool, char*, char*); // Wow.

    ~Encode();

    void run();

    bool next_block();

    void process_block();

    void write_block();

    void write_raw_block();

    void process_read(uint8_t*, size_t, size_t, bool);
};

#endif /* ENCODE_H */
