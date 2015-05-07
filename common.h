#ifndef COMMON_H
#define COMMON_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/resource.h>
#include <sys/time.h>
#include <time.h>

#define MAX_READ_LENGTH 101

extern const uint8_t nucleotide_to_int['T' + 1];

extern const uint8_t int_to_nucleotide[4];

struct Reader {
    uint8_t buffer;
    size_t nbits;
    FILE *fp;
    Reader(const char *path) {
        assert(path);
        buffer = 0;
        nbits = 0;
        fp = fopen(path, "rb");
        assert(fp);
    }
    ~Reader() {
        close();
    }
    void close() {
        if (fp == NULL) return;
        fclose(fp);
        fp = NULL;
    }
    void next() {
        fread(&buffer, 1, 1, fp);
        nbits = 8;
    }
    size_t read() {
        if (nbits == 0) next();
        size_t bit = buffer & 1;
        buffer >>= 1;
        nbits--;
        return bit;
    }
    size_t read(size_t bits) {
        assert(bits > 0);
        size_t value = 0;
        size_t pos = 0;
        while (pos < bits) {
            value |= (read() << pos);
            pos++;
        }
        return value;
    }
};

struct Writer {
    uint8_t buffer;
    size_t nbits;
    FILE *fp;
    Writer(const char *path) {
        assert(path);
        buffer = 0;
        nbits = 0;
        fp = fopen(path, "wb");
        assert(fp);
    }
    ~Writer() {
        close();
    }
    void close() {
        if (fp == NULL) return;
        flush();
        fclose(fp);
        fp = NULL;
    }
    void flush() {
        if (nbits == 0) return;
        fwrite(&buffer, 1, 1, fp);
        buffer = 0;
        nbits = 0;
    }
    void write(size_t bit) {
        assert(bit < 2);
        if (nbits == 8) flush();
        buffer |= (bit << nbits);
        nbits++;
    }
    void write(size_t value, size_t bits) {
        assert(bits > 0);
        do {
            write(value & 1);
            value >>= 1;
            bits--;
        } while (bits > 0);
    }
};

// see: http://www.boost.org/doc/libs/1_44_0/boost/progress.hpp
class Progress {
private:
    FILE *stream;
    uint64_t expected_count;
    uint64_t count;

    uint64_t next_tick_count;
    uint64_t ticks;

    static constexpr double tick_width = 50.0;

    void update() {
        uint64_t ticks_needed;

        ticks_needed = ((double)count / expected_count) / tick_width;

        do {
            fprintf(stream, "*");
            fflush(stream);
        } while (++ticks < ticks_needed);

        next_tick_count = (ticks / tick_width) * expected_count;

        if (count == expected_count && ticks > tick_width + 1)
            fprintf(stream, "*\n");
    }

public:
    Progress(uint64_t _expected_count, FILE *_stream = stdout) :
        stream(_stream),
        expected_count(_expected_count),
        count(0),
        next_tick_count(0),
        ticks(0) {
        if (expected_count == 0)
            expected_count = 1; // Prevent division by zero.

        fprintf(stream,
                "0%%   10   20   30   40   50   60   70   80   90   100%%\n"
                "|----|----|----|----|----|----|----|----|----|----|\n");

        fflush(stream);
    }

    void done() {
        while (ticks++ < tick_width)
            fprintf(stream, "*");

        fprintf(stream, "\n");
    }

    void tick(uint64_t increment) {
        count += increment;

        if (count >= next_tick_count)
            update();
    }
};

#define USECS_PER_SEC 1000000

#define NSECS_PER_SEC 1000000000ULL

class Stopwatch
{
private:
    FILE *stream;

    double wall;

    double rusage;

    double _wall()
    {
        double time;

#if defined _POSIX_SOURCE && _POSIX_C_SOURCE >= 199309L // clock_gettime (-lrt)
#if !defined CLOCK_MONOTONIC_RAW /* since Linux 2.6.28 */
# if !defined CLOCK_MONOTONIC
#  error "Missing CLOCK_MONOTONIC. Upgrade your kernel or comment this out."
# endif
# define CLOCK_MONOTONIC_RAW CLOCK_MONOTONIC
#endif
        struct timespec timespec;

        // If the NTP daemon updates the system clock during a timing run
        // gettimeofday will return the adjusted time.

        // CLOCK_MONOTONIC_RAW provides access to a raw hardware-based
        // time that is not subject to NTP adjustments.

        clock_gettime(CLOCK_MONOTONIC_RAW, &timespec);

        time =
            (double)timespec.tv_sec +
            (double)timespec.tv_nsec / NSECS_PER_SEC;
#else
        struct timeval timeval;

        gettimeofday(&timeval, NULL);

        time =
            (double)timeval.tv_sec +
            (double)timeval.tv_usec / USECS_PER_SEC;
#endif

        return time;
    }

    double _rusage()
    {
        struct rusage usage;
        double usertime;
        double systime;

        getrusage(RUSAGE_SELF, &usage);

        usertime =
            (double) usage.ru_utime.tv_sec +
            (double) usage.ru_utime.tv_usec / USECS_PER_SEC;

        systime =
            (double) usage.ru_stime.tv_sec +
            (double) usage.ru_stime.tv_usec / USECS_PER_SEC;

        return usertime + systime;
    }

    bool running;

public:
    Stopwatch(const char *context = NULL, FILE *_stream = stdout) : stream(_stream) {
        start(context);
    }

    ~Stopwatch() {}

    void start(const char *context = NULL)
    {
        if (context && strlen(context) > 0)
            fprintf(stream, "%s\n", context);

        running = true;

        wall = _wall();
        rusage = _rusage();
    }

    void stop(const char *context = NULL)
    {
        if (running == false) {
            fprintf(stdout, "Error: Stopwatch is not running\n");
            return;
        }

        wall = _wall() - wall;
        rusage = _rusage() - rusage;

        if (context && strlen(context) > 0)
            fprintf(stream, "%s\n", context);

        fprintf(stream,
                "wall: %.8f seconds\n"
                "usr+sys: %.8f seconds\n",
                wall,
                rusage);

        running = false;
    }
};

size_t ceil_log2(size_t);

void rcomplement(uint8_t*, size_t);

void vmpeak(FILE*);

#endif /* COMMON_H */
