#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "common.h"

const uint8_t nucleotide_to_int['T' + 1] = {
    9,9,9,9,9,9,9,9,9,9,
    9,9,9,9,9,9,9,9,9,9,
    9,9,9,9,9,9,9,9,9,9,
    9,9,9,9,9,9,9,9,9,9,
    9,9,9,9,9,9,9,9,9,9,
    9,9,9,9,9,9,9,9,9,9,
    9,9,9,9,9,0/*A*/,9,1/*C*/,9,9,
    9,2/*G*/,9,9,9,9,9,9,9,9,
    9,9,9,9,3/*T*/
};

const uint8_t int_to_nucleotide[4] = {'A', 'C', 'G', 'T'};

size_t ceil_log2(size_t v) {
    assert(v > 0);
    return (size_t)ceil(log(v)/log(2));
}

void rcomplement(uint8_t *read, size_t n) {
    assert(read);
    assert(n > 0);

    size_t i;

    for (i = 0; i < n/2; i++) {
        uint8_t tmp = read[i];
        read[i] = read[n - i - 1];
        read[n - i - 1] = tmp;
    }

    for (i = 0; i < n; i++) {
        if (read[i] == 'A')
            read[i] = 'T';
        else if (read[i] == 'T')
            read[i] = 'A';
        else if (read[i] == 'G')
            read[i] = 'C';
        else if (read[i] == 'C')
            read[i] = 'G';
    }
}

void vmpeak(FILE *stream) {
    assert(stream);

    pid_t pid = getpid();

    char buf[BUFSIZ];
    sprintf(buf, "/proc/%d/status", pid);

    FILE *fp = fopen(buf, "r");
    assert(fp);

    fread(buf, 1, BUFSIZ, fp);

    char *start = strstr(buf, "VmHWM");
    assert(start);

    char *end = strstr(start, "\n");
    assert(end);

    *end = '\0';

    fprintf(stream, "%s", start);

    fclose(fp);
}
