CC = g++ -std=c++11 -march=native #-mfpmath=both

CFLAGS = -Wall -Wextra -Werror -Wshadow -Wcast-qual -Wcast-align -Wwrite-strings -O3 #-Ofast

CPPFLAGS = -DNDEBUG

LDFLAGS = -lrt -lpthread #-flto

INCLUDES =

DEBUG = #-ggdb #-pg

TARGETS = faust-encode faust-decode afin-encode afin-decode

all: $(TARGETS)

faust-encode: faust-encode.o common.o divsufsort.o
	$(CC) $(DEBUG) $(LDFLAGS) -o $@ $^

faust-decode: faust-decode.o common.o
	$(CC) $(DEBUG) $(LDFLAGS) -o $@ $^

afin-encode: afin-encode.o common.o
	$(CC) $(DEBUG) $(LDFLAGS) -o $@ $^

afin-decode: afin-decode.o common.o
	$(CC) $(DEBUG) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	$(CC) $(DEBUG) $(CFLAGS) $(CPPFLAGS) $(INCLUDES) -c $^

clean:
	-$(RM) -f *.o $(TARGETS) >/dev/null 2>&1
