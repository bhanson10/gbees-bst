# makefile, https://github.com/bhanson10/gbees/tree/main/examples/PCR3BP
# Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

CC=gcc

CFLAGS=-pedantic -Wall -std=c99 -D_GNU_SOURCE

# remove -g for release version, and add -DRELEASE -fPIC
LDFLAGS=
LDLIBS=-lm

SRCS = $(wildcard *.c) $(wildcard ../../*.c)
OBJS = $(SRCS:%.c=%.o)

TARGET = PCR3BP

.DEFAULT: all

all: link folders

folders:
	mkdir -p results/P0 results/P1 results/P2 results/P3

link: ${OBJS}
	${CC} ${LDFLAGS} -o ${TARGET} ${OBJS} ${LDLIBS}

%.c.o: %.c
	${CC} ${CFLAGS} -c $< -o $@

clean:
	find ../.. -type f -name '*.o' -delete
	rm -f ${TARGET}