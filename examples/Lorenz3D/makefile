# makefile, https://github.com/bhanson10/gbees/tree/main/examples/Lorenz3D
# Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

CC=gcc

# DEBUG
#CFLAGS=-pedantic -Wall -std=c99 -D_GNU_SOURCE -static -g

# RELEASE
CFLAGS=-pedantic -Wall -std=c99 -D_GNU_SOURCE -fPIC

LDFLAGS=
LDLIBS=-lm

SRCS = $(wildcard *.c) $(wildcard ../../*.c)
OBJS = $(SRCS:%.c=%.o)

TARGET = Lorenz3D

.DEFAULT: all

all: link folders

folders:
	mkdir -p results/P0 results/P1	

link: ${OBJS}
	${CC} ${LDFLAGS} -o ${TARGET} ${OBJS} ${LDLIBS}

%.c.o: %.c
	${CC} ${CFLAGS} -c $< -o $@

clean:
	find ../.. -type f -name '*.o' -delete
	rm -f ${TARGET}