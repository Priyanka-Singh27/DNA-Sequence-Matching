CC = gcc
CFLAGS = -O2 -Wall -std=c11 -fPIC
LDFLAGS = -shared -lm

SRCS = dna_utils.c \
       needleman_wunsch.c \
       suffix_tree.c \
       trie_errors.c \
       segment_tree.c \
       interval_tree.c \
       union_find.c

OBJS = $(SRCS:.c=.o)

all: libdna.so

libdna.so: $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: main_test.exe
	./main_test.exe

main_test.exe: $(SRCS) main_test.c
	$(CC) $(CFLAGS:-fPIC=) $^ -o $@ -lm

clean:
	rm -f *.o libdna.so main_test.exe test_st test_te test_uf test_it test_sgt
