CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors -lm -g

spkmeans: spkmeans.o utils.o
	$(CC) -o spkmeans spkmeans.o utils.o $(CFLAGS)
	python3 setup.py build_ext --inplace

spkmeans.o: spkmeans.c
	$(CC) -c spkmeans.c utils.h $(CFLAGS)

utils.o: utils.c
	$(CC) -c utils.c utils.h $(CFLAGS)

clean:
	rm -f *.o *.gch spkmeans *.so
	rm -rf build

.PHONY: test
test:
	python3 test.py

.PHONY: spk1
spk1:
	python3 spkmeans.py spk test/test_batch/test1.txt