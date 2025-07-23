CC = gcc -g -Wall -O2
LIBS = -lm

all: radial-test

RADIAL_TEST = test.o radial.o

radial-test: $(RADIAL_TEST)
	$(CC) $(RADIAL_TEST) $(LIBS) -o radial-test

.c.o:
	$(CC) -c $<

clean:
	rm -f *.o *.c~ *.h~
