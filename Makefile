CC = gcc -g -Wall -Werror-implicit-function-declaration -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations
LIBS = -lm

all: radial-test

RADIAL_TEST = test.o radial.o

radial-test: $(RADIAL_TEST)
	$(CC) $(RADIAL_TEST) $(LIBS) -o radial-test

.c.o:
	$(CC) -c $<

clean:
	rm -f *.o *.c~ *.h~
