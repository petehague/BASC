#
# Makefile for John Skilling's BayeSys3 - works fine on RedHat and Debian GNU/linux.
#                            			David MacKay         Mon 29/3/04
#

tar: ALWAYS
	tar zcvfh BayesSys.tar.gz  \
 bayesys3.c  license.txt  masstoy.c     random.c \
 bayesapp.c    bayesys3.h  manual.ps    masstoy.log   random.h \
 bayestoy.c    hilbert.c   manual.tex   poisstoy.c    userstr.h \
 bayestoy.log  hilbert.h   massapp.c    poisstoy.log  Makefile

ALWAYS:

CC = gcc

all: bayestoy masstoy poisstoy app.o

clean: ALWAYS
	rm *.o

source: bayesys3.o random.o hilbert.o mcmc2.o
	$(CC) bayesys3.o random.o hilbert.o mcmc2.o -lm -o source

bayestoy:     bayesys3.o random.o hilbert.o bayesapp.o  bayestoy.o
	$(CC) bayesys3.o random.o hilbert.o bayesapp.o  bayestoy.o  -lm -o bayestoy

masstoy:      bayesys3.o random.o hilbert.o massapp.o   masstoy.o
	$(CC) bayesys3.o random.o hilbert.o massapp.o   masstoy.o   -lm -o masstoy

poisstoy:     bayesys3.o random.o hilbert.o massapp.o   poisstoy.o
	$(CC) bayesys3.o random.o hilbert.o massapp.o   poisstoy.o  -lm -o poisstoy
