include Makefile.inc

all: ptrans

LFLAGS =
LIBS =

main.o: main.f90 ptrans.o
	$(FC) -c $(FFLAGS) main.f90

ptrans.o: RTMCphon.f90
	$(FC) -c $(FFLAGS) $< -o $@

ptrans: main.o ptrans.o
	$(FC) main.o ptrans.o $(LFLAGS) $(LIBS) -o $@

clean:
	rm -f ptrans *.o *.mod
