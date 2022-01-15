include Makefile.inc

all: RTMCphon

LFLAGS =
LIBS =

main.o: main.f90 RTMCphon.o
	$(FC) -c $(FFLAGS) main.f90

RTMCphon.o: RTMCphon.f90
	$(FC) -c $(FFLAGS) $< -o $@

RTMCphon: main.o RTMCphon.o
	$(FC) main.o RTMCphon.o $(LFLAGS) $(LIBS) -o $@

clean:
	rm -f RTMCphon *.o *.mod
