F90 = ifort
FFLAG = -O2
LDFLAGS= -L/home/hikage/work/software/cfitsio

PROG = mfcalc

mfcalc	: mfcalc.o mf_tools.o
	$(F90) $(FFLAG) $^ -o $@ $(LDFLAGS) -lcfitsio

mfcalc.o : mfcalc.f90 mf_tools.o
	$(F90) $(FFLAG) -c $<

test	: 
	./mfcalc temp.fits mask.fits test.dat 25 -3.75 3.75

mf_tools.o : mf_tools.f90
	$(F90) $(FFLAG) -c $<

clean	:
	-rm -f $(PROG) *.o *mod core