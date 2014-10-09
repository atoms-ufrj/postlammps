FORT    = gfortran
FOPTS   = -O3 -march=native -ffast-math -funroll-loops -fstrict-aliasing -cpp -Wunused

all: post_lammps

clean:
	rm -rf *.o
	rm -rf *.mod
	rm -rf post_lammps

install:
	cp -f ./post_lammps /usr/local/bin/

post_lammps: post_lammps.f90 mData_Proc.o mProp_List.o mProp_List.o mString.o
	$(FORT) $(FOPTS) -o $@ $^

mData_Proc.o: mData_Proc.f90 mProp_List.o
	$(FORT) $(FOPTS) -c -o $@ $<

mProp_List.o: mProp_List.f90 mConstants.o
	$(FORT) $(FOPTS) -c -o $@ $<

mString.o: mString.f90 mConstants.o
	$(FORT) $(FOPTS) -c -o $@ $<

mConstants.o: mConstants.f90
	$(FORT) $(FOPTS) -c -o $@ $<

