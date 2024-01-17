#Macro definition
FC = ifort -traceback
FFLAGS =  -O3 -qopenmp -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
LDFLAGS =  $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
           $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a \
           -Wl,--start-group \
           $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
           $(MKLROOT)/lib/intel64/libmkl_core.a \
           $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
           -Wl,--end-group -lpthread -lm
OBJ = opt.o minpack.o sub.o orient.o geomtrans.o diab.o hddata.o pes.o nadvibs.o
#end of Macro definition

dat.x: $(OBJ) main.o
	$(FC) $(FFLAGS) $(OBJ) main.o -o dat.x $(LDFLAGS)
findcp.x: $(OBJ) opttools.o findcp.o
	$(FC) $(FFLAGS) $(OBJ) opttools.o findcp.o -o findcp.x $(LDFLAGS)
findmex.x: $(OBJ) opttools.o findmex.o
	$(FC) $(FFLAGS) $(OBJ) opttools.o findmex.o -o findmex.x $(LDFLAGS)
nadvibs.x: $(OBJ) nadvibsinput.o
	$(FC) $(FFLAGS) $(OBJ) nadvibsinput.o -o nadvibs.x $(LDFLAGS)
basis.x: $(OBJ) basis.o
	$(FC) $(FFLAGS) $(OBJ) basis.o -o basis.x $(LDFLAGS)
gf.x: $(OBJ) gf.o
	$(FC) $(FFLAGS) $(OBJ) gf.o -o gf.x $(LDFLAGS)
traj.x: $(OBJ) traj.o
	$(FC) $(FFLAGS) $(OBJ) traj.o -o traj.x $(LDFLAGS)

clean:
	rm -f *.o *.mod a.out *.x

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
