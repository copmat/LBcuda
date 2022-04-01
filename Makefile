FC=mpif90 -DSERIAL
FC=mpif90

CUDAFLAGS = -cuda -gpu=cc70 -O0 -g -Mbounds -Mchkptr -Mchkstk
CUDAFLAGS = -cuda -gpu=cc70,keepptx -O0 -g
CUDAFLAGS = -cuda -fast -gpu=cc80 -DMYDIMESION=256 -DTILE1=8 -DTILE2=4 -DTILE3=4
CUDAFLAGS = -cuda -fast -gpu=cc70 -DMYDIMESION=256 -DTILE1=8 -DTILE2=4 -DTILE3=4
CUDAFLAGS = -cuda -fast -gpu=cc70 -DMYDIMESION=256 


lbCUDA: main.o dimensions_m.o profiling_m.o kernels_fluid.o kernels_fluid_CG.o kernels_fluid_PART.o write_output.o setupMPI.o
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -o $@ main.o dimensions_m.o profiling_m.o kernels_fluid.o write_output.o kernels_fluid_CG.o kernels_fluid_PART.o setupMPI.o

main.o: Makefile dimensions_m.mod profiling_m.o kernels_fluid.o kernels_fluid_CG.o kernels_fluid_PART.o write_output.o main.CUF 
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c main.CUF

dimensions_m.o: Makefile defines.h dimensions_m.CUF
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c dimensions_m.CUF

dimensions_m.mod: Makefile defines.h dimensions_m.CUF
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c dimensions_m.CUF

profiling_m.o: Makefile profiling_m.CUF
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c profiling_m.CUF

kernels_fluid_PART.o: Makefile dimensions_m.mod kernels_fluid_PART.CUF
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c kernels_fluid_PART.CUF

kernels_fluid_CG.o: Makefile dimensions_m.mod kernels_fluid_CG.CUF
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c kernels_fluid_CG.CUF

kernels_fluid.o: Makefile dimensions_m.mod kernels_fluid.CUF
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c kernels_fluid.CUF

write_output.o: Makefile dimensions_m.mod write_output.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c write_output.f90

setupMPI.o: Makefile dimensions_m.mod setupMPI.F90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c setupMPI.F90

clean:
	@rm -rf lbCUDA *.o *.mod *.ptx

