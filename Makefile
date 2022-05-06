CC=mpicc 
CFLAGS=-O2 -cpp
FC=mpif90 -DSERIAL
FC=mpif90


CUDAFLAGS = -cuda -gpu=cc70 -O0 -g -Mbounds -Mchkptr -Mchkstk
CUDAFLAGS = -cuda -gpu=cc70,keepptx -O0 -g
CUDAFLAGS = -cuda -fast -gpu=cc80 -DMYDIMESION=128 -DTILE1=128 -DTILE2=1 -DTILE3=1
CUDAFLAGS = -cuda -fast -gpu=cc70 -cpp -DMYDIMESION=128 -DTILE1=128 -DTILE2=1 -DTILE3=1



lbCUDA: main.o get_mem.o dimensions_m.o profiling_m.o kernels_fluid.o kernels_fluid_CG.o kernels_fluid_BGK.o kernels_fluid_SCP.o kernels_fluid_PART.o write_output.o setupMPI.o 
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -o $@ main.o get_mem.o dimensions_m.o profiling_m.o kernels_fluid.o write_output.o kernels_fluid_CG.o kernels_fluid_BGK.o kernels_fluid_SCP.o kernels_fluid_PART.o setupMPI.o

main.o: Makefile get_mem.o dimensions_m.mod profiling_m.o kernels_fluid.o kernels_fluid_CG.o kernels_fluid_BGK.o kernels_fluid_SCP.o kernels_fluid_PART.o write_output.o main.f90 
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c main.f90

get_mem.o:get_mem.c
	$(CC) $(CFLAGS) -c get_mem.c

dimensions_m.mod: Makefile defines.h dimensions_m.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c dimensions_m.f90

profiling_m.o: Makefile profiling_m.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c profiling_m.f90

kernels_fluid_PART.o: Makefile dimensions_m.mod kernels_fluid_PART.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c kernels_fluid_PART.f90

kernels_fluid_BGK.o: Makefile dimensions_m.mod kernels_fluid_BGK.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c kernels_fluid_BGK.f90

kernels_fluid_SCP.o: Makefile dimensions_m.mod kernels_fluid_SCP.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c kernels_fluid_SCP.f90

kernels_fluid_CG.o: Makefile dimensions_m.mod kernels_fluid_CG.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c kernels_fluid_CG.f90

kernels_fluid.o: Makefile dimensions_m.mod kernels_fluid.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c kernels_fluid.f90

write_output.o: Makefile dimensions_m.mod write_output.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c write_output.f90

setupMPI.o: Makefile dimensions_m.mod setupMPI.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c setupMPI.f90

clean:
	@rm -rf lbCUDA *.o *.mod *.ptx

