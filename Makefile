
FC = pgfortran
FIX = -DSERIAL 
ifdef MPI
 FC = mpif90
 FIX = 
endif

ifdef MPI
 FC = mpif90
endif

CUDAFLAGS = -cuda -gpu=cc70 -O0 -g -Mbounds -Mchkptr -Mchkstk
CUDAFLAGS = -cuda -gpu=cc70,keepptx -O0 -g
CUDAFLAGS = -cuda -fast -gpu=cc70,cuda11.0,lineinfo -Minfo=accel
CUDAFLAGS = -cuda -gpu=cc70 -fast


lbCUDA: main.o dimensions_m.o kernels_fluid.o kernels_fluid_bc.o kernels_fluid_subs.o kernels_particles.o write_output.o setupMPI.o
	$(FC) $(FIX) $(CUDAFLAGS) $(F90FLAGS) -o lbCUDA main.o dimensions_m.o kernels_fluid.o write_output.o setupMPI.o \
		kernels_fluid_bc.o kernels_fluid_subs.o kernels_particles.o

main.o: dimensions_m.mod kernels_fluid.o kernels_fluid_bc.o kernels_fluid_subs.o kernels_particles.o write_output.o setupMPI.o main.CUF 
	$(FC) $(FIX) $(CUDAFLAGS) $(F90FLAGS) -c main.CUF

dimensions_m.o: defines.h dimensions_m.CUF
	$(FC) $(FIX) $(CUDAFLAGS) $(F90FLAGS) -c dimensions_m.CUF

dimensions_m.mod: defines.h dimensions_m.CUF
	$(FC) $(FIX) $(CUDAFLAGS) $(F90FLAGS) -c dimensions_m.CUF

setupMPI.o: dimensions_m.mod setupMPI.CUF
	$(FC) $(FIX) $(CUDAFLAGS) $(F90FLAGS) -c setupMPI.CUF

kernels_fluid_bc.o: dimensions_m.mod kernels_fluid_bc.CUF
	$(FC) $(FIX) $(CUDAFLAGS) $(F90FLAGS) -c kernels_fluid_bc.CUF

kernels_fluid_subs.o: dimensions_m.mod kernels_fluid_subs.CUF
	$(FC) $(FIX) $(CUDAFLAGS) $(F90FLAGS) -c kernels_fluid_subs.CUF

kernels_particles.o: dimensions_m.mod kernels_particles.CUF
	$(FC) $(FIX) $(CUDAFLAGS) $(F90FLAGS) -c kernels_particles.CUF

kernels_fluid.o: dimensions_m.mod kernels_fluid.CUF
	$(FC) $(FIX) $(CUDAFLAGS) $(F90FLAGS) -c kernels_fluid.CUF

write_output.o: dimensions_m.mod write_output.f90
	$(FC) $(FIX) $(CUDAFLAGS) $(F90FLAGS) -c write_output.f90

clean:
	@rm -rf lbCUDA *.o *.mod *.ptx

