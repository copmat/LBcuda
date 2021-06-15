CUDAFLAGS = -cuda -gpu=cc70 -O0 -g -Mbounds -Mchkptr -Mchkstk
CUDAFLAGS = -cuda -gpu=cc70,keepptx -O0 -g
CUDAFLAGS = -cuda -fast -gpu=cc70,cuda11.0,lineinfo -Minfo=accel
CUDAFLAGS = -cuda -gpu=cc70 -fast


lbCUDA: main.o dimensions_m.o kernels_fluid.o kernels_fluid_subs.o write_output.o
	pgfortran $(CUDAFLAGS) $(F90FLAGS) -o $@ main.o dimensions_m.o kernels_fluid.o write_output.o kernels_fluid_subs.o

main.o: dimensions_m.mod kernels_fluid.o kernels_fluid_subs.o write_output.o main.CUF 
	pgfortran $(CUDAFLAGS) $(F90FLAGS) -c main.CUF

dimensions_m.o: defines.h dimensions_m.CUF
	pgfortran $(CUDAFLAGS) $(F90FLAGS) -c dimensions_m.CUF

dimensions_m.mod: defines.h dimensions_m.CUF
	pgfortran $(CUDAFLAGS) $(F90FLAGS) -c dimensions_m.CUF

kernels_fluid_subs.o: dimensions_m.mod kernels_fluid_subs.CUF
	pgfortran $(CUDAFLAGS) $(F90FLAGS) -c kernels_fluid_subs.CUF

kernels_fluid.o: dimensions_m.mod kernels_fluid.CUF
	pgfortran $(CUDAFLAGS) $(F90FLAGS) -c kernels_fluid.CUF

write_output.o: dimensions_m.mod write_output.f90
	pgfortran $(CUDAFLAGS) $(F90FLAGS) -c write_output.f90

clean:
	@rm -rf lbCUDA *.o *.mod *.ptx

