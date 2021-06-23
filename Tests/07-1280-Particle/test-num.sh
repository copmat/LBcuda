#!/bin/bash

echo "Case no Particles"
cat <<EOF > input.dat
50            !niter  number of iterations
50               !stat print on screen every n
10000           !vtk print every n
0.003 0.003 0.003	! vx,vy,vz
NO 0		! with PARTICLES
0.0 0.0 0.0	! ext_fxx,ext_fyy,ext_fzz
YES		! with ROTATION
0.0 0.0 0.0	! ext_tqx,ext_tqy,ext_tqz
EOF

./lbCUDA >& log.nopart.txt


for i in 1 100 200 300 400
do
echo "Case $i"
cat <<EOF > input.dat
50            !niter  number of iterations
50               !stat print on screen every n
10000           !vtk print every n
0.003 0.003 0.003	! vx,vy,vz
YES $i		! with PARTICLES
0.0 0.0 0.0	! ext_fxx,ext_fyy,ext_fzz
YES		! with ROTATION
0.0 0.0 0.0	! ext_tqx,ext_tqy,ext_tqz
EOF

./lbCUDA >& log.part-$i.txt
done
