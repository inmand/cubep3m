cd ../utils/mhd_init

rm -f mhd_init

mpif77 -shared-intel -fpp -g -O3 -xhost -DBINARY -mt_mpi mhd_init.f90 -o mhd_init  -lm -ldl 

cd ../../batch

echo "Sourced mhd_init.csh"

