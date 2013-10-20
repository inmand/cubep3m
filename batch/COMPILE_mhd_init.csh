cd ../utils/mhd_init

rm -f mhd_init
rm -f mhd_init_nu

mpif77 -shared-intel -fpp -g -O3 -xhost -DBINARY -mt_mpi mhd_init.f90 -o mhd_init  -lm -ldl 
mpif77 -shared-intel -fpp -g -03 -xhost -DBINARY -mt_mpi mhd_init_nu.f90 -o mhd_init_nu -lm -ldl
cd ../../batch

echo "Sourced mhd_init.csh"

