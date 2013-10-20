cd ../utils/dist_init_ng

rm -f dist_init_ng

mpif77 -shared-intel -fpp -g -O3 -xhost -DBINARY -mt_mpi dist_init_ng.f90 -o dist_init_ng  -L$SCINET_FFTW_LIB -I$SCINET_FFTW_INC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl 

cd ../../batch

echo "Sourced dist_init_ng.csh"

