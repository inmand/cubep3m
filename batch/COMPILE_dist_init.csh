cd ../utils/dist_init

rm -f dist_init
rm -f dist_init_hybrid
rm -f dist_init_nu

# WRITE INITIAL DENSITY FIELD:
mpif77 -shared-intel -fpp -g -O3 -xhost -Dwrite_den -DBINARY -mt_mpi dist_init_dm.f90 -o dist_init_dm -L$SCINET_FFTW_LIB -I$SCINET_FFTW_INC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl 
#mpif77 -shared-intel -fpp -g -O3 -xhost -Dwrite_den -DBINARY -mt_mpi dist_init_hybrid.f90 -o dist_init_hybrid -L$SCINET_FFTW_LIB -I$SCINET_FFTW_INC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl 
mpif77 -shared-intel -fpp -g -O3 -xhost -Dwrite_den -DBINARY -mt_mpi dist_init_nu.f90 -o dist_init_nu -L$SCINET_FFTW_LIB -I$SCINET_FFTW_INC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl 

cd ../../batch

echo "Sourced dist_init.csh"

