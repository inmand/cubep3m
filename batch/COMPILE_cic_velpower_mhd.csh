cd ../utils/cic_velpower

rm -f cic_velpower_mhd
rm -f ngp_velpower_mhd
rm -f ngp_velpower_mhd_init
rm -f cic_velpower_mhd_init

FFTWLIB=/scinet/gpc/lib/fftw/intel-intelmpi-3.2.2-medium/lib
FFTWINC=-I/scinet/gpc/lib/fftw/intel-intelmpi-3.2.2-medium/include

mpif77 -shared-intel -fpp -g -O3 -mcmodel=medium -xhost -DBINARY -DPPINT -DNGP cic_velpower_mhd.f90 -o ngp_velpower_mhd -L$FFTWLIB -I$FFTWINC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl

mpif77 -shared-intel -fpp -g -O3 -mcmodel=medium -xhost -DBINARY -DPPINT -DNGP cic_velpower_mhd_init.f90 -o ngp_velpower_mhd_init -L$FFTWLIB -I$FFTWINC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl

cd ../../batch/

echo "Sourced COMPILE_cic_velpower_mhd.csh"

