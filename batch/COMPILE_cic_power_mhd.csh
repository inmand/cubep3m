cd ../utils/cic_power

rm -f cic_power_mhd
rm -f ngp_power_mhd

FFTWLIB=/scinet/gpc/lib/fftw/intel-intelmpi-3.2.2-medium/lib
FFTWINC=-I/scinet/gpc/lib/fftw/intel-intelmpi-3.2.2-medium/include

mpif77 -shared-intel -fpp -g -O3 -xhost -mcmodel=medium -DBINARY -DPPINT -DNGP cic_power_mhd.f90 -o ngp_power_mhd -L$FFTWLIB -I$FFTWINC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
mpif77 -shared-intel -fpp -g -O3 -xhost -mcmodel=medium -DBINARY -DDEBUG -DNGP cic_init_power_mhd.f90 -o ngp_init_power_mhd -L$FFTWLIB -I$FFTWINC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl

cd ../../batch/

echo "Sourced COMPILE_cic_power_mhd.csh"

