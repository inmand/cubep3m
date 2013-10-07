cd ../utils/cic_power

rm -f cic_power
rm -f ngp_power
rm -f cic_init_power
rm -f ngp_init_power
rm -f ngp_power_mem
rm -f cic_power_mem

FFTWLIB=/scinet/gpc/lib/fftw/intel-intelmpi-3.2.2-medium/lib
FFTWINC=-I/scinet/gpc/lib/fftw/intel-intelmpi-3.2.2-medium/include

mpif77 -shared-intel -fpp -g -O3 -DBINARY -DPPINT -DNGP -mt_mpi cic_power.f90 -o ngp_power  -L$FFTWLIB -I$FFTWINC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl

mpif77 -shared-intel -fpp -g -O3 -DBINARY -DDEBUG -DNGP -mt_mpi cic_init_power.f90 -o ngp_init_power  -L$SCINET_FFTW_LIB -I$SCINET_FFTW_INC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl

mpif77 -shared-intel -fpp -g -O3 -DBINARY -DPPINT -DNGP -mt_mpi cic_power_mem.f90 -o ngp_power_mem -L$FFTWLIB -I$FFTWINC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl

cd ../../batch/

echo "Sourced COMPILE_cic_power.csh"

