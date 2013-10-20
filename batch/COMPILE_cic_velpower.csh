cd ../utils/cic_velpower

rm -f cic_velpower
rm -f ngp_velpower
rm -f ngp_velpower_mem
rm -f cic_velpower_mem
rm -f cic_velpower_cross
rm -f ngp_velpower_cross

FFTWLIB=/scinet/gpc/lib/fftw/intel-intelmpi-3.2.2-medium/lib
FFTWINC=-I/scinet/gpc/lib/fftw/intel-intelmpi-3.2.2-medium/include

mpif77 -shared-intel -fpp -g -O3 -mcmodel=medium -DBINARY -DPPINT -DNGP -mt_mpi cic_velpower.f90 -o ngp_velpower -L$FFTWLIB -I$FFTWINC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
mpif77 -shared-intel -fpp -g -O3 -mcmodel=medium -DBINARY -DPPINT -DNGP -mt_mpi cic_velpower_mem.f90 -o ngp_velpower_mem -L$FFTWLIB -I$FFTWINC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
mpif77 -shared-intel -fpp -g -O3 -mcmodel=medium -DBINARY -DPPINT -DNGP -mt_mpi cic_velpower_cross.f90 -o ngp_velpower_cross -L$FFTWLIB -I$FFTWINC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl

cd ../../batch/

echo "Sourced COMPILE_cic_velpower.csh"

