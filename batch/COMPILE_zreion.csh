cd ../utils/zreion

rm -f ngp_iondensity
rm -f ngp_iondensity_mem

FFTWLIB=/scinet/gpc/lib/fftw/intel-intelmpi-3.2.2-medium/lib
FFTWINC=-I/scinet/gpc/lib/fftw/intel-intelmpi-3.2.2-medium/include

mpif77 -shared-intel -fpp -g -O3 -mcmodel=medium -DBINARY -DPPINT -DNGP -DREDUCE_OUTPUT_RES -mt_mpi iondensityfield.f90 -o ngp_iondensity -L$FFTWLIB -I$FFTWINC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
mpif77 -shared-intel -fpp -g -O3 -mcmodel=medium -DBINARY -DPPINT -DNGP -mt_mpi iondensityfield_mem.f90 -o ngp_iondensity_mem -L$FFTWLIB -I$FFTWINC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl

cd ../../batch/

echo "Sourced COMPILE_zreion.csh" 

