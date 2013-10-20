# Load proper modules

#module clear
#module unload fftw
#module load intel #/intel-10.0b.017                                     
#module load lam #/lam-7.1.3-intel
#module load fftw #/2.1.5-intel10

# Could add -DKAISER or -DNGP

cd ../utils/cic_power

rm -f halo_power

FFTWLIB=/scinet/gpc/lib/fftw/intel-intelmpi-3.2.2-medium/lib
FFTWINC=-I/scinet/gpc/lib/fftw/intel-intelmpi-3.2.2-medium/include


#mpif77 -shared-intel -fpp -g -O3 -xT -DBINARY -DPPINT -mt_mpi cic_power.f90 -o cic_power  -L$SCINET_FFTW_LIB -I$SCINET_FFTW_INC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
mpif77 -shared-intel -fpp -g -O3 -xT -DBINARY -DPPINT -DNGP -mt_mpi halo_power.f90 -o halo_power  -L$FFTWLIB -I$FFTWINC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
#mpif77 -shared-intel -fpp -g -O3 -xT -DBINARY -DPPINT -DNGP -mt_mpi halo_power.f90 -o halo_power  -L$SCINET_FFTW_LIB -I$SCINET_FFTW_INC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
#mpif77 -shared-intel -fpp -g -O3 -xT -DBINARY -DDEBUG  -mt_mpi cic_init_power.f90 -o cic_init_power  -L$SCINET_FFTW_LIB -I$SCINET_FFTW_INC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
#mpif77 -shared-intel -fpp -g -O3 -xT -DBINARY -DDEBUG -DNGP -mt_mpi cic_init_power.f90 -o ngp_init_power  -L$SCINET_FFTW_LIB -I$SCINET_FFTW_INC -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
cd ../../batch/

echo "Ready to run" 
