#!/usr/bin/python
import numpy as np

dir = '/Users/derekinman/Documents/Research/CUBEP3M/NuNBodyTesting'
zstr = '0.0'
dm = '/dmps/'
nu = '/nups/'

copy = '/copy_dm_dm'
nonc = '/noncopy_dm_dm'
neut = '/neut'

omega_nu = 0.001
omega_l = 0.76
omega_m= 1.0 - omega_l

copy_dmps = np.genfromtxt(dir+copy+dm+zstr+'00ngpps.dat')
k = dmps[:,0]
copy_dmps = copy_dmps[:,1]

copy_nups = np.genfromtxt(dir+copy+nu+zstr+'00ngpps.dat')
for i in range(k.size):
    if k[i]!=copy_nups[i,0]:
        print "ERROR"
        
copy_nups = copy_nups[:,1]

nonc_dmps = np.genfromtxt(dir+nonc+dm+zstr+'00ngpps.dat')
nonc_dmps = nonc_dmps[:,1]

nonc_nups = np.genfromtxt(dir+nonc+nu+zstr+'00ngpps.dat')
for i in range(k.size):
    if k[i]!=nonc_nups[i,0]:
        print "ERROR"
        
nups = nups[:,1]

neut_dmps = np.genfromtxt(dir+neut+dm+zstr+'00ngpps.dat')
neut_dmps = neut_dmps[:,1]

neut_nups = np.genfromtxt(dir+neut+nu+zstr+'00ngpps.dat')
for i in range(k.size):
    if k[i]!=neut_nups[i,0]:
        print "ERROR"
        
neut_nups = neut_nups[:,1]


thePlot =   list_plot(zip(k,copy_dmps),plotjoined=true,color='blue',linestyle='-',legend_label='cdm - copy ident') +\
            list_plot(zip(k,copy_nups),plotjoined=true,color='blue',linestyle='--',legend_label=r'$\nu$ - copy ident') +\
            list_plot(zip(k,nonc_dmps),plotjoined=true,color='red',linestyle='-',legend_label='cdm - copy ident') +\
            list_plot(zip(k,nonc_nups),plotjoined=true,color='red',linestyle='--',legend_label=r'$\nu$ - copy ident') +\
            list_plot(zip(k,neut_dmps),plotjoined=true,color='green',linestyle='-',legend_label='cdm - copy ident') +\
            list_plot(zip(k,neut_nups),plotjoined=true,color='green',linestyle='--',legend_label=r'$\nu$ - copy ident')
            
thePlot.axes_labels(['k (Mpc/h)',r'$\Delta(z='+zstr+')$'])
thePlot.show(scale='loglog')

tfn = np.genfromtxt(dir+'/tfn/camb.dat')
tfPlot =    list_plot(zip(tfn[:,0],tfn[:,1]),plotjoined=true,linestyle='-',legend_label='cdm tfn') +\
            list_plot(zip(tfn[:,0],tfn[:,5]),plotjoined=true,linestyle='--',legend_label='nu tfn') 
tfPlot.axes_labels(['k (Mpc/h)',r'TF$'])
tfPlot.show(scale='loglog')