# -*- coding: utf-8 -*-
import sys
import scipy
import scipy.signal
import numpy as np
import matplotlib.pyplot as plt
import pylab
from sys import exit
from scipy.optimize import curve_fit
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import scipy.integrate as integrate
import scipy.special as special
from scipy.integrate import quad
import csv
import astropy
import os
#import sh
#import pexect
from astropy.cosmology import WMAP9 as cosmo
plt.rcParams['figure.figsize'] = (20.0, 6.0)



def MergeBins(x0edge,y0,x1edge):
    
    #print "len(x1edge)",len(x1edge)
    #print "len(x0edge)",len(x0edge),"leny0",len(y0)
    y1 =np.zeros(len(x1edge)-1)
    if len(x0edge) == len(y0)+1:
        x0_s = x0edge[1:]
        x0_i = x0edge[:len(x0edge)-1]
        x0_av = np.array(x0_s + x0_i)*0.5
    elif(len(x0edge)==len(y0)):
        #print "elif"
        dimIN = np.shape(x0edge)
        lenF = 1
        checkV = 0
        for j in range(len(dimIN)):
            if (dimIN[j]>1):
                checkV = checkV+1
            lenF = dimIN[j]*lenF
        if checkV!=1:
            print "MergeBins function: Dimensional ERROR"
            exit()
        x0edge = np.reshape(x0edge,(lenF))
        indexsort = np.argsort(x0edge)
        x0_av = np.sort(x0edge)
        #print x0_av
        #exit()
        y0 = np.reshape(y0,(lenF))
        y0 = y0[indexsort]
    else:
        print "ERROR: inconsistent dimensions of inputs"
        exit()
    x1_s = x1edge[1:]
    x1_i = x1edge[:len(x1edge)-1]
    

    if len(x1edge)>len(x0edge):
        print "ERROR: cannot merge in a larger number of bins"
        exit()
    sum_to_TOT = 0.0
    sum_to_TOT_l = 0.0
    sum_to_TOT_h = 0.0
    ii1 = 0
    n_ii1 = 0.0
    if((x0_av !=np.sort(x0_av)).all() ):
        print "ERROR: bins not sorted"
        print x0edge
        print np.sort(x0edge)
        exit()
    for hhf in range(len(y0)):
        #print "ii1",ii1,"hhf",hhf
        #print "x0_av[hhf]",x0_av[hhf]
        #print "x1_i[ii1]",x1_i[ii1],"x1_s[ii1]",x1_s[ii1]
        #print "x1_i[ii1]",x1_i[ii1],"x0_av[hhf]",x0_av[hhf],"x1_s[ii1]",x1_s[ii1]
        if(x0_av[hhf]<x1_i[0]):
            #print "fuori inf"
            sum_to_TOT_l = sum_to_TOT_l+ y0[hhf]
            continue
        elif(x0_av[hhf]>=x1_i[ii1] and x0_av[hhf]<x1_s[ii1]):
            #print "adding"
            n_ii1 = n_ii1 + y0[hhf]
            #print "x0_av[hhf]",x0_av[hhf]
            #print "x1_s[ii1]",x1_s[ii1]
            #print "x1_i[ii1]",x1_i[ii1]
            #print "ii1",ii1,"len(x1_i)-1",len(x1_i)-1
            #print "hhf",hhf
            #print
        elif(x0_av[hhf]>=x1_s[ii1] and ii1<len(x1_i)-1):
            #print "saving"
            #y1.append(n_ii1)
            y1[ii1] = n_ii1
            ii1=ii1+1
            while(x0_av[hhf]>=x1_s[ii1] and ii1<len(x1_s)-1):
                ii1 = ii1 +1
            
            n_ii1 = y0[hhf]
            #ii1 = ii1+1
            #print
        elif(x0_av[hhf]>=x1_s[ii1] and ii1>=len(x1_i)-1):
            #print "fuori out"
            sum_to_TOT_h = sum_to_TOT_h+ y0[hhf]
            continue
        else:
            print "ERROR: should not happen"
            print "x0_av[hhf]",x0_av[hhf]
            print "x1_s[ii1]",x1_s[ii1]
            print "x1_i[ii1]",x1_i[ii1]
            print "ii1",ii1,"len(x1_i)-1",len(x1_i)-1
            print "hhf",hhf
            exit()
    #y1.append(n_ii1)
    y1[ii1] = n_ii1
    #if (y1[ii1] == 0.0):
    #    y1[ii1] = n_ii1
    #else:
    #    ii1 = ii1 +1
    #    y1[ii1] = n_ii1
    #y1 = np.array(y1)
    sum_to_TOT = sum_to_TOT_l + sum_to_TOT_h
    #print "len y1",len(y1),"len x1",len(x1edge)
    #print "sum_to_TOT_l",sum_to_TOT_l
    #print "sum_to_TOT_h",sum_to_TOT_h
    #check
    onesA = np.ones(len(y0))
    sum_y0 = np.matrix(y0)*np.transpose(np.matrix(onesA))
    
    onesA1 = np.ones(len(y1))
    sum_y1 = np.matrix(y1)*np.transpose(np.matrix(onesA1)) + sum_to_TOT
    
    if np.abs(sum_y0 - sum_y1)>0.000001:
        print "sum_y0",sum_y0,"sum_y1",sum_y1
        print "ERROR: sums must be equal"
        exit()

    #print "len(y1)",len(y1), "y1",y1
    return y1


def ChooseTABLE_z_V(MODnamef):
    OmegaMf = 0.3 # Omega M
    OmegaLf = 0.7 # Omega vac
    OmegaKf = 1. - OmegaMf - OmegaLf
    Hcf = 70.0
    ######### TABLE V - Z ###########
    print "      # DEFAULT COSMOLOGY: OmegaM = 0.3, OmegaL = 0.7, Hubble constant = 70  (km/s)/Mpc"
    print "      # CHOSEN TABLE redshift(z) - Comoving Volume (CV): VOLUME_table_0p004_Om03_H70.txt                         #"
    print "      #    > step in z is: dz = step in z from the RATE                                                          #"
    print "      # Enter 0 if that is ok                                                                                    #"
    print "      # Enter 1 if you want the more refined table VOLUME_table_em6_Om03_H70.txt                                 #"
    print "      # Enter 2 if you want a more refined table to be built                                                     #"
    print "      # Enter 3 if you want to change cosmology and redefined the table                                          #"
    print "      # NB:  OPTION 3 AVAILABLE ONLY FOR User and Mandel rate models for concistency                             #"
    invalid_option = True
    NameFILE_TABf = ""
    while invalid_option:
            OPT = input("      # Enter the number: ")
            if(OPT==0):
                invalid_option = False
                NameFILE_TABf = "VOLUME_table_0p004_Om03_H70.txt"
                TABzVCf=prefix+"z_CV_TABLES/" +NameFILE_TABf
            #sh.cp(TABzVC,"TEMP/.")
            elif(OPT==1):
                invalid_option = False
                NameFILE_TABf = "VOLUME_table_em6_Om03_H70.txt"
                TABzVCf=prefix+"z_CV_TABLES/"+ NameFILE_TABf
            elif(OPT==2):
                invalid_option = False
                dz = input("Enter desired resolution in redshift dz: ")
                
                from astropy.cosmology import FlatLambdaCDM
                Vp = cosmo.comoving_volume(1)
                cosmo_FL = FlatLambdaCDM(H0=70, Om0=0.3)
                
                MINz = 0.0
                MAXz = np.max(zR)
                
                ROUGH_dV_Om03_H70 = []
                ROUGH_Vcum_Om03_H70 =[]
                Z_Vec = []
                
                zi = MINz
                V_old = cosmo_FL.comoving_volume(zi)
                if (zi>0.0000):
                    Z_Vec.append(zi)
                    ROUGH_dV_Om03_H70.append(V_old.value*Mpc**3)
                    ROUGH_Vcum_Om03_H70.append(V_old.value*Mpc**3)
                while zi < MAXz:
                    zi = zi + dz
                    Z_Vec.append(zi)
                    Vi = cosmo_FL.comoving_volume(zi)
                    DV_i = Vi.value - V_old.value
                    ROUGH_dV_Om03_H70.append(DV_i*Mpc**3)
                    ROUGH_Vcum_Om03_H70.append(Vi.value*Mpc**3)
                    V_old = Vi
                
                print len(ROUGH_dV_Om03_H70)
                print len(ROUGH_Vcum_Om03_H70)
                print len(Z_Vec)
                NameFILE_TABf = "VOLUME_table_%.2E_Om03_H70.txt" % dz
                TABzVCf=prefix+"z_CV_TABLES/" +NameFILE_TABf
                print TABzVCf
                with open(TABzVCf, 'w') as f:
                    writer = csv.writer(f, delimiter='\t')
                    writer.writerows(zip(Z_Vec,ROUGH_dV_Om03_H70,ROUGH_Vcum_Om03_H70))
            elif(OPT==3):
                if (MODnamef!="User" and MODnamef!="Mandel2017"):
                    print "      # Invalid option for rate ",MODnamef
                    continue
                invalid_option = False
                invalidCOSM = True
                while invalidCOSM:
                    Ode0 = input ('# Cosmological constant: ')
                    Om0 = input ('# Matter constant: ')
                    OmegaKf = input ('# Curvature : ')
                    H0 = input ('# Hubble constant [km/s/Mpc]: ')
                    if (np.abs(Ode0+Om0+OmegaKf -1.)<0.000001):
                        invalidCOSM = False
                    else:
                        print "      # INVALID COSMSOLOGY, constants not sum to 1"
                if OmegaKf==0:
                    from astropy.cosmology import FlatLambdaCDM

                    Vp = cosmo.comoving_volume(1)
                    cosmo_FL = FlatLambdaCDM(H0, Om0)
                    dz = input('# Enter desired resolution in redshift dz: ')
                    MINz = 0.0
                    MAXz = input('# Enter maximum redshift for table building (please consider the redshifts required by the rate model): ')
                    MAXz = MAXz + dz

                    ROUGH_dV = []
                    ROUGH_Vcum =[]
                    Z_Vec = []
                
                    zi = MINz
                    V_old = cosmo_FL.comoving_volume(zi)
                    if (zi>0.0000):
                        Z_Vec.append(zi)
                        ROUGH_dV.append(V_old.value*Mpc**3)
                        ROUGH_Vcum.append(V_old.value*Mpc**3)
                    while zi < MAXz:
                        zi = zi + dz
                        Z_Vec.append(zi)
                        Vi = cosmo_FL.comoving_volume(zi)
                        DV_i = Vi.value - V_old.value
                        ROUGH_dV.append(DV_i*Mpc**3)
                        ROUGH_Vcum.append(Vi.value*Mpc**3)
                        V_old = Vi
                
                    print len(ROUGH_dV)
                    print len(ROUGH_Vcum)
                    print len(Z_Vec)
                    NameFILE_TABf = "VOLUME_table_%.2E_Ode%.2f_Om%.2f_H%2.1f.txt" % (dz, Ode0, Om0, H0)
                    TABzVCf=prefix+"z_CV_TABLES/" +NameFILE_TABf
                    print TABzVCf
                    with open(TABzVCf, 'w') as f:
                        writer = csv.writer(f, delimiter='\t')
                        writer.writerows(zip(Z_Vec,ROUGH_dV,ROUGH_Vcum))
                else:
                    from astropy.cosmology import LambdaCDM
                    Vp = cosmo.comoving_volume(1)
                    cosmo_LCDM = LambdaCDM(H0, Om0, Ode0)
                    dz = input('# Enter desired resolution in redshift dz: ')
                    MINz = 0.0
                    MAXz = input('# Enter maximum redshift for table building (please consider the redshifts required by the rate model): ')
        
                    ROUGH_dV = []
                    ROUGH_Vcum =[]
                    Z_Vec = []
                    
                    zi = MINz
                    V_old = cosmo_LCDM.comoving_volume(zi)
                    if (zi>0.0000):
                        Z_Vec.append(zi)
                        ROUGH_dV.append(V_old.value*Mpc**3)
                        ROUGH_Vcum.append(V_old.value*Mpc**3)
                    while zi < MAXz:
                        zi = zi + dz
                        Z_Vec.append(np.round(zi,2))
                        Vi = cosmo_LCDM.comoving_volume(zi)
                        DV_i = Vi.value - V_old.value
                        ROUGH_dV.append(DV_i*Mpc**3)
                        ROUGH_Vcum.append(Vi.value*Mpc**3)
                        V_old = Vi
            
                    print len(ROUGH_dV)
                    print len(ROUGH_Vcum)
                    print len(Z_Vec)
                    NameFILE_TABf = "VOLUME_table_%.2E_Om%.2f_H%2.1f.txt" % (dz, Om0, H0)
                    TABzVCf=prefix+"z_CV_TABLES/" +NameFILE_TABf
                    print TABzVCf
                    with open(TABzVCf, 'w') as f:
                        writer = csv.writer(f, delimiter='\t')
                        writer.writerows(zip(Z_Vec,ROUGH_dV,ROUGH_Vcum))
                    Hcf = H0
                    OmegaLf = Ode0
                    OmegaMf = Om0

            else:
                print "Invalid OPTION"
    
    print "      ############################################################################################################"
    outTAB_Qf = [NameFILE_TABf,TABzVCf,Hcf,OmegaMf,OmegaLf,OmegaKf]
    return outTAB_Qf

