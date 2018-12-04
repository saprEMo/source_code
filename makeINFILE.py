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
from astropy.io import fits
import os, sys
from astropy import units as u
#from astroquery.xmatch import XMatch
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
from scipy.optimize import curve_fit
from astropy.wcs import WCS

RUNname=sys.argv[1]
indexSTRING = "abcdefghilmnopqrtsuvz"

# PARAMETERS
Nbin_GalacticLAT = 18*2
NmainGalacticLAYERS = 3 #LMGLs
lowerLIM_NMGLs = [0.0,12.5,37.5] # on half plane
upperLIM_NMGLs = [12.5,37.5,90.0]
zlimGhirlanda = 6
zlimMandel = 21
#costGHIRL_a  = -0.8
#costGHIRL_b = 0.0
#costGHIRL_c = -0.2
costGHIRL_a  = 0.2
costGHIRL_b = 0.0
costGHIRL_c = 0.8


# COSTANTS
pc = 3.08567758130573*10**18 #cm
Mpc = 10**6*pc
c = 299792458 #m / s
c_cm = c*100.
G = 6.67408* 10**(-11)# m3 kg^(-1) s-2
G_cgs = G*(10**2)**3/1000.
h = 4.135667662*10**(-15.) #eV*s
eV = 1.60218*10**(-12)#1eV to erg
h = h*eV #erg*s
KB =  1.38064852*10**(-16) #erg/K
zmax = 1089 #recombination age
zmin = 0.0

# FUNCTIONS
#### TABLE z V
def ChooseTABLE_z_V(MODnamef):
    OmegaMf = 0.3 # Omega M
    OmegaLf = 0.7 # Omega vac
    OmegaKf = 1. - OmegaMf - OmegaLf
    Hcf = 70.0
    ######### TABLE V - Z ###########
    print "      # DEFAULT COSMOLOGY: OmegaM = 0.3, OmegaL = 0.7, Hubble constant = 70  (km/s)/Mpc"
    print "      # CHOSEN TABLE redshift(z) - Comoving Volume (CV): VOLUME_table_5.00E-04_Om03_H70.txt                      #"
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
                #NameFILE_TABf = "VOLUME_table_0p004_Om03_H70.txt"
                NameFILE_TABf = "VOLUME_table_5.00E-04_Om03_H70.txt"
                TABzVCf=prefix+"z_CV_TABLES/" +NameFILE_TABf
            #sh.cp(TABzVC,"TEMP/.")
            elif(OPT==1):
                print "WARNING: this choice might considerably increase the computational time required to process the data"
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
                MAXz = input('# Enter maximum redshift for table building (please consider the redshifts required by the rate model): ')
                
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
            
                NameFILE_TABf = "VOLUME_table_%.2E_Om03_H70.txt" % dz
                TABzVCf=prefix+"z_CV_TABLES/" +NameFILE_TABf
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
                
                    NameFILE_TABf = "VOLUME_table_%.2E_Ode%.2f_Om%.2f_H%2.1f.txt" % (dz, Ode0, Om0, H0)
                    TABzVCf=prefix+"z_CV_TABLES/" +NameFILE_TABf
                    
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
            
                    NameFILE_TABf = "VOLUME_table_%.2E_Om%.2f_H%2.1f.txt" % (dz, Om0, H0)
                    TABzVCf=prefix+"z_CV_TABLES/" +NameFILE_TABf
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

#### z ARRAY for analytic functions #####

def zARRAY_function(DeltaZ_tabf):
    ZTab_VALID = True
    while (ZTab_VALID):
            
            print "# From the table you choose a step in redshift to calculate the comuving volume of ",min(DeltaZ_tabf)
            print "# "
            print "# THE MODEL IS ANALYTIC, CHOOSE THE REDSHIFT ARRAY TO USE (which sets the step of integrations)"
            print "# Availble options: "
            print "# 1. define a constant step in redshift"
            print "# 2. use an existing file with values of desired z"
            print "# 3. enter manually the desired z values "
            print "#"
            ZRATE_mod = input('# Enter the number correspondent to the desire option: ')
            print "#"
            print "# BE AWARE: different models has been tested up to different maximum redshift"
            print "#        Ghirlanda_et_al : according to the model z>",zlimGhirlanda, "will not be considered "
            print "#        Mandel 2017: according to the model z> ",zlimMandel, "will not be considered "
            print "#"
            
            
            if (ZRATE_mod == 1):
                ZTab_VALID = False
                Zstep = input('# Enter the desired step in z: ')
                Zmax_input = input('# Enter the maximum redshift to be considered ')
                lenVEC = np.int(Zmax_input/Zstep)
                zMAX_EFF = lenVEC*Zstep
                if (np.abs(zMAX_EFF + Zstep - Zmax_input)<0.00001):
                    lenVEC =lenVEC +1
                    zMAX_EFF = lenVEC*Zstep
                zRf = np.linspace(0.0,zMAX_EFF,lenVEC+1)
            
            elif(ZRATE_mod == 2):
                print "WARNING: the number considered are going to be on the first column of the file"
                dirZfiles = "../INPUTS/zDATA_file"
                print "In ",dirZfiles,":"
                CommLINE = "ls "+dirZfiles
                os.system(CommLINE)
                nameZFILE = raw_input('# Enter file name (with relative path, if necessary): ')
                dataZ = np.genfromtxt(nameZFILE)
                shape_dataZ = np.shape(dataZ)
                if len(shape_dataZ)>1:
                    zRf = np.array(dataZ[:,0])
                else:
                    zRf = dataZ
                ZTab_VALID = False
            
            elif(ZRATE_mod == 3):
                NumZ_input=input('# Enter the total number of z-values you would like to add: ')
                zRf = np.zeros(NumZ_input)
                for bb in range(NumZ_input):
                    zi = input('# Enter the z value: ')
                    zRf[bb] = zi
                ZTab_VALID = False
            else:
                print "ERROR: invalid option, choose among the possible numbers"
    
    if len(zRf)> 300:
            print "# WARNING: high number of redshifts will require a considerable amount of time"
    if len(zRf)> 10000:
            print "# ERROR: number of redshifts too high"
            exit()
    return(zRf)


#### RATE
def GhirlandaRate(p1,zp,p2,zARR_f):
    ourARR = (1. + p1*np.array(zARR_f))/(1. + (np.array(zARR_f)/zp)**p2)
    return ourARR

def SelectingPAR_Ghirlanda(par_casef):
    VALtest = True
    while (VALtest):
        print "#"
        print "# Which value of the distribution do you want to use?"
        print "# A for AVERAGE"
        print "# M for MODE"
        print "# L for LOWER BOUND of 68% confidence interval"
        print "# U for UPPER BOUND of 68% confidence interval"
                    
        ColumVAL = raw_input('Type the letter corrispondent to the desired value: ')
        
        if (ColumVAL =='A' or ColumVAL=='a'):
            VALtest = False
            VALf ='mean'
            COL_VAL = 1
        elif (ColumVAL =='M' or ColumVAL=='m'):
            VALtest = False
            VALf = 'mode'
            COL_VAL = 2
        elif (ColumVAL =='L' or ColumVAL=='l'):
            VALtest = False
            VALf = 'lowerCI'
            COL_VAL = 3
        elif (ColumVAL =='U' or ColumVAL=='u'):
            VALtest = False
            VALf = 'upperCI'
            COL_VAL = 4
        else:
            print " ERROR: invalid choice, choose among the available options"
    p1_f = par_casef[0,COL_VAL]
    zp_f = par_casef[1,COL_VAL]
    p2_f = par_casef[2,COL_VAL]
    outPAR = [p1_f,zp_f,p2_f,VALf]
    return outPAR

#### ABS
# from http://ads.harvard.edu/cgi-bin/nph-bprint page 200

# for 0.1 KeV<E<0.53 keV
def sigmaISMfunc1(Earrayf):
    print "WARNING: energies will be consedered in keV "
    SIGMAarrf = []
    for Eel in np.nditer(Earrayf):
        if Eel >= 0.1 and Eel< 0.53:
            sigma_el = 0.65/1e22 * Eel**(-3.)
        elif Eel >= 0.53 and Eel< 5.:
            sigma_el = 2.0/1e22 * Eel**(-2.5)
        SIGMAarrf.append(sigma_el)
    return(SIGMAarrf)

# 0.03 to 10 keV from https://youlituo.wordpress.com/page/3/    http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1983ApJ...270..119M&amp;data_type=PDF_HIGH&amp;whole_paper=YES&amp;type=PRINTER&amp;filetype=.pdf
def sigmaISMfuncANALITIC(Earrayf):
    coefAf = []
    E_efff = []
    for Eel in np.nditer(Earrayf):
        c0 = 0.
        c1 = 0.
        c2 = 0.
        if Eel >= 0.03 and Eel< 0.1:
            c0 = 17.3
            c1 = 608.1
            c2 = -2150.
        elif Eel >= 0.1 and Eel< 0.284:
            c0 = 34.6
            c1 = 267.9
            c2 = -476.1
        elif Eel >= 0.284 and Eel< 0.400:
            c0 = 78.1
            c1 = 18.8
            c2 = 4.3
        elif Eel >= 0.400 and Eel< 0.532:
            c0 = 71.4
            c1 = 66.8
            c2 = -51.4
        elif Eel >= 0.532 and Eel< 0.707:
            c0 = 95.5
            c1 = 145.8
            c2 = -61.1
        elif Eel >= 0.707 and Eel< 0.867:
            c0 = 308.9
            c1 = -380.6
            c2 = 294.0
        elif Eel >= 0.867 and Eel< 1.303:
            c0 = 120.6
            c1 = 169.3
            c2 = -47.7
        elif Eel >= 1.303 and Eel< 1.840:
            c0 = 141.3
            c1 = 146.8
            c2 = -31.5
        elif Eel >= 1.840 and Eel< 2.471:
            c0 = 202.7
            c1 = 104.7
            c2 = -17.
        elif Eel >= 2.471 and Eel< 3.21:
            c0 = 342.7
            c1 = 18.7
            c2 = 0.
        elif Eel >= 3.21 and Eel< 4.038:
            c0 = 352.2
            c1 = 18.7
            c2 = -0.
        elif Eel >= 4.038 and Eel< 7.111:
            c0 = 433.9
            c1 = -2.4
            c2 = 0.75
        elif Eel >= 7.111 and Eel< 8.331:
            c0 = 629.
            c1 = 30.9
            c2 = 0.0
        elif Eel >= 8.331 and Eel<= 10.0:
            c0 = 701.2
            c1 = 25.2
            c2 = 0.0
        elif Eel > 10.0:
            # shall I use this http://www.aanda.org/articles/aa/pdf/2009/26/aa11794-09.pdf instead?
            c0 = 0.
            c1 = 0.
            c2 = 0.
        else:
            print "WARNING: "
            continue
        
        
        coeff = (c0 + c1*Eel + c2*Eel**2)*Eel**(-3.)
        coefAf.append(coeff)
        E_efff.append(Eel)
    E_efff = np.array(E_efff)
    coefAf = np.array(coefAf)
    SIGMAarrf =coefAf/1e24 # cm^2
    out = [E_efff,SIGMAarrf]
    return(out)




#os.system("rm -r TEMP")
TEMP_dir ="TEMP/" +RUNname
command0 ="mkdir " + TEMP_dir
os.system(command0)

prefix = "../INPUTS/"
#os.system("rm -r INJ_FILES")
#os.system("mkdir INJ_FILES")
######## MARK: ENTER 1.: RATE MODEL #######################
print "############################################# STEP 1.: RATE MODEL #############################################"
percentage_EV = -1
while (percentage_EV>1 or percentage_EV<=0.):
    percentage_EV=input('Enter the fraction of events that are expected to emit the light-curve model (number in the range 0-1): ')
##### defining table - make it consistent with the rate - No sens having delta z in rate < delta z in table

OmegaK = 0.0
OmegaM = 0.0
OmegaL = 0.0
invalid_input = True
while invalid_input:
    print "# Choose among the available models: "
    #MODELS= os.listdir("RATEmodels")
    modDIR = prefix+"RATEmodels/"
    MODELS= os.listdir(modDIR)
    cRM = 1
    print
    print "# PAPER PRESCRIPTIONS: "
    for rowMOD in MODELS:
        printROW = rowMOD.replace("['","")
        printROW = printROW.replace("['","")
        print "# ",cRM,")",printROW
        cRM = cRM+1
    print
    print "# ANALITIC MODELS: "
    print "# ",cRM,") Constant with Uncertainities - I. Mandel 2017"
    cRM = cRM+1
    print "# ",cRM,") Constant - User defined"
    cRM = cRM+1
    #print MODELS
    Nmod = np.shape(MODELS)
    if Nmod[0]==1:
        MODname=MODELS[0]
    else:
        MODname = raw_input('Enter the model name: ')
    #pathMOD="RATEmodels/"+MODname

    VERS = ""
    SOURCE = ""
    if (MODname=="Dominik_et_al_2013" or MODname=="1" or MODname == "1)" or MODname =="1 )"):
        MODname = "Dominik_et_al_2013"
        pathMOD=modDIR+MODname

        invalid_input = False
        invalid_input1 = True
        while invalid_input1:
            MODversion = raw_input('Enter HIGH for the optimistic metallicity evolution, LOW for the pessimistic metallicity evolution: ')
            if (MODversion == "HIGH" or MODversion == "H" or MODversion == "high" or MODversion == "h"):
                pathMOD_V = pathMOD+"/stan.rest.high/"
                VERS = "high"
                invalid_input1 = False
            elif (MODversion == "LOW" or MODversion == "L" or MODversion == "low" or MODversion == "l"):
                pathMOD_V = pathMOD+"/stan.rest.low/"
                invalid_input1 = False
                VERS = "low"
            else:
                print "The option is not available, please enter HIGH or LOW"

        invalid_input1 = True
        while invalid_input1:
            print "Which type of binary systems are going to be the sources of the emission?"
            SOURCE = raw_input('Write NSNS or NSBH or BHBH: ')
            
            if(SOURCE =="NSNS" or SOURCE =="nsns"):
                RATEfile = pathMOD_V+"stan.nsns.rest." + VERS+".dat"
                if percentage_EV!= 0.5:
                    check_0 = True
                    while check_0:
                        ans_changePERC = raw_input('About the 50% of BNS mergers are thought to emit this light curve, do you want to change it? ')
                        if (ans_changePERC=="yes" or ans_changePERC=="y" or ans_changePERC=="YES" or ans_changePERC=="Y"):
                            check_0= False
                            percentage_EV = 0.5
                        elif(ans_changePERC == "no" or ans_changePERC == "NO" or ans_changePERC == "N" or ans_changePERC == "n"):
                            check_0 = False
                        else:
                            print "# Invalid answer"
                invalid_input1 = False
            elif(SOURCE =="NSBH" or SOURCE =="BHNS" or SOURCE =="bhns" or SOURCE =="nsbh"):
                RATEfile = pathMOD_V+"stan.bhns.rest." + VERS+".dat"
                invalid_input1 = False
            elif(SOURCE =="BHBH" or SOURCE =="bhbh"):
                RATEfile = pathMOD_V+"stan.bhbh.rest."+VERS+".dat"
                invalid_input1 = False
            else:
                print "Unknown source systems"

        #CommLINE = "cp "+RATEfile+" TEMP/."
        #        print 'CommLINE',CommLINE
        #        os.system(CommLINE)

        DataRATE=np.genfromtxt(RATEfile)
        zR= DataRATE[:,1]
        zR = np.array(zR)
        zR = zR.squeeze()

        RATE_VT = DataRATE[:,2] #[Gpc^{-3}yr^{-1}]
        for j in range (3,12):
            RATE_VT = RATE_VT + DataRATE[:,j]
        #RATE_VT = DataRATE[:,2] #[Gpc^{-3}yr^{-1}]
        RATE_VT = RATE_VT/10.**9 #[Mpc^{-3}yr^{-1}]
        RATE_VT = RATE_VT/(356*24.*3600.) #[Mpc^{-3}s^{-1}]
        RATE_VT = RATE_VT/(Mpc**3) #[cm^{-3}s^{-1}]
        RATE_VT = np.array(RATE_VT)
        
        indzR = np.argsort(zR)
        zR = np.sort(zR)
        RATE_VT = RATE_VT[indzR]
        
        
        RateFileName = TEMP_dir+"/RATE_Dominik_et_al_2013_" + VERS + "_" + SOURCE + "_cm3_sm1.txt"
        with open(RateFileName, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerows(zip(zR,RATE_VT))
        
        ######### MARK: COSMOLOGY FROM ADOPTED MODEL ############
        print "      ######### COSMOLOGY FROM ADOPTED MODEL #####################################################################"
        print "      #                                                                                                          #"
        print "      # BE AWARE: to be consistent with the chosen model H_0 = 70km/s/Mpc, \Omega_{\Lambda} = 0.7, \Omega_M = 0.3#"

        outTAB_Q= ChooseTABLE_z_V(MODname)
        NameFILE_TAB  = outTAB_Q[0]
        TABzVC = outTAB_Q[1]
        Hc = outTAB_Q[2]
        OmegaM  = outTAB_Q[3]
        OmegaL = outTAB_Q[4]
        OmegaK = outTAB_Q[5]

        print NameFILE_TAB
        CommLINE = "cp "+TABzVC+" "+TEMP_dir+"/."
        print 'CommLINE',CommLINE
        os.system(CommLINE)
        print "RATEfile: ",RATEfile

    elif(MODname=="Ghirlanda_et_al_2016" or MODname=="2" or MODname == "2)" or MODname == "2 )"):
        MODname= "Ghirlanda_et_al_2016"
        pathMOD=modDIR+MODname
        print "# BE AWARE: the rate describes the inferred rate(z) of SGRBs"
        print "# The rate comes from an analytic function (see Ghirlanda et al 2016)"
        #OmegaM = 0.3
        #OmegaL = 0.7
        #OmegaK = 1. - OmegaM - OmegaL
        #Hc = 70.0
        invalid_input = False
        invalid_input1 = True
        
        ###### defining table - make it consistent with the rate - No sens having delta z in rate < delta z in table
        print "      ######### COSMOLOGY FROM ADOPTED MODEL #####################################################################"
        print "      #                                                                                                          #"
        print "      # BE AWARE: to be consistent with the chosen model H_0 = 70km/s/Mpc, \Omega_{\Lambda} = 0.7, \Omega_M = 0.3#"
        outTAB_Q= ChooseTABLE_z_V(MODname)
        NameFILE_TAB  = outTAB_Q[0]
        TABzVC = outTAB_Q[1]
        Hc = outTAB_Q[2]
        OmegaM  = outTAB_Q[3]
        OmegaL = outTAB_Q[4]
        OmegaK = outTAB_Q[5]
        print NameFILE_TAB
        CommLINE = "cp "+TABzVC+" "+TEMP_dir+"/."
        print 'CommLINE',CommLINE
        os.system(CommLINE)
        
        ##### defining zR (array of redshift) - necessary since the prescription is analytical #####
        dataTAB = np.genfromtxt(TABzVC)
        zARR_tab = dataTAB[:,0]
        DeltaZ_tab = zARR_tab[1:] - zARR_tab[:len(zARR_tab)-1]
        
        zR = zARRAY_function(DeltaZ_tab)
        zR = np.sort(zR)
        zR = zR[zR<zlimGhirlanda]
        #print 'zR', zR
        ### check min rate > min tab####
        while invalid_input1:
            print "# 3 different set of parameters are available to describe this model: "
            print "# a) the ones derived assuming CORRELATIONS (Ep-Liso) and (Ep-Eiso)"
            print "# b) the ones derived assuming CORRELATIONS (Ep-Liso) and (Ep-Eiso) and a MINIMUM LUMINOSITY"
            print "# c) the ones derived WITHOUT assumining (Ep-Liso) and (Ep-Eiso) CORRELATIONS"
            print "# for futher detail look at Ghirlanda et al 2016 https://arxiv.org/abs/1607.07875 "
            MODversion = raw_input('Enter A for the first case, B for the second, and C for the third and last one: ')
            if (MODversion == "A" or MODversion == "a" or MODversion == "a)"):
                invalid_input1 = False
                parFile = pathMOD+"/case_a.txt"
                par_caseA = np.genfromtxt(parFile)
                outCASEa = SelectingPAR_Ghirlanda(par_caseA)
                print 'outCASEa',outCASEa
                p1_a = outCASEa[0]
                zp_a = outCASEa[1]
                p2_a = outCASEa[2]
                VAL = outCASEa[3]
                p1_a = np.round(p1_a,2)
                zp_a = np.round(zp_a,2)
                p2_a = np.round(p2_a,2)
                RATE_VT = GhirlandaRate(np.round(p1_a,2), np.round(zp_a,2), np.round(p2_a,2), zR)*costGHIRL_a
                VERS = "case_a"
            elif (MODversion == "B" or MODversion == "b" or MODversion == "b)"):
                print "# Not available at the moment"
                #invalid_input1 = False
                #parFile = pathMOD +"/case_b.txt"
                #par_caseB = np.genfromtxt(parFile)
                #outCASEb = SelectingPAR_Ghirlanda(par_caseB)
                #print 'outCASEb',outCASEb
                #p1_b = outCASEb[0]
                #zp_b = outCASEb[1]
                #p2_b = outCASEb[2]
                #VAL = outCASEb[3]
                #RATE_VT = GhirlandaRate(p1_b, zp_b, p2_b, zR) + costGHIRL_b
                #VERS = "case_b"
            elif (MODversion == "C" or MODversion == "c" or MODversion == "c)"):
                invalid_input1 = False
                parFile = pathMOD +"/case_c.txt"
                par_caseC = np.genfromtxt(parFile)
                outCASEc = SelectingPAR_Ghirlanda(par_caseC)
                print 'outCASEc',outCASEc
                p1_c = outCASEc[0]
                zp_c = outCASEc[1]
                p2_c = outCASEc[2]
                VAL = outCASEc[3]
                p1_c = np.round(p1_c,2)
                zp_c = np.round(zp_c,2)
                p2_c = np.round(p2_c,2)
                RATE_VT = GhirlandaRate(p1_c, zp_c, p2_c, zR)*costGHIRL_c
                
                VERS = "case_c"
            else:
                print "The option is not available, please enter A or B or C"

        # ACCOUNTING FOR BEAMING ANGLE:
        BEAMangle0 = input('Enter the beaming angle [degree] for SGRBs (mean suggested by the authors in the range of 3-6 degress): ')
        BEAMangle = BEAMangle0
        print "BEAMangle",BEAMangle
        print "BEAMangle/180.",BEAMangle/180.
        print "BEAMangle/180.*np.pi",BEAMangle/180.*np.pi
        BEAMangle = (BEAMangle/180.)*np.pi
        print "np.cos(BEAMangle)",np.cos(BEAMangle)
        RATE_VT = RATE_VT/(1. - np.cos(BEAMangle))
        # from Gpc^-3yr^-1 to cm^-3s^-1
        RATE_VT = RATE_VT/10.**9 #[Mpc^{-3}yr^{-1}]
        RATE_VT = RATE_VT/(356*24.*3600.) #[Mpc^{-3}s^{-1}]
        RATE_VT = RATE_VT/(Mpc**3) #[cm^{-3}s^{-1}]
        RATE_VT = np.array(RATE_VT)
        
        # check positive

        RateFileName = TEMP_dir+"/RATE_Ghirlanda_et_al_2016_" + VERS + "_"+ VAL +"_cm3_sm1.txt"
        with open(RateFileName, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerows(zip(zR,RATE_VT))

    elif(MODname == "Constant with Uncertainities" or MODname =="Constant with Uncertainities - I. Mandel" or MODname=="I. Mandel 2017" or MODname == "I. Mandel" or MODname =="3" or MODname == "3 )" or MODname == "3)"):
        MODname = "Mandel2017"
        invalid_input = False
        # TAB - V - Z
        print "      ######### COSMOLOGY FROM ADOPTED MODEL #####################################################################"
        print "      #                                                                                                          #"
        outTAB_Q= ChooseTABLE_z_V(MODname)
        NameFILE_TAB  = outTAB_Q[0]
        TABzVC = outTAB_Q[1]
        Hc = outTAB_Q[2]
        OmegaM  = outTAB_Q[3]
        OmegaL = outTAB_Q[4]
        OmegaK = outTAB_Q[5]
        print NameFILE_TAB
        CommLINE = "cp "+TABzVC+" "+TEMP_dir+"/."
        print 'CommLINE',CommLINE
        os.system(CommLINE)
        
        dataTAB = np.genfromtxt(TABzVC)
        zARR_tab = dataTAB[:,0]
        DeltaZ_tab = zARR_tab[1:] - zARR_tab[:len(zARR_tab)-1]
        
        print "# The model predict a constant cosmological rate of 100 per Gpc^3 per year with an uncertainity of ~10% in each direction "
        

        zR = zARRAY_function(DeltaZ_tab)
        
        zR = np.sort(zR)
        zR = zR[zR<zlimMandel]
        
        RATE_VT = 1000*np.ones(len(zR)) # per Gpc^3 per yr
        RATE_VT_low = 100*np.ones(len(zR)) # per Gpc^3 per yr
        RATE_VT_high = 10000*np.ones(len(zR)) # per Gpc^3 per yr
        # from Gpc^-3yr^-1 to cm^-3s^-1
        RATE_VT = RATE_VT/10.**9 #[Mpc^{-3}yr^{-1}]
        RATE_VT = RATE_VT/(356*24.*3600.) #[Mpc^{-3}s^{-1}]
        RATE_VT = RATE_VT/(Mpc**3) #[cm^{-3}s^{-1}]
        RATE_VT = np.array(RATE_VT)
        
        RATE_VT_low = RATE_VT_low/10.**9 #[Mpc^{-3}yr^{-1}]
        RATE_VT_low = RATE_VT_low/(356*24.*3600.) #[Mpc^{-3}s^{-1}]
        RATE_VT_low = RATE_VT_low/(Mpc**3) #[cm^{-3}s^{-1}]
        RATE_VT_low = np.array(RATE_VT_low)

        RATE_VT_high = RATE_VT_high/10.**9 #[Mpc^{-3}yr^{-1}]
        RATE_VT_high = RATE_VT_high/(356*24.*3600.) #[Mpc^{-3}s^{-1}]
        RATE_VT_high = RATE_VT_high/(Mpc**3) #[cm^{-3}s^{-1}]
        RATE_VT_high = np.array(RATE_VT_high)

        RateFileName = TEMP_dir+"/RATE_Constant_Mandel2017_cm3_sm1.txt"
        with open(RateFileName, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerows(zip(zR,RATE_VT))

        RateFileName_low = TEMP_dir+"/RATE_Constant_lower_Mandel2017_cm3_sm1.txt"
        with open(RateFileName_low, 'w') as flow:
            writer = csv.writer(flow, delimiter='\t')
            writer.writerows(zip(zR,RATE_VT_low))

        RateFileName_high = TEMP_dir+"/RATE_Constant_upper_Mandel2017_cm3_sm1.txt"
        with open(RateFileName_high, 'w') as fhigh:
            writer = csv.writer(fhigh, delimiter='\t')
            writer.writerows(zip(zR,RATE_VT_high))
                
    elif(MODname == "Constant - User defined" or MODname=="Constant" or MODname=="User defined" or MODname == MODname =="4" or MODname == "4 )" or MODname == "4)"):
        print "# The model predict a constant cosmological rate with an uncertainities"
        MODname = "User"
        invalid_input = False
        # TAB - V - Z
        outTAB_Q= ChooseTABLE_z_V(MODname)
        NameFILE_TAB  = outTAB_Q[0]
        TABzVC = outTAB_Q[1]
        Hc = outTAB_Q[2]
        OmegaM  = outTAB_Q[3]
        OmegaL = outTAB_Q[4]
        OmegaK = outTAB_Q[5]
        print NameFILE_TAB
        CommLINE = "cp "+TABzVC+" "+TEMP_dir+"/."
        print 'CommLINE',CommLINE
        os.system(CommLINE)
        
        dataTAB = np.genfromtxt(TABzVC)
        zARR_tab = dataTAB[:,0]
        DeltaZ_tab = zARR_tab[1:] - zARR_tab[:len(zARR_tab)-1]
        
        zR = zARRAY_function(DeltaZ_tab)

        checkNEGrate = True
        while checkNEGrate:
            RATE_VT = input('Insert the fiducial constant rate per Gpc^3 per year: ')
            if RATE_VT>=0:
                checkNEGrate =False
            else:
                print "ERROR: negative rate"
    
        checkNEGrate = True
        
        while checkNEGrate:
            RATE_VT_low = input('Insert the minimum constant rate per Gpc^3 per year: ')
            if RATE_VT_low>=0:
                checkNEGrate = False
            else:
                print "ERROR: negative minimum rate"
        
        checkNEGrate = True
        while checkNEGrate:
            RATE_VT_high = input('Insert the maximum constant rate per Gpc^3 per year: ')
            if RATE_VT_high>=0:
                checkNEGrate =False
            else:
                print "ERROR: negative maximum rate"
    
        RATE_VT_low = RATE_VT_low*np.ones(len(zR))/10.**9 #[Mpc^{-3}yr^{-1}]
        RATE_VT_low = RATE_VT_low/(356*24.*3600.) #[Mpc^{-3}s^{-1}]
        RATE_VT_low = RATE_VT_low/(Mpc**3) #[cm^{-3}s^{-1}]
        RATE_VT_low = np.array(RATE_VT_low)
    
        RATE_VT = RATE_VT*np.ones(len(zR))/10.**9 #[Mpc^{-3}yr^{-1}]
        RATE_VT = RATE_VT/(356*24.*3600.) #[Mpc^{-3}s^{-1}]
        RATE_VT = RATE_VT/(Mpc**3) #[cm^{-3}s^{-1}]
        RATE_VT = np.array(RATE_VT)
    
        RATE_VT_high = RATE_VT_high*np.ones(len(zR))/10.**9 #[Mpc^{-3}yr^{-1}]
        RATE_VT_high = RATE_VT_high/(356*24.*3600.) #[Mpc^{-3}s^{-1}]
        RATE_VT_high = RATE_VT_high/(Mpc**3) #[cm^{-3}s^{-1}]
        RATE_VT_high = np.array(RATE_VT_high)
    
        RateFileName = TEMP_dir+"/RATE_Constant_User_cm3_sm1.txt"
        with open(RateFileName, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerows(zip(zR,RATE_VT))
        
        RateFileName_low = TEMP_dir+"/RATE_Constant_lower_User_cm3_sm1.txt"
        with open(RateFileName_low, 'w') as flow:
            writer = csv.writer(flow, delimiter='\t')
            writer.writerows(zip(zR,RATE_VT_low))
        
        RateFileName_high = TEMP_dir+"/RATE_Constant_upper_User_cm3_sm1.txt"
        with open(RateFileName_high, 'w') as fhigh:
            writer = csv.writer(fhigh, delimiter='\t')
            writer.writerows(zip(zR,RATE_VT_high))
    
    
    
    else:
        print "Invalide model"

######## MARK: STEP 2.: LIGHT CURVES #######################
print
print "############################################# STEP 2.: LIGHT CURVES #############################################"
invalid_input = True
ENABLE_LC_z = 1
while invalid_input:
    print "# Choose the emission model (light curve) you would like to test from the following list: "
    modLCDir = prefix+"LightCurves/"
    MOD_LCs= os.listdir(modLCDir)
    print MOD_LCs
    Nmod_LC = np.shape(MOD_LCs)
    if Nmod_LC[0]==1:
        MOD_LC=MOD_LCs[0]
    else:
        MOD_LC = raw_input('# Enter the model name: ')
    pathMOD_LC=modLCDir+MOD_LC
    #print new
    #VERS_LC = ""
    SOURCE_LC = ""
    if (MOD_LC=="Siegel_Ciolfi_2016"):
        invalid_input = False
        invalid_source = True
        while invalid_source:
            print "# Which type of binary systems are going to be the sources of the emission?"
            SOURCE_LC = raw_input('Write SNS for stable neutron star or SMNS for supramassive neutron star: ')
            if (SOURCE_LC=="SNS" or SOURCE_LC=="sns" or SOURCE_LC=="stableNS" or SOURCE_LC=="stablens" or SOURCE_LC=="stable neutron star"):
                invalid_source = False
                SOURCE_LC="stableNS"
                PATH_SOURCE = SOURCE_LC
            elif (SOURCE_LC=="SMNS" or SOURCE_LC=="smns" or SOURCE_LC=="supramassive neutron star" or SOURCE_LC=="supramassive NS"):
                invalid_source = False
                SOURCE_LC="SMNS"
                PATH_SOURCE = "supramassiveNS"
            else:
                print "Unknown source"
    #elif():
    #    invalid_input = False
    else:
        print "Invalid emission model"
#print "Do the selected model allow for proper redshift effects (Are the light curve available for energy also outside the instrument band and involved for high z)?"
Inv_ans = True
LimE = -1
while Inv_ans:
    Ans_moreLC =raw_input('Do the selected model allow for proper redshift effects (Are the light curve available for energy also outside the instrument band and involved for high z)? ')

    if (Ans_moreLC=="yes" or Ans_moreLC=="y" or Ans_moreLC=="YES" or Ans_moreLC=="Y"):
        ENABLE_LC_z = 1
        Inv_ans = False
    elif (Ans_moreLC=="No" or Ans_moreLC=="NO" or Ans_moreLC=="n" or Ans_moreLC=="no" or Ans_moreLC=="N"):
        ENABLE_LC_z = 0
        Inv_ans = False
        LimE = input('To extimate the maximum redshift please enter maximum energy (in keV) at which the maximum of the emission is expected ')
    else:
        print "Invalid answer"

#File_LCfuncE =pathMOD_LC+"/"+PATH_SOURCE+"/LC_funcE_"+SOURCE_LC+".txt"
#File_EnXLC =pathMOD_LC+"/"+PATH_SOURCE+"/Energies_"+SOURCE_LC+"_LC.txt"

invalid_inputLC = True
while invalid_inputLC:
    print "# Choose among the available files: "
    DIR = pathMOD_LC+"/"+PATH_SOURCE+"/"
    listName = File_LCfuncE =pathMOD_LC+"/"+PATH_SOURCE+"/list_LCfiles.txt"
    CommLINE = "ls "+DIR+"*.txt > "+listName
    os.system(CommLINE)
    MOD_LCA = []
    with open(listName,'rb') as flist:
        #next(f) #skip heading row in text file (I cannot get csv.Dictreader instead to work as an alternative to this step)
        data_list = csv.reader(flist) #read text file with csv
        count = 0
        for row in data_list:
            print "# n: ",count,". file:", row[0]
            fileName_LCopt =row[0]
            count = count+1
            MOD_LCA.append(row)

        MOD_LCA = np.array(MOD_LCA)

    test_LCfile0 =0
    test_LCfile1 =0
    enter_MOD_LC =input('Write the number correspondent to the desired Light curves: ')
    enter_MOD_LC1 =input('Write the number correspondent to the correspondent energies: ')
    for aa in range(len(MOD_LCA)):

        if enter_MOD_LC == aa:
            fileName_LCopt=MOD_LCA[aa][0]
            #print "MOD_SURVname",
            test_LCfile0 =1
        if enter_MOD_LC1 == aa:
            fileName_ELCopt=MOD_LCA[aa][0]
            #print "MOD_SURVname",
            test_LCfile1 =1
        if test_LCfile1*test_LCfile0 ==1:
            invalid_inputLC = False
            break

    if(test_LCfile0 ==0 or test_LCfile1 ==0):
        print "# Unknown file, choose one from the available list or add new entries"

COMM4 = "rm "+ listName
os.system(COMM4)

File_LCfuncE =fileName_LCopt
print 'File_LCfuncE ',File_LCfuncE

File_EnXLC =fileName_ELCopt
print 'File_EnXLC ',File_EnXLC

LC_j_LC= np.genfromtxt(File_LCfuncE)
SHAPE_LC = np.shape(LC_j_LC)



ELC_j_LC= np.genfromtxt(File_EnXLC)
SHAPE_ELC = np.shape(ELC_j_LC)


if SHAPE_LC[1]-1!=SHAPE_ELC[0]:
    print 'ERROR: numeber of light curves available different from the numbers of energies in files'
    exit()

maxL =np.max(LC_j_LC[:,3])

print " the dimension of the light curve is: ", SHAPE_LC[0]
checkANS = True
ATT = ""
while checkANS:
    NinLC =input('Would you like to reduce it? If YES, print the number at which you would like to reduce, 0 or less otherwise ')
    if NinLC>0:
        Nlim = NinLC
        if NinLC<SHAPE_LC[0]:
            checkANS = False
            t_LCO =LC_j_LC[:,0]
            fac = np.int(len(t_LCO)/Nlim)
            a_range = np.arange(0,Nlim)*fac
            
            t_LCO0=map(t_LCO.__getitem__, a_range)
            #print 'check desire dimension - time: ',len(t_LCO0)
            
            L_LC_Ared = [t_LCO0]
            for i in range(1,SHAPE_LC[1]):
                L_LCO= LC_j_LC[:,i]
                fac = np.int(len(t_LCO)/Nlim)
                a_range = np.arange(0,Nlim)*fac
                L_LCO0=map(L_LCO.__getitem__, a_range)
                L_LC_Ared.append(L_LCO0)

            L_LC_Ared = np.matrix(L_LC_Ared)
            L_LC_Ared = np.transpose(L_LC_Ared)
            ATT = "red_"
            LC_E_Filename=prefix+"LC_funcE_red_"+ATT+SOURCE_LC+".txt"
            print "LC_E_Filename",LC_E_Filename
            np.savetxt(LC_E_Filename, L_LC_Ared, fmt='%.18e', delimiter='\t', newline='\n', header='', footer='', comments='# ')
        else:
            print "WARNING: the only available option consists in decreasing the number of points in the light curve"
    else:
        print "You chose to maintain the original points of the light curve"
        checkANS = False


############### MARK: STEP 2.1 : ADD POSSIBLE ABSORPTION OF THE HOST GALAXY###################################
print
Ans_HGabs =raw_input('Would you like to add the ABSORPTION of a possible HOST GALAXY? ')
if (Ans_HGabs=="yes" or Ans_HGabs=="y" or Ans_HGabs=="YES" or  Ans_HGabs=="Y"):
    # CHECK IF THE FILE ALREADY EXISTS
    En_j_LC= np.genfromtxt(File_EnXLC)
    # CALCULATE nH and sigma as a function of E
    ABS = []
    for i in range(len(En_j_LC)):
        abs_i = En_j_LC[i]**(-8./3.)/En_j_LC[i]**(-8./3.)
        ABS.append(abs_i)
    ABS = np.array(ABS)
    # READING AND MODIFY THE LUMINOSITY
    LC_j_LC= np.genfromtxt(File_LCfuncE)
    TIME_arr = LC_j_LC[:,0]
    #print "LC_j_LC",np.shape(LC_j_LC)
    DIM = np.shape(LC_j_LC)
    D2 = DIM[1]
    L_LC_A = []
    L_LC_A.append(np.array(TIME_arr))
    print "second dimension",D2
    for j in range(1,D2):
        LC_j_new = ABS[j-1]*np.array(LC_j_LC[:,j])
        L_LC_A.append(np.array(LC_j_new))
    L_LC_A = np.matrix(L_LC_A)
    L_LC_A = np.transpose(L_LC_A)

    ATT = ATT +"ABS_"

    # WRITE NEW FILE WITH THIS INFO AND CHANGE PATH NAME

    LC_E_Filename=prefix+"LC_funcE_"+ATT+SOURCE_LC+".txt"
    np.savetxt(LC_E_Filename, L_LC_A, fmt='%.18e', delimiter='\t', newline='\n', header='', footer='', comments='# ')
    File_LCfuncE = LC_E_Filename
    #SOURCE_LC = "ABS_"+SOURCE_LC

else:
    print "With an answer different from 'yes' you chose NO ABSROPTION at host"
#ASK THIS IN THE CONSTRUCTION ON THE LIGHT CURVE FILES######################################################
#ANS_MaxDimLC = raw_input('Write YES if you want to have a maximum dimension of the light curve vectors? ')
#if (ANS_MaxDimLC == "yes" or "YES" or "y" or "Y"):
#    MaxD_LC = input('Write the maximum dimension of the light curve vector you would like to set: ')
#    Data_LC=np.genfromtxt(RATEfile)
#    Times_LC = Data_LC[]
#else:
############################################################################################################

print File_LCfuncE
CommLINE = "cp "+File_LCfuncE+" "+TEMP_dir+"/LC_"+MOD_LC+"_funcE_"+SOURCE_LC+".txt"
print 'CommLINE',CommLINE
os.system(CommLINE)

#print File_EnXLC
CommLINE = "cp "+File_EnXLC+" "+TEMP_dir+"/Energies_"+MOD_LC+"_"+SOURCE_LC+"_LC.txt"
print 'CommLINE',CommLINE
os.system(CommLINE)
print

######################################## MARK: STEP 3.: SURVEY PROPERTIES ##############################################
print "######################################## STEP 3.: SURVEY PROPERTIES ##############################################"
invalid_input = True
while invalid_input:
    print "# Choose among the available survey of data: "
    #MOD_SURV= os.listdir("SURVEYprop/*.txt")
    #print MOD_SURV
    #listName = prefix+"list_SURVEYprop.txt"
    SurvDIR = prefix+"SURVEYprop/"
    listName = prefix+"list_SURVEYprop.txt"
    CommLINE = "ls "+SurvDIR+"*.txt > "+listName
    os.system(CommLINE)
    MODLIST = []
    with open(listName,'rb') as flist:
        #next(f) #skip heading row in text file (I cannot get csv.Dictreader instead to work as an alternative to this step)
        data_list = csv.reader(flist) #read text file with csv
        count = 0
        for row in data_list:
            print "# n: ",count,". model:", row[0]
            count = count+1
            MODLIST.append(row)

    MODLIST = np.array(MODLIST)

    enter_MODSurvey =input('Write the number correspondent to the desired survey: ')
    #MOD_SURVname = raw_input('Enter the model name: ')
    cIN = 0
    while cIN < count:
        if enter_MODSurvey == cIN:
            MOD_SURVname=MODLIST[cIN][0]
            invalid_input = False
            break
        cIN = cIN+1
    if invalid_input==True:
        print "# Unknown survey, choose one from the available list or add new entries"
    """
    if enter_MODSurvey == 0:
        MOD_SURVname=MODLIST[0][0]
        invalid_input = False
    elif enter_MODSurvey == 1:
        MOD_SURVname=MODLIST[1][0]
        invalid_input = False
        #if(MOD_SURVname=="SURVEYprop/INFO_XMM_SLEW.txt"):
        #invalid_input = False
        #elif(MOD_SURVname=="SURVEYprop/INFO_XMM_POINTED_OBS.txt"):
        #       invalid_input = False
    else:
        print "# Unknown survey, choose one from the available list or add new entries"
    """

SURVname = MOD_SURVname
"""
    with open(SURVname, 'r') as SURVinfoF:
        #SURVinfoData = SURVinfoF.read()
        table = []
        for line in SURVinfoF:
            numbers_str = line.split(delimiter='\t')
            print "line",line," numbers_str",numbers_str
            L = []
            for el in numbers_str:
                L.append(el)
            L = np.array(L)
            table.append(L)
        print 'np.shape(table)',np.shape(table)
        print table[0]
        a = table[0]
        print a[1]
        b = 3+ a[1]
        print b
"""


list=[]
with open(SURVname,'rb') as f:
    #next(f) #skip heading row in text file (I cannot get csv.Dictreader instead to work as an alternative to this step)
    data = csv.reader(f,delimiter='\t') #read text file with csv
    for row in data:
        list.append(row)
        print "row", row

nOBS =  0
AvObsT = 0.0
SD_T = 0.0
FoV = 0.0
FracSkyArea =0.0
Emin = -1.
Emax = -1.
SensFile=""

for el in range (0,len(list)):
    print 'len(list): ',len(list),'el: ',el
    line = list[el]
    print 'line: ',line
    dim_line = np.shape(line)
    if (dim_line[0]!=0):
        label = line[0]
    else:
        continue
    if (label=="minE"):
            # transform units in keV
            if (line[len(line)-1] == "keV" or line[len(line)-1] == "KeV"):
                Emin = float(line[1])
            elif (line[len(line)-1] == "eV" ):
                Emin = float(line[1])/1000.
            elif (line[len(line)-1] =="erg"):
                Emin = float(line[1])*624150647.9963236
            elif (line[len(line)-1] =="J" or line[len(line)-1] == "j"):
                Emin = float(line[1])*6.2415096471204*10**15
            elif (line[len(line)-1] =="kJ" or line[len(line)-1] == "kj"):
                Emin = float(line[1])*6.2415096471204*10**12
            else:
                print "Invalid Units for Emin, available units are: Kev, eV, erg, J and kJ"
                print "Correct INFO-file on the interested survey and start again"
                exit()
            print "Emin",Emin
            Emin = np.array(Emin)

    if (label=="maxE"):
            if (line[len(line)-1] == "keV" or line[len(line)-1] == "KeV"):
                Emax = float(line[1])
            elif (line[len(line)-1] == "eV" ):
                Emax = float(line[1])/1000.
            elif (line[len(line)-1] =="erg"):
                Emax = float(line[1])*624150647.9963236
            elif (line[len(line)-1] =="J" or line[len(line)-1] == "j"):
                Emax = float(line[1])*6.2415096471204*10**15
            elif (line[len(line)-1] =="kJ" or line[len(line)-1] == "kj"):
                Emax = float(line[1])*6.2415096471204*10**12
            else:
                print "Invalid Units for Emax, available units are: Kev, eV, erg, J and kJ"
                print "Correct INFO-file on the interested survey and start again"
                exit()
            print "Emax",Emax
            Emax = np.array(Emax)

    if (label =="nOBS" and len(line)>=2):
                nOBS = line[1]
                print "nOBS",nOBS
                nOBS =int(nOBS)
                print 'type(nOBS)',type(nOBS)

    if (label =="mean observing time"):
            if (line[len(line)-1] == "h" or line[len(line)-1] == "hour" or line[len(line)-1] == "HOUR" or line[len(line)-1] == "hours" or line[len(line)-1] == "HOURS"):
                AvObsT = float(line[1])*3600.
            elif (line[len(line)-1] == "s" or line[len(line)-1] == "S" or line[len(line)-1] == "seconds" or line[len(line)-1] == "second"):
                AvObsT = float(line[1])
            elif (line[len(line)-1] == "year" or line[len(line)-1] == "years" or line[len(line)-1] == "yr" or line[len(line)-1] == "yrs" or line[len(line)-1] == "YR" or line[len(line)-1] == "YRS" or line[len(line)-1] == "YEARS" or line[len(line)-1] == "YEAR"):
                AvObsT = float(line[1])*3600*24*365.
            elif (line[len(line)-1] == "min" or line[len(line)-1] == "minute" or line[len(line)-1] == "MIN" or line[len(line)-1] == "minutes" or line[len(line)-1] == "MINUTE" or line[len(line)-1] == "MINUTES"):
                AvObsT = float(line[1])*60.
            else:
                print "Invalid Units for \"mean observing time\", available units are: s, min, h and yr"
                print "Correct INFO-file on the interested survey and start again"
                exit()
            print "AvObsT",AvObsT
            AvObsT = np.array(AvObsT)

    if (label =="sigma observing time"):
            if (line[len(line)-1] == "h" or line[len(line)-1] == "hour" or line[len(line)-1] == "HOUR" or line[len(line)-1] == "hours" or line[len(line)-1] == "HOURS"):
                SD_T = float(line[1])*3600.
            elif (line[len(line)-1] == "s" or line[len(line)-1] == "S" or line[len(line)-1] == "seconds" or line[len(line)-1] == "second"):
                SD_T = float(line[1])
            elif (line[len(line)-1] == "year" or line[len(line)-1] == "years" or line[len(line)-1] == "yr" or line[len(line)-1] == "yrs" or line[len(line)-1] == "YR" or line[len(line)-1] == "YRS" or line[len(line)-1] == "YEARS" or line[len(line)-1] == "YEAR"):
                SD_T = float(line[1])*3600*24*365.
            elif (line[len(line)-1] == "min" or line[len(line)-1] == "minute" or line[len(line)-1] == "MIN" or line[len(line)-1] == "minutes" or line[len(line)-1] == "MINUTE" or line[len(line)-1] == "MINUTES"):
                SD_T = float(line[1])*60.
            else:
                print "Invalid Units for \"sigma observing time\", available units are: s, min, h and yr"
                print "Correct INFO-file on the interested survey and start again"
                exit()
            print "SD_T",SD_T
            SD_T = np.array(SD_T)
    if (label == "FoV"):
        if len(line)>1:
            if (line[len(line)-1] =="deg^2"):
                FoV =  float(line[1])
            elif (line[len(line)-1] =="rad^2" or line[len(line)-1] =="steradians"):
                FoV =   float(line[1])*(180./np.pi)**2
            elif (line[len(line)-1] =="seconds^2" or line[len(line)-1] =="second^2" or line[len(line)-1] =="'^2" or line[len(line)-1] =="sqs^2"):
                FoV =   float(line[1])/(3600**2)
            elif (line[len(line)-1] =="min^2" or line[len(line)-1] =="minutes^2" or line[len(line)-1] =="minute^2" or line[len(line)-1] =="''^2" or line[len(line)-1] =="arcmin^2")or line[len(line)-1] =="arcminutes^2":
                FoV =   float(line[1])/(3600.)
            else:
                print "Invalid Units for \"FoV\", available units are: deg^2, steradians, seconds^2 and min^2"
                print "Correct INFO-file on the interested survey and start again"
                exit()
            print "FoV",FoV
            FoV = float(FoV)
            print 'type(FoV)',type(FoV)
            print 'FoV',FoV
            print "1+FoV",1+FoV
    if (label == "FracSkyArea" and len(line)>1):
        if (line[1]!=""):
            print line[1]
            FracSkyArea = float(line[1])
            print "FracSkyArea",FracSkyArea
            FracSkyArea = np.array(FracSkyArea)
            
    if (label == "SensCurve" and len(line)>1):
        SensFile = line[1]
        print  "SensCurve",SensFile
        
#CHECK IF MISSING INFORMATION
AREAsky = 4.*np.pi*(180./np.pi)**2 #deg^2
if(AvObsT == 0.0):
     print "Missing average time of observations"
     print "Correct INFO-file on the interested survey and start again"
     exit()
                
if(Emin == -1.):
    print "Missing minimum energy detectable Emin"
    print "Correct INFO-file on the interested survey and start again"
    exit()

if(Emax == -1.):
    print "Missing maximum energy detectable Emin"
    print "Correct INFO-file on the interested survey and start again"
    exit()
                  
if (SensFile==""):
    print "Missing path for sensitivity curve file"
    print "Correct INFO-file on the interested survey and start again"
    exit()
else:
    SensFile = prefix+SensFile

if (nOBS*FoV==0.0 and FracSkyArea ==0.0):
    print "Missing information on covered area, fill \"FracSkyArea\" and/or (\"nOBS\" and \"FoV\")"
    print "Correct INFO-file on the interested survey and start again"
    exit()
                  
elif(nOBS*FoV/AREAsky>0.0 and FracSkyArea >0.0):
    if(FracSkyArea>1.):
        FracSkyArea= FracSkyArea*0.01
    if (np.abs(FracSkyArea-nOBS*FoV/AREAsky)<0.001):
        CoveredSky = FracSkyArea
        print "WARNING: observations of the same sky area should considered"
    else:
                  #test= raw_input('Is FracSkyArea representing the area covered by: observations of different sky regions + observations of the same sky ragions?')
                  #if (test =="yes" or test=="YES" or test =="Y" or test=="y"):
                  #CoveredSky = FracSkyArea
                  #else:
        print "Inconsistent data"
        print "Correct INFO-file on the interested survey and start again"
        exit()
                  
elif(nOBS*FoV/AREAsky<0.0 or FracSkyArea <0.0):
    print "Error: nOBS, FracSkyArea and FoV cannot be negative"
    print "Correct INFO-file on the interested survey and start again"
    exit()
                  
else:
    if (FracSkyArea>0.0):
        CoveredSky = FracSkyArea
        if(CoveredSky>1.):
                  CoveredSky= CoveredSky*0.01
        else:
           print " NOTE:  the value of \"FracSkyArea\" should be >0.0 and <1.0"
    else:
        CoveredSky = nOBS*FoV/AREAsky
        print
        print 'CoveredSky',CoveredSky
        print 'nOBS',nOBS
        print 'FoV',FoV
        print 'AREAsky',AREAsky
        print
print "The equivalent area of sky covered by the survey is: ",CoveredSky*100.," %"

print

print "###########################################################################################################################"

print

###################################### MARK: STEP 4.: USER PARAMETERS - ENERGY BANDS OF THE INSTRUMENT #####################################################
print
print "###################################### STEP 4.: USER PARAMETERS - ENERGY BANDS OF THE INSTRUMENT #####################################################"
print "# Select the number of bands in which deviding the instrument energy range.                                                                          #"
En_Sens= np.genfromtxt(SensFile)
dim_sens = np.shape(En_Sens)
print "dim_sens",dim_sens,"len(dim_sens)",len(dim_sens)

EnAv_A = np.zeros(len(dim_sens))
EnDW_A =  np.zeros(len(dim_sens))
EnUP_A = np.zeros(len(dim_sens))
Sens =  np.zeros(len(dim_sens))
Sens_arr = Sens
dim0= 0
if len(dim_sens) ==1:
    EnAv_A = En_Sens[0]
    EnDW_A = En_Sens[2]
    EnUP_A = En_Sens[3]
    Sens = En_Sens[1]
    dim= np.shape(EnAv_A)
    dim0=1
elif len(dim_sens) >1:
    EnAv_A = En_Sens[:,0]
    EnDW_A = En_Sens[:,2]
    EnUP_A = En_Sens[:,3]
    Sens = En_Sens[:,1]
    dim= np.shape(EnAv_A)
    dim0= dim[0]
else:
    print "ERROR: empty sensitivity file"
    exit()

print "EnAv_A",EnAv_A
print "EnDW_A",EnDW_A
print "EnUP_A",EnUP_A
print "Sens",Sens
"""
EnAv_A = En_Sens[:,0]
EnDW_A = En_Sens[:,2]
EnUP_A = En_Sens[:,3]
Sens = En_Sens[:,1]
"""
Sens_arr = Sens
Ebins_sup = EnUP_A
Ebins_inf = EnDW_A
Ebins_av = EnAv_A

EQ_SP = ""
print "# From the sensitivity files it is divided in ",dim0," band(s), with"
print "# Lower energies: ", EnDW_A,
print "# Average energies: ", EnAv_A
print "# Upper energies: ", EnUP_A

invANS = True
while invANS:
    ANS=raw_input('# Would you like to leave it as it is? ')
    if (ANS == "yes" or ANS == "YES" or ANS == "Y" or ANS == "y"):
        invANS = False
        FILE_SENStoUSE = En_Sens
        EQ_SP="YES"
    elif(ANS == "no" or ANS == "NO" or ANS == "N" or ANS == "n"):
        invANS = False
        print " # The energy range of the chosen instrument is: ",Emin," - ",Emax," keV"
        Nband = input('# If equally spaced, enter the number of bands in which you would like to divide the instrumental energy range (0 if the energy bin are not equally spaced) ')
        
        Ebins_sup = []
        Ebins_inf = []
        Ebins_av = []
        Sens_arr = []
        
        if (Nband >0 ):
            EQ_SP="YES"
            DE = Emax - Emin
            DEbin = DE/Nband
            for i in range(Nband):
                Ebins_sup.append(Emin +(i+1)*DEbin)
                Ebins_inf.append(Emin +(i)*DEbin)
                Ebins_av.append(Emin +(i+0.5)*DEbin)
            Ebins_av= np.array(Ebins_av)
            Ebins_inf= np.array(Ebins_inf)
            Ebins_sup= np.array(Ebins_sup)

        elif(Nband==0):
            
            EQ_SP="NO"
            INen = True
            while INen:
                Nband_fin =input('# Print the number of bands in which you would like to divide the instrumental energy range ')
                if (Nband_fin>0):
                    INen= False
                else:
                    print " Invalid number od bands"
            EnR = True
            EnR_ALLcov = 0
            while EnR:
                ANS_ENrange=raw_input("Do you want to cover all the instrument energy range?")
                if (ANS_ENrange== "yes" or ANS_ENrange=="YES" or ANS_ENrange== "Y" or ANS_ENrange=="y"):
                    EnR = False
                    EnR_ALLcov = 1
                elif(ANS_ENrange == "no" or ANS_ENrange == "NO" or ANS_ENrange == "N" or ANS_ENrange == "n"):
                    EnR = False
                else:
                    print "Invalid answer"
            SUM_rangeEbin = 0.0
            Emax_i = -1.0
            for i in range(Nband_fin):
                Emax_i_old = Emax_i
                sanityCHECK = True
                while sanityCHECK:
                    print "i == ",i
                    #Emin_i (i-band)
                    if i ==0:
                        if EnR_ALLcov==0:
                            print "IMPORTANT NOTE: insert increasing values of energy"
                            Emin_i = input('i-th lower limit of energy bin ')
                        else:
                            Emin_i = Emin
                    else:
                        if EnR_ALLcov==0:
                            Emin_i = input('i-th lower limit of energy bin ')
                            if Emin_i<Emax_i_old:
                                print "Error: definition of overlapping bins"
                                print "i-th lower limit of energy bin will be set to ",Emax_i
                                Emin_i=Emax_i_old
                        else:
                            Emin_i = Emax_i
                        #Emax_i (i-band)
                    if (i<Nband_fin - 1 or EnR_ALLcov==0):
                        Emax_i = input('i-th upper limit of energy bin ')
                    else:
                        Emax_i = Emax

                    if (Emin_i<Emax and Emin_i >= Emin and Emax_i> Emin and Emax_i<=Emax and Emax_i>Emin_i):
                        Eav_i = 0.5*(Emax_i+Emin_i)
                        Ebins_inf.append(Emin_i)
                        Ebins_sup.append(Emax_i)
                        Ebins_av.append(Eav_i)
                        SUM_rangeEbin = SUM_rangeEbin + Emax_i - Emin_i
                        sanityCHECK = False
                    else:
                        print " invalid entry "
                        print "Emin ",Emin, "; Emax ",Emax, "; Emin_i ",Emin_i, "; Emax_i ",Emax_i
                            
            print "SUM_rangeEbin",SUM_rangeEbin
            print "(Emax- Emin)",(Emax- Emin)
            if (SUM_rangeEbin< (Emax- Emin)):
                print "WARNING: not all the instrument band have been covered" # IS it going to be a problem later on? or only if they are disjoint or none?
            elif(SUM_rangeEbin> (Emax- Emin)):
                print "ERROR: with non overlapping bins it is covered more tha the available"
                exit()
            Ebins_av= np.array(Ebins_av)
            Ebins_inf= np.array(Ebins_inf)
            Ebins_sup= np.array(Ebins_sup)
            print "you chose as inferior limits of the energy bins: ",Ebins_inf
            print "you chise as superior limits of the energy bins: ", Ebins_sup
            print "   consequntly the average value of the energy bins are: ", Ebins_av
        
        else:
            print "# Error: negative number of bins - the original one is left"
            EQ_SP="YES"
            DE = Emax - Emin
            DEbin = DE/Nband
            for i in range(Nband):
                Ebins_sup.append(Emin +(i+1)*DEbin)
                Ebins_inf.append(Emin +(i)*DEbin)
                Ebins_av.append(Emin +(i+0.5)*DEbin)
            Ebins_av= np.array(Ebins_av)
            Ebins_inf= np.array(Ebins_inf)
            Ebins_sup= np.array(Ebins_sup)


        #Linear interpolation#
        #NOTE!!! Very ROUGHT
        print "NOTE the SENSITIVITY CURVE EXTRACTED BY LINEALY INTERPOLATING THE AVAILABLE ONE: the value at the average point of the bin is calculated by linearly interpolating the two recorded (original) points which enclose it"
        for j in range(len(Ebins_av)):
            E_i = Ebins_av[j]
            Seq = Sens[EnAv_A==E_i]
            dim = np.shape(Seq)
            if dim[0] !=0:
                #print "if"
                Sens_arr.append(Seq[0])
            else:
                print
                Einf = EnAv_A[EnAv_A<=E_i]
                Einf = np.sort(Einf)
                Einf_i = 0
                if len(Einf)>0:
                    Einf_i = Einf[len(Einf)-1] #error
                Esup = EnAv_A[EnAv_A>E_i]
                Esup = np.sort(Esup)
                Esup_i = 0
                if len(Esup)>0:
                    Esup_i = Esup[0]
                #print "Einf_i",Einf_i,"Esup_i",Esup_i
                if (Esup_i >0 and Einf_i >0):
                    Sens_s = Sens[EnAv_A==Esup_i]
                    Sens_i = Sens[EnAv_A==Einf_i]
                    #print "Sens_s",Sens_s,"Sens_i",Sens_i
                    Sens_av = Sens_i+(Sens_s- Sens_i)/(Esup_i - Einf_i)*(E_i - Einf_i)
                    print "Sens_av",Sens_av
                    #print "Eav_i",E_i,"Sens_av",Sens_av
                    Sens_arr.append(Sens_av[0])
                elif(Esup_i*Einf_i >0):
                    print "ERROR: the chosen energy values are out of available range"
                    exit
                else:
                    if(Esup_i== 0):
                        print "Esup_i (should be 0): ",Esup_i
                        if (len(Einf)>1):
                            Einf_ii = Einf[len(Einf)-2]
                            DE = Einf_i - Einf_ii
                            Sens_inf_i = Sens[EnAv_A==Einf_i]
                            Sens_inf_ii = Sens[EnAv_A==Einf_ii]
                            DS = Sens_inf_i - Sens_inf_ii
                            DE_l = E_i - Einf_i
                            Sens_av = Sens_inf_i + DS/DE*DE_l
                            print "all the 3 qunatities which follow should be positive"
                            print DE
                            print DS
                            print DE_l
                        else:
                            Sens_av = Sens[EnAv_A==Einf_i]
                        print "Sens_av",Sens_av
                        Sens_arr.append(Sens_av[0])
                    else:
                        print "Einf_i (should be 0): ",Einf_i

                        if (len(Esup)>1):
                            Esup_ii = Esup[1]
                            DE = Esup_ii - Einf_i
                            Sens_sup_i = Sens[EnAv_A==Esup_i]
                            Sens_sup_ii = Sens[EnAv_A==Esup_ii]
                            DS = Sens_sup_ii - Sens_sup_i
                            DE_l = Esup_i - E_i
                            print "all the 3 quantities which follow should be positive"
                            print DE
                            print DS
                            print DE_l
                            Sens_av = Sens_sup_i - DS/DE*DE_l

                        else:
                            Sens_av = Sens[EnAv_A==Esup_i]
        
                        print "Sens_av",Sens_av
                        Sens_arr.append(Sens_av[0])
        Sens_arr = np.array(Sens_arr)
    # Write file
    Ebins_av1D = np.atleast_1d(Ebins_av)
    SENS_FileName = TEMP_dir+"/Sensitivity.txt"
    with open(SENS_FileName, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        if len(Ebins_av1D)==1:
            print "[Ebins_av,Sens_arr,Ebins_inf,Ebins_sup]",[Ebins_av,Sens_arr,Ebins_inf,Ebins_sup]
            
            writer.writerow([Ebins_av,Sens_arr,Ebins_inf,Ebins_sup])
        else:
            writer.writerows(zip(Ebins_av,Sens_arr,Ebins_inf,Ebins_sup))


###################################### MARK: STEP 5.: ABSORPTION AT THE OBSERVER #####################################################
print
print "###################################### STEP 5.: ABSORPTION AT THE OBSERVER #####################################################"




Inv_ans = True
ABS_FileName = ""

SENS_MW_THESEUS = 32./13.
### SURVname check THESEUS - BEACAUSE IT HAS ALREADY ABSOPTION INCLUDED> FIXME: only 1 band enabled #####
if "THESEUS" in SURVname:
    Inv_ans = False
    print "# THESEUS sensitivity alreadt contains absorption, default configuration is for out the galactic plane with column density NH~ 5x10^20 /cm^2"
    print "# The second option for galactic plane is NH ~ 10^22 /cm^2"
    ABS_MW =raw_input('# Would you like to include ABSORPTION in Milky Way? ')
    if (ABS_MW == "yes" or ABS_MW == "YES" or ABS_MW == "Y" or ABS_MW == "y"):
        fracGW = input('# Write fraction of galactic observations (1 = 100% > only looking at galactic plane, 0 = 0% back to default configuration)? ')
        if fracGW>0:
            dataSENS = np.genfromtxt(SENS_FileName)
            SENS = dataSENS[1]
            SENS = SENS*(1.-fracGW)+SENS*fracGW*SENS_MW_THESEUS
            print "SENS",SENS
            CommLINE = "rm "+SENS_FileName
            os.system(CommLINE)
            with open(SENS_FileName, 'w') as f:
                writer = csv.writer(f, delimiter='\t')
                print "[Ebins_av,Sens_arr,Ebins_inf,Ebins_sup]",[Ebins_av,SENS,Ebins_inf,Ebins_sup]
                writer.writerow([Ebins_av,SENS,Ebins_inf,Ebins_sup])



r5wqABScorr_SENS = SensFile
ABS_TITLE =""
while(Inv_ans):
    
    ABS_MW =raw_input('# Would you like to add any model for ABSORPTION in Milky Way? ')
    if (ABS_MW == "yes" or ABS_MW == "YES" or ABS_MW == "Y" or ABS_MW == "y"):
        En_Sens= np.genfromtxt(SensFile)
        Inv_ans = False
        Inv_mod = True
        EnAv_A = np.zeros(np.shape(En_Sens))
        Sens_A = np.zeros(np.shape(En_Sens))
        if len(dim_sens) ==1:
            EnAv_A = En_Sens[2]
            Sens_A = En_Sens[3]
        elif len(dim_sens) >1:
            EnAv_A = En_Sens[:,2]
            Sens_A = En_Sens[:,3]
        else:
            print "ERROR: no sensitivity data"
            exit()
        
        ABS_MWDir = prefix+"Absorption/MW"
        while (Inv_mod):
            print "# Now choose between analytic prescriptions or data"
            ABS_mod =input('# Which absorption model would you like to use DATA from file (enter 1) or ANALYTIC (enter 0)? ')
            ABS_arr = []
            # ABSORPTION FROM FILE
            if (ABS_mod):
                Inv_mod = False
                invalid_input = True
                while invalid_input:
                    print "# Choose among the available files: "
                    ABS_MWDir = ABS_MWDir+"/FILES/"
                    listName = ABS_MWDir+"list_ABS_MW.txt"
                    CommLINE = "ls "+ABS_MWDir+"*.txt > "+listName
                    os.system(CommLINE)
                    MODLIST0 = []
                    with open(listName,'rb') as flist:
                        data_list = csv.reader(flist) #read text file with csv
                        count = 0
                        for row in data_list:
                            print "# n: ",count,". model:", row[0]
                            count = count+1
                            MODLIST0.append(row)
                    MODLIST0 = np.array(MODLIST0)
                    if count ==1:
                        MOD_ABSname =MODLIST0[0][0]
                        invalid_input = False
                    else:
                        enter_MODABS =input('Write the number correspondent to the desired survey: ')
                        while invalid_input:
                            if enter_MODABS == 0:
                                MOD_ABSname=MODLIST0[0][0]
                                invalid_input = False
                            elif enter_MODABS == 1:
                                MOD_ABSname=MODLIST0[1][0]
                                invalid_input = False
                            else:
                                print "# Unknown survey, choose one from the available list or add new entries"
                    ABSname = MOD_ABSname
                    ABS_TITLE ="FILE: "+ABSname
                    CommLINE = "rm "+listName
                    os.system(CommLINE)
            
                    #print "ABSname",ABSname
        
                # check the consistency of file range of energy with sensitivity
                En_Abs= np.genfromtxt(ABSname)
                EnAv_Abs_A = En_Abs[:,0]
                #print "EnAv_Abs_A ",EnAv_Abs_A
                ABS_MW_data = En_Abs[:,1]
                #print "ABS_MW_data",ABS_MW_data
                
                #Linear interpolation#
                
                #for j in range(len(Ebins_av)):
                #    E_j = Ebins_av[j]
                #    ABS_j = ABS_MW_data[EnAv_Abs_A==E_j]
                #    dim = np.shape(ABS_j)
                #    if dim[0] !=0:
                #        #print "if"
                #        ABS_arr.append(ABS_j[0])
                #    else:
                #        #print
                #        Einf = EnAv_Abs_A[EnAv_Abs_A<E_j]
                #        Einf = np.sort(Einf)
                #        Einf_i = Einf[len(Einf)-1]
                #        Esup = EnAv_Abs_A[EnAv_Abs_A>E_j]
                #        Esup = np.sort(Esup)
                #        Esup_i = Esup[0]
                #       #print "Einf_i",Einf_i,"Esup_i",Esup_i
                #        Abs_s = ABS_MW_data[EnAv_Abs_A==Esup_i]
                #        Abs_i = ABS_MW_data[EnAv_Abs_A==Einf_i]
                #        #print "Abs_s",Abs_s,"Abs_i",Abs_i
                #        #print "coefficient ",(Abs_s- Abs_i)/(Esup_i - Einf_i)
                #        #print "DE ",(E_j - Einf_i)
                #        Abs_av = Abs_i+(Abs_s- Abs_i)/(Esup_i - Einf_i)*(E_j- Einf_i)
                #        #print "Eav_i",E_j,"Abs_av",Abs_av
                #        ABS_arr.append(Abs_av[0])
                
                for j in range(len(Ebins_av)):
                    E_i = Ebins_av[j]
                    ABS_j= ABS_MW_data[EnAv_Abs_A==E_i]
                    dim = np.shape(ABS_j)
                    if dim[0] !=0:
                        #print "if"
                        ABS_arr.append(ABS_j[0])
                    else:
                        print
                        Einf = EnAv_Abs_A[EnAv_Abs_A<=E_i]
                        Einf = np.sort(Einf)
                        Einf_i = 0
                        if len(Einf)>0:
                            Einf_i = Einf[len(Einf)-1] #error
                        
                        Esup = EnAv_Abs_A[EnAv_Abs_A>E_i]
                        Esup = np.sort(Esup)
                        Esup_i = 0
                        if len(Esup)>0:
                            Esup_i = Esup[0]
                            #print "Einf_i",Einf_i,"Esup_i",Esup_i
                        
                        if (Esup_i >0 and Einf_i >0):
                            Abs_s = ABS_MW_data[EnAv_Abs_A==Esup_i]
                            Abs_i = ABS_MW_data[EnAv_Abs_A==Einf_i]
                            #print "Sens_s",Sens_s,"Sens_i",Sens_i
                            Abs_av = Abs_i+(Abs_s- Abs_i)/(Esup_i - Einf_i)*(E_i - Einf_i)
                            print "Abs_av",Abs_av
                            #print "Eav_i",E_i,"Sens_av",Sens_av
                            ABS_arr.append(Abs_av[0])
                        elif(Esup_i*Einf_i >0):
                            print "ERROR: the chosen energy values are out of available range"
                            exit()
                        else:
                            if(Esup_i== 0):
                                print "Esup_i (should be 0): ",Esup_i
                                if (len(Einf)>1):
                                    Einf_ii = Einf[len(Einf)-2]
                                    DE = Einf_i - Einf_ii
                                    Abs_inf_i = ABS_MW_data[EnAv_Abs_A==Einf_i]
                                    Abs_inf_ii = ABS_MW_data[EnAv_Abs_A==Einf_ii]
                                    DS = Abs_inf_i - Abs_inf_ii
                                    DE_l = E_i - Einf_i
                                    Abs_av = Abs_inf_i + DS/DE*DE_l
                                    print "all the 3 qunatities which follow should be positive"
                                    print DE
                                    print DS
                                    print DE_l
                                else:
                                    Abs_av = ABS_MW_data[EnAv_A==Einf_i]
                                print "Abs_av",Abs_av
                                ABS_arr.append(Abs_av[0])
                            else:
                                print "Einf_i (should be 0): ",Einf_i
                        
                                if (len(Esup)>1):
                                    Esup_ii = Esup[1]
                                    DE = Esup_ii - Einf_i
                                    Abs_sup_i = ABS_MW_data[EnAv_Abs_A==Esup_i]
                                    Abs_sup_ii = ABS_MW_data[EnAv_Abs_A==Esup_ii]
                                    DS = Abs_sup_ii - Abs_sup_i
                                    DE_l = Esup_i - E_i
                                    print "all the 3 qunatities which follow should be positive"
                                    print DE
                                    print DS
                                    print DE_l
                                    Abs_av = Abs_sup_i - DS/DE*DE_l
                        
                                else:
                                    Abs_av = ABS_MW_data[EnAv_Abs_A==Esup_i]
                        
                                print "Abs_av",Abs_av
                                ABS_arr.append(Abs_av[0])
                                
                ABS_arr = np.array(ABS_arr)

    
            #MARK: ABSORPTION ANALYTIC PRESCRIPTION: add once available- now it just writes all 1s
            elif (not ABS_mod):
                Inv_mod = False
                ABS_TITLE ="ANALITIC PRESCRIPTION"
                print "# WARNING: AT THE MOMENT AVAILABLE ONLY FOR SOFT X-RAY, accurate for RANGE: 0.030 - 10 keV, no absorption is assumed for higer energies"
                WEIGHTED_AV_NHg = 0.0
                TESTloc = True
                if  Emax< 0.030 or Emin > 10.0:
                    print "WARNING: instrument range of energy out of the available model range"
                    print "WARNING: NO galactic absorption will be assumed"
                else:
                    WRITEabs = True
                    while(WRITEabs):
                        ansWRITEabs =raw_input('# Would you like to type the value of hydrogen column to use? ')
                        if (ansWRITEabs == "YES" or ansWRITEabs == "yes" or ansWRITEabs == "Y" or ansWRITEabs == "y"):
                            WRITEabs = False
                            ABS_given = input('# Write the column density value in 1e20 $cm^{-2}$ ')
                            WEIGHTED_AV_NHg = ABS_given*10**(20)
                        elif(ansWRITEabs == "NO" or ansWRITEabs == "no" or ansWRITEabs == "N" or ansWRITEabs == "n"):
                            WRITEabs = False
                            print "To evaluate the GALACTIC ABSORPTION saprEMo marginalises NH (hydrogen column) over the sky-locations according to the observations in the survey"
                            # CHOOSE - OPEN THE FILE NH
                            #NH file from https://arxiv.org/pdf/astro-ph/0504140.pdf
                            NHprefix= ABS_MWDir+"/NH/"
                            listNameABS = ABS_MWDir+"list_NHprop.txt"
                    
                            CommLINE = "ls "+NHprefix+"> "+listNameABS
                            os.system(CommLINE)
                    
                            with open(listNameABS,'rb') as flistNH:
                                #next(f) #skip heading row in text file (I cannot get csv.Dictreader instead to work as an alternative to this step)
                                data_listNH = csv.reader(flistNH) #read text file with csv
                                count = 0
                                MOD_NH = []
                                for row in data_listNH:
                                    print "# n: ",count,". NH file:", row[0]
                                    count = count+1
                                    MOD_NH.append(row[0])
                
                            CommLINE = "rm "+listNameABS
                            os.system(CommLINE)
                    
                    
                            if count ==1:
                                NHprefixfile = NHprefix + row[0]
                                print NHprefixfile
                                exit()
                            elif count ==0:
                                print "ERROR: no entry files for NH calculation"
                                exit()
                            else:
                                modelMWabs = False
                                while (not modelMWabs):
                                    MWabsN = input('# Select the desired model: ')
                                    if MWabsN == 0:
                                        NHfile = MOD_NH[0]
                                        modelMWabs = True
                                        print "# you choose: ",NHfile
                                    elif MWabsN == 1:
                                        NHfile = MOD_NH[1]
                                        modelMWabs = True
                                        print "# you choose: ",NHfile
                                    else:
                                        print "Invalid entry, choose among: "
                                        num = 0
                                        for names in MOD_NH:
                                            print "# n: ",num,". NH file: ", MOD_NH[0]
                                            num = num+1
                    
                            NameGWABSfile = NHprefix+str(NHfile)
                            print "Name of file for Milky Way absorption: ",NameGWABSfile
                    
                            # elif () :MODLIST.append(row)
                            obs_data = fits.open(NameGWABSfile)
                            Data = obs_data[0].data
                            dim = np.shape(Data)
                            hdu = obs_data[0]
                    
                            # CALCULATING NHI according to http://www.cv.nrao.edu/course/astr534/HILine.html
                    
                            NHI = (hdu.data[0,:,:])*1.823*1e18*1.030571969000*0.5
                            last_i = -1
                            for i in range(1,dim[0]-1):
                                #for i in range(1,2):
                                if (-458605+i*1030.571969<=399000 and -458605+i*1030.571969>=-399000):
                                    NHI=NHI+(hdu.data[i,:,:])*1.823*1e18*1.030571969000 # n_H [1/cm^2] = 1.82*10^18*int(T/[K]*dv/[km/s]) ; dv = 1m/s
                                    last_i = i
                    
                            NHI=NHI-(hdu.data[last_i,:,:])*1.823*1e18*1.030571969000*0.5
                            """
                                # trapezoidal rule:
                                NHI = (hdu.data[0,:,:])*0.5*1.823*1e18*1.030571969000
                    
                                for i in range(1,dim[0]-2):
                                #for i in range(1,2):
                                if (-458605+i*1030.571969<=399000 and -458605+i*1030.571969>=-399000):
                                NHI=NHI+(hdu.data[i,:,:])*1.823*1e18*1.030571969000 # n_H [1/cm^2] = 1.82*10^18*int(T/[K]*dv/[km/s]) ; dv = 1m/s
                                NHI=NHI+(hdu.data[dim[0]-2,:,:])*0.5*1.823*1e18*1.030571969000
                            """
                            # NB: i valori di NH non sono gli stessi del tool, ma pare introvabile come siano stati ottenuti e anche il file che dono di aver messo a disposizione con questi valori non esiste piu - ci siamo con un errore del ~20%
                            #plt.figure()
                            #plt.imshow(np.log(NHI), origin='lower')
                            #plt.show()
                    
                            Glong = 180 - 0.5*np.arange(dim[2])
                            Glat = -90.0 +0.5*np.arange(dim[1])
                    
                            # DEFAULT VALUES FOR NC NH2 NH2MAX and alpha
                    
                            print "Calculating controbution of NH2 and total galactic NHg using equation 5) amd 6) of http://www.star.le.ac.uk/zrw/xabs/willingale_xabsgal.pdf"
                            # http://www.star.le.ac.uk/zrw/xabs/willingale_xabsgal.pdf
                    
                            NH2max = 7.5*1e20 # cm^-2
                            Nc = 2.37*1e21
                            alpha = 2.
                    
                            print "--- AS FROM THE REFERNCE PAPER, the tool adopts the following values: "
                            print "NH2max: ",NH2max," cm^(-2)"
                            print "Nc: ",Nc," cm^(-2)"
                            print "alpha: ",alpha
                    
                            CHANGEval = input('# Type 1 if you like to change any of these values, 0 otherwise: ')
                            while (CHANGEval):
                                NH2max =input('# NH2max in 1e20 [cm^(-2)]: ')
                                Nc =input('# Nc in 1e20 [cm^(-2)]: ')
                                alpha = 2.
                                print "You set the following values: "
                                print "NH2max: ",NH2max," cm^(-2)"
                                print "Nc: ",Nc," cm^(-2)"
                                print "alpha: ",alpha
                                CHANGEval = input('# Type 0 if you are happy with them: ')
                    
                            # CALCULATING MAP of NH tot for the galaxy
                    
                            NH2 = NH2max*(1. - np.exp(-NHI/Nc))**alpha
                            NHg = NHI + 2.*NH2
                    
                            # INTERGRATION ALONG X AXIS
                            AV_Ng = np.zeros(len(Glat))
                            ones_k = np.ones(len(Glong))
                            ones_k =np.transpose(np.matrix(ones_k))
                            for k in range(len(Glat)):
                                arr_k = np.matrix(NHg[k,:])
                                sum_k = arr_k*ones_k
                                AV_Ng[k] = sum_k
                            AV_Ng = AV_Ng/len(Glong)

                            #plt.figure()
                            #plt.plot(Glat,AV_Ng,'.')
                            #plt.show()
                            #exit()

                            while (TESTloc == True):
                                LOC_Q = raw_input('# Would you like this estimate to be based on actual survey location? ')
                                if ( LOC_Q == "yes" or LOC_Q  == "YES" or LOC_Q  == "Y" or LOC_Q == "y"):
                                    TESTloc = False
                            
                                    # OPEN LOCALIZATION DATA
                                    locDIR = SurvDIR+"OBSERVATIONloc/"
                                    print "In ",locDIR,": "
                                    commandLINE ="ls "+SurvDIR+"OBSERVATIONloc/"
                                    os.system(str(commandLINE))
                                    print "# BE AWARE: the input file should be a txt file with a single column: GALACTIC LATITUTES"
                                    LOCATfile = raw_input('# Write file name (and path if necessary) for file containing locations of survey observations ')
                            
                                    dataLOC_survey = np.genfromtxt(LOCATfile)
                            
                                    # BINNING DATA AND NH MAP
                            
                                    Nbin = Nbin_GalacticLAT

                                    binA = np.min(Glat) +(np.max(Glat) - np.min(Glat))/Nbin*np.arange(Nbin+1)
                                    binA[0] = binA[0] - 0.0001
                                    binA[len(binA)-1] = binA[len(binA)-1] + 0.0001
                            
                                    # DATA: building an histogram with the locations to know which are the preferred observed location
                                    H_bgalS =np.histogram(dataLOC_survey,bins = binA)
                                    HISTvalues = H_bgalS[0]
                                    HISTextrema = H_bgalS[1]
                                    ones_NOMR = np.transpose(np.matrix(np.ones(len(HISTvalues))))
                                    NORM = HISTvalues*ones_NOMR
                                    NORM = np.float(NORM)
                                    HISTvalues = HISTvalues/NORM
                        
                                    # NH
                                    av_NHgperBIN = np.zeros(Nbin)

                                    ll = 0
                                    AVtaken_tot  =0
                                    for kk in range(Nbin):
                                        AVtaken = 0
                                        for ll in range(len(AV_Ng)):
                                            #while(ll<len(AV_Ng)):
                                            if Glat[ll]>=HISTextrema[kk] and Glat[ll]<HISTextrema[kk+1]:
                                                av_NHgperBIN[kk] = av_NHgperBIN[kk] + AV_Ng[ll]
                                                AVtaken = AVtaken +1
                                                #else:
                                                #print 'll',ll,'Glat[ll]',Glat[ll],'HISTextrema[kk]',HISTextrema[kk]
                                        AVtaken_tot = AVtaken_tot+ AVtaken
                                        av_NHgperBIN[kk] = av_NHgperBIN[kk]/AVtaken
                                    
                                        
                                    # CALCULATING WEIGHTED AVERAGE OF NH
                                    ones_a = np.matrix(np.ones(len(av_NHgperBIN)))
                                    HISTvalues = np.matrix(HISTvalues)
                                    av_NHgperBIN = np.transpose(np.matrix(av_NHgperBIN))
                                    WEIGHTED_AV_NHg = HISTvalues*av_NHgperBIN
                                    AV =ones_a*av_NHgperBIN/Nbin
                                    AVreal = np.matrix(AV_Ng)*np.transpose(np.matrix(np.ones(len(AV_Ng))))/len(AV_Ng)
                                    print 'WEIGHTED_AV_NHg',WEIGHTED_AV_NHg,'AV',AV,'AVreal',AVreal,'Not the same because dividing a/3 +b/2 != (a+b)/5'
                            
                            
                                elif (LOC_Q== "no" or LOC_Q == "NO" or LOC_Q == "N" or LOC_Q == "n"):
                                    TESTloc = False
                                    # NH binned in the bands:
                                    av_NHgperBIN = np.zeros(NmainGalacticLAYERS)
                                    AVtaken_tot  =0
                                    upperLIM_NMGLs[len(upperLIM_NMGLs)-1] = upperLIM_NMGLs[len(upperLIM_NMGLs)-1] + 0.000001
                                    for kk in range(NmainGalacticLAYERS):
                                        AVtaken = 0
                                        for ll in range(len(AV_Ng)):
                                            if (Glat[ll]>=lowerLIM_NMGLs[kk] and Glat[ll]<upperLIM_NMGLs[kk]) or (Glat[ll]<-lowerLIM_NMGLs[kk] and Glat[ll]>=-upperLIM_NMGLs[kk]):
                                                av_NHgperBIN[kk] = av_NHgperBIN[kk] + AV_Ng[ll]
                                                AVtaken = AVtaken +1
                                            #else:
                                            #    print 'll',ll
                                        AVtaken_tot = AVtaken_tot+ AVtaken
                                        av_NHgperBIN[kk] = av_NHgperBIN[kk]/AVtaken
                                    print 'AVtaken_tot ',AVtaken_tot ,'len(Glat)',len(Glat)
                            
                                    # OBSERVATIONS binned in the bands:
                                    print " saprEMo considers ",NmainGalacticLAYERS," macro bands (degrees - galactic latitude): "
                                    for hh in range(NmainGalacticLAYERS):
                                        print "BAND ",hh,": [+/-",lowerLIM_NMGLs[hh],",+/-",np.round(upperLIM_NMGLs[hh], decimals=1),"]"
                                    print " Write the percentage of observation inside band: "
                                    percA = []
                                    for hh in range(NmainGalacticLAYERS):
                                        print hh, ": "
                                        perc_hh = input('' )
                                        percA.append(perc_hh)
                                    percA = np.array(percA)
                                    SUMcheck = 0.0
                                    for hh in percA:
                                        if hh< 0:
                                            print "ERROR: INVALID PERCENTAGE"
                                            exit()
                                        if hh > 100:
                                            print "ERROR: INVALID PERCENTAGE"
                                        SUMcheck = SUMcheck + hh

                                    if SUMcheck == 100:
                                        percA = percA/100.
                                    elif SUMcheck!=1.:
                                        print "ERROR: sum is not equal to 1"
                                        exit()

                                    # WEIGHTED AVERAGE
                                    WEIGHTED_AV_NHg  =np.matrix(percA)*np.transpose(np.matrix(av_NHgperBIN))
                                    print "WEIGHTED_AV_NHg",WEIGHTED_AV_NHg
                            
                        else:
                            print "Invalid answer"

                    ABS_arr = []
                    Ebins_inf1D = np.atleast_1d(Ebins_inf)
                    Ebins_sup1D = np.atleast_1d(Ebins_sup)
                    for ee in range(len(Ebins_inf1D)):
                        ENarray = np.linspace(Ebins_inf1D[ee],Ebins_sup1D[ee],5000)
                        
                        # CALCULATING ARRAY OF CROSS-SECTIONS
                        outFUNC = sigmaISMfuncANALITIC(ENarray)
                    
                    
                        energy = outFUNC[0]
                        sigmaISM= outFUNC[1]
                    
                        # CALCULATING ABSORPTION
                    
                        ABS_COEFF_ee = np.exp(-np.array(sigmaISM)*np.float(WEIGHTED_AV_NHg))
                        ABS_COEFF_ee = np.array(ABS_COEFF_ee)
                        ABS_COEFF_AV_BINee = np.average(ABS_COEFF_ee)
                        ABS_arr.append( ABS_COEFF_AV_BINee)
                    print "ABS: ", ABS_arr
                    
                #Inv_mod = False
                #for i in range(len(Ebins_av)):
                #    abs_i = Ebins_av[i]**(-8./3.)/Ebins_av[i]**(-8./3.)
                #                    ABS_arr.append(abs_i)
                #                ABS_arr = np.array(ABS_arr)
                #plt.figure()
                #plt.plot(Ebins_av,ABS_arr,'x')
                #plt.show()
                
            else:
                print "# Invalid choice, enter 0 or 1 only"
                print
            
            
            #write new file or copy it, with Emin, Emax, Eav, Sensitivity, Absorption
            
            ABS_FileName = TEMP_dir+"/LOCAL_ABS_MW.txt"
            with open(ABS_FileName, 'w') as f:
                writer = csv.writer(f, delimiter='\t')
                if len(Ebins_inf)==1:
                    ABS_arr = np.array(ABS_arr)
                    ABS_arr = np.atleast_1d(ABS_arr)
                    writer.writerow([Ebins_av,np.around(ABS_arr[0],2),Ebins_inf,Ebins_sup])
                else:
                    writer.writerows(zip(Ebins_av,np.around(ABS_arr,2),Ebins_inf,Ebins_sup))

                #writer.writerows(zip(Ebins_av,ABS_arr,Ebins_inf,Ebins_sup))

    elif (ABS_MW == "no" or ABS_MW == "NO" or ABS_MW == "N" or ABS_MW == "n"):
        Inv_ans = False
        ABS_FileName = ""

    else:
        print "# Invalide answer"

###################################### MARK: STEP 6.: EXPOSURE TIME DISTRIBUTION  #####################################################
print
print "###################################### STEP 6.: EXPOSURE TIME DISTRIBUTION #####################################################"
Q_dur = True
SD_T = np.sqrt(AvObsT) if SD_T ==0.0 else SD_T
Dist_OBSD = ""
while(Q_dur):
    print "# saprEMo allows to model the exposure time distribution with one of the following functions:"
    print "# 1. MaxwellBoltzmann"
    print "# 2. LogNormal "
    Dur_DIST =input('# Digit the number correspondent to the desire model ')
    if Dur_DIST ==1:
        Q_dur = False
        print "# By default the standard deviation, if not present in the survey file, is considered to be the square root of the average time"
        SD_prov = input('# If you would like to change this value type the new value, 0 otherwise ')
        if SD_prov >0:
            SD_T = SD_prov
        Dist_OBSD = "MaxwellBoltzmann"
    elif Dur_DIST ==2:
        Q_dur = False
        print "# By default the standard deviation, if not present in the survey file, is considered to be the square root of the average time"
        SD_prov = input('# If you would like to change this value type the new value, 0 otherwise ')
        if SD_prov >0:
            SD_T = SD_prov
        Dist_OBSD = "LogNormal"
    else:
        print "Invalide entry"

###################################### MARK: WRITE A SUMMARY - think about adding Nh for intergalacatic medium #####################################################

print
print
print "############################ CHOSEN SET UP: ###############################################"
print "# 1. You selected as EVENT RATE MODEL: ",MODname
print "#     ",percentage_EV*100," % of the systems emit the considered light curve"
if (MODname=="Dominik_et_al_2013"):
    print "# VERSION: ",VERS
    print "# SOURCE: ",SOURCE
if (MODname=="Ghirlanda_et_al_2016"):
    print "#      MODEL BASED ON SHORT GRBs"
    print "#          model case: ", VERS
    print "#          assumened opening angle: ",BEAMangle0
print "#       COSMOLOGY: OmegaK ",OmegaK,"OmegaM",OmegaM,"OmegaL",OmegaL
print "#"
print "# 2. You selected as LIGHT CURVE MODEL: ",MOD_LC
if (MOD_LC=="Siegel_Ciolfi_2016"):
    print "# SOURCE of light curve: ",SOURCE_LC
if ENABLE_LC_z==1:
    print "# The availability of more light curves as function of energy bin will allow a correct consideration of drifting band with z"
else:
    print "# Light curves as function of energy are not available: the maximum redshift at which consider the emission will be set differently"
    #if (Ans_moreLC=="yes" or Ans_moreLC== "YES" or Ans_moreLC== "Y" or Ans_moreLC== "y"):
    #print "# The availability of more light curves as function of energy bin will allow a correct consideration of drifting band with z"
    #LC_z = 1
    #else:
    #print "# Light curves as function of energy are not available: the maximum redshift at which consider the emission will be set differently"
    #LC_z= 0
    print "The limiting energy at low frequency is: ",LimE," keV"
print "# "
print "# ABSORPTION OF HOST GALAXY: ",Ans_HGabs
print "#"
print "# 3. You selected as DATA SURVEY: ",MODLIST[enter_MODSurvey]
print "#"
print "# 4. ENERGY BINS - INSTRUMENT RANGE: ",Ebins_av
print "#                number of bands: ", len(Ebins_av1D)
print "#                equally spaced: ",EQ_SP
print "#                the energy range of the chosen instrumen is: ",Emin," - ",Emax," keV"
print "#                sensitvity curve from: ",SensFile
print "#"
print "# 5. LOCAL ABSORPTION AT MILKY WAY: ",ABS_MW
if (ABS_MW == "yes" or ABS_MW == "YES" or ABS_MW == "Y" or ABS_MW == "y"):
    print "#    DATA from ",ABS_TITLE
print "#"
print "# 6. EXPOSURE TIME DISTRIBUTION modelled with ",Dist_OBSD," function with expected value ",AvObsT," and standard deviation ", SD_T
print "#"
print "####################################################################################"

#===============================================

NAME_OUT_SETUP = "INJ_FILES/"+RUNname+"_SetUp.txt"
f = open(NAME_OUT_SETUP, 'w')
print  >> f, "############################ CHOSEN SET UP: ###############################################"
print  >> f, "# 1. You selected as EVENT RATE MODEL: ",MODname
print >> f, "#     ",percentage_EV*100," % of the systems emit the considered light curve"
if (MODname=="Dominik_et_al_2013"):
    print >> f, "# VERSION: ",VERS
    print >> f, "# SOURCE: ",SOURCE
if (MODname=="Ghirlanda_et_al_2016"):
    print >> f, "#      MODEL BASED ON SHORT GRBs"
    print >> f, "#          model case: ", VERS
    print >> f, "#          assumened opening angle: ",BEAMangle0
print >> f, "#       COSMOLOGY: OmegaK ",OmegaK,"OmegaM",OmegaM,"OmegaL",OmegaL
print >> f, "#"
print >> f, "# 2. You selected as LIGHT CURVE MODEL: ",MOD_LC
if (MOD_LC=="Siegel_Ciolfi_2016"):
    print >> f, "# SOURCE of light curve: ",SOURCE_LC
if ENABLE_LC_z==1:
    print >> f, "# The availability of more light curves as function of energy bin will allow a correct consideration of drifting band with z"
else:
    print >> f, "# Light curves as function of energy are not available: the maximum redshift at which consider the emission will be set differently"
    #if (Ans_moreLC=="yes" or Ans_moreLC== "YES" or Ans_moreLC== "Y" or Ans_moreLC== "y"):
    #print "# The availability of more light curves as function of energy bin will allow a correct consideration of drifting band with z"
    #LC_z = 1
    #else:
    #print "# Light curves as function of energy are not available: the maximum redshift at which consider the emission will be set differently"
    #LC_z= 0
    print >> f, "The limiting energy at low frequency is: ",LimE," keV"
print >> f, "# "
print >> f, "# ABSORPTION OF HOST GALAXY: ",Ans_HGabs
print >> f, "#"
print >> f, "# 3. You selected as DATA SURVEY: ",MODLIST[enter_MODSurvey]
print >> f, "#"
print >> f, "# 4. ENERGY BINS - INSTRUMENT RANGE: ",Ebins_av
print >> f, "#                number of bands: ", len(Ebins_av1D)
print >> f, "#                equally spaced: ",EQ_SP
print >> f, "#                the energy range of the chosen instrumen is: ",Emin," - ",Emax," keV"
print >> f, "#                sensitvity curve from: ",SensFile
print >> f, "#"
print >> f, "# 5. LOCAL ABSORPTION AT MILKY WAY: ",ABS_MW
if (ABS_MW == "yes" or ABS_MW == "YES" or ABS_MW == "Y" or ABS_MW == "y"):
    print >> f, "#    DATA from ",ABS_TITLE
print >> f,"#"
print >> f,"# 6. EXPOSURE TIME DISTRIBUTION modelled with ",Dist_OBSD," function with expected value ",AvObsT," and standard deviation ", SD_T
print >> f, "#"
print >> f, "####################################################################################"

f.close()
#===============================================

Deny = True
NAMEprog = "saprEMo"
while Deny:
    ansCONF= raw_input('Would you like to continue and create the injection file? ')
    if (ansCONF== "yes" or ansCONF == "YES" or ansCONF == "Y" or ansCONF == "y"):
        Deny = False
        print "# You choose to proceed, thank you for using ",NAMEprog
        ###################################### MARK: WRITE INPUT FILE
        #NAME_OUT_FILE = raw_input('Enter the desired name of the injection (output) file (without extension): ')
        
        NAME_OUT_FILE = "INJ_FILES/"+RUNname+".txt"
        TABname = TEMP_dir+"/"+TABzVC
        if (ATT==""):
            LCFile = TEMP_dir+"/LC_"+MOD_LC+"_funcE_"+SOURCE_LC+".txt"
        else:
            LCFile = TEMP_dir+"/LC_"+MOD_LC+"_"+ATT+"_funcE_"+SOURCE_LC+".txt"
            CommLINE = "cp "+LC_E_Filename+" "+TEMP_dir+"/LC_"+MOD_LC+"_"+ATT+"_funcE_"+SOURCE_LC+".txt"
            print "CommLINE",CommLINE
            os.system(CommLINE)
        
        EN_LC_File = TEMP_dir+"/Energies_"+MOD_LC+"_"+SOURCE_LC+"_LC.txt"
        with open(NAME_OUT_FILE, 'w') as fout:
            writer = csv.writer(fout)
            colum  =[Hc, OmegaM, OmegaL,OmegaK,Emin, Emax, CoveredSky, AvObsT, percentage_EV, ENABLE_LC_z,LimE,TEMP_dir+"/"+NameFILE_TAB, RateFileName, LCFile, EN_LC_File, SENS_FileName,ABS_FileName,Dist_OBSD,SD_T]
            writer.writerows(zip(colum))
    elif (ansCONF == "no" or ansCONF== "NO" or ansCONF == "N" or ansCONF == "n"):
        Deny = False
        print "# You choose to stop, thank you for using ",NAMEprog
        command_CANC_setup = "rm "+ NAME_OUT_SETUP
        os.system(command_CANC_setup)
        exit()

command_FREE = "rm "+prefix+"*.txt"
os.system(command_FREE)
