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
import matplotlib
from scipy.integrate import quad
import csv
import astropy
import astropy.units as u
import os
#import sh
#import pexect
from astropy.cosmology import WMAP9 as cosmo
sys.path.insert(0, 'LIBRARY')
from functions import MergeBins
plt.rcParams['figure.figsize'] = (9.0, 6.0)
import math
import matplotlib.ticker as mtick
from scipy.stats import lognorm
font = {'family' : 'Times','weight' : 'normal','size':22}#font = {'family' : 'Times','weight' : 'bold','size':22}
import matplotlib.ticker as ticker
import latex
    
matplotlib.rc('font', **font)
matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

#print "fact 5", math.factorial(5),"5*4*3*2",5*4*3*2


# GLOBAL VARIABLEs #
####################
Ntrials = 200
ResCDF = 1000
nbin_loghisto = 50
minDplot = 0.01
LIN=0 #LIN = 1 for linear graphs LIN = 0 for log
EFF_DUR = 1 # equal 0 if instead we just want the duration calculated on the extrema

pc = 3.08567758130573*10**18 #cm
Mpc = 10**6*pc
c_m = 299792458 #m / s
c = c_m*100.
G = 6.67408* 10**(-11)# m3 kg^(-1) s-2
G_cgs = G*(10**2)**3/1000.
h = 4.135667662*10**(-15.) #eV*s
eV = 1.60218*10**(-12)#1eV to erg
h = h*eV #erg*s
KB =  1.38064852*10**(-16) #erg/K

# INPUTS #
##########


output_dir=sys.argv[1]
print "sys.argv[1]",output_dir
input_file=sys.argv[2]
print "sys.argv[2]",input_file
LC_times_file =sys.argv[3]
print "sys.argv[3]",LC_times_file
Dt_vis_file =sys.argv[4]
print "sys.argv[4]",Dt_vis_file
DIST_EXPT =sys.argv[5]
print "sys.argv[5]",DIST_EXPT
AV_EXPT=float(sys.argv[6])
print "sys.argv[6]",AV_EXPT
SD_EXPT=float(sys.argv[7])
print "sys.argv[7]",SD_EXPT
InFile =sys.argv[8]
print "sys.argv[8]",InFile
######## 20/08/2018
print
print "###### from INPUT FILE #########"
f = open(InFile)
rowARR = []
for l in f:
    row = l
    #print row
    rowARR.append(row)

rowARR = np.array(rowARR)
# a. COSMOLOGY
Hc = float(rowARR[0])
t_Hc = Mpc/(Hc*1e5)
#print t_Hc
t_Hc = cosmo.hubble_time
#print t_Hc
t_Hc = t_Hc.item()/u.Gyr
#print t_Hc
t_Hc = t_Hc*60*60*24*365*1e9
#print t_Hc

print "#"
print "# -------- COSMOLOGY -------"
print '#  Hc',Hc, ' Km/s per Mpc' #km/s per Mpc
OmegaM = float(rowARR[1])
OmegaV = float(rowARR[2])
Omegak = float(rowARR[3])
print '# OmegaM: ',OmegaM,' OmegaV: ',OmegaV, ' Omegak: ',Omegak

print "#"
D_H = c/(10**5.*Hc) #Mpc
D_Hcm = D_H*Mpc

###### MARK: FUNCTIONS ########

# create an array with instrument information
def integ1(x, OMf, Okf, OVf):
    out=(OMf*(1.+x)**3+Okf*(1+x)**2 +OVf)**(-0.5)
    return out


def E_z(OM,Ok,OV,zf):
    Ef = np.sqrt(OM*(1.+zf)**3+Ok*(1+zf)**2 +OV)
    return Ef

#line-of-sight comuving distance
def LoS_CD(zf,OM, Ok, OV):
    res = quad(integ1, 0.0,zf, args=(OM, Ok, OV))
    #print 'LoS_CD zf',zf
    CD0 = D_Hcm*res[0]
    return CD0

#Transverse comuving distance
def T_CD(zf,OM, Ok, OV):
    #print 'T_CD zf',zf
    DC = LoS_CD(zf,OM, Ok, OV)
    DM = 0.0
    if (Ok>0.000001):
        DM = D_Hcm/np.sqrt(Ok)*np.sinh(np.sqrt(Ok)*DC/D_Hcm)
    #print 'magg',DM
    elif (Ok<-0.00001):
        DM = D_Hcm/np.sqrt(np.abs(Ok))*np.sin(np.sqrt(np.abs(Ok))*DC/D_Hcm)
    #print 'min',DM
    else:
        DM = DC
    #print '0',DM
    return DM

def ComuvingVol(zf,OM, Ok, OV): # returns Distance in cm
    #print 'CVol zf',zf
    DM = T_CD(zf,OM, Ok, OV)
    
    if (Ok>0.000001):
        Vc = (4.*np.pi*(D_Hcm)**3)/(2.*Ok)*(DM/D_Hcm*np.sqrt(1.-Ok*DM*DM/(D_Hcm*D_Hcm)) - np.arcsinh(np.sqrt(Ok)*DM/D_Hcm)/np.sqrt(Ok))
    elif (Ok<-0.00001):
        Vc = (4.*np.pi*(D_Hcm)**3)/(2.*Ok)*(DM/D_Hcm*np.sqrt(1.-Ok*DM*DM/(D_Hcm*D_Hcm)) - np.arcsin(np.sqrt(np.abs(Ok))*DM/D_Hcm)/np.sqrt(np.abs(Ok)))
    else:
        Vc = 4./3.*np.pi*DM**3
    return Vc


def fromZtoDL(zf,OM, Ok, OV):  # takes Distance in cm
    
    DLf = (1.+zf)*T_CD(zf,OM, Ok, OV)
    return DLf



#######



# GLOBAL VARIABLES FROM INPUTS #
################################

nsigma = 1.+0.8*9./5*np.log10(Ntrials) # number od sigmas for N% prob
#print nsigma
minEXPT = np.max([1.,AV_EXPT-nsigma*SD_EXPT])
maxEXPT = AV_EXPT+nsigma*SD_EXPT

NBIN = 25
dEXPT= 0.5*(maxEXPT - minEXPT)/(NBIN-1)


# FUNCTION FOR DISTRIBUTION OF EXPOSURE TIMES #
###############################################
def cdf_MaxBol(av,x):
    # comulative distribution Maxwell Boltzmann
    a = av*0.5*np.sqrt(np.pi*0.5)
    cdf = scipy.special.erf(x/(np.sqrt(2.)*a)) - np.sqrt(2./np.pi)*x*np.exp(-0.5*x*x/(a*a))/a
    return cdf



def cdf_LogNorm(av,sd,x):
    # comulative distribution lognormal
    sigma = np.sqrt(np.log(sd*sd/(av*av) +1))
    
    mu = np.log(av) - sigma*sigma*0.5
    
    cdf = 0.5+0.5*scipy.special.erf((np.log(x)-mu)/(np.sqrt(2.)*sigma))
    return cdf


def MB_EXPT(M, nt):
    xAf = []
    
    cdf_try = np.zeros(ResCDF)
    xA = np.linspace(1.e-50,maxEXPT,ResCDF)

    cdf_try = cdf_MaxBol(M,xA)
    
    rA= np.random.uniform(0,1,nt)
    
    # sampling from inverse of a cumulative
    for i in range(nt):
        el = min(cdf_try, key=lambda x:abs(x-rA[i])) # take the x value which corresponds to the smaller difference between the random number rA and the cumulative
        ind = list(cdf_try).index(el)
        x_i = xA[ind]
        xAf.append(x_i)

    xAf = np.array(xAf)
    
    return xAf

def MC_EXPT(DIST, M, S, nt):
    nt_ok = 0
    xAf = []
    count_while = 0
    
    cdf_try = np.zeros(ResCDF)
    xA = np.linspace(1e-50,maxEXPT,ResCDF)
    DIST_MB = 'MaxwellBoltzmann'
    DIST_LN = 'LogNormal'
    if DIST == DIST_MB:
        cdf_try = cdf_MaxBol(AV_EXPT,xA)
    
    elif DIST == DIST_LN:
        cdf_try = cdf_LogNorm(AV_EXPT,SD_EXPT,xA)
    else:
        print"ERROR: UNRECOGNISED with DISTRIBUTION STRING: ",repr(DIST),"LG: ",repr(DIST_LN),"MB: ",repr(DIST_MB)
        exit()
    
    rA= np.random.uniform(0,1,nt)
    for i in range(nt):
        el = min(cdf_try, key=lambda x:abs(x-rA[i]))
        ind = list(cdf_try).index(el)
        x_i = xA[ind]
        xAf.append(x_i)

    count_while = count_while +1
    xAf = np.array(xAf)

    return xAf

#try_samp = MC_EXPT(DIST_EXPT,AV_EXPT,SD_EXPT,Ntrials)
####### method 2: Inverse Tranform Sampling #########

def MaxBol(av,x):
    a = av*0.5*np.sqrt(np.pi*0.5)
    dx = (maxEXPT - minEXPT)/len(x)
    x_norm = np.linspace(minEXPT,maxEXPT,len(x))
    p_norm = np.sqrt(2./np.pi)*x_norm*x_norm*np.exp(-x_norm*x_norm*0.5/(a*a))/(a*a*a)
    norm = np.matrix(dx*np.sqrt(2./np.pi)*x_norm*x_norm*np.exp(-x_norm*x_norm*0.5/(a*a))/(a*a*a))*np.transpose(np.matrix(np.ones(len(x))))
    p = np.sqrt(2./np.pi)*x*x*np.exp(-x*x*0.5/(a*a))/(a*a*a)
    p = p*dx/norm
    return p

def LogNorm(av,sd,x):
    
    sigma = np.sqrt(np.log(sd*sd/(av*av) +1))
    mu = np.log(av) - sigma*sigma*0.5
    x = np.linspace(minEXPT,maxEXPT,10000)
    p = 1./(x*sigma*np.sqrt(2.*np.pi))*np.exp(- 0.5*(np.log(x) - mu)**2./(sigma*sigma))
    p = p*dx
    return p

# READING INPUT FILES WITH DURATION #
#####################################
LC_times = np.genfromtxt(LC_times_file)

Dt_vis_zmin = np.genfromtxt(Dt_vis_file)

# change string part _zmin -> _zeff -> _zmax
Dt_vis_file.replace("zmin","zeff")
Dt_vis_zeff = np.genfromtxt(Dt_vis_file)

Dt_vis_file.replace("zeff","zmax")
Dt_vis_zmax = np.genfromtxt(Dt_vis_file)


DATA=np.genfromtxt(input_file)

zmin = DATA[:,0]
Ti_zmin = DATA[:,1]
Tf_zmin = DATA[:,2]
Np_zmin = DATA[:,3]
Nt_zmin = DATA[:,4]

    
zeff = DATA[:,5]
Ti_zeff = DATA[:,6]
Tf_zeff = DATA[:,7]
Np_zeff = DATA[:,8]
Nt_zeff = DATA[:,9]


zmax = DATA[:,10]
Ti_zmax = DATA[:,11]

Tf_zmax = DATA[:,12]
Np_zmax = DATA[:,13]
Nt_zmax = DATA[:,14]


timeLp = DATA[:,15]


# MARK: effective duration
def DurationEFFECTIVE(Tp_LCobs,Dobs, LCtime, FLf_zeff, LCdt_ind,PorT):
    
    LCdt_ind = np.array(LCdt_ind)
    LCdt_ind = LCdt_ind.astype(int)
    shapeDobs= np.shape(Dobs)

    nzbin = shapeDobs[0]
    nt = shapeDobs[1]
    
    # loop over z
    Ds = []
    Fmax =[]
    for i in range(nzbin):
        FL_i = FLf_zeff[i,:]
        LC_i = LCtime[i,:]# LC*(1+z)
    
        LCdt_i= LC_i[1:] - LC_i[:len(LC_i)-1]
        LCtime_zi = LC_i
        LCdt_ind_i = LCdt_ind[i,:]
        # step 1) choose random starting time
        #          for peaks would be between Tp and Tp- EXPOSURE
        LC_ind_i_eff = LCdt_ind_i[LCdt_ind_i>0.0]
        LC_ind_i_eff = LC_ind_i_eff.astype(int)
        if not len(LC_ind_i_eff):
            D_i = np.zeros(nt)
            continue
        D_i = []
        Fmax_i = []
        Dobs_i = Dobs[i]
        for j in range(nt):
            if PorT =="P": # peaks
                Ts_ij = np.random.uniform(Tp_LCobs[i] - Dobs_i[j], Tp_LCobs[i])
            else:
                minStartTime = LCtime_zi[LCdt_ind_i[1]]
                lastIND = LCdt_ind_i[LCdt_ind_i>0]
                lastIND = lastIND[-2]
                maxStartTime = LCtime_zi[lastIND]
                Ts_ij = np.random.uniform(minStartTime- Dobs_i[j],maxStartTime - Dobs_i[j])
                if Ts_ij>Tp_LCobs[i] - Dobs_i[j]:
                   Ts_ij = Ts_ij+Dobs_i[j]
            Tf_ij = Ts_ij+Dobs_i[j]
            
            if len(LCtime_zi[LCtime_zi<Ts_ij])>0:
                LCtime_zi_first_el = LCtime_zi[LCtime_zi<Ts_ij].max()
            else:
                LCtime_zi_first_el =LCtime_zi[0]

            if len(LCtime_zi[LCtime_zi>(Tf_ij)])>0:
                LCtime_zi_last_el = LCtime_zi[LCtime_zi>(Tf_ij)].min()
            else:
                LCtime_zi_last_el = LCtime_zi[-1]

            INDfirst = list(LCtime_zi).index(LCtime_zi_first_el)
            INDlast = list(LCtime_zi).index(LCtime_zi_last_el) #picks index of right extreme of dt where is the final time Ts+Dobs
            INDlast = INDlast -1 # picks the index of the dt which cointains final time Ts+Dobs
            #summing up dt contained in the extrame (expluding the ones containing startinf time and final time)
            IND_arr = np.arange(len(FL_i))
            Fmax_ij = FL_i[IND_arr<=INDlast]
            IND_arr = IND_arr[IND_arr<=INDlast]
            Fmax_ij = Fmax_ij[IND_arr>=INDfirst]
            Fmax_ij = np.max(Fmax_ij)
            sumIND = LC_ind_i_eff[LC_ind_i_eff>INDfirst]
            sumIND = sumIND[sumIND<INDlast]
            sumTIME = 0.0
            if len(sumIND)>0:
            
                LCdt_iSUM= LCdt_i[sumIND]
                sumTIME = np.sum(LCdt_iSUM)

            if INDfirst in list(LC_ind_i_eff):
                sumTIME = sumTIME+ (LCtime_zi[INDfirst +1]-Ts_ij)
            if INDlast in list(LC_ind_i_eff):
                sumTIME = sumTIME +(-LCtime_zi[INDlast]+Tf_ij)

            if INDlast == INDfirst and (INDlast in list(LC_ind_i_eff)):
                sumTIME = Dobs_i[j]
            D_i.append(sumTIME)
            Fmax_i.append(Fmax_ij)

        Ds.append(D_i)
        Fmax.append(Fmax_i)
    Ds = np.array(Ds)
    Fmax = np.array(Fmax)
    return [Ds,Fmax]


def DurationEXTREMA_p(Tif,Tff,Tp_LC,Dobs, nt):
    # loop over z
    Ds = []
    count = 0
    for i in range(len(Tif)):
        # step 1) choose random starting time
        #          for peaks would be between Tp and Tp- EXPOSURE
        Ts_i = []
        for j in range(nt):
        
            Ts_ij = np.random.uniform(Tp_LC[i] - Dobs[j], Tp_LC[i])
            Ts_i.append(Ts_ij)
        
        Ts_i = np.array(Ts_i)
        Ti_i = Tif[i]*np.ones(nt)
        # effective starting of signal
        Ts_ieff = np.maximum(Ts_i, Ti_i)
        
        Tstop_i= Ts_i+Dobs
        
        # effective end of signal
        Tf_i = Tff[i]* np.ones(nt)
        Tf_ieff = np.minimum(Tstop_i, Tf_i)
        Ds_i = Tf_ieff - Ts_ieff
        for j in range(nt):
            
            if Ds_i[j]<=0.0:
                print Ts_ieff[j],Ts_i[j], Ti_i[j],"Tp_LC[i] - Dobs[j]",Tp_LC[i] - Dobs[j],"Tp_LC[i]",Tp_LC[i],"Dobs[j]",Dobs[j]
                print "Tf_ieff ",Tf_ieff[j],"Tstop_i",Tstop_i[j], "Tf_i",Tf_i[j]
                print "Ts_ieff",Ts_ieff[j],"Ts_i",Ts_i[j], "Ti_i",Ti_i[j]
                print "Tf_ieff ",Tf_ieff[j], "Ts_ieff",Ts_ieff[j]
                count = count +1
                exit()
        
        Ds.append(Ds_i)
    Ds = np.array(Ds)
    return Ds

def DurationEXTREMA_t(Tif,Tff,Tp_LC,Dobs, nt):
    # loop over z
    Ds = []
    for i in range(len(Tif)):
        # step 1) choose random starting time
        #          for peaks would be between Tp and Tp- EXPOSURE
        Ts_i = []
        for j in range(nt):
            Ts_ij = np.random.uniform(Tif[i] - Dobs[j], (Tp_LC[i] - Dobs[j]) + (Tff[i] - Tp_LC[i])) # aggiungere altro pezzo
            offset = Dobs[j]
            
            Ts_ij = Ts_ij if Ts_ij<=Tp_LC[i] - Dobs[j] else Ts_ij + offset
            
            Ts_i.append(Ts_ij)
        
        Ts_i= np.array(Ts_i)
        Ti_i = Tif[i]*np.ones(nt)
        
        # effective starting of signal
        Ts_ieff = np.maximum(Ts_i, Ti_i)
        Tstop_i= Ts_i+Dobs
        
        # effective end of signal
        Tf_i = Tff[i]* np.ones(nt)
        Tf_ieff = np.minimum(Tstop_i, Tf_i)
        Ds_i = Tf_ieff - Ts_ieff
        Ds.append(Ds_i)
    
    Ds = np.array(Ds)
    # step 1) choose random starting time
    return Ds

# MAIN #
########
lenz = len(zmin)


DIST_MB = 'MaxwellBoltzmann'
DIST_LN = 'LogNormal'

if DIST_EXPT == DIST_MB:
    EXPOSURE = MB_EXPT(AV_EXPT, Ntrials)

if DIST_EXPT == DIST_LN:
    
    m2 = np.exp(np.log(AV_EXPT)*2.)*np.exp(SD_EXPT*SD_EXPT)
    m = np.power(m2,0.5)
    v2 =m2*(np.exp(SD_EXPT*SD_EXPT)-1.)
    v =np.power(v2,0.5)
    
    mu = np.log(AV_EXPT/(1.+(SD_EXPT*SD_EXPT/AV_EXPT*AV_EXPT))**0.5)
    s =(np.log(1.+(SD_EXPT*SD_EXPT/AV_EXPT*AV_EXPT)))**0.5
    EXPOSURE = np.random.lognormal(mu,s, Ntrials)





Ti_zminOBS = np.multiply(Ti_zmin,(1.+zmin))
Tf_zminOBS = np.multiply(Tf_zmin,(1.+zmin))
timeLpzmin = np.multiply(timeLp,(1.+zmin))

Ti_zeffOBS = np.multiply(Ti_zeff,(1.+zeff))
Tf_zeffOBS = np.multiply(Tf_zeff,(1.+zeff))
timeLpzeff = np.multiply(timeLp,(1.+zeff))

Ti_zmaxOBS = np.multiply(Ti_zmax,(1.+zmax))
Tf_zmaxOBS = np.multiply(Tf_zmax,(1.+zmax))
timeLpzmax = np.multiply(timeLp,(1.+zmax))

D_sf_zmin_p = []
D_sf_zmin_t = []
D_sf_zmax_p = []
D_sf_zmax_t = []
D_sf_zeff_p = []
D_sf_zeff_t = []

FL_zmin_p = []
FL_zmin_t = []
FL_zmax_p = []
FL_zmax_t = []
FL_zeff_p = []
FL_zeff_t = []

#
##### FLUXES
# for peaks
FLUX_file = output_dir + "/txt4ANALYSIS/FL_zeff.txt"
FLUX_data = np.genfromtxt(FLUX_file)
Deff = []
Dmin = []
Dmax = []
zmin[0] = 1e-50
for i in range(len(zeff)):
    Deff.append(fromZtoDL(zeff[i],OmegaM, Omegak, OmegaV))
    Dmin.append(fromZtoDL(zmin[i],OmegaM, Omegak, OmegaV))
    Dmax.append(fromZtoDL(zmax[i],OmegaM, Omegak, OmegaV))

Deff = np.array(Deff)
Dmin = np.array(Dmin)
Dmax = np.array(Dmax)

Dmin2 = np.multiply(Dmin,Dmin)
Deff2 = np.multiply(Deff,Deff)
Dmax2 = np.multiply(Dmax,Dmax)

fact_zmin = np.multiply(Deff2,1./Dmin2)
fact_zmax = np.multiply(Deff2,1./Dmax2)
#print "Deff[0], np.shape(Deff), np.shape(Deff2), np.shape(fact_zmin)"
#print Deff[0], np.shape(Deff), np.shape(Deff2), np.shape(fact_zmin)

FLUX_data_zmin = FLUX_data * fact_zmin[:, np.newaxis]
FLUX_data_zmax = FLUX_data * fact_zmax[:, np.newaxis]
#FLUX_data_zmin  = np.multiply(FLUX_data, fact_zmin)
#FLUX_data_zmax = np.multiply(FLUX_data, fact_zmax)

#print "Deff[0], np.shape(Deff), np.shape(Deff2), np.shape(fact_zmin), np.shape(FLUX_data_zmin)"
#print Deff[0], np.shape(Deff), np.shape(Deff2), np.shape(fact_zmin), np.shape(FLUX_data_zmin)


FL_peaks = np.amax(FLUX_data,1)


if EFF_DUR ==1:
    EXPOSURE_eff = MC_EXPT(DIST_EXPT,AV_EXPT, SD_EXPT, Ntrials*len(zeff))
    EXPOSURE_eff = np.reshape(EXPOSURE_eff,(len(zeff),Ntrials))
    
    LCtime_zeff= np.array(np.transpose(np.matrix(1.+zeff))*np.matrix(LC_times))
    LCtime_zmin= np.array(np.transpose(np.matrix(1.+zmin))*np.matrix(LC_times))
    LCtime_zmax= np.array(np.transpose(np.matrix(1.+zmax))*np.matrix(LC_times))
    
    OUT_zmin_p = DurationEFFECTIVE(timeLpzmin,EXPOSURE_eff, LCtime_zmin, FLUX_data_zmin, Dt_vis_zmin,"P")
    D_sf_zmin_p = OUT_zmin_p[0]
    FL_zmin_p = OUT_zmin_p[1]
    
    OUT_zmin_t = DurationEFFECTIVE(timeLpzmin,EXPOSURE_eff, LCtime_zmin, FLUX_data_zmin, Dt_vis_zmin,"T")
    D_sf_zmin_t = OUT_zmin_t[0]
    FL_zmin_t = OUT_zmin_t[1]
    
    OUT_zeff_p = DurationEFFECTIVE(timeLpzeff,EXPOSURE_eff, LCtime_zeff, FLUX_data, Dt_vis_zeff,"P")
    D_sf_zeff_p = OUT_zeff_p[0]

    FL_zeff_p = OUT_zeff_p[1]
    OUT_zeff_t = DurationEFFECTIVE(timeLpzeff,EXPOSURE_eff, LCtime_zeff, FLUX_data, Dt_vis_zeff,"T")
    D_sf_zeff_t = OUT_zeff_t[0]
    FL_zeff_t = OUT_zeff_t[1]
    
    OUT_zmax_p = DurationEFFECTIVE(timeLpzmax,EXPOSURE_eff, LCtime_zmax, FLUX_data_zmax, Dt_vis_zmax,"P")
    D_sf_zmax_p = OUT_zmax_p[0]
    FL_zmax_p = OUT_zmax_p[1]
    
    OUT_zmax_t = DurationEFFECTIVE(timeLpzmax,EXPOSURE_eff, LCtime_zmax, FLUX_data_zmax, Dt_vis_zmax,"T")
    D_sf_zmax_t = OUT_zmax_t[0]
    FL_zmax_t = OUT_zmax_t[1]
    
    #print AV_EXPT , SD_EXPT, np.average(EXPOSURE_eff), np.std(EXPOSURE_eff), "av D_zeff p ",np.average(D_sf_zeff_p), "av D_zeff t ", np.average(D_sf_zeff_t)
    
else:
    EXPOSURE = MC_EXPT(DIST_EXPT,AV_EXPT, SD_EXPT, Ntrials)

    D_sf_zmin_p = DurationEXTREMA_p(Ti_zminOBS,Tf_zminOBS,timeLpzmin,EXPOSURE, Ntrials)
    D_sf_zmin_t = DurationEXTREMA_t(Ti_zminOBS,Tf_zminOBS,timeLpzmin,EXPOSURE, Ntrials)

    D_sf_zeff_p = DurationEXTREMA_p(Ti_zeffOBS,Tf_zeffOBS,timeLpzeff,EXPOSURE, Ntrials)
    D_sf_zeff_t = DurationEXTREMA_t(Ti_zeffOBS,Tf_zeffOBS,timeLpzeff,EXPOSURE, Ntrials)

    D_sf_zmax_p = DurationEXTREMA_p(Ti_zmaxOBS,Tf_zmaxOBS,timeLpzmax,EXPOSURE, Ntrials)
    D_sf_zmax_t = DurationEXTREMA_t(Ti_zmaxOBS,Tf_zmaxOBS,timeLpzmax,EXPOSURE, Ntrials)


                                        
######### zmin ########### Ti_zmin,Tf_zmin,timeLp,Dt_vis_zmin,Ntrials
Dfinal_zmin_p = D_sf_zmin_p #OUTzmin_p
Ffinal_zmin_p = FL_zmin_p

Nrel_p_zmin = np.transpose(np.matrix(Np_zmin/Ntrials))*np.ones(Ntrials)
Nrel_p_zmin = np.array(Nrel_p_zmin)

Nrel_t_zmin = np.transpose(np.matrix(Nt_zmin/Ntrials))*np.ones(Ntrials)
Nrel_t_zmin = np.array(Nrel_t_zmin)

Dfinal_zmin_t = D_sf_zmin_t# OUTzmin_t
Ffinal_zmin_t = FL_zmin_t
Dfinal_zmin_tot = np.concatenate([Dfinal_zmin_p,Dfinal_zmin_t])
Ffinal_zmin_tot = np.concatenate([Ffinal_zmin_p,Ffinal_zmin_t])

Nrel_tot_zmin =np.concatenate([Nrel_p_zmin,Nrel_t_zmin])

######### zmax ###########

Dfinal_zmax_p = D_sf_zmax_p#OUTzmax_p
Ffinal_zmax_p = FL_zmax_p

Dfinal_zmax_t = D_sf_zmax_t#OUTzmax_t
Ffinal_zmax_t = FL_zmax_t

Nrel_p_zmax = np.transpose(np.matrix(Np_zmax/Ntrials))*np.ones(Ntrials)
Nrel_p_zmax = np.array(Nrel_p_zmax)

Nrel_t_zmax = np.transpose(np.matrix(Nt_zmax/Ntrials))*np.ones(Ntrials)
Nrel_t_zmax = np.array(Nrel_t_zmax)

Dfinal_zmax_tot = np.concatenate([Dfinal_zmax_p,Dfinal_zmax_t])
Ffinal_zmax_tot = np.concatenate([Ffinal_zmax_p,Ffinal_zmax_t])

Nrel_tot_zmax =np.concatenate([Nrel_p_zmax,Nrel_t_zmax])

######### zeff ###########

Dfinal_zeff_p = D_sf_zeff_p#OUTzeff_p
Ffinal_zeff_p = FL_zeff_p

Dfinal_zeff_t = D_sf_zeff_t#OUTzeff_t
Ffinal_zeff_t  = FL_zeff_t

#print np.shape(Np_zeff)


Nrel_p_zeff = np.transpose(np.matrix(Np_zeff/Ntrials))*np.ones(Ntrials)
Nrel_p_zeff = np.array(Nrel_p_zeff)


Nrel_t_zeff = np.transpose(np.matrix(Nt_zeff/Ntrials))*np.ones(Ntrials)
sum0 = np.sum(Nrel_t_zeff) + np.sum(Nrel_p_zeff)
#print np.shape(Nrel_t_zeff), sum0


Nrel_t_zeff = np.array(Nrel_t_zeff)

Dfinal_zeff_tot = np.concatenate([Dfinal_zeff_p,Dfinal_zeff_t])
Ffinal_zeff_tot = np.concatenate([Ffinal_zeff_p,Ffinal_zeff_t])

Nrel_tot_zeff = np.concatenate([Nrel_p_zeff,Nrel_t_zeff]) #Nrel_p_zeff#

# MARK: FUTURE WORK - ADD EFFECTIVE TIME OVER THRESHOLD #
#########################################################


#building apposite directory
output_dir = output_dir+"/DURATION_DISTRIB"
if not(os.path.exists(output_dir)):
    CommLINE = "mkdir "+output_dir
    os.system(CommLINE)


# FIGURES #
###########


bin_zmin_i = np.zeros(20)
Hsum_zmin =  np.zeros(20)

bin_zeff_i = np.arange(math.floor(np.log10(np.min(Dfinal_zmax_tot[Dfinal_zmax_tot>0.]))),  math.ceil(np.log10(np.max(Dfinal_zmin_tot))))
bin_zeff_i = 10**bin_zeff_i

Hsum_zeff =  np.zeros(20)
Hsum_zeff_p =  np.zeros(20)
Hsum_zeff_t =  np.zeros(20)

bin_zmin_i = np.zeros(20)
Hsum_zmax =  np.zeros(20)


#print "np.shape(Dfinal_zmin_tot)",np.shape(Dfinal_zmin_tot),"np.shape(Nrel_tot_zmin)",np.shape(Nrel_tot_zmin),"len(zmin)",len(zmin)
#print "np.shape(Dfinal_zmax_tot)",np.shape(Dfinal_zmax_tot),"np.shape(Nrel_tot_zmax)",np.shape(Nrel_tot_zmax),"len(zmin)",len(zmin)
#print "np.shape(Dfinal_zeff_tot)",np.shape(Dfinal_zeff_tot),"np.shape(Nrel_tot_zeff)",np.shape(Nrel_tot_zeff),"len(zmin)",len(zmin)
#print "np.shape(Dfinal_zeff_p)",np.shape(Dfinal_zeff_p),"np.shape(Nrel_zeff_p)",np.shape(Nrel_p_zeff),"len(zmin)",len(zmin)
#print "np.shape(Dfinal_zeff_t)",np.shape(Dfinal_zeff_t),"np.shape(Nrel_zeff_t)",np.shape(Nrel_t_zeff),"len(zmin)",len(zmin)

dlog = 0.1
dlogF = 0.5


Dfinal_zeff_tot_dlog =Dfinal_zeff_tot

Nrel_tot_zeff_dlogD = Nrel_tot_zeff/dlog
Nrel_tot_zeff_dlogF = Nrel_tot_zeff/dlogF

Ffinal_zeff_tot_dlog =Ffinal_zeff_tot

valuesHist_zeff_tot, bins = np.histogram(Dfinal_zeff_tot_dlog, weights=Nrel_tot_zeff_dlogD, normed=False)
valuesHist_zeff_totF, binsF = np.histogram(Ffinal_zeff_tot_dlog, weights=Nrel_tot_zeff_dlogF, normed=False)


mindDdlog = math.ceil(np.log10(np.min(bins[bins>0.]))) -1 #Dfinal_zeff_tot_dlog[Dfinal_zeff_tot_dlog>0.])))
mindFdlog = math.ceil(np.log10(np.min(binsF[binsF>0.])))

maxdDdlog = math.floor(np.log10(np.max(bins))) +1 #Dfinal_zeff_tot_dlog)))
maxdFdlog = math.floor(np.log10(np.max(binsF)))

#mindDdlog = 3.
#maxdDdlog = 6.

#maxdFdlog = -6.
#mindFdlog =  -15.

nbin = (maxdDdlog - mindDdlog)/dlog
nbinF =(maxdFdlog - mindFdlog)/dlogF

#print "nbin", nbin, " minDdlog ", mindDdlog, " maxDdlog ", maxdDdlog
#print "nbinF", nbinF, " minFdlog ", mindFdlog, " maxFdlog ", maxdFdlog

binsValues = np.linspace(mindDdlog, maxdDdlog,nbin+1)
binsValuesF = np.linspace(mindFdlog, maxdFdlog,nbinF+1)

binsValues = 10.**binsValues
binsValuesF = 10.**binsValuesF



###### MARK: figure zeff tot-peak-tail distrib

# zeff - tot
values_zeff_tot = Dfinal_zeff_tot_dlog
weights_zeff_tot = Nrel_tot_zeff_dlogD

shapeVal  = np.shape(values_zeff_tot)
values_zeff_tot = values_zeff_tot.reshape((shapeVal[0]*shapeVal[1]))
weights_zeff_tot  = weights_zeff_tot.reshape((shapeVal[0]*shapeVal[1]))

valuesHist_zeff_tot, bins = np.histogram(values_zeff_tot, weights=weights_zeff_tot, bins=binsValues, normed=False)

# zeff - peaks
Dfinal_zeff_p_dlog =Dfinal_zeff_p#/dlog
values_zeff_p = Dfinal_zeff_p_dlog
weights_zeff_p = Nrel_p_zeff/dlog

shapeVal_p  = np.shape(values_zeff_p)
values_zeff_p = values_zeff_p.reshape((shapeVal_p[0]*shapeVal_p[1]))
weights_zeff_p  = weights_zeff_p.reshape((shapeVal_p[0]*shapeVal_p[1]))

valuesHist_zeff_p, bins = np.histogram(values_zeff_p, weights=weights_zeff_p, bins=binsValues, normed=False)

# zeff - tails
Dfinal_zeff_t_dlog =Dfinal_zeff_t#/dlog
values_zeff_t = Dfinal_zeff_t_dlog
weights_zeff_t = Nrel_t_zeff/dlog

shapeVal_t  = np.shape(values_zeff_t)
values_zeff_t = values_zeff_t.reshape((shapeVal_t[0]*shapeVal_t[1]))
weights_zeff_t  = weights_zeff_t.reshape((shapeVal_t[0]*shapeVal_t[1]))

valuesHist_zeff_t, bins = np.histogram(values_zeff_t, weights=weights_zeff_t, bins=binsValues, normed=False)

centerBin = (bins[1:]+bins[:-1])/2.

sbin = bins[:-1]
ebin = bins[1:]

fig, axes = plt.subplots(1,1)

oldH = 0.0
for ii in range(len(ebin)):
    x = [sbin[ii],ebin[ii]]
    y2 = [valuesHist_zeff_tot[ii],valuesHist_zeff_tot[ii]]
    y1 = [valuesHist_zeff_p[ii],valuesHist_zeff_p[ii]]
    y0 = [0.0,0.0]
    axes.semilogx([sbin[ii],ebin[ii]],[valuesHist_zeff_tot[ii],valuesHist_zeff_tot[ii]],color="black",linewidth = 2.5)
    if ii >0:
        axes.fill_between(x, y0, y1, color = '#4A838B')
        axes.fill_between(x, y1, y2, color = '#82BC92')
        axes.semilogx([sbin[ii],sbin[ii]],[oldH,valuesHist_zeff_tot[ii]],color="darkgray")
        axes.semilogx([ebin[ii],ebin[ii]],[valuesHist_zeff_tot[ii],oldH],color="darkgray")
        axes.semilogx([sbin[ii],sbin[ii]],[valuesHist_zeff_tot[ii-1],valuesHist_zeff_tot[ii]],color="black",linewidth = 2.5)
    else:
        axes.fill_between(x, y0, y1, color = '#4A838B', label = "PEAKS" )
        axes.fill_between(x, y1, y2, color = '#82BC92',  label = "TAILS")
        axes.semilogx([sbin[ii],sbin[ii]],[oldH,valuesHist_zeff_tot[ii]],color="black",linewidth = 2.5)


axes.semilogx([ebin[-1],ebin[-1]],[valuesHist_zeff_tot[-1],0.0],color="black",label='TOTAL',linewidth = 2.5)


plt.legend(fontsize = 15)
plt.ticklabel_format(style='sci', axis='y')
plt.xlabel(r'$\mathbf{Duration~~D~}s[s]$')
plt.ylabel(r'$\mathbf{dN/dlogD}$')
axes.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
axes.yaxis.set_major_formatter(mtick.FormatStrFormatter('%'))
axes.set_yticklabels(axes.get_yticks(), {'weight': 'normal'})
#plt.show()
name_zeff = output_dir+"/SIGNALderivative_DURATION_IN_OBS_zeff_new_EFF_DUR%i_Ntrials%i_flux.pdf" %(EFF_DUR,Ntrials)
fig.savefig(name_zeff ,dpi=100,bbox_inches='tight')

#####FLUX
# zeff - tot
values_zeff_tot = Ffinal_zeff_tot_dlog
weights_zeff_tot = Nrel_tot_zeff_dlogF

shapeVal  = np.shape(values_zeff_tot)
values_zeff_tot = values_zeff_tot.reshape((shapeVal[0]*shapeVal[1]))
weights_zeff_tot  = weights_zeff_tot.reshape((shapeVal[0]*shapeVal[1]))

valuesHist_zeff_tot, binsF = np.histogram(values_zeff_tot, weights=weights_zeff_tot, bins=binsValuesF, normed=False)

# zeff - peaks
Ffinal_zeff_p_dlog =Ffinal_zeff_p#/dlogF
values_zeff_p = Ffinal_zeff_p_dlog
weights_zeff_p = Nrel_p_zeff/dlogF

shapeVal_p  = np.shape(values_zeff_p)
values_zeff_p = values_zeff_p.reshape((shapeVal_p[0]*shapeVal_p[1]))
weights_zeff_p  = weights_zeff_p.reshape((shapeVal_p[0]*shapeVal_p[1]))

valuesHist_zeff_p, binsF = np.histogram(values_zeff_p, weights=weights_zeff_p, bins=binsValuesF, normed=False)

# zeff - tails
Ffinal_zeff_t_dlog =Ffinal_zeff_t#/dlogF
values_zeff_t = Ffinal_zeff_t_dlog
weights_zeff_t = Nrel_t_zeff/dlogF

shapeVal_t  = np.shape(values_zeff_t)
values_zeff_t = values_zeff_t.reshape((shapeVal_t[0]*shapeVal_t[1]))
weights_zeff_t  = weights_zeff_t.reshape((shapeVal_t[0]*shapeVal_t[1]))

valuesHist_zeff_t, binsF = np.histogram(values_zeff_t, weights=weights_zeff_t, bins=binsValuesF, normed=False)

centerBinF = (binsF[1:]+binsF[:-1])/2.

sbin = binsF[:-1]
ebin = binsF[1:]

figF, axesF = plt.subplots(1,1)

oldH = 0.0

for ii in range(len(ebin)):
    x = [sbin[ii],ebin[ii]]
    y2 = [valuesHist_zeff_tot[ii],valuesHist_zeff_tot[ii]]
    y1 = [valuesHist_zeff_p[ii],valuesHist_zeff_p[ii]]
    y0 = [0.0,0.0]
    axesF.semilogx([sbin[ii],ebin[ii]],[valuesHist_zeff_tot[ii],valuesHist_zeff_tot[ii]],color="black",linewidth = 2.5)
    if ii >0:
        axesF.fill_between(x, y0, y1, color = '#4A838B')
        axesF.fill_between(x, y1, y2, color = '#82BC92')
        axesF.semilogx([sbin[ii],sbin[ii]],[oldH,valuesHist_zeff_tot[ii]],color="darkgray")
        axesF.semilogx([ebin[ii],ebin[ii]],[valuesHist_zeff_tot[ii],oldH],color="darkgray")
        axesF.semilogx([sbin[ii],sbin[ii]],[valuesHist_zeff_tot[ii-1],valuesHist_zeff_tot[ii]],color="black",linewidth = 2.5)
    else:
        axesF.fill_between(x, y0, y1, color = '#4A838B', label = "PEAKS" )
        axesF.fill_between(x, y1, y2, color = '#82BC92',  label = "TAILS")
        axesF.semilogx([sbin[ii],sbin[ii]],[oldH,valuesHist_zeff_tot[ii]],color="black",linewidth = 2.5)


axesF.semilogx([ebin[-1],ebin[-1]],[valuesHist_zeff_tot[-1],0.0],color="black",label='TOTAL',linewidth = 2.5)


plt.legend(fontsize = 15)
plt.ticklabel_format(style='sci', axis='y')
plt.xlabel(r'$\mathbf{Flux~~F~~}[erg~cm^{-2}~s^{-1}]$')
plt.ylabel(r'$\mathbf{dN/dlogF}$')
axesF.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
axesF.yaxis.set_major_formatter(mtick.FormatStrFormatter('%'))
axesF.set_yticklabels(axesF.get_yticks(), {'weight': 'normal'})
#plt.show()
name_zeff = output_dir+"/SIGNALderivative_FLUX_IN_OBS_zeff_new_EFF_DUR%i_Ntrials%i_flux.pdf"%(EFF_DUR,Ntrials)
figF.savefig(name_zeff ,dpi=100,bbox_inches='tight')

####### SAVE RESULTS IN FILES


#nameFile_zmin_OUT = output_dir+"/RESULTS_DURATIONS_ZMINtxt"
#with open(nameFile_zmin_OUT, 'w') as f:
#    writer = csv.writer(f, delimiter='\t')
#    writer.writerows(zip(Dfinal_zmin_p_dlog,Nrel_p_zmin,Dfinal_zmin_t_dlog,Nrel_t_zmin))

shape_p = np.shape(Dfinal_zeff_p)
xp = Dfinal_zeff_p.reshape((shape_p[0]*shape_p[1]))
yp = Nrel_p_zeff/dlog
yp = yp.reshape((shape_p[0]*shape_p[1]))

shape_t = np.shape(Dfinal_zeff_t)
xt = Dfinal_zeff_t.reshape((shape_t[0]*shape_t[1]))
yt = Nrel_t_zeff/dlog
yt = yt.reshape((shape_t[0]*shape_t[1]))


nameFile_zeff_OUT = output_dir+"/RESULTS_DURATIONS_ZEFF.txt"
with open(nameFile_zeff_OUT, 'w') as f1:
    writer = csv.writer(f1, delimiter='\t')
    writer.writerows(zip(xp,yp,xt,yt))

cump_zeff = np.cumsum(Np_zeff)
cum_p = np.cumsum(Nrel_p_zeff)

#nameFile_zmax_OUT = output_dir+"/RESULTS_DURATIONS_ZMAX.txt"
#with open(nameFile_zmax_OUT, 'w') as f2:
#    writer = csv.writer(f2, delimiter='\t')
#    writer.writerows(zip(Dfinal_zmax_p_dlog,Nrel_p_zmax,Dfinal_zmax_t_dlog,Nrel_t_zmax))

shape_p = np.shape(Ffinal_zeff_p)
xp = Ffinal_zeff_p.reshape((shape_p[0]*shape_p[1]))
yp = Nrel_p_zeff/dlogF
yp = yp.reshape((shape_p[0]*shape_p[1]))

shape_t = np.shape(Ffinal_zeff_t_dlog)
xt = Ffinal_zeff_t_dlog.reshape((shape_t[0]*shape_t[1]))
yt = Nrel_t_zeff/dlogF
yt = yt.reshape((shape_t[0]*shape_t[1]))


nameFile_FLUXzeff_OUT = output_dir+"/RESULTS_FLUX_ZEFF.txt"
with open(nameFile_FLUXzeff_OUT, 'w') as f3:
    writer = csv.writer(f3, delimiter='\t')
    writer.writerows(zip(xp,yp,xt,yt))
data = np.genfromtxt(nameFile_FLUXzeff_OUT)

print
print "The analysis is done"
print ""
print "Thank you for using saprEMo"

