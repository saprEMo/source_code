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
import astropy.units as u
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib
import latex

#import basic_functions as bf
plt.rcParams['figure.figsize'] = (9.0, 6.0)

# 
InFile=sys.argv[1]
#InFile_check = InFile[-3:]
if (InFile[-4:]!= '.txt'):
    InFile=InFile + '.txt'
dir_OUT=InFile.replace(".txt","")

listName = "listINJ.txt"
CommLINE = "ls INJ_FILES> "+listName
os.system(CommLINE)
invalid_in = True
INJ_LIST = []
with open(listName,'rb') as flist:
    data_list = csv.reader(flist) #read text file with csv
    for row in data_list:
        if ((row[0]==InFile) and ('SetUp' not in row[0])):
            invalid_in = False
        if 'SetUp' not in row[0]:
            INJ_LIST.append(row)
if invalid_in==True:
    print 'Invalid input, choose among:',INJ_LIST
    exit()
COMM = "rm "+listName
os.system(COMM)
InFile= "INJ_FILES/"+InFile


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
zmax =6.0  #1089 #recombination age
zmin = 0.0

NN = 50
nn =500

# ---- Defining variables ----
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

# b. SURVEY
print "# -------- SURVEY -------"
minI = float(rowARR[4])
print '# minimum energy considered: ',minI,'eV'
maxI = float(rowARR[5])
print '# maximum energy considered: ',maxI,' keV'

CoveredSky = float(rowARR[6])
print '# CoveredSky ',CoveredSky*100,' %'
t_AV = float(rowARR[7])
print '# averaged time for observation: ',t_AV,' s'
#SD_T = float(rowARR[len(rowARR)-1])
#print '# standard deviation time for observation: ',SD_T,' s'
print "#"

print "# ---- OTHER INPUTS ---- "
# c. PECENTAGE OF EMITTING SYSTEMS
f_NS = float(rowARR[8])
print "# Percentace of emitting systems: ",f_NS*100,"%"
# d. AVAILABILITY OF DIFFERENT LC as function of energy
LC_z = int(rowARR[9])
if LC_z ==1:
    print "# More light curves in function of energy are available for studing the emission at higher redshift"
elif LC_z == 0:
    print "# A rought estimate it's going to be taken into account since light curves as function of energy have not being provided"
else:
    print
    print "ERROR: invalid entry for light curve usage"
    exit()

Emax_LC = float(rowARR[10])
# e. VOLUME(z)
NameTAB4 = rowARR[11]
#NameTAB4 = NameTAB4.replace('\r','')
NameTAB4 = NameTAB4[:-2]
print "# ComuvingV(z)-table: ",NameTAB4
#CommLINE =  "ls "+NameTAB4
#os.system(CommLINE)
#print "NameTAB4:",repr(NameTAB4)

#### saving data ####
Tab4=np.genfromtxt(NameTAB4)
z4Tab = Tab4[:,0]
dV4Tab = Tab4[:,1]
Vc4Tab = Tab4[:,2]
z4Tab = np.array(z4Tab)
dV4Tab = np.array(dV4Tab)
Vc4Tab = np.array(Vc4Tab)
print "#"

# f. RATE MODEL
print "# ---- RATE MODEL ---- "
RateFileName = rowARR[12]
#RateFileName = RateFileName.replace('\r','')
RateFileName = RateFileName[:-2]
print "# File name for RATE-MODEL: ",RateFileName
listFILErate = "FileRate.txt"
CommLINE =  "ls TEMP/"+dir_OUT+"/RATE*cm3_sm1.txt>> "+listFILErate
print CommLINE
os.system(CommLINE)

#### saving data ####

DR_r1=np.genfromtxt(RateFileName)

Rshape= np.shape(DR_r1)
print "#        Shape Rate data: ",Rshape," # col: ", Rshape[1]," # rows: ",Rshape[0]
if Rshape[1]<2:
    print
    print "ERROR: one column for REDSHIFT and one column for RATE MUST BE PROVIDED"
    exit()
if Rshape[1]>2:
    print "# WARNING: the rate file might contain more columns, for the ANALYSIS the SECOND is considered "
zR_r1 = DR_r1[:,0]
RSum1 = DR_r1[:,1]

RateFileName_low = ""
RateFileName_high = ""
with open(listFILErate,'rb') as flistrate:
    #next(f) #skip heading row in text file (I cannot get csv.Dictreader instead to work as an alternative to this step)
    rate_list = csv.reader(flistrate) #read text file with csv
    count = 0
    for row in rate_list:
        #print row
        if "_lower_" in row[0]:
            RateFileName_low = row[0]
            print "# RateFileName_low:", RateFileName_low
        if "_upper_" in row[0]:
            RateFileName_high = row[0]
            print "# RateFileName_high:", RateFileName_high
    #print
CommLINE =  "rm "+listFILErate
#print CommLINE
os.system(CommLINE)
#print "RateFileName_low!=",RateFileName_low
#print "RateFileName_high!=",RateFileName_high


##########################################################

print "#"
print "# ---- LIGHT CURVE MODEL ---- "
# g. LIGHT CURVEs
LCFile = rowARR[13]
LCFile = LCFile[:-2]
print "# File name for LIGHT CURVE -MODEL: ",LCFile
print "# WARNING: bin equally spaced in energy (write to serena.vinciguerra89@gmail.com to ask for changes)"
#CommLINE =  "ls "+LCFile
#os.system(CommLINE)
LC_data=np.genfromtxt(LCFile)
dim = np.shape(LC_data)
print '#        Shape Light Curve data:',dim," # col: ",dim[1], " # rows: ", dim[0]

T_forLC = LC_data[:,0]

dT_forLC = T_forLC[1:] - T_forLC[:len(T_forLC)-1]
#print dT_forLC
#exit()
# h. ENRGIES for LIGHT CURVEs
LC_E_File = rowARR[14]
LC_E_File = LC_E_File[:-2]
print "# File name for LIGHT CURVE ENERGY: ", LC_E_File
E_LC_data = np.genfromtxt(LC_E_File)
E_LC_bin = np.array(E_LC_data)
DE_LC = (E_LC_bin[1] - E_LC_bin[0])
hDE_LC = 0.5*DE_LC
print '#        Shape Light Curve Energies: ',np.shape(E_LC_data)

Diff_Ebin = E_LC_bin[1:] - E_LC_bin[:(len(E_LC_bin)-1)]
Diff_Ebin  = np.array(Diff_Ebin)/DE_LC
oneARR = np.ones(len(Diff_Ebin))
if (Diff_Ebin.all == oneARR.all):
    print " ERROR: bin in energy for light curve not equally spaced"
    exit()

if LC_z == 1:
    Emax_LC = E_LC_bin[len(E_LC_bin)-1]+hDE_LC


LC_EA = []
LC_tot = np.zeros(len(T_forLC))
#print "dim[1]",dim[1]
#print "E_LC_bin",np.shape(E_LC_bin)


for a in range(1,dim[1]):
    LC_Ei=LC_data[:,a]
    LC_Ei =np.array(LC_Ei)
    LC_EA.append(LC_Ei)
    #print a-1, "len(E_LC_bin)",len(E_LC_bin)
    if (E_LC_bin[a-1] - hDE_LC <=maxI and E_LC_bin[a-1] - hDE_LC >minI) or (E_LC_bin[a-1] + hDE_LC >=minI and E_LC_bin[a-1] + hDE_LC <maxI):
        #if a <10:
        #print 'a-1',a-1
        minELCi = E_LC_bin[a-1]-hDE_LC
        maxELCi = E_LC_bin[a-1]+hDE_LC
        minEi = minI
        maxEi =maxI
        if minELCi>minEi:
            minEi = minELCi
        if maxELCi< maxEi:
            maxEi = maxELCi
        DE_Ei = maxEi - minEi
        LC_tot = LC_tot+ LC_Ei*DE_Ei
    elif E_LC_bin[a-1] - hDE_LC <=maxI:
        break
    else:
        continue

LC_tot = LC_tot/DE_LC

LC_max = LC_tot

#print 'np.max(LC_max)',np.max(LC_max)

LC_EA = np.array(LC_EA)
print '#         Shape Light Curve Array:', np.shape(LC_EA)


if dim[1]-1-len(E_LC_bin)!=0:
    print "ERROR: number of likelihood available not consistent with number of energies from file"
    exit()

print "#"
print "# ---- SURVEY SENSITIVITY ---- "

# i. SESITIVITY
SENS_FileName = rowARR[15]
SENS_FileName =SENS_FileName[:-2]
print "# File name for SURVEY SENSITIVITY: ", SENS_FileName
SENS_data =np.genfromtxt(SENS_FileName)
print "# SENS_data",SENS_data

shapeSENS=np.shape(SENS_data)

Ebin_AvD = np.zeros(len(shapeSENS))
Ebin_iD= np.zeros(len(shapeSENS))
Ebin_sD = np.zeros(len(shapeSENS))
SensD = np.zeros(len(shapeSENS))

dim0= 0
if len(shapeSENS) ==1:
    Ebin_AvD = SENS_data[0]
    Ebin_iD= SENS_data[2]
    Ebin_sD = SENS_data[3]
    SensD = SENS_data[1]
    dim0=1
elif len(shapeSENS) >1:
    Ebin_AvD = SENS_data[:,0]
    Ebin_iD= SENS_data[:,2]
    Ebin_sD = SENS_data[:,3]
    SensD = SENS_data[:,1]
    dim= np.shape(Ebin_AvD)
    dim0= dim[0]
else:
    print "ERROR: empty sensitivity file"
    exit()

print "# Ebin_AvD",Ebin_AvD
print "# Ebin_iD",Ebin_iD
print "# Ebin_sD",Ebin_sD
print "# SensD",SensD


print "#"
print "# ---- ABSORPTION ---- "

# L. GALACTIC ABSORPTION
ABS_FileName = rowARR[16]
ABS_FileName=ABS_FileName[:-2]
Abs = np.zeros(10)
#print repr(ABS_FileName)
if ABS_FileName=='""':
    print "# No galactic absorption has been selected"
    Abs= np.ones(dim0)
else:
    print "# File name for ABSORPTION: ", ABS_FileName
    ABS_data = np.genfromtxt(ABS_FileName)
    
    En_abs_av = np.zeros(len(shapeSENS))
    Abs= np.zeros(len(shapeSENS))
    
    if len(shapeSENS) ==1:
        En_abs_av = ABS_data[0]
        Abs= ABS_data[1]
    elif len(shapeSENS) >1:
        En_abs_av = ABS_data[:,0]
        Abs= ABS_data[:,1]
    else:
        print "ERROR: empty sensitivity file"
        exit()
    

    #if (np.array(En_abs_av)!=np.array(Ebin_AvD)):
    if (not (np.array(En_abs_av)==np.array(Ebin_AvD)).all()):
        print "ERROR: energies used for defining absorption are different from the ones used to define the sensitivity"
        print "Energies from sensitivities: ",Ebin_AvD
        print "Energies from absorption: ", En_abs_av

    #print "# En_abs_av ",En_abs_av

print "# Trasmission coefficients:",Abs
print "#"
print "# ---- EXPOSURE TIME DISTRIBUTION ----"
DIST_EXPT = rowARR[17]
DIST_EXPT=DIST_EXPT[:-2]
SD_EXPT = float(rowARR[18])
print "# Represented by a ",DIST_EXPT," distribution with standard deviation ",SD_EXPT
print "#"
print "######################################################"
print
# MARK: END DEFINITION

#[Hc, OmegaM, OmegaL,OmegaK,Emin, Emax, CoveredSky, AvObsT, percentage_EV, LC_z, LimE, TABname, RateFileName, LCFile, EN_LC_File, SENS_FileName,ABS_FileName, Dist_OBSD, SD_T]


#INPUTS= np.genfromtxt(InFile)
#print INPUTS

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

def fromDLtoZ0(DLf,OM, Ok, OV,Err):  # returns Distance in cm
    
    z_max = zmax
    z_min = zmin
    z_m = (zmax+zmin)*0.5
    zr_i =z_min
    zr_s =z_max
    Ntrial= 1000
    DLmaxF = fromZtoDL(zmax,OM, Ok, OV)
    if (DLmaxF<DLf):
        print 'z correspondent to DL is bigger than zmax',zmax
        return zmax
    for j in range(Ntrial):
        DLF_i = fromZtoDL(z_m,OM, Ok, OV)
        Diff = DLf-DLF_i
        if Diff<0.:
            #print 'min'
            zr_s = z_m
            z_m = (zr_i + zr_s)*0.5
        else:
            #print 'mag'
            zr_i = z_m
            z_m = (zr_i+ zr_s)*0.5
        if (np.abs(Diff)<=Err):
            continue
    print 'Error',Diff
    return z_m

def fromDLtoZ(DLf,OM, Ok, OV): # from distance in [cm] to z
    z_max = zmax
    z_min = zmin
    z_m = (zmax+zmin)*0.5
    zr_i =z_min
    zr_s =z_max
    Ntrial= 1000
    ErrTh = 0.1
    DLF_i= fromZtoDL(z_m,OM, Ok, OV)
    DLmaxF = fromZtoDL(zmax,OM, Ok, OV)
    if (DLmaxF<DLf):
        #print 'DL is larger than the one correspondent to zmax'
        return zmax
    ErrRel = np.abs(DLF_i - DLf)/DLf
    count = 0
    Ntrial = 1000
    while (count<Ntrial or ErrRel>ErrTh):
        DLF_i = fromZtoDL(z_m,OM, Ok, OV)
        Diff = DLf-DLF_i
        if Diff<0.:
            zr_s = z_m
            z_m = (zr_i + zr_s)*0.5
        else:
            zr_i = z_m
            z_m = (zr_i+ zr_s)*0.5
        ErrRel = np.abs(Diff)/DLf
        count = count+1
    if ErrRel>0.0001:
        print 'Relative error: ',ErrRel
    return z_m


def fromDLtoZErr(DLf,OM, Ok, OV,Err):
    z_max = zmax
    z_min = zmin
    z_m = (zmax+zmin)*0.5
    zr_i =z_min
    zr_s =z_max
    Ntrial= 1000
    DLmaxF = fromZtoDL(zmax,OM, Ok, OV)
    if (DLmaxF<DLf):
        return zmax
    for j in range(Ntrial):
        DLF_i = fromZtoDL(z_m,OM, Ok, OV)
        Diff = DLf-DLF_i
        if Diff<0.:
            zr_s = z_m
            z_m = (zr_i + zr_s)*0.5
        else:
            zr_i = z_m
            z_m = (zr_i+ zr_s)*0.5
        if (np.abs(Diff)<=Err):
            continue
    return Diff



def OBS_array_func(Flimf, CoveredSkyf, t_AVf, SD_Tf,minIf, maxIf):
    OBSERV = [Flimf, CoveredSkyf, t_AVf, SD_Tf,minIf,maxIf]
    return OBSERV


def returnIND(zTab, z_try):
    zTAR = np.min(np.abs(zTab - z_try))
    zIND = np.where(np.abs(zTab - z_try)==zTAR)[0][0]
    return zIND

def z_effF(z1,z2,nf):
    zslice = (z2-z1)/nf
    DVtot = 0
    z_eff0 = 0
    DVf = ComuvingVol(z2,OmegaM, Omegak, OmegaV)-ComuvingVol(z1,OmegaM, Omegak, OmegaV)
    for j in range(1,nf):
        zjh = z1 + zslice*(j-0.5)
        zj = z1 + zslice*j
        zjm1 = z1 + zslice*(j-1)
        dVc = ComuvingVol(zj,OmegaM, Omegak, OmegaV)-ComuvingVol(zjm1,OmegaM, Omegak, OmegaV)
        z_eff0=z_eff0+zjh*dVc
        DVtot = DVtot + dVc
    
    z_eff0= z_eff0/DVtot
    return z_eff0

def CV_T(z1,z2,nf, zTab, VcTab):
    dzf = (z2-z1)/nf
    #print 'dzf',dzf
    zim1f = z1
    zp1dVf = 0.0
    if zim1f == 0.0:
        vim1f = 0.0
    else:
        ind_if = returnIND(zTab,zim1f)
        vim1f = VcTab[ind_if]
    
    for xx in range(nf):
            zif =(xx+1)*dzf + z1
            ind_if = returnIND(zTab,zif)
            vif = VcTab[ind_if]
            dvif = vif- vim1f
            zp1dVf = zp1dVf+dvif*(1./(1.+zif) +1./(1.+zim1f))*0.5
            zim1f = zif
            vim1f = vif
    
    return zp1dVf

######################### LIMIT IN REDSHIFT ###############################


DIM_LCE = np.shape(E_LC_data)
DE_LC = 0.
DE_LC_h = 0.

if (DIM_LCE<2):
    print "WARNING: since only one light curve has been provided it will be considered to cover ALL the energy band of the instrument until ", Emax_LC, " keV are out of the instrument range"
    DE_LC = maxI - minI
    DE_LC_h = DE_LC*0.5
    LC_max = LC_EA[0]
else:
    DE_LC =E_LC_data[1] - E_LC_data[0]
    DE_LC_h = DE_LC*0.5

#Rough estimate of the total limit in flux [erg/(s cm^2)]
sumAv = 0.0
SensD = np.atleast_1d(SensD)
for a in range(len(SensD)):
    sumAv = sumAv+SensD[a]

F_lim =sumAv

#Rough estimate of the maximum z achievable
# there are 3 majour limits:
#       1) from maximum LC (just a rough estimate since usure at which Energy it's)
#       2) model of the RATE-> unphysical but real
#       3) from shift of the maximum peak or end of light curves
# NB: 1) and 3) are related to the LC -> if LC_z >0 both the estimation are not necessary


def maxZ(Emax_LCf,LCf,zmaxf,zR_r1f, minIf):
    # 1)
    maxDf = np.sqrt(np.amax(LCf)/(4.*np.pi*np.amin(SensD)))#to be more conservative-> the selection is anyway done in the calculation(4.*np.pi*F_lim))
    zmaxLf = fromDLtoZ(maxDf,OmegaM,Omegak,OmegaV)
    z_maxf = zmaxLf
    # 2)
    zmaxRATEf = zR_r1f[len(zR_r1f)-1] #+ (zR_r1[len(zR_r1)-2] - zR_r1[len(zR_r1)-1])*0.5
    z_maxf = zmaxRATEf if z_maxf > zmaxRATEf else z_maxf
    # 3)
    zmax_specLCf = Emax_LCf/minIf - 1.
    z_maxf = zmax_specLCf if z_maxf > zmax_specLCf else z_maxf
    zall = [zmaxLf,zmaxRATEf,zmax_specLCf]
    out = [z_maxf,zall]
    
    return out

def MergeBins(x0edge,y0,x1edge):
    y1 =np.zeros(len(x1edge)-1)
    if len(x0edge) == len(y0)+1:
        x0_s = x0edge[1:]
        x0_i = x0edge[:len(x0edge)-1]
        x0_av = np.array(x0_s + x0_i)*0.5
    elif(len(x0edge)==len(y0)):
        indexsort = np.argsort(x0edge)
        x0_av = np.sort(x0edge)
        y0 = y0[indexsort]
    
    else:
        print "ERROR: inconsistent dimensions of inputs"
        exit()
    x1_s = x1edge[1:]
    x1_i = x1edge[:len(x1edge)-1]
    

    if len(x1edge)>len(x0edge):
        print "ERROR: cannot merge in a larger number of bins"
        return y0
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
    sum_to_TOT = sum_to_TOT_l + sum_to_TOT_h
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


zmaxOUT = maxZ(Emax_LC,LC_max,zmax,zR_r1,minI)
zmax = zmaxOUT[0]
zmaxARR = zmaxOUT[1]
zmaxL = zmaxARR[0]
zmaxRATE  = zmaxARR[1]
zmax_specLC = zmaxARR[2]


###############################################

###### CODE CORE #######

def ErrFUN(Np,pre_fac,tav,sd,tLC,Reff_if,zeffA,zRf,Rf,OM,Ok,OV,zmax):

    ########### Numerical error REDSHIFT
    dz = 0.001
    f_old = Rf[0]/(1.+zRf[0])*(ComuvingVol(zRf[0]+dz,OM, Ok, OV)-ComuvingVol(np.max([zRf[0]-dz,0.001]),OM, Ok, OV))/(2*dz)
    z_old = zRf[0]
    
    f_i = Rf[1]/(1.+zRf[1])*(ComuvingVol(zRf[1]+dz,OM, Ok, OV)-ComuvingVol(np.max([zRf[1]-dz,0.001]),OM, Ok, OV))/(2*dz)
    
    fp_i = (f_i-f_old)/(zRf[1] - z_old)
    zp_i = (z_old+zRf[1])*0.5
    
    Err = 0.0
    for i in range(1,len(zRf)-1):
        
        f_ii = Rf[i+1]/(1.+zRf[i+1])*(ComuvingVol(zRf[i+1]+dz,OM, Ok, OV)-ComuvingVol(np.max([zRf[i+1]-dz,0.001]),OM, Ok, OV))/(2*dz)
        
        fp_ii = (f_ii-f_i)/(zRf[i+1] - zRf[i])
        zp_ii = (zRf[i]+zRf[i+1])*0.5
        
        fs_i = (fp_ii - fp_i)/(zp_ii - zp_i)
        zs_i = zRf[i]

        Err= Err+np.abs(fs_i)/12.*(zp_ii - z_old)
    
        z_old = zp_ii
        zp_i = zp_ii
        fp_i = fp_ii
    
        if zRf[i]>zmax:
            break
    
    Err = Err*tav*pre_fac
    ##############################

    Err_Tobs = Np/tav*sd
    
    #############################

    ErrP = np.sqrt(Err**2+Err_Tobs**2)
    #print "PEAK ERRORS ",ErrP,"Err num z",Err,"Err_Tobs",Err_Tobs
    #############################
    #############################
    
    dt_max = np.max(tLC[1:] - tLC[:len(tLC)-1])
    
    # both contribution are over estimates
    #Err_TLC = Np/tav*dt_max*(np.matrix((1+zeffA))*np.transpose(np.matrix(np.ones(len(zeffA)))))
    Err_TLC = 0.0

    for i in range(len(Reff_if)):
        Err_TLC = Err_TLC+(pre_fac*Reff_if[i]*dt_max*(1+zeffA[i]))**2

    Err_TLC = np.sqrt(Err_TLC)

    ErrT = np.sqrt(Err**2+Err_TLC**2)
    #print "TAILS ERRORS ",ErrT,"Err num z",Err,"Err_TLC",Err_TLC
    
    ErrTOT = np.sqrt(ErrT**2+ErrP**2)
    
    outErr=[ErrTOT,ErrT,ErrP]
    return outErr




###### 1. CORE FUNCTIONS #######

def PEAKS_ALL(fImin,fImax,SensDf,Ebin_iDf,Ebin_sDf,Absf,tAf,E_LC_binf, L_LCAf,z_Tab,Vc_Tab,zDAf,RsDVAf,OM,Ok,OV,zMM):
    #print "ASSUMPTION: detection is made if at least in one of the sensitivity band the FLUX > FLUX_lim"
    
    #from smaller to greater values of z to have positive DV
    
    s_indf = np.argsort(zDAf)
    zDAf = np.sort(zDAf)
    RsDVAf = RsDVAf[s_indf]
        
    zM = zMM
    
    #zmin_abs = RsDVAf[0] if z_Tab[0]<RsDVAf[0] else z_Tab[0]
    

        
    NbinEI = len(Ebin_iDf)
        
    #--ordering energetic bins in the instrument band
    ind0 = np.argsort(Ebin_iDf)
    Ebin_iDf = np.sort(Ebin_iDf)
    Ebin_sDf = np.sort(Ebin_sDf)
    SensDf = np.array(SensDf)
    SensDf = SensDf[ind0]
        
    #--ordeing energetic bins in the LCs
    indLC = np.argsort(E_LC_binf)
    E_LC_binf = np.sort(E_LC_binf)
    L_LCAf = np.array(L_LCAf)
    L_LCAf = L_LCAf[indLC,:]
    F_zeff = np.zeros([len(zDAf),len(L_LCAf[0,:])])
        
    DeltaE_LC = E_LC_binf[1] - E_LC_binf[0]
    hDeltaE_LC = DeltaE_LC*0.5
    
    
    R_s = 0.0
    DVtot0 = 0.0
    Vcum =[]
    DV_A = []
    DV_zw =[]
    Rei=[]
    Rec=[]
    ziA =[]
    zim1A=[]
    zeffA = []
    Flux_zeffA = []
    Flux_zmaxA = []
    Flux_zminA = []
    

    # distinguishing the cases of more bind in z - rate are used or only one
    #zMbin1 = zDAf[0] if zM > zDAf[0] else zM
    
    #zmin_abs: the biggest between first z in rate and in table
    zmin_abs = zDAf[0] if z_Tab[0]<zDAf[0] else z_Tab[0]
    zMbin1 = zmin_abs if zM > zmin_abs else zM

    z_old = zMbin1
    ind_j = returnIND(z_Tab,zMbin1)
    VC_old = Vc_Tab[ind_j] # Volume in the table correspondent to zDAf[0]
            
    # Calculating dV/(z+1) from 0 to first step -> NB: NON HO GIA UNA FUNZIONE CHE LO FA?
        
        
    zp1dV0 = CV_T(0.0, z_old, NN, z_Tab, Vc_Tab)
    z0_eff = z_effF(0.,zMbin1,NN) # most probable z, considering sources prop to volume
    
    iLC=0
    iLCzmin = 0
    iLCzmax = 0
    z_seq =[]
    
    indLIST_zeff = []
    indLIST_zmin = []
    indLIST_zmax = []
    
    D_A_zeff = []
    D_A_zmin = []
    D_A_zmax = []
    
    D_A_zeff_BP =[]
    D_A_zmin_BP =[]
    D_A_zmax_BP =[]
    
    ####### ADDED 18 SEPT 2017 FOR DURATIONS #######
    
    Ti_zmin = []
    Tf_zmin = []
    
    Ti_zeff = []
    Tf_zeff = []
    
    Ti_zmax = []
    Tf_zmax = []
    
    timeLpArr = []
    #################################################
    
    
    Dt_A = np.array(tAf[1:]) - np.array(tAf[:len(tAf)-1])

    index_tot = np.arange(len(tAf))
    """loop over the instrument Energy bin"""
    L_LC_totf = np.zeros(len(tAf)-1)
    LC_totf = np.zeros(len(tAf))
    count = 0
    ##### MARK: ENERGY LOOP #####
    for k in range(NbinEI):
                DEcov = 0.0
                
                Av_LC_tot = np.zeros(len(tAf)-1)
                LC_tot = np.zeros(len(tAf))
                
                #THINK ABOUT IT!
                zj = z_old
                
                
                Emin_k = Ebin_iDf[k]*(1.+z0_eff) if LC_z==1 else Ebin_iDf[k]
                Emax_k = Ebin_sDf[k]*(1.+z0_eff) if LC_z==1 else Ebin_sDf[k]
                ############ USEFUL FOR ERROR CALCULATION: NOT USED AT THE MOMENT
                #
                #Emin_k_zmin= Ebin_iDf[k]
                #Emax_k_zmin = Ebin_sDf[k]
                #Emin_k_zmax = Ebin_iDf[k]*(1.+zMbin1) if LC_z==1 else  Ebin_iDf[k]
                #Emax_k_zmax = Ebin_sDf[k]*(1.+zMbin1) if LC_z==1 else Ebin_sDf[k]
                #
                #################################################################

                countW0 = 0
                
                while (E_LC_binf[iLC]+hDeltaE_LC<=Emin_k and iLC< (len(E_LC_binf)-1)):
                    iLC = iLC+1
                    #print 'iLC',iLC,'(len(E_LC_binf)',(len(E_LC_binf))
                    countW0 = countW0 +1
                    if countW0==len(E_LC_binf)-1:
                        print'coutW0==len(E_LC_binf)-1:', countW0
                        exit()

                # NO ERROR CALCULATION # -- 24/03/2017
                if (E_LC_binf[iLC]-hDeltaE_LC>Emax_k):
                    continue
                LCi_Emax = E_LC_binf[iLC]+hDeltaE_LC
                LCi_Emin = E_LC_binf[iLC]-hDeltaE_LC
                
                while (iLC< (len(E_LC_binf))and LCi_Emin< Emax_k and LCi_Emax > Emin_k):
                    #print 'iLC',iLC
                    LCi_Emax = E_LC_binf[iLC]+hDeltaE_LC
                    LCi_Emin = E_LC_binf[iLC]-hDeltaE_LC
                
                    MINe = LCi_Emin if LCi_Emin > Emin_k else Emin_k
                    MAXe = LCi_Emax if LCi_Emax < Emax_k else Emax_k
                
                    vL = np.array(L_LCAf[iLC,:])
                    
                    Av_LC=vL[1:] + vL[:len(vL)-1]
                    Av_LC = np.array(Av_LC)
                    Av_LC = Av_LC*0.5
                    W = (MAXe - MINe)/DeltaE_LC
                    #print "len(Av_LC_tot)",len(Av_LC_tot),"len(Av_LC)",len(Av_LC)
                    Av_LC_tot = Av_LC_tot+ Av_LC* W
                    
                    LC_tot = LC_tot+ vL* W
                    
                    DEcov = DEcov + (MAXe - MINe)
                    
                    # condition to stop the "while"
                    #if (LCi_Emax>= Emax_k or LCi_Emax-Emax_k<0.000001):
                    if (LCi_Emax>= Emax_k):
                        #print "dentro if: break, count ",count,"LCi_Emax",LCi_Emax,"Emax_k",Emax_k
                        break
                    
                    iLC = iLC+1
                    #print 'DEcov', DEcov


                if(DEcov-(Emax_k - Emin_k)>0.0000000001):
                    print 'CODE BUG: the energy band considered for the light curve is wider than the one energey band of the instrument!'
                    print 'DEcov>(Emax_k - Emin_k)'
                    print 'DEcov',DEcov
                    print 'Emax_k - Emin_k',Emax_k - Emin_k
                    exit()
                    
                Av_LC_tot = Av_LC_tot*Absf[k]
                #print Absf
                
                LC_tot = LC_tot*Absf[k]
                
                #--- DEFINING LIGHT CURVES AT z_eff, z_min, z_max ---#
                #--- zeff ---#

                dist_z0eff = fromZtoDL(z0_eff, OM, Ok, OV)
                
                
                Lmax_k = 4.*np.pi*SensDf[k]*dist_z0eff*dist_z0eff
                # saving LC(z_min) which at zmax have a flux above the sensitivity threshold
                LthAf_zEj = Av_LC_tot[Av_LC_tot>Lmax_k]
                indLIST_zeff_i = index_tot[LC_tot>Lmax_k] #index_tot[LC_totf>Lmax_k]
                indLIST_zeff.extend(indLIST_zeff_i)

                #--- zmin ---#
                dist_zmin = fromZtoDL(0.0, OM, Ok, OV)
                Lmax_k_zmin = 4.*np.pi*SensDf[k]*dist_zmin*dist_zmin
                # saving LC(z_eff) which at zmax have a flux above the sensitivity threshold
                LthAf_zmin = Av_LC_tot[Av_LC_tot>Lmax_k_zmin]
                indLIST_zmin_i = index_tot[LC_tot>Lmax_k_zmin] #index_tot[LC_totf>Lmax_k_zmin]
                #print "len(indLIST_zmin_i)",len(indLIST_zmin_i)
                indLIST_zmin.extend(indLIST_zmin_i)

                #--- zmax ---#
                dist_zmax = fromZtoDL(zMbin1, OM, Ok, OV)
                Lmax_k_zmax = 4.*np.pi*SensDf[k]*dist_zmin*dist_zmax
                # saving LC(z_max) which at zmax have a flux above the sensitivity threshold
                LthAf_zmax = Av_LC_tot[Av_LC_tot>Lmax_k_zmax]
                indLIST_zmax_i = index_tot[LC_tot>Lmax_k_zmax] #index_tot[LC_totf>Lmax_k_zmax]
                indLIST_zmax.extend(indLIST_zmax_i)
                
                list_AVLC_k = list(Av_LC_tot)
                
                ##### TEST 1 18/09/2017
                #print "len(Lmax_k_zmax)",np.shape(LthAf_zmax),"len(Av_LC_tot)",np.shape(Av_LC_tot)
                #print "len(LthAf_zmin)",np.shape(LthAf_zmin),"len(Av_LC_tot)",np.shape(Av_LC_tot)
                #exit()
                # RES: LENGTHS ARE EQUALS
                #######

                L_LC_totf = L_LC_totf + Av_LC_tot
                LC_totf = LC_totf +LC_tot

                count = count +1
                #Calculating maximum distance reachable
                DL_max_jk = np.sqrt(np.max(Av_LC_tot)/(4.*np.pi*SensDf[k]))
                z_max_kj=fromDLtoZ(DL_max_jk,OM,Ok,OV)
                #print 'z_max_kj',z_max_kj
                z_seq.append(z_max_kj)
                if iLC == len(E_LC_binf):
                    #print "at the end of LIGHT CURVE BIN before the completing the spectrum range of the instrument"
                    break
                        
    
    "---- end of loop over energy bands"
    
    
    z_seq = np.array(z_seq)
            
    if(np.max(z_seq)<zMbin1):
            print "the maximum redshift at which the light curve is visible is less than the first z saved for rate -> no accuarte estimate"
            return(-1)

    maxLC = np.max(LC_totf)
    LC_totfL = list(LC_totf)
    indexPEAK = LC_totfL.index(maxLC)
    F_tot = LC_totf/(4.*np.pi*dist_z0eff*dist_z0eff)
    F_zeff[0,:] = F_tot

    R_ot_i = RsDVAf[0]*zp1dV0
    j =1
    Vcum.append(VC_old)
    DV_A.append(VC_old)
    DV_zw.append(zp1dV0)
    Rei.append(R_ot_i)
    Rec.append(R_ot_i)
    ziA.append(zMbin1)
    zim1A.append(0.0)


    Dz0_eff = dist_z0eff #fromZtoDL(z0_eff,OM, Ok, OV)
    Dz0_min= dist_zmin #fromZtoDL(0.0,OM, Ok, OV)
    Dz0_max = dist_zmax #fromZtoDL(zMbin1,OM, Ok, OV)
    
    # to avoid inf, we set z =0.0 at a distance of 0.5Mpc which roughly correspond to z = 0.0001
    Dz0_min = 0.5*Mpc
    FL_0 = np.amax(L_LC_totf)/(4*np.pi*Dz0_eff*Dz0_eff)
    FL_0min = np.amax(L_LC_totf)/(4*np.pi*Dz0_min*Dz0_min)
    FL_0max = np.amax(L_LC_totf)/(4*np.pi*Dz0_max*Dz0_max)
    
    
    zeffA.append(z0_eff)
    
    Flux_zeffA.append(FL_0)
    Flux_zminA.append(FL_0min)
    Flux_zmaxA.append(FL_0max)
    
    
    ########## ADDED 18 SEPT FOR DURATION PLOTS ########
    
    if not indLIST_zmin:
        Ti_zmin.append(indLIST_zmin)
        Tf_zmin.append(indLIST_zmin)
    else:
        Ti_zmin.append(tAf[np.min(np.atleast_1d(indLIST_zmin))])
        Tf_zmin.append(tAf[np.max(indLIST_zmin)])
    
    if not indLIST_zeff:
        Ti_zeff.append(indLIST_zeff)
        Tf_zeff.append(indLIST_zeff)
    else:
        Ti_zeff.append(tAf[np.min(indLIST_zeff)])
        Tf_zeff.append(tAf[np.max(indLIST_zeff)])

    if not indLIST_zmax:
        Ti_zmax.append(indLIST_zmax)
        Tf_zmax.append(indLIST_zmax)
    else:
        Ti_zmax.append(tAf[np.min(indLIST_zmax)])
        Tf_zmax.append(tAf[np.max(indLIST_zmax)])
    
    timeLpArr.append(tAf[indexPEAK])
    
    ##################


    del indLIST_zeff[:]
    del indLIST_zmin[:]
    del indLIST_zmax[:]
    del list_AVLC_k[:]
    
    """loop over redshifts"""
    while (zDAf[j] - zM<0.0000000000001 and j<len(zDAf)):
        
                indLIST_zeff = []
                indLIST_zmin = []
                indLIST_zmax = []
        
                L_LC_totf = np.zeros(len(tAf)-1)
                LC_totf = np.zeros(len(tAf))
                ind_zj = returnIND(z_Tab,zDAf[j])
                VC_j = Vc_Tab[ind_zj]
                DV = VC_j - VC_old
                
                zEj = z_effF(z_old,zDAf[j],NN)
                
                #########################
                # maximum z of the bin to be considered for the volume -> since there is a degree of recursion, w define a new variable whose value can be changed accordingly to the TEST value - see later -
                zBINmax = zDAf[j]
                #########################
                
                iLC =0
                z_seq =[]

                '''loop over energy bin of the instrument'''
                count = 0
                ##### MARK: ENERGY LOOP #####
                for k in range(NbinEI):
                    #print 'k',k
                    DEcov = 0.0
        
                    Av_LC_tot = np.zeros(len(tAf)-1)
                    
                    ######################################################################################################################
                    # NB: EVEN HERE THERE IS AN ERROR DUE TO THE REDSHIFT BIN - INFLUENCE ON THE ENERGY BAND CONSIDERED
                    # NB: despite the broadening if the band the spectral bins remain disjoint - if I consider a single redshift
                    #
                    #if instead I use zBINmax and zBINmin this is no more true - however It wont be correct to consider a light curve
                    # given by the sum of the light curve which includes such a large z bin! I would sum more spectral band than the true
                    # what I should do is to weight the light curve calculated at each z from z_old to zDA[j] with the volume
                    # I think considering zEj is a good approximation
                    #########################################################################################################################
                    
                    Emin_k = Ebin_iDf[k]*(1.+zEj) if LC_z==1 else Ebin_iDf[k]
                    Emax_k = Ebin_sDf[k]*(1.+zEj) if LC_z==1 else Ebin_sDf[k]
                    
                    ############ USEFUL FOR ERROR CALCULATION: NOT USED AT THE MOMENT
                    #
                    #Emin_k_zmin= Ebin_iDf[k]*(1.+z_old) if LC_z==1 else Ebin_iDf[k]
                    #Emax_k_zmin = Ebin_sDf[k]*(1.+z_old) if LC_z==1 else Ebin_sDf[k]
                    #Emin_k_zmax = Ebin_iDf[k]*(1.+zDAf[j]) if LC_z==1 else  Ebin_iDf[k]
                    #Emax_k_zmax = Ebin_sDf[k]*(1.+zDAf[j]) if LC_z==1 else Ebin_sDf[k]
                    #
                    ###################################################################
                
                
                    countW0 = 0
                    while (E_LC_binf[iLC]+hDeltaE_LC<=Emin_k and iLC< (len(E_LC_binf)-1)):
                        iLC = iLC+1
                        countW0 = countW0 +1
                        if countW0> 200:
                            print'coutW0>200', countW0
                            exit()
                    
                    
                    if (E_LC_binf[iLC]-hDeltaE_LC>Emax_k):
                        continue
        
                    LCi_Emax = E_LC_binf[iLC]+hDeltaE_LC
                    LCi_Emin = E_LC_binf[iLC]-hDeltaE_LC
                    
                    while (iLC< (len(E_LC_binf))and LCi_Emin< Emax_k and LCi_Emax > Emin_k):

                        LCi_Emax = E_LC_binf[iLC]+hDeltaE_LC
                        LCi_Emin = E_LC_binf[iLC]-hDeltaE_LC
                    
                        MINe = LCi_Emin if LCi_Emin > Emin_k else Emin_k
                        MAXe = LCi_Emax if LCi_Emax < Emax_k else Emax_k
                        
                        vL = np.array(L_LCAf[iLC,:])
                        
                        Av_LC=vL[1:] + vL[:len(vL)-1]
                        Av_LC = np.array(Av_LC)
                        Av_LC = Av_LC*0.5
                        
                        W = (MAXe - MINe)/DeltaE_LC
                        
                        Av_LC_tot = Av_LC_tot+ Av_LC* W
                        LC_tot = LC_tot+ vL* W
                        
                        DEcov = DEcov + (MAXe - MINe)
                        
                        
                        if (LCi_Emax>= Emax_k):
                            break
                    
                        iLC = iLC+1

                    #------ while closed -----
                
                    if(DEcov - (Emax_k - Emin_k)>0.00000000001):
                        print 'CODE BUG: the energy band considered for the light curve is wider than the one energey band of the instrument!'
                        print 'DEcov>(Emax_k - Emin_k)'
                        print 'DEcov',DEcov
                        print 'Emax_k - Emin_k',Emax_k - Emin_k
                        exit()
                    
                    ##################################

                    Av_LC_tot = Av_LC_tot*Absf[k]
                    LC_tot = LC_tot*Absf[k]

                    L_LC_totf = L_LC_totf + Av_LC_tot
                    LC_totf = LC_totf + LC_tot
                    
                    dist_zeff = fromZtoDL(zEj, OM, Ok, OV)
                    Lmax_k = 4.*np.pi*SensDf[k]*dist_zeff*dist_zeff
                    LthAf_zEj = Av_LC_tot[Av_LC_tot>Lmax_k]
                    indLIST_zeff_i = index_tot[LC_tot>Lmax_k] #index_tot[LC_totf>Lmax_k]
                    indLIST_zeff.extend(indLIST_zeff_i)
        
                    dist_zmin = fromZtoDL(z_old, OM, Ok, OV)
                    Lmax_k_zmin = 4.*np.pi*SensDf[k]*dist_zmin*dist_zmin
                    LthAf_zmin = Av_LC_tot[Av_LC_tot>Lmax_k_zmin]
                    indLIST_zmin_i = index_tot[LC_tot>Lmax_k_zmin] #index_tot[LC_totf>Lmax_k_zmin]
                    indLIST_zmin.extend(indLIST_zmin_i)
                
                    dist_zmax = fromZtoDL(zDAf[j], OM, Ok, OV)
                    Lmax_k_zmax = 4.*np.pi*SensDf[k]*dist_zmin*dist_zmax
                    LthAf_zmax = Av_LC_tot[Av_LC_tot>Lmax_k_zmax]
                    indLIST_zmax_i = index_tot[LC_tot>Lmax_k_zmax] #index_tot[LC_totf>Lmax_k_zmax]
                    indLIST_zmax.extend(indLIST_zmax_i)

                    list_AVLC_k = list(Av_LC_tot)

                    ##############################

                    
                    #Calculating the maximum redshift achievable in the k band
                    DL_max_jk = np.sqrt(np.max(Av_LC_tot)/(4.*np.pi*SensDf[k]))
                    z_max_kj=fromDLtoZ(DL_max_jk,OM,Ok,OV)
                    z_seq.append(z_max_kj)
                    
                    if iLC == len(E_LC_binf):
                        #print "at the end of LIGHT CURVE BIN before the completing the spectrum range of the instrument"
                        break
    
                z_seq = np.array(z_seq)
                
                # TEST on REDSHIFT - New Version
                #
                #   1) IF z_lim := np.max(z_seq) < zold -> rejected
                #   2) IF z_lim < zDAf[j] and z_lim > zold -> count for volume until zlim
                #   3) ELSE take all the volume

                #if(np.max(z_seq)<zDAf[j]):
                
                # 1)
                if(np.max(z_seq)<z_old):
                    print 'the maximum redshift at which the light curve is visible is less than the considered z saved for rate'
                
                    #############################################
                    # It seems more correct to have "continue" than "break" since for some models the peak of luminosity can happen at higher energies
                    # and at higher energies means higher redshift and it can actually be that the L increase faster in redhift than how
                    # the luminosity distance squared increases - until a certain point, at least in principle
                    #
                    #VC_old = VC_j
                    #z_old =  zDAf[j]
                    #continue
                    #
                    # However "break" is easier to deal with - it is a bit confusing how would it goes if the Flux doesn't decrease with redhift
                    # but it doesn't happen for GRBs so it would be difficult
                    ################################################
                    
                    print " !!! NB: ASSUMPTION: with increasing z -> the flux decreases (it doesn't have to be so though, since the light curve could peak at higher energies, and the energy bin included by the instrument is larger with redshift) !!!"
                    break

                    ###############################################
                    # I SHOULD DIVIDE IN TWO DIFFERENT CASES: ACCORDINGLY TO THE LIGHT CURVE PEAKS
                    ####################################################

                # 2) changing the redshift at the superior limit of the bin
                elif (np.max(z_seq)<zDAf[j] and np.max(z_seq)>z_old):
                
                   zBINmax = np.max(z_seq)
                
                # 3) I don't need to change anything in this case
                #===========================================================================
                #Calculating the total number until redshift z_j
                #THINK ABOUT IT: I'M NOTE SURE 0.5 IS THE RIGHT WEIGHT TO IT

                
                zp1dV_i = CV_T(z_old,zBINmax,NN, z_Tab, Vc_Tab)
                
                
                ########### MARK: TEST 18 SEPT 2017: if I switch from trapezoidal to rectangular at zj it changes at the cent point the peak value -> of about 3-4%
                #
                #R_ot_i =R_ot_i+ RsDVAf[j]*zp1dV_i
                R_ot_i =R_ot_i+ (RsDVAf[j-1]*0.5 + RsDVAf[j]*0.5)*zp1dV_i # ASSUMPTION: Rate chances proportionally with redshift - and not with Volume
                # up to zM = 6 N from peaks 0.664296826243 Rect N from peaks 0.670686206875
                
                # Calculating the expected flux from the event

                Dzj_eff = fromZtoDL(zEj,OM, Ok, OV) # z most probable considering source rate proportional to volume
                Dzj_min = fromZtoDL(z_old, OM, Ok, OV)
                Dzj_max = fromZtoDL(zDAf[j], OM, Ok, OV)
                FL_j = np.amax(L_LC_totf)/(4*np.pi*Dzj_eff*Dzj_eff)  # from the most probable z in the z-bin consiredered by the rate model
                FL_jmin = np.amax(L_LC_totf)/(4*np.pi*Dzj_min*Dzj_min) # from inferior border z in the z-bin consiredered by the rate model
                FL_jmax = np.amax(L_LC_totf)/(4*np.pi*Dzj_max*Dzj_max) # from superior border z in the z-bin consiredered by the rate model
        
                zeffA.append(zEj)
                
                Flux_zeffA.append(FL_j)
                Flux_zminA.append(FL_jmin)
                Flux_zmaxA.append(FL_jmax)
                #  Having all these fluxes will help in putting ERROR ON THE LUMINOSITY DISTRIBUTION coming from the RATE bin volume

                DV_A.append(DV)
                DV_zw.append(zp1dV_i)


                maxLC = np.max(LC_totf)
                LC_totfL = list(LC_totf)
                indexPEAK = LC_totfL.index(maxLC)
                ######### AGAIN TRAPEZOIDAL RULE -> LINEAR IN REDSHIFT ##########
                #calculating single contribution to the total number of events
                Rei.append((RsDVAf[j-1]*0.5 + RsDVAf[j]*0.5)*zp1dV_i)
                Rec.append(R_ot_i)
                ziA.append(zDAf[j])
                zim1A.append(z_old)
                Vcum.append(VC_j)
                VC_old = VC_j
                z_old =  zDAf[j]
                
                ########## ADDED 18 SEPT FOR DURATION PLOTS ########


                Ti_zmin.append(tAf[np.min(indLIST_zmin)])
                Tf_zmin.append(tAf[np.max(indLIST_zmin)])
    
                Ti_zeff.append(tAf[np.min(indLIST_zeff)])
                Tf_zeff.append(tAf[np.max(indLIST_zeff)])
    
                Ti_zmax.append(tAf[np.min(indLIST_zmax)])
                Tf_zmax.append(tAf[np.max(indLIST_zmax)])
    
                timeLpArr.append(tAf[indexPEAK])
                
                F_tot = LC_totf/(4.*np.pi*dist_zeff*dist_zeff)
                F_zeff[j,:] = F_tot

        
                ##################

                del indLIST_zeff[:]
                del indLIST_zmin[:]
                del indLIST_zmax[:]
                del list_AVLC_k[:]
        
                if (j == len(zDAf) - 1):
                    break

                j = j+1
                
    ##################### ADDING LAST PART BETWEEN Z-BINS ###############################
    #
    if (np.max(z_seq)>z_old and zM -z_old>0.0000000000001):
        
            indLIST_zeff = []
            indLIST_zmin = []
            indLIST_zmax = []
            
            L_LC_totf =  np.zeros(len(tAf)-1)
            LC_totf = np.zeros(len(tAf))
            
            z_seq = []
            #### loop over energy bin####
            count = 0
            
            iLC = 0
            
            zE_last = z_effF(z_old, zM, NN)
            
            ##### MARK: ENERGY LOOP #####
            for k in range(NbinEI):
                    #print 'k',k
                    DEcov = 0.0
                    
                    Av_LC_tot = np.zeros(len(tAf)-1)
                    LC_tot = np.zeros(len(tAf))
                    ######################################################################################################################
                    # NB: EVEN HERE THERE IS AN ERROR DUE TO THE REDSHIFT BIN - INFLUENCE ON THE ENERGY BAND CONSIDERED
                    # NB: despite the broadening if the band the spectral bins remain disjoint - if I consider a single redshift
                    #
                    #if instead I use zBINmax and zBINmin this is no more true - however It wont be correct to consider a light curve
                    # given by the sum of the light curve which includes such a large z bin! I would sum more spectral band than the true
                    # what I should do is to weight the light curve calculated at each z from z_old to zDA[j] with the volume
                    # I think considering zEj is a good approximation
                    #########################################################################################################################
                    
                    
                    Emin_k = Ebin_iDf[k]*(1.+zE_last) if LC_z==1 else Ebin_iDf[k]
                    Emax_k = Ebin_sDf[k]*(1.+zE_last) if LC_z==1 else Ebin_sDf[k]
                    
                    ############ USEFUL FOR ERROR CALCULATION: NOT USED AT THE MOMENT
                    #
                    #Emin_k_zmin= Ebin_iDf[k]*(1.+z_old) if LC_z==1 else Ebin_iDf[k]
                    #Emax_k_zmin = Ebin_sDf[k]*(1.+z_old) if LC_z==1 else Ebin_sDf[k]
                    #Emin_k_zmax = Ebin_iDf[k]*(1.+zM) if LC_z==1 else  Ebin_iDf[k]
                    #Emax_k_zmax = Ebin_sDf[k]*(1.+zM) if LC_z==1 else Ebin_sDf[k]
                    #
                    ##################################################################
                    
                    countW0 = 0
                    
                    while (E_LC_binf[iLC]+hDeltaE_LC<=Emin_k and iLC< (len(E_LC_binf)-1)):
                        iLC = iLC+1
                        countW0 = countW0 +1
                        if countW0> 200:
                            print'coutW0>200', countW0
                            exit()
                    
                    
                    if (E_LC_binf[iLC]-hDeltaE_LC>Emax_k):
                        continue
                
                    LCi_Emax = E_LC_binf[iLC]+hDeltaE_LC
                    LCi_Emin = E_LC_binf[iLC]-hDeltaE_LC
                    
                    while (iLC< (len(E_LC_binf))and LCi_Emin< Emax_k and LCi_Emax > Emin_k):
                        
                        LCi_Emax = E_LC_binf[iLC]+hDeltaE_LC
                        LCi_Emin = E_LC_binf[iLC]-hDeltaE_LC
                        
                        MINe = LCi_Emin if LCi_Emin > Emin_k else Emin_k
                        MAXe = LCi_Emax if LCi_Emax < Emax_k else Emax_k
                        
                        vL = np.array(L_LCAf[iLC,:])
                        
                        Av_LC=vL[1:] + vL[:len(vL)-1]
                        Av_LC = np.array(Av_LC)
                        Av_LC = Av_LC*0.5
                        W = (MAXe - MINe)/DeltaE_LC
                        
                        Av_LC_tot = Av_LC_tot+ Av_LC* W
                        LC_tot = LC_tot+ vL* W
                        
                        DEcov = DEcov + (MAXe - MINe)
                        
                        if (LCi_Emax>= Emax_k):
                            break
                    

                        iLC = iLC+1
    
    
                    #------ while closed -----
                    if(DEcov - (Emax_k - Emin_k)>0.00000000001):
                        print 'CODE BUG: the energy band considered for the light curve is wider than the one energey band of the instrument!'
                        print 'DEcov>(Emax_k - Emin_k)'
                        print 'DEcov',DEcov
                        print 'Emax_k - Emin_k',Emax_k - Emin_k
                        exit()
                                                        
                    ##################################
                                                        
                    Av_LC_tot = Av_LC_tot*Absf[k]
                    LC_tot = LC_tot*Absf[k]
                    
                    L_LC_totf = L_LC_totf + Av_LC_tot
                    LC_totf = LC_totf + LC_tot
                                                                
                    dist_zeff = fromZtoDL(zE_last, OM, Ok, OV)
                    Lmax_k = 4.*np.pi*SensDf[k]*dist_zeff*dist_zeff
                    LthAf_zEj = Av_LC_tot[Av_LC_tot>Lmax_k]
                    indLIST_zeff_i = index_tot[LC_tot>Lmax_k]
                    indLIST_zeff.extend(indLIST_zeff_i)
                                                                            
                    dist_zmin = fromZtoDL(z_old, OM, Ok, OV)
                    Lmax_k_zmin = 4.*np.pi*SensDf[k]*dist_zmin*dist_zmin
                    LthAf_zmin = Av_LC_tot[Av_LC_tot>Lmax_k_zmin]
                    indLIST_zmin_i = index_tot[LC_totf>Lmax_k_zmin]
                    indLIST_zmin.extend(indLIST_zmin_i)
                    
                    dist_zmax = fromZtoDL(zM, OM, Ok, OV)
                    max_k_zmax = 4.*np.pi*SensDf[k]*dist_zmin*dist_zmax
                    LthAf_zmax = Av_LC_tot[Av_LC_tot>Lmax_k_zmax]
                    indLIST_zmax_i = index_tot[LC_totf>Lmax_k_zmax]
                    indLIST_zmax.extend(indLIST_zmax_i)
                    
                    
                    ##############################
                                                                                                                                                    
                                                                                                                                                    
                    #Calculating the maximum redshift achievable in the k band
                    DL_max_jk = np.sqrt(np.max(Av_LC_tot)/(4.*np.pi*SensDf[k]))
                    z_max_kj=fromDLtoZ(DL_max_jk,OM,Ok,OV)
                    z_seq.append(z_max_kj)
                    if iLC == len(E_LC_binf):
                        #print "at the end of LIGHT CURVE BIN before the completing the spectrum range of the instrument"
                        break

            
            
            #### end loop ######
            
            #if(np.max(z_seq)<z_old):
            if(np.max(z_seq)>=z_old):
                zMf = np.max(z_seq)
                if (zMf> zM):
                    zMf = zM
                print 'Imposing maximum redshift considered equal to maximum visible redshift (consistent with the rate)'
                print 'Calculating contribution of the last part of volume, from last z in rate model < maximum redshift (zM) to zM'

                ind_zM = returnIND(z_Tab,zMf)
                VC_j = Vc_Tab[ind_zM]
                
                DV = VC_j - VC_old
                
                zp1dV =DV*(0.5/(1.+z_old)+0.5/(1.+zMf))
                
                zp1dV_i = CV_T(z_old, zMf, NN, z_Tab, Vc_Tab)
        
                Reff_zM = RsDVAf[j-1]
                
                if j<len(zDAf):
                    ind_ip1 = returnIND(z_Tab,zDAf[j]) # it is the next redshift in the RATEtable
                    vjp1 = Vc_Tab[ind_ip1]
            
            
            
                    #CALCULATING RATE by
                    # 1) estimating the rate at the superior z of the bin
                    # 2) averaging the rate between first and last redshift
            
                    # 1)
                    #R_EFF(zM) calculated as an average weighted on the volume difference from VM - which is supposed to be< V(z_j) and > V(z_(j-1))
                    #linear interpolation
            
                    Reff_zM = (RsDVAf[j-1]) +  ((RsDVAf[j]-RsDVAf[j-1])/(zDAf[j] - z_old))*(zMf- z_old)
                    # weighted by volume, but it is not what I am doing for the other part
                    #Reff_zM = (RsDVAf[j-1]/(VC_j - VC_old) + RsDVAf[j]/(vjp1 - VC_j))/(1./(VC_j - VC_old) +1./(vjp1 - VC_j))
                    #R_ot_i =R_ot_i+ (RsDVAf[j-1]+Reff_zM)*0.5*zp1dV_i
            
                    # 2)
                    #*** THINK ABOUT THE 0.5 ****

                R_ot_i =R_ot_i+ (RsDVAf[j-1]+Reff_zM)*0.5*zp1dV_i
            
                # Calculating relevant Fluxes
                #
                #
                
            
                ######################################## SAVING DURATION
            
                maxLC = np.max(LC_totf)
                LC_totfL = list(LC_totf)
                indexPEAK = LC_totfL.index(maxLC)
            
            
                #########################################################
            
                Dzj_effL = dist_zeff# z most probable considering source rate proportional to volume
                Dzj_minL = dist_zmin
                Dzj_maxL = dist_zmax
            
                #################################################################
                # NB: It is possible to leave the light curve calculatd before (on the last j on redshifts) here - since in the if it is actually required it max(z_seq)>zold
                # another approach would be to redo exactly the same thing as done before instead and recalculating the light curve
                # NB: GOOD ANALITICAL CHECK, if at a z_i the test is not satified, it is possible to have it satified at larger redshift?
                #       YES: if the L increases with broader band of energy and at higher energies MORE than the flux decreases with z--- not a common case though, I would guess
                # CAN BE SOMETHING USEFUL ALSO FOR OTIMAZION THE REST OF THE CODE (put a break somewhere)
                ##################################################################
            
                FL_jL = np.amax(L_LC_totf)/(4*np.pi*Dzj_effL*Dzj_effL)  # from the most probable z in the z-bin consiredered by the rate model
                FL_jminL = np.amax(L_LC_totf)/(4*np.pi*Dzj_minL*Dzj_minL) # from inferior border z in the z-bin consiredered by the rate model
                FL_jmaxL = np.amax(L_LC_totf)/(4*np.pi*Dzj_maxL*Dzj_maxL) # from superior border z in the z-bin consiredered by the rate model
            
                Flux_zeffA.append(FL_jL)
                Flux_zminA.append(FL_jminL)
                Flux_zmaxA.append(FL_jmaxL)
            
                zeffA.append(zE_last)
            
                DV_A.append(DV)
                DV_zw.append(zp1dV)
                #print "zi",zMf,"n_zi",(RsDVAf[j-1]+Reff_zM)*0.5*zp1dV_i
                Rei.append((RsDVAf[j-1]+Reff_zM)*0.5*zp1dV_i)
                Rec.append(R_ot_i)
                ziA.append(zMf)
                zim1A.append(z_old)
                Vcum.append(VC_j)
                    
                ########## MARK: ADDED 18 SEPT FOR DURATION PLOTS ########
                if not indLIST_zmin:
                    Ti_zmin.append(indLIST_zmin)
                    Tf_zmin.append(indLIST_zmin)
                else:
                    Ti_zmin.append(tAf[np.min(np.atleast_1d(indLIST_zmin))])
                    Tf_zmin.append(tAf[np.max(indLIST_zmin)])

                if not indLIST_zeff:
                    Ti_zeff.append(indLIST_zeff)
                    Tf_zeff.append(indLIST_zeff)
                else:
                    Ti_zeff.append(tAf[np.min(indLIST_zeff)])
                    Tf_zeff.append(tAf[np.max(indLIST_zeff)])

                if not indLIST_zmax:
                    Ti_zmax.append(indLIST_zmax)
                    Tf_zmax.append(indLIST_zmax)
                else:
                    Ti_zmax.append(tAf[np.min(indLIST_zmax)])
                    Tf_zmax.append(tAf[np.max(indLIST_zmax)])
                timeLpArr.append(tAf[indexPEAK])
            
                F_tot = LC_totf/(4.*np.pi*Dzj_effL*Dzj_effL)
                F_zeff[j,:] = F_tot
                ##################
            
            else:
                print 'the maximum redshift at which the light curve is visible is less than the considered z saved for rate'
                    
                    #############################################
                    # It seems more correct to have "continue" than "break" since for some models the peak of luminosity can happen at higher energies
                    # and at higher energies means higher redshift and it can actually be that the L increase faster in redhift than how
                    # the luminosity distance squared increases - until a certain point
                    #
                    #VC_old = VC_j
                    #z_old =  zDAf[j]
                    #continue
                    #
                    # However "break" is easier to deal with - it is a bit confusing how would it goes if the Flux doesn't decrease with redhift
                    # but it doesn't happen for GRBs so it would be difficult
                    ################################################
                    
                print " !!! NB: ASSUMPTION: with increasing z -> the flux decreases (it doesn't have to be so though, since the light curve could peak at higher energies, and the energy bin included by the instrument is larger with redshift) !!!"


    Flux_zminA=np.array(Flux_zminA)
    Flux_zmaxA= np.array(Flux_zmaxA)
    Flux_zeffA = np.array(Flux_zeffA)
    zeffA = np.array(zeffA)

    out_arr=[zim1A,ziA,Rei,Rec,DV_zw,DV_A,Vcum,Flux_zminA,Flux_zmaxA,Flux_zeffA,zeffA, Ti_zmin, Tf_zmin,Ti_zeff,Tf_zeff,Ti_zmax,Tf_zmax,timeLpArr]

    out_arr2=[out_arr,np.array(F_zeff)]
    return out_arr2



######## TEST 1 PEAKS ############

SensFake = F_lim/(len(SensD))*np.ones(len(SensD))

#############################################

def InvINT_TAILS(fImin,fImax,SensDf,Ebin_iDf,Ebin_sDf,Absf,tAf,E_LC_binf, L_LCAf,zEA,zsA_f,ReiV,FluxBINSi,FluxBINSs,OM,Ok,OV):
    #print "ASSUMPTION: the duration considered for detection is the Union (in the esamble theory sense) of the times visible in each band"
    NbinEI = len(SensDf)
        
    if not(len(Ebin_iDf) == len(Ebin_sDf)and len(Ebin_sDf) == len(SensDf) and len(SensDf)== len(Absf)):
            print "Ebin_iDf",Ebin_iDf
            print "Ebin_sDf",Ebin_sDf
            print "SensDf",SensDf
            print "Absf",Absf
            print "TAILS: some energy dimensions are inconsistent"
            exit()
    
    #---calculating the difference in energy between two consecutive light curves
    ##### assuming that the difference is always the same
    
    DeltaE_LC = E_LC_binf[1] - E_LC_binf[0]
    hDE_LC = 0.5*DeltaE_LC
        
    #--calculating time bin width in LCs
        
    Dt_A = np.array(tAf[1:]) - np.array(tAf[:len(tAf)-1])
        
        
    #--ordering ergetic bins in the instrument band
    ind0 = np.argsort(Ebin_iDf)
    Ebin_iDf = np.sort(Ebin_iDf)
    Ebin_sDf = np.sort(Ebin_sDf)
    SensDf = np.array(SensDf)
    SensDf = SensDf[ind0]
        
    #--ordeing energetic bins in the LCs
    indLC = np.argsort(E_LC_binf)
    E_LC_binf = np.sort(E_LC_binf)
    L_LCAf = np.array(L_LCAf)
    L_LCAf = L_LCAf[indLC,:]
        
    n_perFlBin = np.zeros(len(FluxBINSi))
    n_perFlBin_zmax = np.zeros(len(FluxBINSi))
    n_perFlBin_zmin = np.zeros(len(FluxBINSi))
    
    NBinF = len(FluxBINSi)
    
    index_LC = np.arange(len(tAf)-1)
    
    dn_z = []
    dn_cumz = []
    
    dn_z_zmin = []
    dn_cumz_zmin = []
    
    dn_z_zmax = []
    dn_cumz_zmax = []

    dnCumj = 0.0
    dnCumj_zmin = 0.0
    dnCumj_zmax = 0.0

    z_old = 0.0
    DL_zold = 0.5*Mpc
    SUM_dzdt =0.0
    
    # to check all the number of events not in the bin are not there only because of the flux limits
    SUM_rem_zmax = 0.0
    SUM_rem_zmin = 0.0
    SUM_rem_zeff = 0.0

    OBSdur_zeff = []
    OBSdur_zmin = []
    OBSdur_zmax = []
    
    OBSdur_zeff_noz = []

    
    INDEX_Dt_zeff = []#Dth_jk_indA_Uniq_zE
    INDEX_Dt_zmin = []
    INDEX_Dt_zmax = []

    '''loop over redshift'''
    #-------------------------------------------
    for j in range(len(zsA_f)):
            ##### upper limit zBIN_j ####
            #print 'j',j,'len(zsA_f)',len(zsA_f)
            z_j= zsA_f[j]
            #print 'z: ',round(z_j,2)
            # conservative definition of distance (Distance at the hupper limit of redshif bin)
            DL_j = fromZtoDL(z_j, OM,OV,Ok)
            
            
            ##### effective z - (z in the zBIN_j most probably acccording to the volume)
            # definition of distance defined at the equavalent zEA - correspondent to $sum _bin_range z/(1+z)?
            DL_zEj = fromZtoDL(zEA[j], OM,OV,Ok)
            
            
            # definition of energy band of instrument inside the z_bin(j)
            # Consistently with the limit in energy for PEAK function
            
            
            #! NB: to be fare to consider the whole bin it should be
            EminI = (z_old+1.)*fImin if LC_z ==1 else fImin
            EmaxI = (z_j+1.)*fImax if LC_z ==1 else fImax
            
            
            
            ''' # CALCLATING ERRORS ?#
            '''
            
            iLC = 0
            maxiLC = len(E_LC_binf)-1
            dnCumjk = 0.0
            
            #If light curves start at energy higher then the maximum energy seen by the instrument - EXIT
            if (E_LC_binf[iLC]-hDE_LC>EmaxI):
                print 'ERROR: light curves are provided for energies higher then the maximum energy seen by the instrument'
                exit()
            
            #If light curves are provided up to energies lower then the minimum energy seen by the instrument - EXIT
            if (E_LC_binf[maxiLC]+hDE_LC<EminI):
                print 'ERROR: light curves are provided up to energies lower then the minimum energy seen by the instrument'
                exit()
            
            
            Dth_j = 0.0
            Av_LC_tot_AllE = np.zeros(len(tAf)-1)
            
            
            DEcov = 0.0
            Dth_jA =[]
            Dth_jk_zE = []
            Dth_jk_zold = []
            
            
            Dth_jk_indA =[]
            Dth_jk_indA_zE =[]
            Dth_jk_indA_zold =[]
            
            indD_j =[]
            indD_j_zE =[]
            indD_j_zold = []
            
            Flux_j =[]
            Flux_zmin = []
            Flux_zmax = []
            
            '''loop over energy bin'''
            #### MARK: ENERGY LOOP ####
            for k in range(NbinEI):
                
                indD_jk=[]
                indD_jk_zE = []
                indD_jk_zold = []

                DEcov = 0.0
                
                Av_LC_tot = np.zeros(len(tAf)-1)
                #! NB: to be fare to consider the whole bin it should be
                Emin_k = (z_old+1.)*Ebin_iDf[k] if LC_z ==1 else Ebin_iDf[k]
                Emax_k = (z_j+1.)*Ebin_sDf[k] if LC_z ==1 else Ebin_sDf[k]
                
                countW0 = 0
                while (E_LC_binf[iLC]+hDE_LC - Emin_k< 0.000001 and iLC< (len(E_LC_binf)-1)):
                    iLC = iLC+1
                    countW0 = countW0 +1
                    
                    if countW0> len(E_LC_binf):
                        print'cout>len(E_LC_binf)', countW0
                        exit()
            
                while (iLC< (len(E_LC_binf)-1) and E_LC_binf[iLC]+hDE_LC>Emin_k and E_LC_binf[iLC]-hDE_LC<Emax_k):
                
                    countW0 = countW0 + 1
                    
            
                    LCi_Emax = E_LC_binf[iLC]+hDE_LC
                    LCi_Emin = E_LC_binf[iLC]-hDE_LC
                

                    MINe = LCi_Emin if LCi_Emin > Emin_k else Emin_k
                    MAXe = LCi_Emax if LCi_Emax < Emax_k else Emax_k
                    
                    DEcov = DEcov + (MAXe - MINe)
                    
                    #03/04/2017: substitute this if with the one at line 1367
                    vL = np.array(L_LCAf[iLC,:])
                    Av_LC=vL[1:] + vL[:len(vL)-1]
                    Av_LC = np.array(Av_LC)

                    #Av_LC_tot = Av_LC*0.5*(MAXe - MINe)/DeltaE_LC
                    # NB: I think it should be so
                    Av_LC_tot = Av_LC_tot + Av_LC*0.5*(MAXe - MINe)/DeltaE_LC
                    
                    if (LCi_Emax - Emax_k>0.00000000001):
                        break
                    
                    iLC = iLC+1
                
                 #------ while closed -----
                 #print "end of while"
                if(DEcov - (Emax_k - Emin_k)> 0.00000000001):
                    
                    print 'CODE BUG: the energy band considered for the light curve is wider than the one energey band of the instrument!'
                    print 'DEcov>(Emax_k - Emin_k)'
                    print 'DEcov',DEcov
                    print 'Emax_k - Emin_k',Emax_k - Emin_k
                    exit()
                
                Av_LC_tot = Av_LC_tot*Absf[k]
                
                Dth_jk = 0.0
                Dth_jk_zE = 0.0
                Dth_jk_zold = 0.0
                
                Lth_j = SensDf[k]*4.*np.pi*DL_j*DL_j # conservative way: this is the minimum L which can be seen along
                
                # trying to add a bit of error:
                Lth_zEj = SensDf[k]*4.*np.pi*DL_zEj*DL_zEj #considering most probable z in the zBIN (Volume weighted)
                Lth_zold = SensDf[k]*4.*np.pi*DL_zold*DL_zold
                
                # selection of part of the light curve which is above the threshold
                LthAf_zj = Av_LC_tot[Av_LC_tot>Lth_j]
                
                LthAf_zEj = Av_LC_tot[Av_LC_tot>Lth_zEj]
                
                LthAf_zold = Av_LC_tot[Av_LC_tot>Lth_zold]

                list_AvL_A = list(Av_LC_tot)
                
                dn_fixZj_fixEk_dt =[]
                dL_fixZj_fixEk_dt =[]
                dF_fixZj_fixEk_dt =[]
                
                ################# ALTERNATIVE CODE ##############
                #
                tryIndL = index_LC[Av_LC_tot>Lth_j]
                Dth_jk_indA.extend(tryIndL)
                Flux_jk_zmax = LthAf_zj/(4.*np.pi*DL_j*DL_j)
                
                tryIndL_zE = index_LC[Av_LC_tot>Lth_zEj]
                Dth_jk_indA_zE.extend(tryIndL_zE)
                Flux_jk =LthAf_zEj/(4.*np.pi*DL_zEj*DL_zEj)
                
                tryIndL_zold = index_LC[Av_LC_tot>Lth_zold]
                Dth_jk_indA_zold.extend(tryIndL_zold)
                Flux_jk_zmin = LthAf_zold/(4.*np.pi*DL_zold*DL_zold)
                
            
                #For checking: calculate the total light curve in the band of the intrument
                Av_LC_tot_AllE = Av_LC_tot_AllE + Av_LC_tot
                
                Flux_jk = np.array(Flux_jk)
                Flux_jk_zmin = np.array(Flux_jk_zmin)
                Flux_jk_zmax = np.array(Flux_jk_zmax)
                
                indD_jk = np.array(indD_jk)
                indD_jk_zE = np.array(indD_jk_zE)
                indD_jk_zold = np.array(indD_jk_zold)
                
                #Saving the array of fluxes correspondent to each light-curve bin, whose L> L_th
                Flux_j.append(Flux_jk)
                Flux_zmin.append(Flux_jk_zmin)
                Flux_zmax.append(Flux_jk_zmax)
                #Saving the array of LC indeces correspondent to each light-curve bin, whose L> L_th
                indD_j.append(tryIndL)
                indD_j_zE.append(tryIndL_zE)
                indD_j_zold.append(tryIndL_zold)
                
            "----fine loop over energy bin"
            Dth_jk_indA = np.array(Dth_jk_indA) # zmax
            Dth_jk_indA_zE = np.array(Dth_jk_indA_zE)
            Dth_jk_indA_zold = np.array(Dth_jk_indA_zold)
    
            #Selecting the time that have been observed by the whole intrument (including the contribution of each band k)
            Dth_jk_indA_Uniq = np.unique(Dth_jk_indA)
            Dth_jk_indA_Uniq_zE = np.unique(Dth_jk_indA_zE)
            Dth_jk_indA_Uniq_zold = np.unique(Dth_jk_indA_zold)
        
            INDEX_Dt_zmax.append(np.concatenate([Dth_jk_indA_Uniq,np.zeros(len(Av_LC_tot) - len(Dth_jk_indA_Uniq))]))
            INDEX_Dt_zeff.append(np.concatenate([Dth_jk_indA_Uniq_zE,np.zeros(len(Av_LC_tot) - len(Dth_jk_indA_Uniq_zE))]))
            INDEX_Dt_zmin.append(np.concatenate([Dth_jk_indA_Uniq_zold,np.zeros(len(Av_LC_tot) - len(Dth_jk_indA_Uniq_zold))]))
        
            # === THIS PART IS REPEATED 3 TIMES TO CALCULATE HOW IT IS VARYING WITHIN THE Z BIN, taking different z as reference
            # 1) zmax
            
            if(len(Dth_jk_indA_Uniq)!=0):
                Flux_uniq_zmax = np.zeros(len(Dth_jk_indA_Uniq))
                
                Dt_A = np.array(Dt_A)
                Dth_jk_indA_Uniq = np.array(Dth_jk_indA_Uniq)
                Dth_jAm = np.matrix(Dt_A[Dth_jk_indA_Uniq])
                Dth_jA = Dt_A[Dth_jk_indA_Uniq]#/(1.+ z_j)
                #print "1.+ zj",1.+ z_j
                OnesA = np.transpose(np.matrix(np.ones(len(Dth_jk_indA_Uniq))))
                
                Dth_j_zmax = Dth_jAm*OnesA
                #print 'len(Dth_jk_indA_Uniq)',len(Dth_jk_indA_Uniq),'len(Flux_uniq_zmax)',len(Flux_uniq_zmax)
                #print "len(Dth_jk_indA_Uniq)",len(Dth_jk_indA_Uniq)
                ################################
                ####### TRY TO UNDERSTAND THE SECOND LOOP OVER C AND IF indD_j ARE IMPORTANT BECAUSE IT SEEMS EMPTY
                "loop over the visible light-curve bins"
                for b in range(len(Dth_jk_indA_Uniq)):
                    "loop over the energy bin"
                    for c in range(len(indD_j)):
                        Flux_jc=Flux_zmax[c]
                        ind_loop_c = list(indD_j[c])
                        if (Dth_jk_indA_Uniq[b] in ind_loop_c):
                            ind_cb = ind_loop_c.index(Dth_jk_indA_Uniq[b])
                            Flux_uniq_zmax[b] = Flux_uniq_zmax[b] + Flux_jc[ind_cb]
                    "----------end of energy loop"
                "-------end of loop over light-curve bins"
                "loop over total bins of fluxes in visible LC"
                
                CtotIN_IF = 0.0
                for e in range(len(Flux_uniq_zmax)):
                    "loop over bins in fluxes for graphs"
                    #print "I FLUX loop: ",e
                    FlBINcount = 0
                    for d in range(NBinF):
                        #print "II FLUX loop: ",d
                        if(Flux_uniq_zmax[e]>FluxBINSi[d] and Flux_uniq_zmax[e]<FluxBINSs[d]):
                            FlBINcount = FlBINcount +1
                            CtotIN_IF = CtotIN_IF +1
                            n_perFlBin_zmax[d] = n_perFlBin_zmax[d]+ ReiV[j]*Dth_jA[e]*(1.+ z_j)
                    if FlBINcount ==0:
                        SUM_rem_zmax = SUM_rem_zmax + ReiV[j]*Dth_jA[e]*(1.+ z_j)
                            
                "------- end of both the loops"
            else:
                Dth_j_zmax = 0.0
            Dth_j_zmax =Dth_j_zmax*(1.+ z_j)
            OBSdur_zmax.append(Dth_j_zmax)
            ###########################################
            # 2)

            if(len(Dth_jk_indA_Uniq_zE)!=0):
                Flux_uniq= np.zeros(len(Dth_jk_indA_Uniq_zE))
                Dt_A = np.array(Dt_A)#/(1.+ zEA[j])
                #print "1.+ zEA[j]",1.+ zEA[j]
                Dth_jk_indA_Uniq_zE = np.array(Dth_jk_indA_Uniq_zE)
                #Calculating the duration of the visible LC considering all the energy bands of the instrument
                Dth_jAm = np.matrix(Dt_A[Dth_jk_indA_Uniq_zE])
                Dth_jA = Dt_A[Dth_jk_indA_Uniq_zE]
                OnesA = np.transpose(np.matrix(np.ones(len(Dth_jk_indA_Uniq_zE))))
                
                Dth_j_zE = Dth_jAm*OnesA
                
                Dth_j_zE_0 = Dth_j_zE[0]
                Dth_j_zE_00 = Dth_j_zE_0[0]
                #print
                #print "TAILS"
                #print 'len(Dth_jk_indA_Uniq_zE)',len(Dth_jk_indA_Uniq_zE),'len(Flux_uniq)',len(Flux_uniq)
                "loop over the visible light-curve bins"
                for b in range(len(Dth_jk_indA_Uniq_zE)):
                    "loop over the energy bin"
                    for c in range(len(indD_j_zE)):
                        Flux_jc=Flux_j[c]
                        ind_loop_c = list(indD_j_zE[c])
                        if (Dth_jk_indA_Uniq_zE[b] in ind_loop_c):
                            ind_cb = ind_loop_c.index(Dth_jk_indA_Uniq_zE[b])
                            Flux_uniq[b] = Flux_uniq[b] + Flux_jc[ind_cb]
                    "----------end of energy loop"
                "-------end of loop over light-curve bins"
                
                "loop over total bins of fluxes in visible LC"
                for e in range(len(Flux_uniq)):
                    "loop over bins in fluxes for graphs"
                    FlBINcount = 0
                    for d in range(NBinF):
                        #print "Flux_uniq[e]",Flux_uniq[e]
                        if(Flux_uniq[e]>FluxBINSi[d] and Flux_uniq[e]<FluxBINSs[d]):
                            FlBINcount = FlBINcount+1
                            n_perFlBin[d] = n_perFlBin[d]+ ReiV[j]*Dth_jA[e]*(1.+zEA[j])
                    if FlBINcount ==0:
                        SUM_rem_zeff = SUM_rem_zeff + ReiV[j]*Dth_jA[e]*(1.+zEA[j])
                "------- end of both the loops"
            else:
                Dth_j_zE = 0.0
            OBSdur_zeff_noz.append(Dth_j_zE)
            Dth_j_zE=Dth_j_zE*(1.+zEA[j])
            Dth_j = Dth_j_zE
            OBSdur_zeff.append(Dth_j)

            #######################################
            # 3) zmin
            if(len(Dth_jk_indA_Uniq_zold)!=0):
                Flux_uniq_zmin= np.zeros(len(Dth_jk_indA_Uniq_zold))
                Dt_A = np.array(Dt_A)
                #/(1.+ z_old)
                #print "(1.+ z_old)",(1.+ z_old)
                Dth_jk_indA_Uniq = np.array(Dth_jk_indA_Uniq_zold)
                Dth_jAm = np.matrix(Dt_A[Dth_jk_indA_Uniq_zold])
                Dth_jA = Dt_A[Dth_jk_indA_Uniq_zold]
                OnesA = np.transpose(np.matrix(np.ones(len(Dth_jk_indA_Uniq_zold))))
                
                Dth_j_zmin = Dth_jAm*OnesA
                "loop over the visible light-curve bins"
                for b in range(len(Dth_jk_indA_Uniq_zold)):
                    "loop over the energy bin"
                    for c in range(len(indD_j_zold)):
                        Flux_jc=Flux_zmin[c]
                        ind_loop_c = list(indD_j_zold[c])
                        if (Dth_jk_indA_Uniq_zold[b] in ind_loop_c):
                            ind_cb = ind_loop_c.index(Dth_jk_indA_Uniq_zold[b])
                            Flux_uniq_zmin[b] = Flux_uniq_zmin[b] + Flux_jc[ind_cb]
                    "----------end of energy loop"
                "-------end of loop over light-curve bins"
                
                "loop over total bins of fluxes in visible LC"
                for e in range(len(Flux_uniq_zmin)):
                    "loop over bins in fluxes for graphs"
                    FlBINcount =0
                    for d in range(NBinF):
                        if(Flux_uniq_zmin[e]>FluxBINSi[d] and Flux_uniq_zmin[e]<FluxBINSs[d]):
                            FlBINcount =FlBINcount+1
                            n_perFlBin_zmin[d] = n_perFlBin_zmin[d]+ ReiV[j]*Dth_jA[e]*(1.+ z_old)
                        "------- end of both the loops"
                    if FlBINcount ==0:
                        SUM_rem_zmin = SUM_rem_zmin + ReiV[j]*Dth_jA[e]*(1.+ z_old)
            else:
                Dth_j_zmin = 0.0
            
            Dth_j_zmin = Dth_j_zmin*(1.+ z_old)
            OBSdur_zmin.append(Dth_j_zmin)

            # calculating the dn contribution from z_j: dV(z_j)/(zj+1)*R_v(z_j)*DeltaT
            dnjk = ReiV[j]*Dth_j_zE

            #print 'dnik',dnjk,'z_j',z_j,'z_old',z_old
    
            dnjk_zmin = ReiV[j]*Dth_j_zmin
            dnjk_zmax = ReiV[j]*Dth_j_zmax
            
            #summing to obtain the total number of events
            dnCumj = dnCumj + dnjk
            
            dnCumj_zmin = dnCumj_zmin + dnjk_zmin
            dnCumj_zmax = dnCumj_zmax + dnjk_zmax
            
            #print 'dnCumj',dnCumj
            dn_z.append(dnjk)
            dn_cumz.append(dnCumj)
            
            dn_z_zmin.append(dnjk_zmin)
            dn_cumz_zmin.append(dnCumj_zmin)
            
            dn_z_zmax.append(dnjk_zmax)
            dn_cumz_zmax.append(dnCumj_zmax)
            
            z_old =  z_j
            DL_zold = DL_j


            ##########################################

            '''fine loop over redshift'''

    OBSdur_zeff = np.array(OBSdur_zeff)
    OBSdur_zmax = np.array(OBSdur_zmax)
    OBSdur_zmin = np.array(OBSdur_zmin)
    #print
    SUM_zmin = 0.0
    for a in range(len(n_perFlBin_zmin)):
        SUM_zmin = SUM_zmin + n_perFlBin_zmin[a]
    SUM_zeff = 0.0
    for a in range(len(n_perFlBin)):
        SUM_zeff = SUM_zeff + n_perFlBin[a]
    SUM_zmax = 0.0
    for a in range(len(n_perFlBin_zmax)):
        SUM_zmax = SUM_zmax + n_perFlBin_zmax[a]
    OBSdur_zeff = np.squeeze(np.array(OBSdur_zeff))
    OBSdur_zmin = np.squeeze(np.array(OBSdur_zmin))
    OBSdur_zmax = np.squeeze(np.array(OBSdur_zmax))

    ind_zeff = [zEA,OBSdur_zmin,OBSdur_zeff_noz]


    out_arr=[n_perFlBin,dn_z,dnCumj, n_perFlBin_zmin,dn_z_zmin,dnCumj_zmin,n_perFlBin_zmax,dn_z_zmax,dnCumj_zmax,OBSdur_zeff,OBSdur_zmin,OBSdur_zmax,INDEX_Dt_zeff,INDEX_Dt_zmin,INDEX_Dt_zmax]
    return out_arr


###############################################

def COMPUTE_RATECUM(ztab,Vctab,zmodf,Rmodf,zmax_f):
    iif = 1
    z_oldf = zmodf[0]
    R_oldf = Rmodf[0]
    RCum_out = []
    zsCum_out = []
    ziCum_out = []
    RCum_i= 0.0
    
    
    
    while(zmodf[iif]<=zmax_f and iif<len(zmodf)-1):
        
        R_presf = Rmodf[iif]
        z_presf = zmodf[iif]
        
        ziCum_out.append(z_oldf)
        zsCum_out.append(z_presf)

        zdV_i = CV_T(z_oldf,z_presf,NN, ztab, Vctab)
        RCum_i = RCum_i+(R_oldf*0.5+R_presf*0.5)*zdV_i
        RCum_out.append(RCum_i)
        z_oldf = z_presf
        R_oldf = R_presf
        iif = iif+1

    if(z_oldf<zmax_f):
        zdV_i = CV_T(z_oldf,zmax_f,NN, ztab, Vctab)
        if zmax_f> zmodf[len(zmodf)-1]: # if (iif == len(Rmodf)-1):
            fac = Rmodf[len(Rmodf)-1]
    
        else:
            fac = Rmodf[iif-1]+ (Rmodf[iif] - Rmodf[iif-1])/(zmodf[iif]-z_oldf)*(zmax_f-z_oldf)
        RCum_i = RCum_i+ fac*zdV_i
        RCum_out.append(RCum_i)
        ziCum_out.append(z_oldf)
        zsCum_out.append(zmax_f)

    RCum_out = np.array(RCum_out)
    ziCum_out = np.array(ziCum_out)
    zsCum_out = np.array(zsCum_out)
    OUT_COMPUTE_RATECUM = [ziCum_out,zsCum_out,RCum_out]
    return OUT_COMPUTE_RATECUM



def modelOUT(SensFakef,Ebin_iDf,Ebin_sDf,AbsFakef,T_forLCf,E_LC_binf, LC_EAFakef,z_Tab,Vc_Tab,OBS_arr, zmaxf,zRf,Rf,OM,Ok,OV,nameRATEopt):
    

    # creating dir for files useful for flux and duration distribution
    CommLINEmkdirTXTfiles = "mkdir "+"OUTPUTS/"+dir_OUT+nameRATEopt+"/txt4ANALYSIS"
    os.system(CommLINEmkdirTXTfiles)
    
    Flim = OBS_arr[0]
    covSA = OBS_arr[1]
    t_of = OBS_arr[2]
    sd_to = OBS_arr[3]
    minI_f = OBS_arr[4]
    maxI_f = OBS_arr[5]
    SensFakef = np.atleast_1d(SensFakef)
    Ebin_iDf= np.atleast_1d(Ebin_iDf)
    Ebin_sDf= np.atleast_1d(Ebin_sDf)
    AbsFakef= np.atleast_1d(AbsFakef)

    minLlog = np.int(np.log10(np.min(SensFakef)))-2
    Dmin = fromZtoDL(zRf[1],OM,Ok,OV)
    maxLlog = np.int(np.log10(np.max(LC_max)/(4.*np.pi*Dmin*Dmin)))
    
    
    ZMAX= zmaxf
    print  "ZMAX", ZMAX
    
    
    ind_zf = np.argsort(z_Tab)
    z_Tab = np.sort(z_Tab)
    Vc_Tab = np.array(Vc_Tab)
    Vc_Tab = Vc_Tab[ind_zf]

    dim_L=np.shape(LC_EAFakef)
    if (len(T_forLCf) != dim_L[1]):
        print "ERROR TIME AND LUMINOSITY HAVE DIFFERENT DIMENSIONS IN LIGHT CURVES ",len(T_forLCf)," ",np.shape(LC_EAFakef)
        exit()
    #MARK: peak contribution-------------------------------
    outARRf2 =PEAKS_ALL(minI_f,maxI_f,SensFakef,Ebin_iDf,Ebin_sDf,AbsFakef,T_forLCf,E_LC_binf, LC_EAFakef, z_Tab, Vc_Tab,zRf,Rf,OM,Ok,OV,ZMAX)
    outARRf = np.array(outARRf2[0])
    FL_zeff = outARRf2[1]
    nameFileFLzeff = "OUTPUTS/"+dir_OUT+nameRATEopt+"/txt4ANALYSIS/FL_zeff.txt"
    datafile_id0 = open(nameFileFLzeff,'w')
    data =FL_zeff
    np.savetxt(datafile_id0,data)
    datafile_id0.close()

    ziAmodf = outARRf[0]
    zsAmodf = outARRf[1]
    Reff_if = outARRf[2]
    rateCumf = outARRf[3]

    
    ziEDGE = list(ziAmodf)
    ziEDGE.append(zsAmodf[len(zsAmodf)-1])
    ziEDGE = np.array(ziEDGE)

    LumP_zmin = outARRf[7]
    LumP_zmax = outARRf[8]
    LumP_zeff = outARRf[9]
    zE_pff = outARRf[10]

    
    Ti_zmin0 =outARRf[11]
    Tf_zmin0 =outARRf[12]
    Ti_zeff0 =outARRf[13]
    Tf_zeff0 =outARRf[14]
    Ti_zmax0 =outARRf[15]
    Tf_zmax0 =outARRf[16]

    timeLp = outARRf[17]


    R_equiv0_1f = rateCumf[len(rateCumf)-1]
    N_Lmax0_1f = f_NS*covSA*t_of*R_equiv0_1f
    print 'N from peaks',N_Lmax0_1f
    
    W = np.array(zsAmodf) -np.array(ziAmodf)
    
    ###### MARK: CALCULATING RATE OF FAVOURABLE SOURCES - independently on the survey properties #######
    #
    outRATE_UNIV = COMPUTE_RATECUM(z_Tab, Vc_Tab, zRf,Rf,zsAmodf[len(zsAmodf)-1])
    
    zi_cum=outRATE_UNIV[0]
    zs_cum=outRATE_UNIV[1]
    W_or = zs_cum - zi_cum
    rate_cum=f_NS*covSA*t_of*outRATE_UNIV[2]
    #
    #
    ##############################################################################################
    
    zAmodED = list(zsAmodf)
    zAmodED.append(ziAmodf[0])
    zAmodED = np.array(zAmodED)
    zAmodED = np.sort(zAmodED)
    
    dn_pf = f_NS*covSA*t_of*np.array(outARRf[2])

    dn_cum_pf = f_NS*covSA*t_of*np.array(outARRf[3])

    zEF_Af = zE_pff


    zbin= np.linspace(-3.,3*np.log(ZMAX),100)
    bins= np.exp(zbin)
    
    zbin_o= np.linspace(-3.,np.log(ZMAX),10)
    bins_o= np.exp(zbin_o)
    
    
    binsED = list(bins)
    binsED.append(0.0)
    binsED = np.sort(binsED)
    binsED = np.array(binsED)
    
    bins_oED = list(bins_o)
    bins_oED.append(0.0)
    bins_oED = np.sort(bins_oED)
    bins_oED = np.array(bins_oED)
    
    LC_max_av = (LC_max[1:] + LC_max[:len(LC_max)-1])*0.5
    
    DLbins = []
    Lbin =[]
    Lbin_i = [0.0]
    
    for jj in range(len(bins)):
        DL_jj = fromZtoDL(bins[jj],OM,Ok,OV)
        
        DLbins.append(DL_jj)
        Lbin_jj = np.max(LC_max_av)/(4.*np.pi*(DL_jj**2))
        Lbin.append(np.array(Lbin_jj))
        
        if (jj>0):
            Lbin_i.append(Lbin_jj)

    DLbins = np.array(DLbins)
    DLbins = np.array(DLbins)
    Lbin = np.array(Lbin)
    Lbin = np.sort(Lbin)
    Lbin_s =  Lbin
    
    Lbin_i = np.array(Lbin_i)
    Lbin_i = np.sort(Lbin_i)
    
    lin = -np.linspace(np.abs(maxLlog),np.abs(minLlog),(np.abs(minLlog)-np.abs(maxLlog)+1))
    lin = np.sort(lin)
    Lbin_o = np.power(10,lin)
    
    Lbin_o = np.array(Lbin_o)
    Lbin_o = np.sort(Lbin_o)
    Lbin_s_o =  Lbin_o[1:]
    
    Lbin_i_o = np.array(Lbin_o[:len(Lbin_o)-1])
    Lbin_i_o = np.sort(Lbin_i_o)
    
    
    LbinED =list(Lbin)
    LbinED.append(Lbin_i[0])
    LbinED = np.sort(LbinED)
    LbinED = np.array(LbinED)
    Lbin_oED = Lbin_o
    Lbin_o = Lbin_s_o
    Lbin_oED = np.array(Lbin_oED)
    Lbin_oED =np.sort(Lbin_oED)

    
    out04= InvINT_TAILS(minI_f,maxI_f,SensFakef,Ebin_iDf,Ebin_sDf,AbsFakef,T_forLCf,E_LC_binf, LC_EAFakef,zEF_Af,zsAmodf,Reff_if,Lbin_i,Lbin_s,OM,Ok,OV)

    #=============================
    #==== Zeff ====
    NDtail04 = covSA*f_NS*out04[2]
    
    import matplotlib.ticker as mtick
    # MARK: interesting plot N vs z
    plt.rcParams['figure.figsize'] = (18.0, 6.0)
    dimNdiff = np.shape(np.array(out04[1]))
    diffN_tails = np.reshape(np.array(out04[1]),(dimNdiff[0]))
    Rf0 = np.reshape(np.array(Rf[:dimNdiff[0]]),(dimNdiff[0]))
    Reff_if0 = np.reshape(np.array(Reff_if),(dimNdiff[0]))
    diffN_tails = covSA*f_NS*diffN_tails

    font = {'family' : 'Times','weight' : 'bold','size':22}

    matplotlib.rc('font', **font)
    matplotlib.rc('xtick', labelsize=18)
    matplotlib.rc('ytick', labelsize=18)
    matplotlib.rc('text', usetex=True)
    matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

    #MARK: figure dN dz
    fig = plt.figure()
    
    host = host_subplot(111, axes_class=AA.Axes)
    plt.subplots_adjust(right=0.75)
    
    par1 = host.twinx()
    par2 = host.twinx()
    
    DzEF = np.array(zsAmodf) - np.array(ziAmodf)
    
    nameFiledNT = "OUTPUTS/"+dir_OUT+nameRATEopt+"/DNdzTails.txt"
    datafile_id0 = open(nameFiledNT,'w')
    data = [zEF_Af,diffN_tails/DzEF]
    np.savetxt(datafile_id0,data)
    datafile_id0.close()

    #nameFiledNT = "OUTPUTS/"+dir_OUT+nameRATEopt+"/integral.txt"
    #datafile_id0 = open(nameFiledNT,'w')
    #data = [zEF_Af,covSA*f_NS*Reff_if0/DzEF]
    #np.savetxt(datafile_id0,data)
    #datafile_id0.close()

    nameFiledNT = "OUTPUTS/"+dir_OUT+nameRATEopt+"/dNdzPeaks.txt"
    datafile_id0 = open(nameFiledNT,'w')
    data = [zEF_Af,covSA*f_NS*Reff_if0*t_of/DzEF]
    np.savetxt(datafile_id0,data)
    datafile_id0.close()
    
    #nameFiledNT = "OUTPUTS/"+dir_OUT+nameRATEopt+"/duratfig.txt"
    #datafile_id0 = open(nameFiledNT,'w')
    #data = [zEF_Af,out04[9]]
    #np.savetxt(datafile_id0,data)
    #datafile_id0.close()
    
    #nameFiledNT = "OUTPUTS/"+dir_OUT+nameRATEopt+"/peaks_funcz.txt"
    #datafile_id0 = open(nameFiledNT,'w')
    #data = [zEF_Af,covSA*f_NS*Reff_if0*t_of/DzEF]
    #np.savetxt(datafile_id0,data)
    #datafile_id0.close()

    offset = 100
    new_fixed_axis = par2.get_grid_helper().new_fixed_axis
    par2.axis["right"] = new_fixed_axis(loc="right", axes=par2, offset=(offset, 0))
    par1.tick_params(length=6)
    par2.axis["right"].toggle(all=True)
    par1.axis["right"].toggle(all=True)
    host.set_xlim(0, np.max(zEF_Af))

    host.set_xlabel(r'$\mathbf{z}$')
    host.set_ylabel(r'$\mathbf{dN/dz\, of\, tails}')
    par1.set_ylabel(r'$\mathbf{R_V(z)~dV_C/dz~(1+z)^{-1}\, [s^{-1}]}$', fontname='Times',color = 'darkmagenta')
    par2.set_ylabel(r'$\mathbf{Light\, Curve\, Span\, LCS\, [s]}$', fontname='Times',color = 'Darkblue')
    
    
    zEL = [0.0]
    
    
    plt.hist(zEF_Af,weights = diffN_tails/DzEF, bins = ziEDGE,histtype='stepfilled', stacked=True,color = 'LightGray',rwidth=1.5, edgecolor='white')
    
    plt.hist(zEF_Af,weights = diffN_tails/DzEF, bins = ziEDGE,histtype='bar', color = 'white', edgecolor='black',alpha = 0.1,linewidth=0.3)

    
    p1, = host.plot(zEF_Af,diffN_tails/DzEF, color = 'black', label = 'TAILS', linestyle = '--',linewidth=3)
    
    p2, = par1.plot(zEF_Af,Reff_if0/DzEF, dashes=[8, 4, 2, 4, 2, 4], color = 'darkmagenta', linewidth=3)
    
    LCS = out04[9]
    p3, = par2.plot(zEF_Af[LCS>0.0],LCS[LCS>0.0], dashes=[8, 4, 2, 4, 2, 4], color ='Darkblue', linewidth=3)

    
    plt.plot(zEF_Af,covSA*f_NS*Reff_if0*t_of/DzEF, color = 'black',label = 'PEAKS',linewidth=2)
    
    host.set_ylim(0, np.max(diffN_tails/DzEF)*1.1)
    par1.set_ylim(0, np.max(Reff_if0/DzEF)*1.1) # covSA*f_NS*
    par2.set_ylim(0,np.max(out04[9])*1.1)
    
    host.axis["left"].label.set_color(p1.get_color())
    par1.axis["right"].label.set_color(p2.get_color())
    par2.axis["right"].label.set_color(p3.get_color())
    
    
    par1.get_yaxis().set_major_formatter(plt.LogFormatter(10,  labelOnlyBase=False))
    par2.get_yaxis().set_major_formatter(plt.LogFormatter(10,  labelOnlyBase=False))

    plt.rc('grid', linestyle='dashed', color='Silver')
    plt.grid(True)
    plt.xlabel('z')
    plt.ylabel('dN/dz') #tail
    plt.legend(fontsize = 15)
    
    par2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    par1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    
    nomeNvsZ = "OUTPUTS/"+dir_OUT+nameRATEopt+"/Ntails_vs_z.pdf"
    fig.savefig(nomeNvsZ)
    
    plt.rcParams['figure.figsize'] = (9.0, 6.0)
    errors= ErrFUN(N_Lmax0_1f,f_NS*covSA,t_of,sd_to,T_forLCf,Reff_if,zEF_Af,zRf,Rf,OM,Ok,OV,ZMAX)
    
    dn_L = np.array(out04[0])
    dn_L = covSA*f_NS*dn_L
    dn_z = np.array(out04[1])
    dn_z = covSA*f_NS*dn_z
    
    
    Dout_t_zeff =out04[9]
    
    dn_cumT = []
    dn_cumT_i = 0.0
    for jjf in range(len(dn_z)):
        dn_cumT_i = dn_cumT_i + dn_z[jjf]
        dn_cumT.append(dn_cumT_i)

    dn_cumT = np.array(dn_cumT)
    dn_cumT = np.squeeze(dn_cumT)
    dn_cum_tot = np.array(dn_cumT) + np.array(dn_cum_pf)

    #MARK: Cumulative_z_distributions.pdf
    """
    fig,axtot=plt.subplots(1,1,sharex=True)
    axtot = plt.subplot(1, 1, 1)
    plt.plot(ziAmodf,np.array(dn_cum_pf)+np.array(dn_cumT), color = 'black',label = 'TOTAL',linewidth=2)
    plt.plot(zi_cum,rate_cum, color = 'Darkgray', label = 'N of BNS MERGERS in {\it S}',linewidth=2)
    plt.hist(zEF_Af,weights = dn_cum_pf+dn_cumT, bins = ziEDGE,histtype='bar', color = 'white', alpha = 0.3, edgecolor='black',linewidth=0.3)
    oldY = np.zeros(len(dn_cum_pf)+2)
    ziPLOT = [ziEDGE[0]]
    ziPLOT.extend(zEF_Af)
    ziPLOT.append(ziEDGE[-1])
    newY = [dn_cum_pf[0]]
    newY.extend(dn_cum_pf)
    newY.append(dn_cum_pf[-1])
    newY = np.array(newY)
    axtot.fill_between(ziPLOT,oldY,  newY, color='#1780FB',lw=0.0, label = 'PEAKS')
    oldY = newY
    newY = [dn_cum_pf[0]+dn_cumT[0]]
    newY.extend(dn_cum_pf+dn_cumT)
    newY.append(dn_cum_pf[-1]+dn_cumT[-1])
    newY = np.array(newY)
    axtot.fill_between(ziPLOT,oldY,  newY, color='#400080',lw=0.0,label = 'TAILS')
    plt.xlabel('z')
    plt.ylabel('Number of observations')
    plt.xlim([0,ZMAX])
    plt.legend(loc="upper left",fontsize=15)
    nameCUM="OUTPUTS/"+dir_OUT+nameRATEopt+"/Cumulative_z_distributions.pdf"
    fig.savefig(nameCUM,dpi=100,bbox_inches='tight')
    """
    #==============================
    #==== Zmin ====
    
    NDtail04_zmin = covSA*f_NS*out04[5]
    dn_z_zmin = np.array(out04[4])
    dn_z_zmin = covSA*f_NS*dn_z_zmin
    dn_L_zmin = np.array(out04[3])
    dn_L_zmin = covSA*f_NS*dn_L_zmin

    Dout_t_zmin =out04[10]

    print
    print 'N from tails zmin',NDtail04_zmin
    

    #================================
    #==== Zmax ====

    NDtail04_zmax = covSA*f_NS*out04[8]
    dn_z_zmax = np.array(out04[7])
    dn_z_zmax = covSA*f_NS*dn_z_zmax
    dn_L_zmax = np.array(out04[6])
    dn_L_zmax = covSA*f_NS*dn_L_zmax

    Dout_t_zmax =out04[11]
    
    print
    print 'N from tails zmax',NDtail04_zmax

    

    aggregate_dn = np.zeros(len(bins))
    aggregate_dnL = np.zeros(len(bins))
    aggregate_dnL_zmin = np.zeros(len(bins))
    aggregate_dnL_zmax = np.zeros(len(bins))
    
    #tails
    aggregate_dnt = np.zeros(len(bins))
    aggregate_dnLt = np.zeros(len(bins))
    aggregate_dnLt_zmin = np.zeros(len(bins))
    aggregate_dnLt_zmax = np.zeros(len(bins))
    
    #peaks
    aggregate_dnp = np.zeros(len(bins))
    aggregate_dnLpold = np.zeros(len(bins))
    aggregate_dnLp = np.zeros(len(bins))
    aggregate_dnLp_zmin = np.zeros(len(bins))
    aggregate_dnLp_zmax = np.zeros(len(bins))
    #-----------
    
    #------------
    
    a = zEF_Af
    b = dn_z
    for tt in range(len(a)):
        oldBIN = 0.0
        for k in range(len(bins)):
            if(a[tt]<bins[k] and a[tt]>oldBIN):
                aggregate_dn[k] = aggregate_dn[k] + b[tt]
                aggregate_dnt[k] = aggregate_dnt[k] + b[tt]
                oldBIN = bins[k]


    for k in range(len(bins)):
        aggregate_dnL[k] = aggregate_dnL[k] + dn_L[k]
        #print "aggregate_dnL[k]",aggregate_dnL[k],"dn_L[k]",dn_L[k]
        aggregate_dnLt[k] = aggregate_dnLt[k] + dn_L[k]
        aggregate_dnLt_zmin[k] = aggregate_dnLt_zmin[k] + dn_L_zmin[k]
        aggregate_dnLt_zmax[k] = aggregate_dnLt_zmax[k] + dn_L_zmax[k]
        aggregate_dnL_zmax[k] = aggregate_dnL_zmax[k] + dn_L_zmax[k]
        aggregate_dnL_zmin[k] = aggregate_dnL_zmin[k] + dn_L_zmin[k]

    #print "aggregate_dnLt",aggregate_dnL
    #print "dn_L",dn_L
    
    SUM_dn_pfL = 0.0
    SUM_dn_pfL_zmin = 0.0
    SUM_dn_pfL_zmax = 0.0
    for i in range(len(zE_pff)):
        oldBIN = 0.0
        oldBINL = 0.0
        #I DON'T SEE ANY POINT IN LOOKING IF THE MIN or THE MAX Z IN THE RATE MODEL BIN IS IN THE GRAPH BIN
        for k in range(len(bins)):
            if(zE_pff[i]<bins[k] and zE_pff[i]>oldBIN):
                aggregate_dn[k] = aggregate_dn[k] + dn_pf[i]
                aggregate_dnp[k] = aggregate_dnp[k] + dn_pf[i]
                oldBIN = bins[k]
            if(LumP_zeff[i]<Lbin[k] and LumP_zeff[i]>oldBINL):
                #print 'if 1) i',i,'k',k
                aggregate_dnL[k] = aggregate_dnL[k] + dn_pf[i]
                SUM_dn_pfL = SUM_dn_pfL + dn_pf[i]
                aggregate_dnLp[k] = aggregate_dnLp[k] + dn_pf[i]
            if(LumP_zmin[i]<Lbin[k] and LumP_zmin[i]>oldBINL):
                aggregate_dnLp_zmin[k] = aggregate_dnLp_zmin[k] + dn_pf[i]
                aggregate_dnL_zmin[k]= aggregate_dnL_zmin[k] + dn_pf[i]
                SUM_dn_pfL_zmin = SUM_dn_pfL_zmin +dn_pf[i]
            if(LumP_zmax[i]<Lbin[k] and LumP_zmax[i]>oldBINL):
                aggregate_dnLp_zmax[k] = aggregate_dnLp_zmax[k] + dn_pf[i]
                aggregate_dnL_zmax[k]= aggregate_dnL_zmax[k] + dn_pf[i]
                
                SUM_dn_pfL_zmax = SUM_dn_pfL_zmax +dn_pf[i]
            oldBINL = Lbin[k]

    #producing images-------------------------------------------



    #MERGE_BIN_z_tot = MergeBins(binsED,aggregate_dn,bins_oED)
    dn_z_zmin = np.squeeze(dn_z_zmin)
    dn_z_zmax = np.squeeze(dn_z_zmax)
    dn_z = np.squeeze(dn_z)
    dn_pf = np.squeeze(dn_pf)
    
    dn_ztot =dn_z+dn_pf
    dn_ztot_zmax =dn_z_zmax+dn_pf
    dn_ztot_zmin =dn_z_zmin+dn_pf

    #=====================================
    # TOTAL DURATION


    ############### MARK: OUTPUT FILES FOR DURATIONS
    nameFileLC_time = "OUTPUTS/"+dir_OUT+nameRATEopt+"/txt4ANALYSIS/LC_times.txt"
    data_LC = T_forLCf
    datafile_id0 = open(nameFileLC_time,'w')
    np.savetxt(datafile_id0, data_LC)
    datafile_id0.close()

    ind_zeff = out04[12]
    ind_zmin = out04[13]
    ind_zmax = out04[14]

    nameFileDt_vis_zmin = "OUTPUTS/"+dir_OUT+nameRATEopt+"/txt4ANALYSIS/dt_vis_index_zmin.txt"
    datafile_id1_zmin = open(nameFileDt_vis_zmin,'wb')
    np.savetxt(datafile_id1_zmin, ind_zmin,fmt='%s')
    datafile_id1_zmin.close()
    
    nameFileDt_vis_zeff = "OUTPUTS/"+dir_OUT+nameRATEopt+"/txt4ANALYSIS/dt_vis_index_zeff.txt"
    datafile_id1_zeff = open(nameFileDt_vis_zeff,'wb')
    np.savetxt(datafile_id1_zeff, ind_zeff,fmt='%s')
    datafile_id1_zeff.close()

    nameFileDt_vis_zmax = "OUTPUTS/"+dir_OUT+nameRATEopt+"/txt4ANALYSIS/dt_vis_index_zmax.txt"
    datafile_id1_zmax = open(nameFileDt_vis_zmax,'wb')
    np.savetxt(datafile_id1_zmax, ind_zmax,fmt='%s')
    datafile_id1_zmax.close()
    
    nameFilePOSTPROC = "OUTPUTS/"+dir_OUT+nameRATEopt+"/txt4ANALYSIS/DURATIONS_z.txt"
    data_post = np.array([ziAmodf,Ti_zmin0,Tf_zmin0,dn_pf,dn_z_zmin,zE_pff,Ti_zeff0,Tf_zeff0,dn_pf,dn_z, zsAmodf,Ti_zmax0,Tf_zmax0,dn_pf, dn_z_zmax,timeLp])
    datafile_id = open(nameFilePOSTPROC,'w')
    np.savetxt(datafile_id, data_post.T)
    
    datafile_id.close()
    
    
    dn_ztot_merge =list(dn_z)+list(dn_pf)
    dn_ztot_zmax_merge =list(dn_z_zmax)+list(dn_pf)
    dn_ztot_zmin_merge =list(dn_z_zmin)+list(dn_pf)
    
    dn_ztot_merge = np.array(dn_ztot_merge)
    dn_ztot_zmax_merge = np.array(dn_ztot_zmax_merge)
    dn_ztot_zmin_merge = np.array(dn_ztot_zmin_merge)
    
    # DURATION PLOTS
    """#1) three peaks -> do not make sense"""
    #2) three tails
    """#3) three total -> do not make sense"""
    """#4) contibution peaks, tails to total (zeff case) -> DO NOT MAKE SENSE AGAIN"""

    margin = np.min(Dout_t_zmax[Dout_t_zmax>0.0])*0.5#np.min(Dout_tot_zmax)*0.5
    min_bin = np.min(Dout_t_zmax[Dout_t_zmax>0.0]) - margin #np.min(Dout_tot_zmax) - margin
    max_bin = np.max(Dout_t_zmin) + margin
    exp_bin = np.linspace(np.log10(min_bin.item()),np.log10(max_bin.item()),16)

    
    bin_AED = np.power(10.,exp_bin)
    Wbar = bin_AED[1:] - bin_AED[:len(bin_AED)-1]
    bin_A_av = (bin_AED[1:] + bin_AED[:len(bin_AED)-1])*0.5


    #========== GRAPHS ===============================
    MERGE_BIN_z_p = MergeBins(zAmodED,dn_pf,bins_oED)
    
    MERGE_BIN_z_tot = MergeBins(zAmodED,dn_ztot,bins_oED)
    MERGE_BIN_z_tot_zmin = MergeBins(zAmodED,dn_ztot,bins_oED)
    MERGE_BIN_z_tot_zmax = MergeBins(zAmodED,dn_ztot,bins_oED)
    
    MERGE_BIN_z_tail = MergeBins(zAmodED,dn_z,bins_oED)
    MERGE_BIN_z_tail_zmin = MergeBins(zAmodED,dn_z_zmin,bins_oED)
    MERGE_BIN_z_tail_zmax = MergeBins(zAmodED,dn_z_zmax,bins_oED)
    
    
    MERGE_BIN_F_p = MergeBins(LumP_zeff,dn_pf,Lbin_oED)
    MERGE_BIN_F_p_zmax = MergeBins(LumP_zmax,dn_pf,Lbin_oED)
    MERGE_BIN_F_p_zmin = MergeBins(LumP_zmin,dn_pf,Lbin_oED)
    
    MERGE_BIN_F_tails = MergeBins(LbinED,aggregate_dnLt,Lbin_oED)
    MERGE_BIN_F_tails_zmin = MergeBins(LbinED,aggregate_dnLt_zmin,Lbin_oED)
    MERGE_BIN_F_tails_zmax = MergeBins(LbinED,aggregate_dnLt_zmax,Lbin_oED)
    
    MERGE_BIN_F_tot = MERGE_BIN_F_tails+MERGE_BIN_F_p
    MERGE_BIN_F_tot_zmin = MERGE_BIN_F_tails_zmin+MERGE_BIN_F_p_zmin
    MERGE_BIN_F_tot_zmax = MERGE_BIN_F_tails_zmax+MERGE_BIN_F_p_zmax

    HISTcum_L_tail_zmax = []
    HISTcum_L_tail_zeff = []
    HISTcum_L_tail_zmin = []
    HISTcum_L_peak= []
    HISTcum_L_peak_zmin = []
    HISTcum_L_peak_zmax = []

    HISTcum_L_tot = []
    HISTcum_L_tot_zmin = []
    HISTcum_L_tot_zmax = []

    HISTcum_L_peaks_i = 0.0
    HISTcum_L_peaks_zmin_i = 0.0
    HISTcum_L_peaks_zmax_i = 0.0

    HISTcum_L_tail_zeff_i = 0.0
    HISTcum_L_tail_zmin_i = 0.0
    HISTcum_L_tail_zmax_i = 0.0
    
    HISTcum_L_tot_i = 0.0
    HISTcum_L_tot_zmin_i = 0.0
    HISTcum_L_tot_zmax_i = 0.0

    for kk in range(len(Lbin)):
        HISTcum_L_tail_zeff_i  = HISTcum_L_tail_zeff_i+ aggregate_dnLt[len(Lbin) -1-kk]
        HISTcum_L_tail_zeff.append(HISTcum_L_tail_zeff_i)
        HISTcum_L_tail_zmax_i  = HISTcum_L_tail_zmax_i+ aggregate_dnLt_zmax[len(Lbin) -1-kk]
        HISTcum_L_tail_zmax.append(HISTcum_L_tail_zmax_i)
        HISTcum_L_tail_zmin_i  = HISTcum_L_tail_zmin_i+ aggregate_dnLt_zmin[len(Lbin) -1-kk]
        HISTcum_L_tail_zmin.append(HISTcum_L_tail_zmin_i)
        
        HISTcum_L_peaks_i  = HISTcum_L_peaks_i+aggregate_dnLp[len(Lbin)- 1 -kk]
        HISTcum_L_peak.append(HISTcum_L_peaks_i)
        HISTcum_L_peaks_zmin_i  = HISTcum_L_peaks_zmin_i+aggregate_dnLp_zmin[len(Lbin) -1-kk]
        HISTcum_L_peak_zmin.append(HISTcum_L_peaks_zmin_i)
        HISTcum_L_peaks_zmax_i  = HISTcum_L_peaks_zmax_i+aggregate_dnLp_zmax[len(Lbin) -1-kk]
        HISTcum_L_peak_zmax.append(HISTcum_L_peaks_zmax_i)
        
        HISTcum_L_tot_i  = HISTcum_L_tot_i+ aggregate_dnL[len(Lbin) -1 -kk]
        HISTcum_L_tot.append(HISTcum_L_tot_i)
        HISTcum_L_tot_zmin_i  = HISTcum_L_tot_zmin_i+ aggregate_dnL_zmin[len(Lbin) -1-kk]
        HISTcum_L_tot_zmin.append(HISTcum_L_tot_zmin_i)
        HISTcum_L_tot_zmax_i  = HISTcum_L_tot_zmax_i+ aggregate_dnL_zmax[len(Lbin) -1-kk]
        HISTcum_L_tot_zmax.append(HISTcum_L_tot_zmax_i)

    #####################################################
    #####################################################

    Nout = [NDtail04 +N_Lmax0_1f,NDtail04,NDtail04_zmin, NDtail04_zmax,N_Lmax0_1f,rateCumf, zEF_Af,zsAmodf,Reff_if,Lbin_s,Lbin_i] #Ntot, N from Tail, N from preaks
    return Nout

OBSinfo = OBS_array_func(F_lim, CoveredSky,t_AV,SD_EXPT, minI,maxI)
name_fileOUT = "OUTPUTS/"+dir_OUT+"/results.txt"


RATEmod_label = ""
COMM2 = "mkdir OUTPUTS/"+dir_OUT +RATEmod_label
os.system(COMM2)

outTRY = modelOUT(SensD,Ebin_iD,Ebin_sD,Abs,T_forLC,E_LC_bin, LC_EA,z4Tab,Vc4Tab,OBSinfo,zmax,zR_r1,RSum1,OmegaM,Omegak,OmegaV,RATEmod_label)
Ntot_out = outTRY[0]
Ntails_out = outTRY[1]
Ntails_zmin_out = outTRY[2]
Ntails_zmax_out = outTRY[3]
Npeaks_out = outTRY[4]

print "WARNING: at the moment observations are drawn from a MAXWELL-BOLTZMANN/LOG-NORMAL DISTRIBUTION with expected value the average of the oberved time average"

TobsSTRING  = "%10.2f " %(t_AV)
SD_EXPTSTRING = "%10.2f " %(SD_EXPT)


TAV_STR = "%10.2f " %t_AV
COMM_OBSdur = "python OBSERVED_D_Lmax.py "+ "OUTPUTS/" +dir_OUT + RATEmod_label+" "+"OUTPUTS/" +dir_OUT + RATEmod_label+"/txt4ANALYSIS/DURATIONS_z.txt "+ "OUTPUTS/"+dir_OUT+RATEmod_label+"/txt4ANALYSIS/LC_times.txt "+"OUTPUTS/"+dir_OUT+RATEmod_label+"/txt4ANALYSIS/dt_vis_index_zmin.txt "+DIST_EXPT+" "+TAV_STR+" "+ SD_EXPTSTRING + " " +InFile


line1 = "N tot: %3.4f" % Ntot_out
line2 = "N tail: %3.4f" % Ntails_out
line3 = "N tail zmin: %3.4f" % Ntails_zmin_out
line4 = "N tail zmax: %3.4f" % Ntails_zmax_out
line5 = "N peaks: %3.4f \n" % Npeaks_out

with open(name_fileOUT,'a+') as f:
    f.writelines([line1,line2,line3,line4,line5])
f.close()
print
print "IF YOU WANT TO DISTRIBUTIONS OF DURATION AND FLUXES, run:"
print COMM_OBSdur

