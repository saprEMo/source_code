# SaprEMo: a simplified algorithm for predicting electromagnetic observations

Here we present *saprEMo*.
*saprEMo* is a Python algorithm designed to predict the number of electromagnetic signals characterised by a specific light curve and spectrum, expected in a particular sky survey.  
By looking at past surveys, saprEMo allows us to constrain models of electromagnetic emission or event rates. Applying saprEMo to proposed astronomical missions/observing campaigns provides a perspective on their scientific impact and tests the effect of adopting different observational strategies. <br>
More details can be found in the relative paper at [arXiv:1809.08641](https://arxiv.org/pdf/1809.08641.pdf).

Contacts:<br>
*saprEMo.tool@gmail.com*<br>
*serena.vinciguerra@aei.mpg.de*

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

## Prerequisites

What things you need to install the software:
 - ```Python2.7 ```: to install ```Python2.7 ``` follow the instructions [here](https://www.python.org/downloads/release)<br>
 
required standard packges:<br>
 - numpy, scipy, csv, sys, os <br>
  
  
required *other* packages:<br>
 - astropy<br>
 
To install packages with pip run the following ```pip install name_package```<br>
For further information click [here](https://packaging.python.org/tutorials/installing-packages/)
## Abbreviations
z : redshift
LC : Light Curve
BNS : binary neutron stars

## Instructions

A step by step series of examples that tell you how to get a development env running.<br>
The main running file "runTOOL.py" requires an input set file "*INJ_FILES/name_run.txt*" and the input data stored in the directory "*TEMP/name_run/*".
You can create both of them using the file "makeINFILE.py".

### BUILDING ENVIRONMENT
```
saprEMo_PATH="/your_PATH" # change "your_PATH" with the path where you want to run saprEMo
cd $saprEMo_PATH
mkdir saprEMo
cd saprEMo
mkdir RUN
cd RUN
mkdir INJ_FILES TEMP LIBRARY OUTOUTS
cd ..
```
The following lines are needed only if the input file is created with "makeINFILE.py"
```########### needed to create input file with ** , but not for running
mkdir INPUTS 
cd INPUTS
```
if you want to add your own input data and still be able to run "makeINFILE.py", make sure to have the following structure:
```
mkdir Absorption	RATEmodels	SURVEYprop LightCurves	zDATA_file #1 directory for each input
mkdir z_CV_TABLES # tables redshift z - Comoving volume [cm^{3}] - having premade tables speed up the analysis
mkdir zDATA_file # put here non uniform pre-made z steps (e.g. from rate models)
cd ..

```

**WHAT SHOULD BE IN ALL THESE DIRECTORIES?**<br>
*Absorption*<br>
This directory should contain 2 sub-directories *MW* and *HOST_GAL*, the first for collecting information on the absorption in the Milky Way, the second for the absorption in the host galaxy. 
We only tested the Galactic absorption. In *MW*, there can be 2 further sub-directories *FILES* and *NH*. In the directorory *FILES*, add directly txt files with 2 columns: 
- 1st: center of the energy bin [keV];
- 2nd: the trasmission coefficient. <br>
Alternatively the "makeINFILE.py" can create the same file for us given the Hydrogen Column Density, NH and cross sections. In the latter case, files (fit format) containing maps of Galactic of NH should be present in the directory *Absrption/MW/NH*. The uploaded files 'lab.fit'	and 'labh.fit', contain these map. The first contain measured data, the second windowed, more details can be found in [Kalberla et al 2005](http://adsabs.harvard.edu/abs/2005A%26A...440..775K).<br>

*RATEmodels*<br>
This directory should contain the rate models of interest, if data are necessary to build Rate(z) (if analytic formulas already implemented are used, this is not necessary). Specifically it should contains sub-directories indicating the name of the relative model. According to the complexity of the considered model more sub-directories can be present, to allow the use of different sub-cases. If you add your own model and still want to create the input file with *makeINFILE.py*, you should modify the appropriate section in *makeINFILE.py*. 
The implemented cases, based on Dominik et al 2013 and Ghirlanda et al 2016, follow the instructions appearing while running ```python makeINFILE.py```.

*SURVEYprop*<br>
This directory should contain all the relevant information of the survey of interest. 
In particular should contain a sub-directory *SensCURVES*, where txt files on the survey sensitivity should be found.
This txt files should contain 4 columns:
- average energy bin, on which is based the flux limit of the second coulumn [keV]
- sensity in terms of flux limit [erg/s/cm^2]
- lower limit of the energy bin [keV]
- upper limit of the energy bin [keV]
Directly in the *SURVEYprop* dir, there should be a txt file for each survey of interest, containing other general info. 
In particular these files should contain 3 columns:
- 1st: labels
- 2nd: quantity #
- 3rd: unit of measure.<br>
Here is an example on how it should look like (all the information here reported are necessary and sufficient, unless specified, for the analysis, after # you find a mini description, which should not be present in the file):
```
minE    0.2     keV #lower limit of the entire energy band at which the survey is sensitive
maxE    12.0    keV #upper limit of the entire energy band at which the survey is sensitive
nOBS # number of observations (*)
mean observing time     7.1     s #average time for each exposure
sigma observing time     2.4    s #standard deviation for each exposure - not necessary
FoV # Filed of View (*)
FracSkyArea     0.84 # fraction of sky area (*)
SensCurve       SURVEYprop/SensCURVES/XMM_SLEW_AvE_Sens_Ehard_Esoft.txt # path to sensitivity curve
# (*): assuming no overlaps, nOBS x FoV should be the same as FracSkyArea, so you can decide which of the 2 parameter to fill (more details on the reference paper Vinciguerra et al 2018)
```
*LightCurves*<br>
This directory should contain the light curve models of interest. Specifically it should contains sub-directories indicating the name of the relative model. According to the complexity of the considered model more sub-directories can be present, to allow the use of different sub-cases. If you add your own model and still want to create the input file with *makeINFILE.py*, you should modify the appropriate section in *makeINFILE.py*. 
The final data that should be present are 2 txt files: (i)
*zDATA_file*<br>
This directory should contain 
Alternatively you can download the entire *INPUTS* directory from the git.

```
cd
mv Downloads/INPUTS $saprEMo_PATH/.
###########
```

To build the *input file* and *input dir* from the data, download the file *makeINFILE.py*, move it to the *RUN* directory, run it and follow the instructions.

```
cd
mv Downloads/makeINFILE.py $saprEMo_PATH/RUN/.
cd $saprEMo_PATH/RUN/
python makeINFILE.py name_run # substitute "name_run" with the name you like

```
**MODIFYING INPUT FILE and DIRECTLY RUN**
```
mkdir saprEMo_RUN
cd saprEMo_RUN
mkdir INJ_FILES TEMP LIBRARY OUTOUTS # 
```
When you run the "runTOOL.py" script this is what you need to have:
```
INJ_FILES/name_run.txt
TEMP/name_run/Energies_*LC_model*_LC.txt	# e.g "TEMP/name_run/Energies_Siegel_Ciolfi_2016_stableNS_LC.txt" in this txt there should be the average energy correspondent to the available light curves. E.g. if you have a single light curve obtained integrating the emission between 1keV and 5keV, the file should contain only the number 3.0 
TEMP/name_run/LC_*LC_model*.txt	 # e.g. "TEMP/name_run/LC_Siegel_Ciolfi_2016_funcE_stableNS.txt" in this txt there should be the light curves correpondent to the energy file described above. If interested in considering absortion at the HOST GALAXY, the LCs here should already be multiplied by the transmission coefficient (Done while creating input file, if done by using "makeINFILE.py").
TEMP/name_run/RATE_*rate_model*_cm3_sm1.txt # e.g "TEMP/name_run/RATE_Dominik_et_al_2013_high_nsns_cm3_sm1.txt" where there should be 2 columns: 1 one with z, the second one with rate in [cm^{-3}s^{-1}] 
TEMP/name_run/Sensitivity.txt # in this txt there should be 3 coulumns: 1st with average energy bin (tested only in [keV]); 2nd with limiting fluxes (tested only with [erg/s/cm^2]); 3rd with lower limit energy bin (tested only in [keV]); 4th with upper limit energy bin (tested only in [keV]).
TEMP/name_run/LOCAL_ABS_MW.txt	# NOT MANDATORY - this txt should contain the absorption for each sensitivity-instrument energy-band; the file should contain 4 coloumns; the 1st, 3rd and 4th as in the Sensitivity.txt file; the second should contain the transmission coefficient derived from the absorption model for each band: exp^{-N*sigma}(E) where $N$ is the effective H column denisty in [1/cm^2] and sigma the cross section [cm^2]. See paper for more details. 
TEMP/name_run/VOLUME_table_*step_z*_*cosmology*.txt # e.g. "TEMP/name_run/VOLUME_table_5.00E-04_Om03_H70.txt", here the z-step is of 0.0005; the cosmology chosen is flat with matter density 0.3 with Hubble constant 70 km/s/Mpc 
```
### Running test

Explain how to run the automated tests for this system

## Contributing

<!---
Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.
--->

## Authors

* **Serena Vinciguerra**  
* **Marica Branchesi** 
* **Riccardo Ciolfi** 
* **Ilya Mandel** 
* **Coenraad Neijssel**
* **Giulia Stratta**
<!---
(*Initial work* - [PurpleBooth](https://github.com/PurpleBooth))
--->
<!---
See also the list of [contributors](https://github.com/contributors) who participated in this project.
--->

## License and Citation

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details.
We encourage use of these data in derivative works. If you use the material provided here, please cite the paper using the reference:

```
@article{Vinciguerra:2018,
      author         = "Serena Vinciguerra, Marica Branchesi, Riccardo Ciolfi, Ilya Mandel, Coeanraad J. Neijssel, Giulia Stratta",
      title          = "{saprEMo: a simplified algorithm for predicting detections of electromagnetic transients in surveys}",
      year           = "2018",
      eprint         = "1809.08641",
      archivePrefix  = "arXiv",
      primaryClass   = "astro-ph"
}
```

## Acknowledgments
The authors thank L. Amati for his assistance with the case of THESEUS and for the useful comments.
We thank A. Belfiore, A. De Luca, M.Marelli, D. Salvietti and  A. Tiengo for the help in understanding and interpreting XMM-Newton data. We thank R. Salvaterra for the useful suggestions and discussions. 
The research leading to these results has received funding from the People Programme (Marie Curie Actions) of the European Union's Seventh Framework Programme FP7/2007-2013/ (PEOPLE-2013-ITN) under REA grant agreement no.\~[606176]. This paper reflects only the authors' view and the European Union is not liable for any use that may be made of the information contained therein.
G.S. acknowledges EGO support through a VESF fellowship (EGO-DIR-133-2015).
This research has made use of data obtained from XMMSL2, the Second XMM-Newton Slew Survey Catalogue, produced by members of the XMM SOC, the EPIC consortium, and using work carried out in the context of the EXTraS project ("Exploring the X-ray Transient and variable Sky", funded from the EU's Seventh Framework Programme under grant agreement no.\~[607452]).
This research has made use of data obtained from the 3XMM XMM-Newton serendipitous source catalogue compiled by the 10 institutes of the XMM-Newton Survey Science Centre selected by ESA.
