# Z-scan-PARAMETER-EXTRACTION-PROGRAM
Data processing code for nonlinear materials parameter extraction.

This python-based open-source software is made public along with the conference paper *"Nonlinear photonics in undergraduate curriculum: hands-on training to meet the demands of a qualified workforce"*, associated with the work presented at the Optical Engineering + Applications Conference of the SPIE event Optics+Photonics 2022 under the same title. This paper contains the information regarding the theoretical framework, experimental procedures, and numerical strategies required for a successful Z-scan nonlinear characterization. This is the software that is used to extract the nonlinear and experimental setup parameters from the Z-scan measurements explained in the proceeding. Please refer to this document for in-depth explanation of the optics and data processing methodologies that support this software. In the following, the instructions to operate the software are presented, and some common errors are addressed.

The software may be downloaded from this repository, in the file `Z_Scan_Data_Processing.py`, and is also available for online execution at [*Google Colab*](https://colab.research.google.com/drive/19VCOF2bhVFJW9PjAbUNMA_QtJN85g9i0?usp=sharing).

## Requirements for measurements and input data preparation

The software has strict requirements on the input data format, in order to run properly. So, unless it is modified by a third party to make it more flexible or adjusted to external specific needs, the following conditions must be met to ensure software stability and functionality.

\
**For the measurements**
- A complete set of measurements consists on a series of pairs of open-aperture and closed-aperture Z-scans. Every pair of Z-scans within this series (power sweep) is made for a different incident power level (input optical power).
- The power sweep should comprise at least two different incident power levels. Performing a more comprehensive power sweep (more Z-scans pairs) is suggested, to increase the acuracy of the results (by mitigating random experimetal error). Between 5 and 10 Z-scans pairs might provide enough information for the material nonlinear parameter extraction.
- Both open-aperture and closed-aperture measurements must be performed for every incident power level in the sweep.
- **All** Z-scans must share the same z-axis sampling: i.e., the same amount of samples and the same z-coordinates.
- The same aperture should be used for all the closed-aperture measurements.

*Naturally, previous conditions apply for the Z-scans subjected to the same analysis (parameter extraction process). Different data sets, for different analyses, might differ -for example- in aperture size or z-axis sampling.

\
**For the input data file**
- The input file with the information to be passed to the software, should be a comma-separated values files (with .csv extension). Although, it may be created with any program in capacity of saving in this format; e.g., Microsoft Excel.
- The structure of the file is organized by columns, and every column has a header.
- The headers for all columns are arbitrary. The user may employ them to label the data according to the followed experimental protocol.
- The file contains only information about z-coordinates and optical power values. All columns associated to z-coordinates should have values in units of milimeters (mm); and all columns associated to optical power should contain values in watts (W). For example, if a power measurement of 800 microwatts is made for a sample located at 1 inch in the z-axis, the values corresponding to this measurement should be saved as 0.0008 and 25.4.
- The sampling in the z-axis is arbitrary, and the reference value for distance along the z-axis is arbitrary. However, it is recommended that the z-coordinates selected for the measurements allow the detection of power transmission under linear conditions at the beginning and at the end of the scan (i.e., the z-coordinates of the Z-scans would be such that the first and last z-coordinates are far enough from the focus as to allow the transmission to return to the linear regime).
- The first column of the file should (only) contain the values of all the incident power levels for which a pair of open- and closed-aperture Z-scans was performed. For example, if the whole experiment comprised a power sweep of 100 mW, 200 mW and 300 mW, the first column of the file should only contain the header and the values 0.1, 0.2, and 0.3.
- The following columns always go in sets of three. Each set corresponds to a different incident power level; so, there should be as many sets as incident power levels were registered in the first column of the file. Each of these sets corresponds to a pair of open- and closed-aperture Z-scans. The first column of each set contains the information of the Z-coordinates of the Z-scans, the second column of the set contains the information of the output optical power detected during the open-aperture Z-scan, and the third column of each set contains the information of the output optical power detected during the closed-aperture Z-scan. To have a consistent file, all the z-coordinates columns (of all sets) should be identical (only varying in their header).
- These sets of columns should be registered in ascending order of incident power level.

*If the file is created from a worksheet document, make sure that the final CSV file does not include fully idle/empty columns or rows (e.g., several lines just with `,,,,,,`). This may prevent the proper reading of the file.

**To avoid non-detected errors derived from indequeate file structure, the software performs a preliminary check of the file to ensure that the basic guidelines are met.

\
To facilitate the understanding of the required structure for the input file, the following example table is included. It follows the aforementioned guidelines. The headers are arbitrary strings, all other values are numbers. `P_in1`<`P_in2`<`P_in3`. Columns labelled with headers "*Z_input_power_1*", "*OA_input_power_1*", and "*CA_input_power_1*" correspond to the incident power level `P_in1`; and so forth.

| Power_levels | Z_input_power_1 | OA_input_power_1 | CA_input_power_1 | Z_input_power_2 | OA_input_power_2 | CA_input_power_2 | Z_input_power_3 | OA_input_power_3 | CA_input_power_3 |
| ----------- | ----------- | ----------- | ----------- | ----------- | ----------- | ----------- | ----------- | ----------- | ----------- |
|P_in1 | z1 | OA1_output_P1 | CA1_output_P1 | z1 | OA2_output_P1 | CA2_output_P1 | z1 | OA3_output_P1 | CA3_output_P1|
|P_in2 | z2 | OA1_output_P2 | CA1_output_P2 | z2 | OA2_output_P2 | CA2_output_P2 | z2 | OA3_output_P2 | CA3_output_P2|
|P_in3 | z3 | OA1_output_P3 | CA1_output_P3 | z3 | OA2_output_P3 | CA2_output_P3 | z3 | OA3_output_P3 | CA3_output_P3|
| | z4 | OA1_output_P4 | CA1_output_P4 | z4 | OA2_output_P4 | CA2_output_P4 | z4 | OA3_output_P4 | CA3_output_P4|
| | z5 | OA1_output_P5 | CA1_output_P5 | z5 | OA2_output_P5 | CA2_output_P5 | z5 | OA3_output_P5 | CA3_output_P5|
| | z6 | OA1_output_P6 | CA1_output_P6 | z6 | OA2_output_P6 | CA2_output_P6 | z6 | OA3_output_P6 | CA3_output_P6|
| | z7 | OA1_output_P7 | CA1_output_P7 | z7 | OA2_output_P7 | CA2_output_P7 | z7 | OA3_output_P7 | CA3_output_P7|
| | z8 | OA1_output_P8 | CA1_output_P8 | z8 | OA2_output_P8 | CA2_output_P8 | z8 | OA3_output_P8 | CA3_output_P8|
| | z9 | OA1_output_P9 | CA1_output_P9 | z9 | OA2_output_P9 | CA2_output_P9 | z9 | OA3_output_P9 | CA3_output_P9|
| | z10 | OA1_output_P10 | CA1_output_P10 | z10 | OA2_output_P10 | CA2_output_P10 | z10 | OA3_output_P10 | CA3_output_P10|

\
The CSV file associated to the previous table (or worksheet) should look like like the following example:

```
Power_levels,Z_input_power_1,OA_input_power_1,CA_input_power_1,Z_input_power_2,OA_input_power_2,CA_input_power_2,Z_input_power_3,OA_input_power_3,CA_input_power_3
P_in1,z1,OA1_output_P1,CA1_output_P1,z1,OA2_output_P1,CA2_output_P1,z1,OA3_output_P1,CA3_output_P1
P_in2,z2,OA1_output_P2,CA1_output_P2,z2,OA2_output_P2,CA2_output_P2,z2,OA3_output_P2,CA3_output_P2
P_in3,z3,OA1_output_P3,CA1_output_P3,z3,OA2_output_P3,CA2_output_P3,z3,OA3_output_P3,CA3_output_P3
,z4,OA1_output_P4,CA1_output_P4,z4,OA2_output_P4,CA2_output_P4,z4,OA3_output_P4,CA3_output_P4
,z5,OA1_output_P5,CA1_output_P5,z5,OA2_output_P5,CA2_output_P5,z5,OA3_output_P5,CA3_output_P5
,z6,OA1_output_P6,CA1_output_P6,z6,OA2_output_P6,CA2_output_P6,z6,OA3_output_P6,CA3_output_P6
,z7,OA1_output_P7,CA1_output_P7,z7,OA2_output_P7,CA2_output_P7,z7,OA3_output_P7,CA3_output_P7
,z8,OA1_output_P8,CA1_output_P8,z8,OA2_output_P8,CA2_output_P8,z8,OA3_output_P8,CA3_output_P8
,z9,OA1_output_P9,CA1_output_P9,z9,OA2_output_P9,CA2_output_P9,z9,OA3_output_P9,CA3_output_P9
,z10,OA1_output_P10,CA1_output_P10,z10,OA2_output_P10,CA2_output_P10,z10,OA3_output_P10,CA3_output_P10
```

\
Additionally, an [example input data file](https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/blob/main/Z_Scan_Data.csv) is included in the repository (i.e., `Z_Scan_Data.csv`), which contains the data of real measurements carried out for the nonlinear characterization of a silicon sample. This file may be used as template.


## Software initial settings

After carrying out the measurements and having them organized in a file according to the previous guidelines, open the python script ´Z_Scan_Data_Processing.py´. It must be stored in the same directory as the CSV input data file. As python interpreter, *Spyder* from [*Anaconda*](https://www.anaconda.com/) is recommended for offline execution, and the [*Google Colab* notebook](https://colab.research.google.com/drive/19VCOF2bhVFJW9PjAbUNMA_QtJN85g9i0?usp=sharing) is recommended for online execution. The interpreter must support the python libraries `datetime`, `numpy`, `pandas`, `matplotlib`, `matplotlib.pyplot`, `bokeh`, `scipy`, `lmfit` and `sklearn`. But there is no problem if they are not installed, the software can perform a library installation/check.

Before running the program, users must manually configure the settings at the beginning of the source code, according to their setup specifications and analysis preferences. The user-configured section is explicitly indicated; delimited by a box, as shown below:

```
# 📝 Modify user settings within the box 📝
#╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
#│          

data_filename = "Z_Scan_Data2.csv" #W vs mm [Transmitted power for open aperture and closed apperture scans]
                                             # *Check templates or instructions for data file format.
                                             
                                                       ·
                                                       ·
                                                       ·
                                                       
N_obj = 1007  # [Resolution of beam cross-section for discrete Fourier transform propagation]
                # *Slightly adjust this number, if matrix shape incompatibility errors arise.

#│                                                                                                                                │
#╘════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╛
```

The name of the variables within user settings sections must not be modified, only their values, and preserving the type of variable. Variables are presented in the original code with default values. The settings that must be adjusted for each execution are the following:

- `data_filename`: It is a string variable, ending in `.csv`, indicating the name of the input data file. Distance and optical power stored values must be in units of milimeters (mm) and watts (W).
- `install_libs`: It is a boolean variable (i.e., admits `True` or `False` values only). If true, the software performs an external library installation/check. It is recommended to be set `True` for the first execution in an environment. Disable afterwards to avoid launch delay.
- `data_res_survey`: It is a boolean variable (i.e., admits `True` or `False` values only). If true, enables the examination and diagnostics mode of the software (detailed in a dedicated section). The software executes all the instructions of the ordinary mode, but also prints console notifications and shows graphics/plots that may be used to survey the data and results. It may be used to check if the program is running properly, and if input data show the expected trends and is read correctly by the software. If false, the ordinary mode is executed, only printing and showing the messages and plots required for the normal use of the software.
- `save_to_file`: It is a boolean variable (i.e., admits `True` or `False` values only). If true, creates a results folder in which the parameters extracted are externally saved. If false, results are only showed in console.
- `k_ext`: It is a float variable. User must provide the extinction coefficient *k* of the material of the sample. This is the imaginary part of the linear  complex refractive index of the material. It is adimensional.
- `L_sample`: It is a float variable. User must provide the thickness (in m) of the sample exposed to the Z-scans. 
- `tau`: It is a float variable. User must provide the duration (in s) of the light source pulses. If the source is not pulsed, but continuous-wave, an arbitrary pulse duration must be set, equating the repetition period (i.e., tau=1/f_rep).
- `f_rep`: It is a float variable. User must provide the repetition rate (in Hz) of the pulsed light source. If the source is not pulsed, but continuous-wave, an arbitrary repetition rate must be set, such that the repetition period equates the pulse duration (i.e., tau=1/f_rep).
- `lda_0`: It is a float variable. User must provide the central vacuum wavelength (in m) of the light pulses if source is pulsed, or the predominant wavelength (in m) of the light spectral content if the source is continuous-wave.
- `L_setup`: It is a float variable. User must provide the length (in m) of the optical setup, measured from lens to aperture along optical axis.
- `f_lens`: It is a float variable. User must provide the focal length (in m) of the lens used to focus the beam, evaluated for the wavelength given by lda_0.
- `n_prop_media`: It is a float variable. User must provide the linear (real part of the) refractive index of the external medium. It is adimensional.
- `a`: It is a float variable. User can provide the radius (in m) of the circular aperture used in the closed-aperture Z-scans. If this value is not known, software can also estimate it from data. The same aperture must be used for all closed-aperture measurements.
- `n_2_search`: It is a list of two float variables. They are the lower and upper limits passed by user to the software to perform the search of the n_2 value. It is used to reduce the time it takes to the optimization process to converge, but the software is capable of finding the value outside of the given interval if the best estimation is not within the passed range. If the n_2 is expected to be positive, **both** values must be positive, and viceversa. If the sign of the n_2 is not known *a priori*, it is possible to override this values during software execution, after viewing the normalized closed-aperture transmission trend. These values must be entered with the units m^2/W.
- `N_obj`:It is an integer variable. It is a parameter that determines the amount of samples considered for the square cross-section evaluation of the laser beam transverse profile, used for numerical propagations. The total amount of 'pixels' of each cross-section of the beam would be N_obj^2; thereby, this parameter determines the resolution of the of beam cross-section usded in discrete Fourier transform propagation. It is adimensional.

**Important.** Values for all float variables must be entered in SI units, without prefixes (i.e., meters, seconds, hertz, W, etc).

## Software execution in Ordinary Mode and console-based settings

After completing the manual adjustment of the initial settings, run the program. If `install_libs` was set as `True`, the library installation/check will be performed first; this may take several seconds.

If the library import is successful, the software will show the start message:

```
═══════════════════════════════════════╯ |╱
                                         ⚪⚟
═══════════════════════════════════════╮ |╲

  ▗▗▗        ▞▚   ▞▚  ▞▚▞  ▚▞▚
     ▞  ▗▗   ▚    ▗   ▗  ▗   ▗  ▗
   ▞            ▚   ▘   ▘  ▘  ▗  ▗
  ▝▝▝       ▚▞    ▚▞  ▚▞▚  ▘  ▘

       PARAMETER EXTRACTION PROGRAM

                    Juan José Arango. 2022 ┃
          Universidad Nacional de Colombia ┃
              Bridgewater State University ┃
```

Next, before starting with the calculations, the software will perform the input daya file check. If the file structure does not follow the [aforementioned guidelines](https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/edit/main/README.md#requirements-for-measurements-and-input-data-preparation), program execution will stop and an error message is printed, along with information of the read data. If this occurs please check input data file format and structure; check the [example input data file](https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/blob/main/Z_Scan_Data.csv) included in this repository.

```
⚠️ Amount of reported power values is NOT consistent with data provided.
⚠️ Please verify that the format of the CSV file meets the requirements.
⚠️ Please make sure that your CSV file does not contain idle columns (fully blank).
⚠️ Your data:
     Power_values(W)  Z_0.075W  OA_0.075W  ...  Z_0.250W  OA_0.250W  CA_0.250W
0              0.075       0.0   0.037787  ...       0.0    0.12265   0.050551
1              0.100       0.1   0.037787  ...       0.1    0.12266   0.050575
2              0.125       0.2   0.037787  ...       0.2    0.12268   0.050597
3              0.150       0.3   0.037789  ...       0.3    0.12268   0.050447
4              0.175       0.4   0.037790  ...       0.4    0.12269   0.050579
..               ...       ...        ...  ...       ...        ...        ...
166              NaN      16.6   0.037948  ...      16.6    0.12285   0.050200
167              NaN      16.7   0.037950  ...      16.7    0.12286   0.050242
168              NaN      16.8   0.037952  ...      16.8    0.12286   0.050255
169              NaN      16.9   0.037954  ...      16.9    0.12287   0.050296
170              NaN      17.0   0.037957  ...      17.0    0.12287   0.050288

[171 rows x 25 columns]
⚠️ All columns of your data:
['Power_values(W)', 'Z_0.075W', 'OA_0.075W', 'CA_0.075W', 'Z_0.100W', 'OA_0.100W', 'CA_0.100W', 'Z_0.125W', 'OA_0.125W', 'CA_0.125W', 'Z_0.150W', 'OA_0.150W', 'CA_0.150W', 'Z_0.175W', 'OA_0.175W', 'CA_0.175W', 'Z_0.200W', 'OA_0.200W', 'CA_0.200W', 'Z_0.225W', 'OA_0.225W', 'CA_0.225W', 'Z_0.250W', 'OA_0.250W', 'CA_0.250W']
Traceback (most recent call last):

  File "C:\Users\xxx\.conda\envs\spyder\lib\site-packages\spyder_kernels\py3compat.py", line 356, in compat_exec
    exec(code, globals, locals)

  File "c:\users\xxx\z_scan_data_processing.py", line 345, in <module>
    raise TypeError("Inadequate data format")

TypeError: Inadequate data format
```

If the file structure check is successful, the program keeps its execution. The program will show the input data as *bokeh* plots (in browser or in notebook).

 <img src="https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/blob/main/Reference%20figures/bokeh_plot%20OA%20raw%20data.png" width="400" height="250"><img src="https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/blob/main/Reference%20figures/bokeh_plot%20CA%20raw%20data.png" width="400" height="250">

User must review measurement data, and determine which are the z-coordinates ranges for which the sample behavior may be considered linear; i.e., intervals of z-coordinates for which the output transmission level did not experience notorious alterations. User would typically identify two ranges, one at the beginning and other at the end of the Z-scans, as the beam focus is expected to be located near to the center of the z-coordinates selection. The program requires that the user enter the amount of identified linear ranges. If the transmission did not return to the linear regime for one of the Z-scan extreme regions, enter 1. Next, the program ask the user to enter the upper and lower limits of each range (in mm). Users can use *bokeh* plots to zoom in, to better determine the limits.

```
Enter number of linear ranges (typically 2):
>>: 2

Enter lower limit for range 1 (in mm):
>>: 0

Enter upper limit for range 1 (in mm):
>>: 0.5

Enter lower limit for range 2 (in mm):
>>: 16.5

Enter upper limit for range 2 (in mm):
>>: 17
```

Software uses linear ranges entered to determine a reference asymptote line for each Z-scan. This reference line is used to determine effective incident power levels, and if plots have linear deviations (e.g., introduced by alignment problems) they are 'straightened', taking the slope of these lines as reference of the curve 'tilt'. Next, data is prepared to start with the calculations.

The normalized closed-aperture curves are computed, and showed as *bokeh* plots.

<img src="https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/blob/main/Reference%20figures/bokeh_plot%20CAn%20data.png" width="500" height="300">

Users must check the trend of the closed-aperture curves, in order to verify that the measurements have the expected form, and to determine the expected sign of the nonlinear index n_2. User must select between "positive" and "negative" by entering `0` or `1`, respectively.

```
According to the plotted data and the following coordinate axis convention:

 ╭╮                           ▉
 ││            ║              ▉
 ││            ║              ▉    Z
┉││┉┉┉┉┉┉┉┉┉┉┉┉║┉┉┉┉┉┉┉┉┉┉┉┉┉┉┉┉┉┉┉┉>
 ││            ║              ▉
 ││            ║              ▉
 ╰╯                           ▉
Lens        Sample         Aperture

Select the expected sign of the nonlinear refractive index (n_2):

(0) Positive.
                      χ3, n2 > 0
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃                   ◉◉          ┃
┃                 ◉    ◉        ┃
┃               ◉        ◉      ┃
┃◉ ◉ ◉       ◉          ◉ ◉ ◉┃
┃       ◉    ◉                  ┃
┃         ◉◉                    ┃
┃                       ┉┉┉┉┉> Z ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

(1) Negative.
                      χ3, n2 < 0
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃                       ┉┉┉┉┉> Z ┃
┃           ◉◉                  ┃  
┃         ◉   ◉                 ┃
┃       ◉       ◉               ┃
┃◉ ◉ ◉          ◉       ◉ ◉ ◉┃
┃                   ◉    ◉      ┃
┃                     ◉◉        ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

>>: 0
```

Next, users are allowed to select the fitting method to be used in the determination of w_0 and beta_TPA. Users can select the exact logarithmic fitting method, by entering `1`, or a Taylor expansion-based method, by entering `0`. *For more information on the methods, check the SPIE proceeding.

```
Select fitting method for beta_TPA estimation:
(0) Use Taylor expansion approximation for logarithmic expression.
(1) Perform optimization with exact logarithmic expression.
>>: 0
```

If users select the "Taylor expansion approximation" option, they must select among the two Taylor expansion-based methods, by entering `0` or `1`. *For more information on the methods, check the SPIE proceeding.

```
Select fitting method:
(0) [z->P] First normalized lorentzian fitting, then average across power levels.
(1) [P->z] First linear slopes, then global lorentzian fitting.
>>: 1
```

Next, users are allowed to select the method for determining circular aperture radius. If this value is known beforehand, users can prefer to use the value initially entered in the [source code intial settings](https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/edit/main/README.md#software-initial-settings), by entering `0`. If the value is not known, or users prefer a numerical verification, a optimization-based estimation routine to determine the radius of the aperture can be executed; in this case, enter `1`. *For more information on the aperture estimation method, check the SPIE proceeding.

```
Select method for aperture setting:
(0) Use value in source code, provided manually.
(1) Execute estimation from data.
>>: 1
```

If the "estimation from data" option is selected, the program will ask users to provide a maximum threshold (in m) for the search of the radius value estimation. This information helps the optimization routine to converge faster.

```
Provide a maximun value for the estimation (in m):
>>: 0.005
```

Finally, users are required to determine the search interval for the n_2 value estimation; which is also made with an optimization routine. Users can prefer to use the values initially entered in the [source code intial settings](https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/edit/main/README.md#software-initial-settings), by entering `0`. If users prefer to override these values, and enter new ones the must enter `1`.

```
Select the method for n_2 search initialization:
(0) Use search interval limits in source code, provided manually.
(1) Enter search interval limits via console.
>>: 1
```

If the "interval limits via console" option is selected, users must enter the lower and the upper limits for the n_2 search interval. The sign of the two limits must be the same, and it was user-determined previosly during the program execution. Instructions for this step are printed in the console, if these instructions are not followed, the n_2 estimation routine may converge in a wrong value, or may not even converge.

```
**Search interval limits must be positive.
 It is possible to use scientific notation;
 e.g., enter 1.5e-20, which stands for 1.5 × 10^(-20)
Enter lower limit for n_2 search interval (must be less than next one):
>>: 1e-19
Enter upper limit for n_2 search interval (must be greater than previous one):
>>: 1e-17
```

Finally, the program executes the numerical fitting and optimization routines, and the results of the execution are presented in the console. If the `save_to_file` configuration was set `True`, the results are also saved in a CSV file stored in a folder, that is saved in the directory where the python script is stored. 

```
╓─── ▾▾▾
║ Parameters extracted:
║ r_a: 0.0013786008230452674 m
║ w_0: 1.9298252947742058e-05 m
║ n_2: 5.593013060924584e-18 m^2/W
║ beta_TPA: 2.4652329525255596e-11 m/W
║ FOM_TPA: 0.1454331031463361
╙─── ▴▴▴
```

The results report consists on the estimations for the aperture radius r_a (if required by the user), beam's radius at focus w_0, nonlinear refractive index n_2, two-photon absorption coefficient beta_TPA, and the figure of merit of the nonlinear process FOM_TPA. Being FOM_TPA defined as n_2/(lda_0·beta_TPA).

## Software execution in Examination and Diagnostics Mode: Data and results survey

If the `data_res_survey` variable is set `True` in the source code intial settings, the program will execute the same instructions as for the Ordiary Mode operation, but it also print additional notifications through the console and shows complementary plots and graphics. This information may be useful to survey the input data, and troubleshoot potential software issues. This mode is aimed to provide a 'transparent' operation alternative for the software execution, to extend user insights about the input data or to help solving errors/unexpected results.

\
In this mode, the software prints in the console the following additional messages:

- A report after performing the input data file check:

```
📢 Amount of power values detected in power sweep: 12.

📢 Amount of reported power values is consistent with data provided.
```

- A report of the maximum effective incident power level used in the execution of the program (calculated from the input data):

```
📢 Take into account that input power levels reported in the first column
    of user-provided data are not used for calculations, just for the
    identification of the number of power levels that the power sweep comprises.
    The incident power values that are used in calculations are inferred from
    output data (in linear regime of the scan), and estimations of linear losses (given by alpha).

📢 Max. average power value used for program execution:0.12267366927752621 W

📢 Max. peak power value used for program execution:13630.407697502911 W
```

\
In this mode, the software shows the following additional *bokeh* plots and graphics:

- Unprocessed input data curves, with the corresponding reference asymptote lines, used for the calculation of effective incident power levels and curves preparation.

<img src="https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/blob/main/Reference%20figures/bokeh_plot%20OA%20raw%20data%20with%20ref%20lines.png" width="400" height="250"><img src="https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/blob/main/Reference%20figures/bokeh_plot%20CA%20raw%20data%20with%20ref%20lines.png" width="400" height="250">

- Representative plots illustrating the fitting process made with the open-aperture curves, for w_0 and beta_TPA estimation. Plotted figures vary according to the fitting method selected.

<img src="https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/blob/main/Reference%20figures/B_TPA%20first%20fitting.png" width="400" height="250"><img src="https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/blob/main/Reference%20figures/B_TPA%20second%20fitting.png" width="400" height="250">

- Beam irradiance profile at aperture plane, with and without the aperture-induced power extinction effect. The aperture size is set according to the numerically optimized radius or the one provided by the user.

<img src="https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/blob/main/Reference%20figures/beam%20at%20aperture%20plane%20OA.png" width="300" height="270"><img src="https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/blob/main/Reference%20figures/beam%20at%20aperture%20plane%20CA.png" width="300" height="270">

- Representative plot illustrating the calculation of the experimental slope of the peak-valley transmission difference vs incident peak power curve. The slope obtained is used as target for the n_2 search process.

<img src="https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/blob/main/Reference%20figures/CAn%20data%20processing.png" width="500" height="313">

- After a value for n_2 estimation has been found, the beam irradiance profile evaluated at aperture plane is shown, for sample placed at -1.7 z_0, 0, and 1.7 z_0. Where z_0 is the beam Rayleigh length.

<img src="https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/blob/main/Reference%20figures/beam%20after%20NL%20sample1.png" width="250" height="225"><img src="https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/blob/main/Reference%20figures/beam%20after%20NL%20sample2.png" width="250" height="225"><img src="https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/blob/main/Reference%20figures/beam%20after%20NL%20sample3.png" width="250" height="225">

- After a value for n_2 estimation has been found, an illustrative figure is prepared, showing the simulated closed-aperture normalized transmission, for different input power values. It is possible to observe asymmetries in this plot if there is a mismatch between the value of beta_TPA and the one found for n_2.

<img src="https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/blob/main/Reference%20figures/n_2%20power%20sweep%20CA%20z-scan%20simulation.png" width="700" height="400">

## Potential errors

The following are common errors that may arise from the execution of the software. The way of solving these issues is described below. For any physical inconsistency or software error, either derived from the code performance or the quality of the experimental data, make use of the [Examination and Diagnostics Mode](https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/edit/main/README.md#software-execution-in-examination-and-diagnostics-mode-data-and-results-survey) to troubleshoot the problem.

> **Just after running the code, the program raised the error "*Inadequate data format*".**
> 
> It occurs because the input data file provided does not meet the [file structure guidelines](https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/edit/main/README.md#requirements-for-measurements-and-input-data-preparation). Check that the file has the CSV format, the first column of the file includes the poweer sweep values only, the following columns are organized in sets of three (one set per incident power level), and there are 1+3N columns in the file (being N the number of power values reported in the first column). Check the [template file](https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/blob/main/Z_Scan_Data.csv) for more information about the proper file structure.

> **The program found some estimations for the nonlinear parameters, but they have unreasonable orders of magnitude.**
> 
> Please make sure that the data provided in the input file has units of milimeters (mm) and watts (W) for all length and power measurements. Make sure that the values provided in the [initial settings section](https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/edit/main/README.md#software-initial-settings) are in SI units, without prefixes (i.e., m, s, Hz, etc).

> **The program execution does not converge (approx. more than an hour of execution) or the sign of the estimated parameter is not correct.**
> 
> Please make sure that you are following the instructions presented by the console, regarding the sign of the entered values, their arithmetic relations (>,<,=), and their units. Please also use the [Examination and Diagnostics Mode](https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/edit/main/README.md#software-execution-in-examination-and-diagnostics-mode-data-and-results-survey) to verify that experimental data have the expected functional form and dependence with input power.

> **The simulated power detection curves provided by the software execution in 'Examination and Diagnostics Mode' exhibit a very strange functional form, not related at all with expected Z-scan typical curves.**
> 
> This may be a numerical error, related with the sampling used for the beam's field profile square cross-sections, used for the propagation calculations. It may be necessary to increase the resolution of this sampling. Try increasing the value of `N_obj` in the [initial settings section](https://github.com/juanjosearango/Z-scan-PARAMETER-EXTRACTION-PROGRAM/edit/main/README.md#software-initial-settings). This will increase the execution time and the computational resources required for the execution.

> **The beam irradiance profiles shown by the program in 'Examination and Diagnostics Mode' illustrating the simulated aperture and the simulated final detection at the aperture plane for different sample locations, appear to be extremely small within the used cross-sections, or overflow the cross-sections borders and strange interference-like patterns appear within the irradiance profile.**
> 
> This may be a numerical error, related with the size of the beam's field profile square cross-sections, used for the propagation calculations. It may occur when the experimental optical setup have atypical metric parameters (e.g., very short lens focal length and very large setup length, or viceversa). In this rare case, it may be necessary to adjust the source code to solve the problem. Look for the following lines in the code:
> ```
> # 🔎🕳KERR EFFECT PARAMETER EXTRACTION🔎🕳
> L_obj = 30*(2*10*w_0)
> ```
> `L_obj` determines the physical (not computational) size of the square cross-sections, used for the propagation calculations. This value takes as reference the radius of the beam at focues, w_0. If necessary, reduce or increase the scale factor (30, by default) until the visualization of the irradiance profiles shown by the software in 'Examination and Diagnostics Mode' becomes reasonable.
