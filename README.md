# Z-scan-PARAMETER-EXTRACTION-PROGRAM
Data processing code for nonlinear materials parameter extraction.

This open-source software is made public along with the conference paper *"Nonlinear photonics in undergraduate curriculum: hands-on training to meet the demands of a qualified workforce"*, associated with the work presented at the Optical Engineering + Applications Conference of the SPIE event Optics+Photonics 2022 under the same title. This paper contains the information regarding the theoretical framework, experimental procedures, and numerical strategies required for a successful Z-scan nonlinear characterization. This is the software that is used to extract the nonlinear and experimental setup parameters from the Z-scan measurements explained in the proceeding. Please refer to this document for in-depth explanation of the optics and data processing methodologies that support this software. In the following, the instructions to operate the software are presented, and some common errors are addressed.

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
To facilitate the understanding of the required structure for the input file, the following example table is included. It follows the aforementioned guidelines. The headers are arbitrary strings, all other values are numbers. P_in1<P_in2<P_in3. Columns labelled with headers "Z_input_power_1", "OA_input_power_1", and "CA_input_power_1" correspond to the incident power level P_in1; and so forth.

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
Additionally, an example file is included in the repository (i.e., `Z_Scan_Data.csv`), which contains the data of real measurements carried out for the nonlinear characterization of a silicon sample. This file may be used as template.


## Software initial settings

```
# 📝 Modify user settings within the box 📝
#╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
#│          

data_filename = "Z_Scan_Data2.csv" #W vs mm [Transmitted power for open aperture and closed apperture scans]
                                             # *Check templates or instructions for data file format.

install_libs = False #boolean [Enables the installation of the external libraries required by the software]
                              # *Disable after first installation to avoid launch delay.

data_res_survey = False #boolean [Enables external plotting of input data and obtained results, and functionality-related messages]
                                 # All plots appear in browser, to allow simultaneous execution of the program.

save_to_file = True #boolean [Enables automatic CSV file generation, with parameters estimations]
                              # Results folder with the file is saved in the same directory as this python script.
                                    
k_ext = 3.70187e-16 #adim [Extinction coefficient] (( k_ext( Si ) @1.56 um: 3.70187e-16 ))

L_sample = 350e-6 #m [Sample width]

tau = 90e-15 #s [Laser pulse duration] *For CW lasers set any tau and f_ref such that tau = 1/f_rep.

f_rep = 100e6 #Hz [Laser repetition rate] *For CW lasers set any tau and f_ref such that tau = 1/f_rep.

lda_0 = 1.56e-6 #m [(Central) Wavelength in vacuum]

L_setup = 152.4e-3 #m [Setup length] *Measured from lens to aperture along optical axis.

f_lens = 50e-3 #m [Lens focal length @lda_0]

n_prop_media = 1.0 #adim [Linear refractive index of external medium]

a = 1.42e-3 #m [Aperture radius] *Software can also estimate this from data, if it was generated preserving the same aperture.

n_2_search = [3e-19, 2e-17] #m^2/W [Range to initialize n_2 estimation method] *It can be also provided manually during execution.

N_obj = 1007  # [Resolution of beam cross-section for discrete Fourier transform propagation]
                # *Slightly adjust this number, if matrix shape incompatibility errors arise.

#│                                                                                                                                │
#╘════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╛
```


## Software execution and console-based settings

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

Select fitting method for beta_TPA estimation:
(0) Use Taylor expansion approximation for logarithmic expression.
(1) Perform optimization with exact logarithmic expression.
>>: 1

Select method for aperture setting:
(0) Use value in source code, provided manually.
(1) Execute estimation from data.
>>: 1
Provide a maximun value for the estimation (in m):
>>: 0.005

Select the method for n_2 search initialization:
(0) Use search interval limits in source code, provided manually.
(1) Enter search interval limits via console.
>>: 1

**Search interval limits must be positive.
 It is possible to use scientific notation;
 e.g., enter 1.5e-20, which stands for 1.5 × 10^(-20)
Enter lower limit for n_2 search interval (must be less than next one):
>>: 1e-19
Enter upper limit for n_2 search interval (must be greater than previous one):
>>: 1e-17
```


## Results

```
╓─── ▾▾▾
║ Parameters extracted:
║ r_a: 0.0014403292181069957 m
║ w_0: 1.796569577250944e-05 m
║ n_2: 4.6305771694311395e-18 m^2/W
║ beta_TPA: 2.2542349961316023e-11 m/W
║ FOM_TPA: 0.13167742952474348
╙─── ▴▴▴
```


## Examination and diagnostics mode: Data and results survey
