# Z-scan-PARAMETER-EXTRACTION-PROGRAM
Data processing code for nonlinear materials parameter extraction.

This open-source software is made public along with the conference paper *"Nonlinear photonics in undergraduate curriculum: hands-on training to meet the demands of a qualified workforce"*, associated with the work presented at the Optical Engineering + Applications Conference of the SPIE event Optics+Photonics 2022 under the same title. This paper contains the information regarding the theoretical framework, experimental procedures, and numerical strategies required for a successful Z-scan nonlinear characterization. This is the software that is used to extract the nonlinear and experimental setup parameters from the Z-scan measurements explained in the proceeding. Please refer to this document for in-depth explanation of the optics and data processing methodologies that support this software. In the following, the instructions to operate the software are presented, and some common errors are addressed.

## Requirements for measurements and input data preparation

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
