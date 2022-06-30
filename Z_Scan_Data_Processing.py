"""
Z-scan: PARAMETER EXTRACTION PROGRAM

Data processing code for nonlinear materials parameter extraction

@author: Juan José Arango
Universidad Nacional de Colombia
Bridgewater State University
June 2022.
"""

# 📝 Modify user settings within the box 📝
#╒════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
#│          

data_filename = "Z_Scan_Data.csv" #W vs mm [Transmitted power for open aperture and closed apperture scans]
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

if(install_libs):
    import subprocess
    import sys
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'datetime'])
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'pandas'])
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'numpy'])
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'matplotlib'])
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'bokeh'])
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'scipy'])
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'lmfit'])
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'sklearn'])
import os
import datetime
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from bokeh.plotting import figure as bokeh_fig
from bokeh.plotting import show as bokeh_show
from scipy import stats
from scipy.optimize import curve_fit
from lmfit.models import LorentzianModel
from sklearn.linear_model import LinearRegression

def check_out_of_lin_range(x, lin_ranges):
    """
    Determines whether a number is out of a set of ranges.
    
        N ranges are provided in a Nx2 list. Each row states the lower and upper limits of the corresponding range.
    """
    aux=1
    for ii in np.arange(len(lin_ranges)):
        if(x>=lin_ranges[ii,0]):
            if(x<=lin_ranges[ii,1]):
                aux*=0
    return bool(aux)


def error_percentage(x, x_ref):
    """
    Returns conventional error percentage.
    
    It compares the value of x with respecto to the reference value x_ref.
    """
    return abs((x-x_ref)/x_ref)


def linear_diff(a, x0, x): 
    """Correction function used to 'straighten' a curve with certain slope a."""
    return a*(x0-x)


def linear_f(a, b, x):
    """Conventional linear function."""
    return a*x+b


def lorentzian(A, mu, sigma, x):
    """Conventional lorentzian function."""
    return (A/np.pi)*(sigma/((x-mu)**2+sigma**2))


def log_exp(P_in,F):
    """Logarithmic function used for the exact TPA functional form fitting of P_out vs P_in curves."""
    return (1/F)*np.log(1+(F*P_in))


def T_diff(X):
    """
    Returns the transmission peak-valley difference.
    
    Used for the analysis of renormalized closed aperture z-scans.
    """
    return max(X)-min(X)


def beam(X,Y, lda_0, P_PEAK, w_0, z_0, n_prop_media, z, t): #[L],[z],[w0]=mm #[lda0]=um #[t]=s #[E_norm]=N/C
    """
    Returns the electric field's complex amplitude transverse profile of a conventional gaussian beam.
    
        Parameters:
            X,Y: rectangular coordinate mesh grid used as transverse spatial coordinates for field evaluation.
            lda_0: vacuum wavelenght.
            P_PEAK: power carried by the beam.
            w_0: beam's waist at focus.
            z_0: beam's Rayleigh length.
            n_prop_media: refractive index of the external medium.
            z: longitudinal spatial coordinate for the field evaluation.
            t: temporal coordinate for the field evaluation.
    """
    k=n_prop_media*(2*pi/lda_0) #1/m
    I_0=2*P_PEAK/(pi*(w_0**2)) #W/(m^2)
    E_0_norm=np.sqrt(I_0/(2*c_light*epsilon_0*n_prop_media)) #N/C
    omega=c_light*n_prop_media*k #1/s
    w_z=w_0*np.sqrt(1+(z/z_0)**2) #m
    R_z=z*(1+(z_0/z)**2) #m
    r=np.sqrt(X**2+Y**2) #m
    E_zrt=E_0_norm*(w_0/w_z)*np.exp(-((r/w_z)**2)+(I*k*(r**2)/(2*R_z)))*np.exp(I*((n_prop_media*k*z)-(omega*t))) #N/C
    return E_zrt #N/C


def beam_NL(L_obj, N_obj, lda_0, P_PEAK, w_0, z_0, n_prop_media, z, t, L_sample, alpha, L_eff, n_2, beta_TPA):
    """"
    Returns the electric field's complex amplitude transverse profile of a gaussian beam that has propagated through a sample with nonlinear optical properties.
    
    It uses the function beam to evaluate the field profile, before adding the propagation effects.
        Parameters:
           L_obj: side length of the beam square cross-section.
           N_obj: resolution (number of divisions per dimension) of the beam square cross-section.
           lda_0, P_PEAK, w_0, z_0, n_prop_media, z, t: beam function parameters.
           L_sample: lenght of the sample along the propagation axis.
           alpha: linear attenuation coefficient of sample material.
           L_eff: effective length of the sample.
           n_2: nonlinear refractive index of sample material.
           beta_TPA: two-photon absorption coefficient of sample material.
    """
    dx = L_obj/N_obj #m #Coordinate system sampling minimum division
    x = np.arange(-N_obj*dx/2,N_obj*dx/2,dx) #m   # Coordinate space construction
    X, Y = np.meshgrid(x,x) #m # 2-dimensional transverse coordinate space
    k_0=2*pi/lda_0 #1/m
    E_zrt=beam(X,Y, lda_0, P_PEAK, w_0, z_0, n_prop_media, z, t) #N/C
    I_zrt=2*c_light*epsilon_0*n_prop_media*(abs(E_zrt)**2) #W/m^2 #Irradiance calculation
    factor_1=np.exp(-0.5*alpha*L_sample)*E_zrt #N/C
    factor_2=(1+(beta_TPA*L_eff*I_zrt))**(-0.5+(I*k_0*n_2/beta_TPA)) #adim 
    E_after_sample=np.multiply(factor_1,factor_2) #N/C #Element-wise multiplication
    return E_after_sample #N/C


def propaFT(u1, lda_0, dist, dx):
    """"
    Returns the electric field's complex amplitude transverse profile of a beam after propagation, taking the initial field profile as input.
    
    This function was built using the code provided in Ref. 3.
        Parameters:
            u1: initial beam electric field transverse profile.
            lda_0: beam vacuum wavelength.
            dist: propagation distance.
            dx: minimum division length of the mesh used for the discrete evaluation of the field profile.
    """
    k = 2*np.pi/lda_0
    tam = len(u1)
    df = 1/(tam*dx)
    f = np.arange(-tam*df/2,tam*df/2,df)
    Fx , Fy = np.meshgrid(f,f)
    H = np.exp(-1j*np.pi*lda_0*dist*(Fx**2 + Fy**2))
    U1 = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(u1)))
    return np.exp(1j*k*dist)*np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(U1*H)))


def iris(L_obj, N_obj, u1,a):
    """"
    Returns the electric field's complex amplitude transverse profile of a beam after propagation thorugh a iris/diaphragm.
    
    No diffraction effects are considered. Field gets nulled for all transverse coordinates outside of a circular region, with the given radius of the aperture.
        Parameters:
            L_obj: side length of the beam square cross-section.
            N_obj: resolution (number of divisions per dimension) of the beam square cross-section.
            u1: incident beam electric field transverse profile.
            a: circular aperture radius.
    """
    dx = L_obj/N_obj #Coordinate system sampling minimum division
    x2 = np.arange(-N_obj*dx/2,N_obj*dx/2,dx)   #Coordinate space construction
    X2, Y2 = np.meshgrid(x2,x2)  #2-dimensional transverse coordinate space
    r=np.sqrt(X2**2+Y2**2) #polar coordinate
    aperture=np.heaviside(-(r-a),1) #polar heaviside-like function
    return np.multiply(aperture,u1)


def detect(L_obj, N_obj, n_prop_media, u1):
    """"
    Returns the discrete integral of a two-dimensional function. Used to calculate power from the transverse irradiance profile of a light beam.
    
        Parameters:
            L_obj: side length of the beam square cross-section.
            N_obj: resolution (number of divisions per dimension) of the beam square cross-section.
            n_prop_media: refractive index of the external medium.
            u1: beam electric field transverse profile.
    """
    dx = L_obj/N_obj #Coordinate system sampling minimum division
    u_aux=2*n_prop_media*c_light*epsilon_0*np.abs(u1)**2 #Irradiance distribution calculation
    integration=np.sum(u_aux)*(dx**2) #Discrete integral calculation
    return integration


def Tpv_P_slope(n_2, beta_TPA, alpha, L_sample, L_eff, P_peak_max, lda_0, w_0, z_0, f_lens, L_setup, aa, n_prop_media, L_obj, N_obj,  t):
    """"
    Calculates the renormalized closed aperture Z-scan transmission for different power levels, and returns the slope of the linear relation between peak-valley transmission difference and incident power.
    
        Parameters:
            n_2: nonlinear refractive index of sample material.
            beta_TPA: two-photon absorption coefficient of sample material.
            alpha: linear attenuation coefficient of sample material.
            L_sample: lenght of the sample along the propagation axis.
            L_eff: effective length of the sample.
            P_peak_max: maximum power level of power sweep.
            lda_0: vacuum wavelength.
            w_0: beam's waist at focus.
            z_0: beam's Rayleigh length.
            f_lens: lens focal length.
            L_setup: setup length measured from lens to aperture.
            aa: aperture radius.
            n_prop_media: refractive index of the external medium.
            L_obj: side length of the beam square cross-section used for discrete Fourier transform propagation.
            N_obj: resolution (number of divisions per dimension) of the beam square cross-section used for discrete Fourier transform propagation.
            t: temporal coordinate for the field evaluation.
    """
    dx = L_obj/N_obj #Coordinate system sampling minimum division
    Tpv_vs_Power=[] #Peak-valley transmission difference, obtained for different incident power levels
    for Power in np.linspace(0.001*P_peak_max, P_peak_max, 3): #Incident power sweep
        PD_oa=[]
        PD_ca=[]
        z_sample_pos=np.linspace(-2.01*z_0, 2*z_0, 15) #z coordinates, for different sample positions
        for Z in z_sample_pos:            
            BEAM_sampled=beam_NL(L_obj, N_obj, lda_0, Power, w_0, z_0, n_prop_media, Z, t, L_sample, alpha, L_eff, n_2, beta_TPA)
            BEAM_out=propaFT(BEAM_sampled,lda_0,(L_setup-f_lens)-Z,dx)    
            BEAM_irised=iris(L_obj, N_obj, BEAM_out, aa)    
            PD_oa.append(detect(L_obj, N_obj, n_prop_media, BEAM_out)) #Open aperture z-scan
            PD_ca.append(detect(L_obj, N_obj, n_prop_media, BEAM_irised)) #Closed aperture z-scan
        BEAM_sampled=beam_NL(L_obj, N_obj, lda_0, Power, w_0, z_0, n_prop_media, 0.99*(L_setup-f_lens), t, L_sample, alpha, L_eff, n_2, beta_TPA)
        BEAM_out=propaFT(BEAM_sampled,lda_0,0.01*(L_setup-f_lens),dx)    
        BEAM_irised=iris(L_obj, N_obj, BEAM_out, aa)    
        PD_REF_oa=detect(L_obj, N_obj, n_prop_media, BEAM_out) #Open aperture reference measurement (in linear regime)
        PD_REF_ca=detect(L_obj, N_obj, n_prop_media, BEAM_irised) #Closed aperture reference measurement (in linear regime)
        PD_oa=PD_oa/PD_REF_oa #Open aperture normalized z-scan
        PD_ca=PD_ca/PD_REF_ca #Closed aperture normalized z-scan
        PD_renorm=np.divide(np.array(PD_ca), np.array(PD_oa)) #Closed aperture renormalized z-scan
        Tpv_vs_Power.append(T_diff(PD_renorm)) #Calculation of peak-valley transmission difference for closed aperture renormalized z-scan
    lin_regress = LinearRegression(fit_intercept=False).fit((np.linspace(0.001*P_peak_max, P_peak_max, 3)).reshape((-1, 1)), np.array(Tpv_vs_Power))
    slope=lin_regress.coef_[0]
    return slope


I=0+1j #adim #imaginary unit

pi=np.pi #adim #pi number

t = 0 #s #instant for field evaluation

c_light=299792458 #m/s #speed of light in vacuum

epsilon_0=8.854188e-12  #C^2/(N*m^2) #vacuum electrical permittivity

raw_data=pd.read_csv(data_filename) #W vs mm [Transmitted power for open-aperture and closed-aperture scans]

alpha=4*pi*k_ext/lda_0 #1/m #linear attenuation/absorption coefficient

k_0=2*pi/lda_0 #1/m #vacuum wavenumber

L_eff=(1-np.exp(-alpha*L_sample))/alpha #m #effective length (see Ref. 1)

colormap =cm.get_cmap("magma")

bokehpalette = [mpl.colors.rgb2hex(m) for m in colormap(np.arange(colormap.N))]

print("")
print("  ═══════════════════════════════════════╯ |╱")
print("                                         ⚪⚟")
print("  ═══════════════════════════════════════╮ |╲")
print("")
print("  ▗▗▗        ▞▚   ▞▚  ▞▚▞  ▚▞▚")
print("     ▞  ▗▗   ▚    ▗   ▗  ▗   ▗  ▗")
print("   ▞            ▚   ▘   ▘  ▘  ▗  ▗")
print("  ▝▝▝       ▚▞    ▚▞  ▚▞▚  ▘  ▘")
print("")
print("       PARAMETER EXTRACTION PROGRAM")
print("")
print("                    Juan José Arango. 2022 ┃")
print("          Universidad Nacional de Colombia ┃")
print("              Bridgewater State University ┃")
print("")

# 📜📜DATA PREPARATION📜📜

lin_data=raw_data.copy()
data=raw_data.copy()
number_of_power_sweeps=raw_data.iloc[:,0].count()

if(data_res_survey):
    print("\n📢 Amount of power values detected in power sweep: "+str(number_of_power_sweeps)+".")
if(float(number_of_power_sweeps)==float((len(list(raw_data))-1)/3)):
    if(data_res_survey):
        print("\n📢 Amount of reported power values is consistent with data provided.")
else:
    print("⚠️ Amount of reported power values is NOT consistent with data provided.")
    print("⚠️ Please verify that the format of the CSV file meets the requirements.")
    print("⚠️ Please make sure that your CSV file does not contain idle columns (fully blank).")
    print("⚠️ Your data:")
    print(raw_data)
    print("⚠️ All columns of your data:")
    print(list(data))
    raise TypeError("Inadequate data format") 

#Raw data is plotted in browser, so users can check linear ranges for data adjustment
fig= bokeh_fig(plot_width=1000, plot_height=600)
for ii in np.arange(number_of_power_sweeps):    
    dataX=np.array(raw_data.iloc[:,1+3*ii])
    dataY=np.array(raw_data.iloc[:,2+3*ii])
    fig.line(dataX, [float(i) for i in dataY], line_width=2, line_dash='solid', color='#1f77b4')
    fig.title = "Open aperture measurements"
    fig.xaxis.axis_label = "Sample position (mm)"
    fig.xaxis.axis_label_text_font_style='normal'
    fig.yaxis.axis_label = "Transmittance (W)"
    fig.yaxis.axis_label_text_font_style='normal'
bokeh_show(fig)

fig2= bokeh_fig(plot_width=1000, plot_height=600)
for ii in np.arange(number_of_power_sweeps):    
    dataX=np.array(raw_data.iloc[:,1+3*ii])
    dataY=np.array(raw_data.iloc[:,3+3*ii])
    fig2.line(dataX, [float(i) for i in dataY], line_width=2, line_dash='solid', color='#52a736')
    fig2.title = "Closed aperture measurements"
    fig2.xaxis.axis_label = "Sample position (mm)"
    fig2.xaxis.axis_label_text_font_style='normal'
    fig2.yaxis.axis_label = "Transmittance (W)"
    fig2.yaxis.axis_label_text_font_style='normal'
bokeh_show(fig2)

lin_range_num=int(input("\nEnter number of linear ranges (typically 2):\n>>: "))
lin_ranges=np.zeros((lin_range_num,2))
for ii in np.arange(len(lin_ranges)):
    lin_ranges[ii,0]=float(input("\nEnter lower limit for range "+str(ii+1)+" (in mm):\n>>: "))
    lin_ranges[ii,1]=float(input("\nEnter upper limit for range "+str(ii+1)+" (in mm):\n>>: "))
    
for ii in np.arange(number_of_power_sweeps):
    Z_header=list(lin_data)[1+3*ii]
    P_header=list(lin_data)[2+3*ii]
    P2_header=list(lin_data)[3+3*ii]
    for jj in np.arange(len(lin_data)):
        if(check_out_of_lin_range(lin_data[Z_header][jj],lin_ranges)):
            lin_data[P_header][jj]=float('nan')
            lin_data[P2_header][jj]=float('nan')

#Length coordinate information is adjusted, in order to perform all calculations with length data in m, rather than mm.
for ii in np.arange(number_of_power_sweeps):
    Z_header=list(lin_data)[1+3*ii]
    lin_data[Z_header]=1e-3*lin_data[Z_header]
    data[Z_header]=1e-3*data[Z_header]

#Data linear deviations are calculated (open aperture measurements)
lin_params=np.zeros((number_of_power_sweeps,6))
for ii in np.arange(number_of_power_sweeps):
    dataX=np.array(lin_data.iloc[:,1+3*ii])
    dataY=np.array(lin_data.iloc[:,2+3*ii])    
    mask = ~np.isnan(dataX) & ~np.isnan(dataY)
    #slope, intercept, r, p, std_err
    lin_params[ii,0], lin_params[ii,1], lin_params[ii,2], lin_params[ii,3], lin_params[ii,4] = stats.linregress(dataX[mask], dataY[mask])

#Calculation of the maximum peak power that the sampled was exposed to.
P_peak_max=lin_params[-1,1]/(tau*f_rep) #W
if(data_res_survey):
    print("\n📢 Take into account that input power levels reported in the first column\n    of user-provided data are not used for calculations, just for the\n    identification of the number of power levels that the power sweep comprises.\n    The incident power values that are used in calculations are inferred from\n    output data (in linear regime of the scan), and estimations of linear losses (given by alpha).")
    print("\n📢 Max. average power value used for program execution:"+str(lin_params[-1,1])+" W")
    print("\n📢 Max. peak power value used for program execution:"+str(P_peak_max)+" W")

#Data linear deviations are calculated (closed aperture measurements)
lin_params2=np.zeros((number_of_power_sweeps,6))
for ii in np.arange(number_of_power_sweeps):
    dataX=np.array(lin_data.iloc[:,1+3*ii])
    dataY=np.array(lin_data.iloc[:,3+3*ii])    
    mask = ~np.isnan(dataX) & ~np.isnan(dataY)
    #slope, intercept, r, p, std_err
    lin_params2[ii,0], lin_params2[ii,1], lin_params2[ii,2], lin_params2[ii,3], lin_params2[ii,4] = stats.linregress(dataX[mask], dataY[mask])
    Z_header=list(data)[1+3*ii]
    P_header=list(data)[2+3*ii]
    P2_header=list(data)[3+3*ii]
    data[P2_header]=data[P2_header]/linear_f(lin_params2[ii,0],lin_params2[ii,1],data[Z_header])
    data[P2_header]=data[P2_header]/(data[P_header]/linear_f(lin_params[ii,0],lin_params[ii,1],data[Z_header]))

#Plots for input data survey, if required by user
if(data_res_survey):
    fig4= bokeh_fig(plot_width=1000, plot_height=600)
    for ii in np.arange(number_of_power_sweeps):    
        dataX=np.array(data.iloc[:,1+3*ii])
        dataY=np.array(data.iloc[:,2+3*ii])
        dataY2=lin_params[ii,0]*dataX+lin_params[ii,1]
        fig4.line(dataX, [float(i) for i in dataY], line_width=2, line_dash='solid', color='#1f77b4')
        fig4.line(dataX, [float(i) for i in dataY2], line_width=2, line_dash='dashed', color='black')
        fig4.title = "Open aperture measurements (solid) and linear correction lines (dashed)"
        fig4.xaxis.axis_label = "Sample position (m)"
        fig4.xaxis.axis_label_text_font_style='normal'
        fig4.yaxis.axis_label = "Transmittance (W)"
        fig4.yaxis.axis_label_text_font_style='normal'
    bokeh_show(fig4)
    
    fig6= bokeh_fig(plot_width=1000, plot_height=600)
    for ii in np.arange(number_of_power_sweeps):    
        dataX=np.array(data.iloc[:,1+3*ii])
        dataY=np.array(raw_data.iloc[:,3+3*ii])
        dataY2=lin_params2[ii,0]*dataX+lin_params2[ii,1]
        fig6.line(dataX, [float(i) for i in dataY], line_width=2, line_dash='solid', color='#52a736')
        fig6.line(dataX, [float(i) for i in dataY2], line_width=2, line_dash='dashed', color='black')
        fig6.title = "Closed aperture measurements (solid) and linear correction lines (dashed)"
        fig6.xaxis.axis_label = "Sample position (m)"
        fig6.xaxis.axis_label_text_font_style='normal'
        fig6.yaxis.axis_label = "Transmittance (W)"
        fig6.yaxis.axis_label_text_font_style='normal'
    bokeh_show(fig6)

#Renormalized closed aperture transmission data is plotted in browser, so users can check power trend and determine the expected sign of the third-order susceptibility
fig3= bokeh_fig(plot_width=1000, plot_height=600)
for ii in np.arange(number_of_power_sweeps):    
    dataX=np.array(data.iloc[:,1+3*ii])
    dataY=np.array(data.iloc[:,3+3*ii])
    fig3.line(dataX, [float(i) for i in dataY], line_width=2, line_dash='solid', color='#64b066')
    fig3.title = "Closed aperture RENORMALIZED measurements"
    fig3.xaxis.axis_label = "Sample position (m)"
    fig3.xaxis.axis_label_text_font_style='normal'
    fig3.yaxis.axis_label = "Renormalized transmittance (adim)"
    fig3.yaxis.axis_label_text_font_style='normal'
bokeh_show(fig3)

negative_n_2=bool(int(input("\nAccording to the plotted data and the following coordinate axis convention:\n\n ╭╮                           ▉\n ││            ║              ▉\n ││            ║              ▉    Z\n┉││┉┉┉┉┉┉┉┉┉┉┉┉║┉┉┉┉┉┉┉┉┉┉┉┉┉┉┉┉┉┉┉┉>\n ││            ║              ▉\n ││            ║              ▉\n ╰╯                           ▉\nLens        Sample         Aperture\n\nSelect the expected sign of the nonlinear refractive index (n_2):\n\n(0) Positive.\n                      χ3, n2 > 0\n┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓\n┃                   ◉◉          ┃\n┃                 ◉    ◉        ┃\n┃               ◉        ◉      ┃\n┃◉ ◉ ◉         ◉           ◉ ◉ ◉┃\n┃       ◉    ◉                  ┃\n┃         ◉◉                    ┃\n┃                      ┉┉┉┉┉> Z ┃\n┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛\n\n(1) Negative.\n                      χ3, n2 < 0\n┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓\n┃                      ┉┉┉┉┉> Z ┃\n┃           ◉◉                  ┃  \n┃         ◉   ◉                 ┃\n┃       ◉       ◉               ┃\n┃◉ ◉ ◉            ◉        ◉ ◉ ◉┃\n┃                   ◉    ◉      ┃\n┃                     ◉◉        ┃\n┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛\n\n>>: ")))

# 🔴🔴TWO-PHOTON ABSORPTION PARAMETER EXTRACTION👻👻

#Fitting method is selected
use_log_fit=bool(int(input("\nSelect fitting method for beta_TPA estimation:\n(0) Use Taylor expansion approximation for logarithmic expression.\n(1) Perform optimization with exact logarithmic expression.\n>>: ")))
if(use_log_fit):
    #Data is standardized
    for ii in np.arange(number_of_power_sweeps):
        Z_header=list(data)[1+3*ii]
        P_header=list(data)[2+3*ii]
        data[P_header]=np.exp(alpha*L_eff)*(data[P_header]+linear_diff(lin_params[ii,0],0,data[Z_header]))
    
    #Exact fitting is performed along P_in-axis, for every sample position z
    F=np.zeros(len(data))
    P_in=lin_params[:,1]*np.exp(alpha*L_sample) #W #Estimation of input power within sample, from linear output power
    for jj in np.arange(len(data)):
        P_out=np.array(data.iloc[jj,2:3*number_of_power_sweeps:3])
        popt, pcov = curve_fit(log_exp, P_in, P_out, bounds=(0,np.inf))
        F[jj]=popt[0]  
    
    #Plots for input data survey, if required by user
    if(data_res_survey):
        fig12= bokeh_fig(plot_width=1000, plot_height=600)
        dataX=np.array(P_in)
        dataY=np.array(data.iloc[int(np.floor(len(data)/2)), 2:3*number_of_power_sweeps:3])
        dataY2=log_exp(dataX,F[int(np.floor(len(data)/2))])
        fig12.scatter(dataX, dataY, color='#1374d1', size=10)
        fig12.line(dataX, [float(i) for i in dataY2], line_width=2, line_dash='dashed', color='#A80000')
        fig12.title = "Variation of the open aperture transmission as average input power is swept\nObtained from experimental open aperture normalized curves at z="+str(data.iloc[int(np.floor(len(data)/2)), 1])+".\nExperimental results are scattered dots.\nThe dashed curve corresponds to the parametrized logarithmic function fitted to extract the parameter that describes data trend."
        fig12.xaxis.axis_label = "Average input power (W)"
        fig12.xaxis.axis_label_text_font_style='normal'
        fig12.yaxis.axis_label = "NON-Normalized transmission (adim)"
        fig12.yaxis.axis_label_text_font_style='normal'
        bokeh_show(fig12)
    
    #The function F previously obtained is fitted as a lorentzian
    Z_positions=data.iloc[:,1] #m
    mod = LorentzianModel()
    pars = mod.guess(F, x=Z_positions)
    out = mod.fit(F, pars, x=Z_positions)
    
    #Plots for results survey, if required by user
    if(data_res_survey):
        fig7= bokeh_fig(plot_width=1000, plot_height=600)
        for ii in np.arange(number_of_power_sweeps):    
            dataX=np.array(Z_positions)
            dataY=np.array(F)
            dataY2=lorentzian(0.85*pars.get('amplitude').value, pars.get('center').value, pars.get('sigma').value, Z_positions)
            fig7.line(dataX, [float(i) for i in dataY], line_width=2, line_dash='solid', color='#1374d1')
            fig7.line(dataX, [float(i) for i in dataY2], line_width=2, line_dash='dashed', color='#A80000')
            fig7.title = "TPA fitting: Exact logarithmic method results. \nStandardized experimental data (solid) and fitting (dashed)."
            fig7.xaxis.axis_label = "Sample position (m)"
            fig7.xaxis.axis_label_text_font_style='normal'
            fig7.yaxis.axis_label = "Log. parameter value F (1 / W)"
            fig7.yaxis.axis_label_text_font_style='normal'
        bokeh_show(fig7)
    
    #TPA absorption is obtained from the amplitude of the lorentzian. The Rayleigh length and beam waist at focus are obtained from the lorentzian FWHM.
    beta_TPA=0.85*pars.get('amplitude').value*(lda_0*tau*f_rep)/(2*pi*L_eff)
    z_0=pars.get('sigma').value
    w_0=np.sqrt(lda_0*z_0/pi)
    
else:
    use_method_1=bool(int(input("\nSelect fitting method:\n(0) [z->P] First normalized lorentzian fitting, then average across power levels.\n(1) [P->z] First linear slopes, then global lorentzian fitting.\n>>: ")))
    if(use_method_1):
        #Data is standardized
        for ii in np.arange(number_of_power_sweeps):
            Z_header=list(data)[1+3*ii]
            P_header=list(data)[2+3*ii]
            #Mid-value of linear power transmission
            lin_params[ii,5]=(lin_params[ii,0]*data[Z_header][int(np.floor(len(data[Z_header])/2))])+lin_params[ii,1]
            #data[P_header]=(data[P_header]+linear_diff(lin_params[ii,0],lin_params[ii,5],data[Z_header]))/lin_params[ii,1]
            data[P_header]=data[P_header]/linear_f(lin_params[ii,0],lin_params[ii,1],data[Z_header])
        
        #Linear approximation fitting is performed along P_in-axis, for every sample position z
        A=np.zeros(len(data))
        B=np.zeros(len(data))
        F=np.zeros(len(data))        
        P_in=lin_params[:,1]*np.exp(alpha*L_sample) #W #Estimation of input power within sample, from linear output power
        for jj in np.arange(len(data)):
            P_ratio=np.array(data.iloc[jj,2:3*number_of_power_sweeps:3])
            A[jj], B[jj], aux_r, aux_p, aux_std_err = stats.linregress(P_in, P_ratio)
        F=-2*A
        
        #Plots for input data survey, if required by user
        if(data_res_survey):
            fig12= bokeh_fig(plot_width=1000, plot_height=600)
            dataX=np.array(P_in)
            dataY=np.array(data.iloc[int(np.floor(len(data)/2)), 2:3*number_of_power_sweeps:3])
            A_plot, B_plot, aux_r, aux_p, aux_std_err = stats.linregress(P_in, np.array(data.iloc[int(np.floor(len(data)/2)), 2:3*number_of_power_sweeps:3]))
            dataY2=A_plot*dataX+B_plot
            fig12.scatter(dataX, dataY, color='#1374d1', size=10)
            fig12.line(dataX, [float(i) for i in dataY2], line_width=2, line_dash='dashed', color='#A80000')
            fig12.title = "Variation of the open aperture transmission as average input power is swept\nObtained from experimental open aperture normalized curves at z="+str(data.iloc[int(np.floor(len(data)/2)), 1])+".\nExperimental results are scattered dots, and the line corresponds to the linear regression used to extract the slope."
            fig12.xaxis.axis_label = "Average input power (W)"
            fig12.xaxis.axis_label_text_font_style='normal'
            fig12.yaxis.axis_label = "Normalized transmission (adim)"
            fig12.yaxis.axis_label_text_font_style='normal'
            bokeh_show(fig12)
        
        #The function F previously obtained is fitted as a lorentzian
        Z_positions=data.iloc[:,1] #m
        mod = LorentzianModel()
        pars = mod.guess(F, x=Z_positions)
        out = mod.fit(F, pars, x=Z_positions)
        
        #Plots for results survey, if required by user
        if(data_res_survey):
            fig7= bokeh_fig(plot_width=1000, plot_height=600)
            for ii in np.arange(number_of_power_sweeps):    
                dataX=np.array(Z_positions)
                dataY=np.array(F)
                dataY2=lorentzian(0.85*pars.get('amplitude').value, pars.get('center').value, pars.get('sigma').value, Z_positions)
                fig7.line(dataX, [float(i) for i in dataY], line_width=2, line_dash='solid', color='#1374d1')
                fig7.line(dataX, [float(i) for i in dataY2], line_width=2, line_dash='dashed', color='#A80000')
                fig7.title = "TPA fitting: Taylor expansion approximated method :: First linear slopes, then global lorentzian fitting. \nStandardized experimental data (solid) and fitting (dashed)."
                fig7.xaxis.axis_label = "Sample position (m)"
                fig7.xaxis.axis_label_text_font_style='normal'
                fig7.yaxis.axis_label = "Log. parameter value F (1 / W)"
                fig7.yaxis.axis_label_text_font_style='normal'
            bokeh_show(fig7)
        
        #TPA absorption is obtained from the amplitude of the lorentzian. The Rayleigh length and beam waist at focus are obtained from the lorentzian FWHM.
        beta_TPA=0.85*pars.get('amplitude').value*(lda_0*tau*f_rep)/(2*pi*L_eff)
        z_0=pars.get('sigma').value
        w_0=np.sqrt(lda_0*z_0/pi)
        
    else:
        #The open aperture measurements are fitted as lorentzians, for each power level
        Amp=np.zeros(number_of_power_sweeps)
        sig=np.zeros(number_of_power_sweeps)
        cent=np.zeros(number_of_power_sweeps)
        for ii in np.arange(number_of_power_sweeps):
            Z_header=list(data)[1+3*ii]
            P_header=list(data)[2+3*ii]
            #Mid-value of linear power transmission
            lin_params[ii,5]=(lin_params[ii,0]*data[Z_header][int(np.floor(len(data[Z_header])/2))])+lin_params[ii,1]
            #Data is standardized
            data[P_header]=2*(lin_params[ii,1]-(data[P_header]+linear_diff(lin_params[ii,0],0,data[Z_header])))/(np.exp(alpha*L_sample)*(lin_params[ii,1])**2)
            Z_positions=data[Z_header] #m
            mod = LorentzianModel()
            pars = mod.guess(data[P_header], x=Z_positions)
            out = mod.fit(data[P_header], pars, x=Z_positions)
            Amp[ii]=0.85*pars.get('amplitude').value
            sig[ii]=pars.get('sigma').value
            cent[ii]=pars.get('center').value
            
            #Plots for results survey, if required by user
            if(data_res_survey):
                fig7= bokeh_fig(plot_width=1000, plot_height=600)
                dataX=np.array(Z_positions)
                dataY=np.array(data[P_header])
                dataY2=lorentzian(0.85*pars.get('amplitude').value, pars.get('center').value, pars.get('sigma').value, Z_positions)
                fig7.line(dataX, [float(i) for i in dataY], line_width=2, line_dash='solid', color='#1374d1')
                fig7.line(dataX, [float(i) for i in dataY2], line_width=2, line_dash='dashed', color='#A80000')
                fig7.title = "TPA fitting: Taylor expansion approximated method :: First normalized lorentzian fitting, then average across power levels. \nStandardized experimental data (solid) and fitting (dashed).\n Incident power: "+str(lin_params[ii,1])+" W"
                fig7.xaxis.axis_label = "Sample position (m)"
                fig7.xaxis.axis_label_text_font_style='normal'
                fig7.yaxis.axis_label = "Log. parameter value F (1 / W)"
                fig7.yaxis.axis_label_text_font_style='normal'
                bokeh_show(fig7)
        
        #TPA absorption is obtained from the mean amplitude of all obtained lorentzian fittings. The Rayleigh length and beam waist at focus are obtained from the mean FWHM of all obtained lorentzian fittings.
        beta_TPA=np.mean(Amp)*(lda_0*tau*f_rep)/(2*pi*L_eff)
        z_0=np.mean(sig)
        w_0=np.sqrt(lda_0*z_0/pi)


# 🔎🕳KERR EFFECT PARAMETER EXTRACTION🔎🕳

L_obj = 30*(2*10*w_0) #m #Side length of the beam square cross-section used for discrete Fourier transform propagation.
dx = L_obj/N_obj #Coordinate system sampling minimum division

#If it is preferred by the user, the aperture radius is estimated by comparing the power measurements with open and closed aperture (with sample far from focus -linear regime)
estimate_aperture=bool(int(input('\nSelect method for aperture setting:\n(0) Use value in source code, provided manually.\n(1) Execute estimation from data.\n>>: ')))
if(estimate_aperture):
    #Ratio between open and closed aperture transmission in linear regime is calculated
    extinction_objective_value=np.mean(np.array(lin_params2[:,1])/np.array(lin_params[:,1]))
    
    #An iterative optimization process is carried out, to obtaine an estimation for the aperture that corresponds to the OA/CA power ratio experimentally obtained.
    dx = L_obj/N_obj #Coordinate system sampling
    x = np.arange(-N_obj*dx/2,N_obj*dx/2,dx)    # Coordinate space construction
    X, Y = np.meshgrid(x,x)   # 2-dimensional transverse coordinate space  
    max_guess_val=float(input('Provide a maximun value for the estimation (in m):\n>>: '))
    lower_lim_guess=0
    upper_lim_guess=max_guess_val
    iteration=0
    while(iteration<=10):
        iteration+=1
        err_iter_list=[]
        for a_guess in np.linspace(lower_lim_guess,upper_lim_guess,10):
            val=detect(L_obj, N_obj, n_prop_media, iris(L_obj, N_obj, beam(X,Y, lda_0, 1, w_0, z_0, n_prop_media, L_setup-f_lens, 0),a_guess))/detect(L_obj, N_obj, n_prop_media, beam(X,Y, lda_0, 1, w_0, z_0, n_prop_media, L_setup-f_lens, 0))
            err_iter_list.append(error_percentage(val,extinction_objective_value))
        step_guess=(upper_lim_guess-lower_lim_guess)/9
        best_guess=np.linspace(lower_lim_guess,upper_lim_guess,10)[err_iter_list.index(min(err_iter_list))]
        lower_lim_guess=best_guess-step_guess
        upper_lim_guess=best_guess+step_guess
        if(min(err_iter_list)<=0.01):
            break
    aa=best_guess
else:
    aa=a

#Graphics for results survey, if required by user
if(data_res_survey):
    x3 = np.arange(-N_obj*dx/2,N_obj*dx/2,dx)   #Coordinate space construction
    X3, Y3 = np.meshgrid(x3,x3)  #2-dimensional transverse coordinate space
    beam_to_plot=abs(beam(X3,Y3, lda_0, 1, w_0, z_0, n_prop_media, L_setup-f_lens, 0))**2
    fig8= bokeh_fig(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")])
    fig8.x_range.range_padding = fig8.y_range.range_padding = 0
    fig8.image(image=[beam_to_plot], x=-L_obj/2, y=-L_obj/2, dw=L_obj, dh=L_obj, palette=bokehpalette, level="image")
    fig8.grid.grid_line_width = 0.5
    fig8.title = "Calculated laser beam transverse profile at aperture plane.\nOpen aperture. Without sample exposition."
    fig8.xaxis.axis_label = "x (m)"
    fig8.xaxis.axis_label_text_font_style='normal'
    fig8.yaxis.axis_label = "y (m)"
    fig8.yaxis.axis_label_text_font_style='normal'
    fig8.xgrid.grid_line_color = None
    fig8.ygrid.grid_line_color = None
    bokeh_show(fig8)
    
    
    beam_to_plot=abs(iris(L_obj, N_obj, beam(X3,Y3, lda_0, 1, w_0, z_0, n_prop_media, L_setup-f_lens, 0),aa))**2
    fig8= bokeh_fig(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")])
    fig8.x_range.range_padding = fig8.y_range.range_padding = 0
    fig8.image(image=[beam_to_plot], x=-L_obj/2, y=-L_obj/2, dw=L_obj, dh=L_obj, palette=bokehpalette, level="image")
    fig8.grid.grid_line_width = 0.5
    fig8.title = "Calculated laser beam transverse profile at aperture plane.\nClosed aperture. Without sample exposition."
    fig8.xaxis.axis_label = "x (m)"
    fig8.xaxis.axis_label_text_font_style='normal'
    fig8.yaxis.axis_label = "y (m)"
    fig8.yaxis.axis_label_text_font_style='normal'
    fig8.xgrid.grid_line_color = None
    fig8.ygrid.grid_line_color = None
    bokeh_show(fig8)


#Calculation of the slope of the linear function that describes the relation between peak-valley differences of renormalized closed aperture transmissions and the corresponding incident power levels. This is performed with experimental data, and the slope is used as objective value for further numerical optimization.
peak_valley_diff=[]
for ii in np.arange(number_of_power_sweeps):
    peak_valley_diff.append(T_diff(np.array(data.iloc[:,3+3*ii])))
lin_regress = LinearRegression(fit_intercept=False).fit((np.array(lin_params[:,1])/(tau*f_rep)).reshape((-1, 1)),np.array(peak_valley_diff))
pvd_slope_objective=lin_regress.coef_[0]

#Plots for input data survey, if required by user
if(data_res_survey):
    fig11= bokeh_fig(plot_width=1000, plot_height=600)
    dataX=np.array(lin_params[:,1])/(tau*f_rep)
    dataY=np.array(peak_valley_diff)
    dataY2=pvd_slope_objective*dataX
    fig11.scatter(dataX, dataY, color='#64b066', size=10)
    fig11.line(dataX, [float(i) for i in dataY2], line_width=2, line_dash='dashed', color='#6d068f')
    fig11.title = "Variation of the peak-valley transmission difference as average input power is swept\nObtained from experimental closed aperture RENORMALIZED curves.\nExperimental results are scattered dots, and the line corresponds to the linear regression used to extract the slope."
    fig11.xaxis.axis_label = "Peak input power (W)"
    fig11.xaxis.axis_label_text_font_style='normal'
    fig11.yaxis.axis_label = "Peak-valley transmission difference (adim)"
    fig11.yaxis.axis_label_text_font_style='normal'
    bokeh_show(fig11)

#The n_2 search is initialized by providing a range to initialize the search
if(bool(int(input("\nSelect the method for n_2 search initialization:\n(0) Use search interval limits in source code, provided manually.\n(1) Enter search interval limits via console.\n>>: ")))):
    n_2_search=[]
    if(negative_n_2):
        print("\n**Search interval limits must be negative\n It is possible to use scientific notation;\n e.g., enter -1.5e-20, which stands for -1.5 × 10^(-20)")
    else:
        print("\n**Search interval limits must be positive.\n It is possible to use scientific notation;\n e.g., enter 1.5e-20, which stands for 1.5 × 10^(-20)")
    n_2_search.append(float(input("Enter lower limit for n_2 search interval (must be less than next one):\n>>: ")))
    n_2_search.append(float(input("Enter upper limit for n_2 search interval (must be greater than previous one):\n>>: ")))

#The n_2 estimation is obtained by performing multiple numerical simulations of the Z-scan, with different n_2 values; until obtaining the expected relation between peak-valley transmission difference and incident power.
#Optimization process: first stage
preliminary_calculated_slopes=[]
for N_2 in np.linspace(n_2_search[0],n_2_search[1],3):
    preliminary_calculated_slopes.append(Tpv_P_slope(N_2, beta_TPA, alpha, L_sample, L_eff, P_peak_max, lda_0, w_0, z_0, f_lens, L_setup, aa, n_prop_media, L_obj, N_obj,  t))
lin_regress2 = LinearRegression(fit_intercept=False).fit((np.linspace(n_2_search[0],n_2_search[1],3)).reshape((-1, 1)), np.array(preliminary_calculated_slopes))
slope_vs_n2_relation=lin_regress2.coef_[0]
n_2_first_prediction=(1/slope_vs_n2_relation)*pvd_slope_objective #Preliminary guess for the n_2 estimation

#Optimization process: second stage
if(negative_n_2):
    lower_lim_search=10*n_2_first_prediction
    upper_lim_search=0.1*n_2_first_prediction
else:
    lower_lim_search=0.1*n_2_first_prediction
    upper_lim_search=10*n_2_first_prediction

iteration=0
while(iteration<=10):
    iteration+=1
    err_iter_list=[]
    for N_2 in np.linspace(lower_lim_search, upper_lim_search, 5):
        val=Tpv_P_slope(N_2, beta_TPA, alpha, L_sample, L_eff, P_peak_max, lda_0, w_0, z_0, f_lens, L_setup, aa, n_prop_media, L_obj, N_obj,  t)
        err_iter_list.append(error_percentage(val,pvd_slope_objective))
    step_search=abs(upper_lim_search-lower_lim_search)/4 #Set integer accordingly
    best_val=np.linspace(lower_lim_search, upper_lim_search, 5)[err_iter_list.index(min(err_iter_list))] #Set integer accordingly
    
    if(negative_n_2):
        lower_lim_search=best_val-step_search
        if((best_val+step_search)<0):
            upper_lim_search=best_val+step_search
        else:
            upper_lim_search=0
    else:
        upper_lim_search=best_val+step_search
        if((best_val-step_search)>0):
            lower_lim_search=best_val-step_search
        else:
            lower_lim_search=0
    
    if(min(err_iter_list)<=0.01): #Error percentage tolerance threshold is set to 1%
        break
    
    if(iteration==11): #If the last iteration is reached, the estimation process finishes and a warning message is printed; since the "Error percentage tolerance threshold" condition has not been fulfilled.
        print("\n**Maximun number of iterations for n_2 optimization was reached: iter."+str(iteration)+"/"+str(iteration)+"**")
        print("Best estimation attained for n_2 might exceed the error limit.")
        print("Error percentage associated to the optim. process: "+str(100*min(err_iter_list))+"%")
    
n_2=best_val #n_2 estimation is determined as the value that produced the most accurate relation between Tpv and P_in

#Graphics and plots for results survey, if required by user
if(data_res_survey):
    x3 = np.arange(-N_obj*dx/2,N_obj*dx/2,dx)   #Coordinate space construction
    X3, Y3 = np.meshgrid(x3,x3)  #2-dimensional transverse coordinate space
    BEAM_sampled_check=beam_NL(L_obj, N_obj, lda_0, P_peak_max, w_0, z_0, n_prop_media, -1.7*z_0, t, L_sample, alpha, L_eff, n_2, beta_TPA) 
    beam_to_plot=abs(propaFT(BEAM_sampled_check,lda_0,(L_setup-f_lens)-(-1.7*z_0),dx))**2
    fig9= bokeh_fig(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")])
    fig9.x_range.range_padding = fig8.y_range.range_padding = 0
    fig9.image(image=[beam_to_plot], x=-L_obj/2, y=-L_obj/2, dw=L_obj, dh=L_obj, palette=bokehpalette, level="image")
    fig9.grid.grid_line_width = 0.5
    fig9.title = "Calculated laser beam transverse profile at aperture plane.\nOpen aperture. With sample exposition at z = -1.7 z_0."
    fig9.xaxis.axis_label = "x (m)"
    fig9.xaxis.axis_label_text_font_style='normal'
    fig9.yaxis.axis_label = "y (m)"
    fig9.yaxis.axis_label_text_font_style='normal'
    fig9.xgrid.grid_line_color = None
    fig9.ygrid.grid_line_color = None
    bokeh_show(fig9)
    
    BEAM_sampled_check=beam_NL(L_obj, N_obj, lda_0, P_peak_max, w_0, z_0, n_prop_media, 1e-25*z_0, t, L_sample, alpha, L_eff, n_2, beta_TPA) 
    beam_to_plot=abs(propaFT(BEAM_sampled_check,lda_0,(L_setup-f_lens)-(1e-25*z_0),dx))**2
    fig9= bokeh_fig(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")])
    fig9.x_range.range_padding = fig8.y_range.range_padding = 0
    fig9.image(image=[beam_to_plot], x=-L_obj/2, y=-L_obj/2, dw=L_obj, dh=L_obj, palette=bokehpalette, level="image")
    fig9.grid.grid_line_width = 0.5
    fig9.title = "Calculated laser beam transverse profile at aperture plane.\nOpen aperture. With sample exposition at z = 0."
    fig9.xaxis.axis_label = "x (m)"
    fig9.xaxis.axis_label_text_font_style='normal'
    fig9.yaxis.axis_label = "y (m)"
    fig9.yaxis.axis_label_text_font_style='normal'
    fig9.xgrid.grid_line_color = None
    fig9.ygrid.grid_line_color = None
    bokeh_show(fig9)
    
    BEAM_sampled_check=beam_NL(L_obj, N_obj, lda_0, P_peak_max, w_0, z_0, n_prop_media, 1.7*z_0, t, L_sample, alpha, L_eff, n_2, beta_TPA) 
    beam_to_plot=abs(propaFT(BEAM_sampled_check,lda_0,(L_setup-f_lens)-(1.7*z_0),dx))**2
    fig9= bokeh_fig(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")])
    fig9.x_range.range_padding = fig8.y_range.range_padding = 0
    fig9.image(image=[beam_to_plot], x=-L_obj/2, y=-L_obj/2, dw=L_obj, dh=L_obj, palette=bokehpalette, level="image")
    fig9.grid.grid_line_width = 0.5
    fig9.title = "Calculated laser beam transverse profile at aperture plane.\nOpen aperture. With sample exposition at z = 1.7 z_0."
    fig9.xaxis.axis_label = "x (m)"
    fig9.xaxis.axis_label_text_font_style='normal'
    fig9.yaxis.axis_label = "y (m)"
    fig9.yaxis.axis_label_text_font_style='normal'
    fig9.xgrid.grid_line_color = None
    fig9.ygrid.grid_line_color = None
    bokeh_show(fig9)
    
    prop_cycle = plt.rcParams['axes.prop_cycle']
    color_set=prop_cycle.by_key()['color']
    color_index=0
    fig10= bokeh_fig(plot_width=1000, plot_height=600)
    for Power_level_check in np.linspace(0.001*P_peak_max, P_peak_max, 5):
        PD_oa=[]
        PD_ca=[]
        z_sample_pos_check=np.linspace(-5.01*z_0, 5*z_0, 80) #z coordinates, for different sample positions
        for Z in z_sample_pos_check:            
            BEAM_sampled=beam_NL(L_obj, N_obj, lda_0, Power_level_check, w_0, z_0, n_prop_media, Z, t, L_sample, alpha, L_eff, n_2, beta_TPA)
            BEAM_out=propaFT(BEAM_sampled,lda_0,(L_setup-f_lens)-Z,dx)    
            BEAM_irised=iris(L_obj, N_obj, BEAM_out, aa)    
            PD_oa.append(detect(L_obj, N_obj, n_prop_media, BEAM_out)) #Open aperture z-scan
            PD_ca.append(detect(L_obj, N_obj, n_prop_media, BEAM_irised)) #Closed aperture z-scan
        BEAM_sampled=beam_NL(L_obj, N_obj, lda_0, Power_level_check, w_0, z_0, n_prop_media, 0.99*(L_setup-f_lens), t, L_sample, alpha, L_eff, n_2, beta_TPA)
        BEAM_out=propaFT(BEAM_sampled,lda_0,0.01*(L_setup-f_lens),dx)    
        BEAM_irised=iris(L_obj, N_obj, BEAM_out, aa)    
        PD_REF_oa=detect(L_obj, N_obj, n_prop_media, BEAM_out) #Open aperture reference measurement (in linear regime)
        PD_REF_ca=detect(L_obj, N_obj, n_prop_media, BEAM_irised) #Closed aperture reference measurement (in linear regime)
        PD_oa=PD_oa/PD_REF_oa #Open aperture normalized z-scan
        PD_ca=PD_ca/PD_REF_ca #Closed aperture normalized z-scan
        PD_renorm=np.divide(np.array(PD_ca), np.array(PD_oa)) #Closed aperture renormalized z-scan
        dataX=np.array(z_sample_pos_check)
        dataY=np.array(PD_renorm)
        fig10.line(dataX, [float(i) for i in dataY], line_width=2, line_dash='solid', color=color_set[color_index])
        color_index+=1
    fig10.title = "Simulated renormalized closed aperture z-scan.\nFor n_2="+str(n_2)+" with incident peak power swept from 0 to "+str(P_peak_max)+" W"
    fig10.xaxis.axis_label = "Sample position (m)"
    fig10.xaxis.axis_label_text_font_style='normal'
    fig10.yaxis.axis_label = "Renormalized transmittance (adim)"
    fig10.yaxis.axis_label_text_font_style='normal'
    bokeh_show(fig10)

# 📰💾RESULTS PRESENTATION AND STORING📰💾

if(save_to_file):
    date=datetime.datetime.now()
    date_now=str(date.year)+'-'+str(date.month)+'-'+str(date.day)+'; '+str(date.hour)+'-'+str(date.minute)+'-'+str(date.second)
    if (not os.path.exists(date_now)):
        os.makedirs(date_now)
    params_file=open(date_now+'/'+str(data_filename[0:-4])+'--Nonlinear_extracted_parameters.csv','w')
    if(estimate_aperture):
        params_file.write('r_a, '+str(aa)+',')
        params_file.write('\n')
    params_file.write('w_0, '+str(w_0)+',')
    params_file.write('\n')
    params_file.write('n_2, '+str(n_2)+',')
    params_file.write('\n')
    params_file.write('beta_TPA, '+str(beta_TPA)+',')
    params_file.write('\n')
    params_file.write('FOM_TPA, '+str(n_2/(lda_0*beta_TPA))+',')
    params_file.write('\n')
    params_file.close()

print("")
print("╓─── ▾▾▾")
print("║ Parameters extracted:")
if(estimate_aperture):
    print('║ r_a: '+str(aa)+' m')
print('║ w_0: '+str(w_0)+' m')
print("║ n_2: "+str(n_2)+" m^2/W")
print('║ beta_TPA: '+str(beta_TPA)+' m/W')
print("║ FOM_TPA: "+str(n_2/(lda_0*beta_TPA)))
print("╙─── ▴▴▴")


"""
References:
    [1] M. Sheik-Bahae, A.A. Said, T.H. Wei, D.J. Hagan, and E.W. Van Stryland. Sensitive Measurement of Optical Nonlinearities Using a Single Beam. IEEE Journal of Quantum Electronics, 26(4). 1990.
    [2] S. Serna. Design and characterization of Silicon Photonic structures for third order nonlinear effects. Ph.D. thesis, Université Paris-Saclay. 2016.
    [3] E. Rueda. Simulación de propagación de la luz en el espacio libre usando la transformada de Fourier Rápida: propagación básica. Video tutorial and open-source python notebook. 2020.
    [4] E. Rueda. Simulación de propagación de la luz en el espacio libre usando la transformada de Fourier Rápida: transformada de Fourier Óptica. Video tutorial and open-source python notebook. 2020.
"""
