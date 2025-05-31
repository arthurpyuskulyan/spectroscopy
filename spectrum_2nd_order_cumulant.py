import numpy as np
import matplotlib.pyplot as plt
vees= np.loadtxt('testvee_traj1.dat')
time= np.loadtxt('testtimefs.txt')
time=time*41.341374575751 #fs to au
energy=vees[:,0]
t=time[1:]   #time starts at 0. This is taking the positive part
t2=-t #to get negative time
t3=t2[::-1] #reverse to put in order
t3 = np.append(t3, time)  #attach times to have negative and positive time
#2
#convert vee from eV to au
energy=energy/27.211324570273
#3 calculate avg vee and get vee fluctuations
mean=np.mean(energy)
fluctuation=energy-mean
len(fluctuation)
acf=np.correlate(fluctuation,fluctuation,mode='full')/len(fluctuation) # compute correlation function with np.correlate
print(acf)
#damping function
damp = np.exp(-np.abs(t3)/12402.412372725299)   #damp=-exp(|t|/Tau) 12402.412372725299 is for Tau=300fs in au
plt.plot(t3,damp)
plt.xlim(-25000,25000)
plt.title('300fs damping function')
plt.show()
plt.plot(t3,acf,label='C(t)')
plt.plot(t3,damp*acf, label='300fs damped C(t)')
plt.xlabel('time')
plt.ylabel('C(t)')
plt.title('Correlation function')
plt.legend()
plt.xlim(-20000,20000)
plt.show()
#define beta Boltzmann constant
k = 1.380649e-23  #in Joules
T = 300  #kelvin
b = 1/(k*T)
b = b*4.359744e-18 #convert to atomic unit
w = np.arange( -000,4000, 1) #in wavenumber define omega grid
w=w*0.0000046 #in au
def fourier_transform_trapezoid(t3, w, time_func):
    # calculate time step
    dt = t3[1] - t3[0]
    # initialize the result array for the Fourier transform
    fw = np.zeros(len(w), dtype=complex)
    # loop over the frequencies
    for i in range(len(w)):
        # Define the integrand as e^(i * w * t) multiplied by the time function
        f = np.exp(complex(0, 1) * w[i] * t3) * time_func
        # trapezoidal integration
        fw[i] = np.trapz(f, dx=dt)
    return fw
integrand = acf*damp
fw_with_time_func = fourier_transform_trapezoid(t3, w, integrand)
print(fw_with_time_func)
Jwnew=(b*w/2)*fw_with_time_func
print(Jwnew)
plt.plot(w*27.211396*8065.5,Jwnew*27.211396) #divide by 8065.56 to convert to wavenumber, multiply by 27.2113 to go to eV
plt.xlabel('Spectral density')
plt.ylabel('Wavenumber (cm$^-$$^1$)')
plt.legend()
plt.xlim(-1000,2000)
plt.show()
dw=w[1]-w[0]
# Define new time grid from 0 to 600 fs with 0.3 fs intervals
time_fs = np.arange( 0, 600, .03)
time = time_fs * 41.341374575751  # Convert to atomic units (au)

# Calculate g(t) using the new time grid
gt_real = np.zeros(len(time))
gt_imag = np.zeros(len(time))
for t_idx in range(len(time)):
    t_val = time[t_idx]
    real = Jwnew * (1/np.pi) * (1 - np.cos(t_val * w)) / (w**2 * np.tanh(b * w / 2))
    real[0] = 0  #Put analytical time 0 expression
    imag = Jwnew * (1/np.pi) * (np.sin(t_val * w) - w * t_val) / (w**2)
    imag[0] = 0
    # Use np.trapz to integrate the real and imaginary parts
    gt_real[t_idx] = np.trapz(real, dx=dw)
    gt_imag[t_idx] = np.trapz(imag, dx=dw)
# Combine real and imaginary parts into g(t)
gt = gt_real + complex(0, 1) * gt_imag
response = np.zeros(len(time), dtype=complex)
#response=np.exp(complex(0,1)*mean*time)*np.exp(-gt)
# Plot the real and imaginary parts of g(t)
plt.plot(time, gt_real, label='Real part')
plt.plot(time, gt_imag, label='Imaginary part')
plt.ylabel('g$_2$(t)')
plt.xlabel('Time (fs)')
plt.legend()
plt.title('g(t) Real and Imaginary Parts')
plt.show()
response=np.exp(-complex(0,1)*mean*time)*np.exp(-gt)
#forpositive=np.exp(-complex(0,1)*mean*time)*np.exp(-gt)
#responsecc=np.exp(complex(0,1)*mean*time)*np.exp(-gt)
#response = np.where(time < 0, responsecc, responseforpositive)
w1=np.arange(2.4,3.6,0.01)/27.2114
#trapezoidial absorbance
dt=time[1]-time[0]
f2=np.zeros(len(time))
absw=np.zeros(len(w1))
for i in range(len(w1)):
    f2=response*np.exp(complex(0, 1)*w1[i]*time)
    absw[i] = np.trapz(f2, dx=dt)
abswnew=((4*np.pi)/3)*w1*absw
print(mean)
#plt.plot(xsdh,Ih/np.max(Ih),label='Molspeck')
plt.plot(w1*27.2114,abswnew/np.max(abswnew),label='300 fs') #divide by 0.0000046 to wavenumber multiply by 27.2113 to go to eV
plt.xlabel('Energy (eV)')
plt.ylabel('Intensity norm')
plt.legend()
plt.show()