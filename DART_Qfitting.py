# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 17:13:15 2019

@author: teb3
"""

# DART Q and f fitting

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


csv = pd.read_csv('trackx2_DART_tracking.csv')
csv.columns = (['time','Amp','Amp2','Phase','Phase2','Freq','Freq0','Q'])
Freq2=np.array(csv.Freq+60000) #60 kHz Dual AC mode frequency width
csv.Phase = csv.Phase/180*np.pi # degress to rad
csv.Phase2 = (csv.Phase2)/180*np.pi

# Initialize arrays
Omega = np.ndarray(len(csv))
Phi = np.ndarray(len(csv))
X1 = np.ndarray(len(csv))
X2 = np.ndarray(len(csv))
f0 = np.ndarray(len(csv))
Q = np.ndarray(len(csv))
Adrive = np.ndarray(len(csv))
PhaseDrive = np.ndarray(len(csv))

# Calculate driven damped harmonic oscillator parameters (Gannepalli et al., 2011)
for t in range(len(csv)):
    Omega[t] = csv.Freq[t]*csv.Amp[t]/(Freq2[t]*csv.Amp2[t])
    Phi[t] = np.tan(csv.Phase2[t]-csv.Phase[t])
    X1[t] = -(1-np.sign(Phi[t])*np.sqrt(1+Phi[t]**2)/Omega[t])/Phi[t] # Equations from the corrigendum
    X2[t] = (1-np.sign(Phi[t])*Omega[t]*np.sqrt(1+Phi[t]**2))/Phi[t]
    
    f0[t] = np.sqrt(csv.Freq[t]*Freq2[t]) # np.sqrt(csv.Freq[t]*Freq2[t]*(Freq2[t]*X1[t]-csv.Freq[t]*X2[t])/(csv.Freq[t]*X1[t]-Freq2[t]*X2[t]))
    Q[t] = np.sqrt(csv.Freq[t]*Freq2[t]*(Freq2[t]*X1[t]-csv.Freq[t]*X2[t])*(csv.Freq[t]*X1[t]-Freq2[t]*X2[t]))/(Freq2[t]**2-csv.Freq[t]**2)
    
    Adrive[t] = csv.Amp[t]*np.sqrt((f0[t]**2-csv.Freq[t]**2)**2+(f0[t]*csv.Freq[t]/Q[t]))/f0[t]**2
    PhaseDrive[t] = csv.Phase[t]-np.arctan(f0[t]*csv.Freq[t]/(Q[t]*(f0[t]**2-csv.Freq[t]**2)))


csv2 = pd.read_csv('trackx2_before_after.csv')
csv2.columns = (['Freq','AmpAfter','Freq','AmpBefore'])
plt.plot(csv2.Freq,csv2.AmpBefore,label='before') 
plt.plot(csv2.Freq,csv2.AmpAfter,label='after')
for t2 in range(0,int(len(csv)),20000):
    f = np.linspace(1e6,1.5e6,100)
    A = np.ndarray(100)
    A = f0[t2]**2*Adrive[t2]/np.sqrt((f0[t2]**2-f**2)**2+(f0[t2]*f/csv.Q[t2])**2)
    plt.semilogy(f,A)
#plt.legend()
plt.xlim(7e5,1.8e6)
plt.ylim(.0005,.3)
plt.ylabel('Amp [V]')
plt.xlabel('Frequency [Hz]')
plt.savefig('./figures/trackx2.png',dpi=300)
plt.show() 

    
plt.plot(csv.time,csv.Amp, csv.time, csv.Amp2)
plt.legend(['Amp1','Amp2'])
plt.xlabel('Time [s]')
plt.ylabel('Amp [V]')
plt.savefig('./figures/trackx2_Amp.png',dpi=300)
plt.show()

plt.plot(csv.time, csv.Phase,csv.time, csv.Phase2)
plt.legend(['Phase1','Phase2'])
plt.xlabel('Time [s]')
plt.ylabel('Phase')
plt.savefig('./figures/trackx2_Phase.png',dpi=300)
plt.show()

plt.plot(csv.time,csv.Freq0, csv.time, f0, ':', csv.time, csv.Freq, csv.time, Freq2)
plt.ylim(1200000, 1300000)
plt.legend(['Freq0_Igor','Freq0_Calc','f1','f2'])
plt.xlabel('Time [s]')
plt.ylabel('Frequency [Hz]')    
plt.savefig('./figures/trackx2_Frequency.png',dpi=300)
plt.show()

plt.plot(csv.time, csv.Q,csv.time, Q, ':')
plt.legend(['Q_Igor','Q_Calc'])
plt.xlabel('Time [s]')
plt.ylabel('Q')
plt.savefig('./figures/trackx2_Q.png',dpi=300)
plt.show()