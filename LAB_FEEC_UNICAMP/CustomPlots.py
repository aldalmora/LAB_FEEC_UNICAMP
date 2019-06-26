import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from .Utils import *

#Plot all the microphones in time
def Plot_Mics_Time(arr):
    colors = ['r','g','b','c','m','k','','','','','']
    for i in range(0 ,len(arr.Microphones)):
        mic_Time = arr.Microphones[i].signal.time()
        plt.plot(np.linspace(0,len(mic_Time)/arr.fs,len(mic_Time)),mic_Time,color=colors[i],linewidth=0.25)

#Given the Spectrum, plot the signal in time
def plotTime(fs,spectrum,cor='',linewidth=0.25):
    t2 = toTime(spectrum)
    plt.plot(np.linspace(0,len(t2)/fs,len(t2)),t2,cor,linewidth=linewidth)
    plt.xlabel('Tempo(s)')
    plt.ylabel('Amplitude')
    
#Plot the Spectrum Magnitude
def plotFreq(fs,spectrum,cor='',linewidth=0.25):
    spc = np.absolute(spectrum)
    plt.semilogx(np.linspace(0,fs/2 -1,len(spc)),20 * np.log10(spc),cor,linewidth=linewidth)
    plt.xlabel('Frequência(Hz)')
    plt.ylabel('Amplitude')

#Plot the Spectrum Phase
def plotPhase(fs,spectrum,cor='',linewidth=0.25):
    ang = np.angle(spectrum)*180/np.pi
    plt.semilogx(np.linspace(0,fs/2 -1,len(ang)),ang,cor,linewidth=linewidth)
