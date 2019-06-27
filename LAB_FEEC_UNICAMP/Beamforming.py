import cmath
import math
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
from sklearn.linear_model import Ridge

from .CustomPlots import plotFreq, plotTime
from .Signal import Signal
from .Utils import toFreq, toTime

#Inversion of Matrix using Moore-Penrose
def inv(x):
    return sp.linalg.pinv(x)

#Calculates the weights per frequency for a time delay
def calcWeightDelay(fs,lenSpec,s):
    samples = fs*s
    ret = np.array([cmath.exp(-complex(0,math.pi*m*samples/lenSpec)) for m in range(0,lenSpec)])
    return ret

#Calculates the array manifold vector
#For both ways: from source to sample and sample to source
def calcManifoldVector(arr, src_to_smpl):
    lenSpec = len(arr.Microphones[0].signal.Spectrum)
    cnt_mics = len(arr.Microphones)
    pos_last = arr.Microphones[-1].dist_source
    pos_first = arr.Microphones[0].dist_source
    v = []

    for m in range(0,len(arr.Microphones)):
        if (src_to_smpl):
            v.append(calcWeightDelay(arr.fs,lenSpec,-(pos_first - arr.Microphones[m].dist_source)/arr.c))
        else:
            v.append(calcWeightDelay(arr.fs,lenSpec,(pos_first - arr.Microphones[m].dist_source)/arr.c))
        
    return v

#Calculates and applies DAS weights and return the result
def applyDAS(arr, src_to_smpl):
    v = calcManifoldVector(arr, src_to_smpl)
    
    for i in range(0,len(arr.Microphones)):
        arr.Microphones[i].signal.weight(np.conj(v[i]))

    DAS_Spectrum = np.sum([m.signal.weighted() for m in arr.Microphones],axis=0)/len(arr.Microphones)

    return DAS_Spectrum

#Calculates and applies MVDR weights and return the result
#Given the noise covariance matrix
def applyMVDR(arr , src_to_smpl :bool , Rxx_inv):
    v_mvdr = np.matrix(calcManifoldVector(arr,src_to_smpl), dtype=np.complex_)
    w_mvdr = np.zeros(v_mvdr.shape, dtype=np.complex_)

    for i in range(0,Rxx_inv.shape[-1]):
        v_freq = v_mvdr[:,i]

        Rxx_inv_freq = np.matrix(Rxx_inv[...,i])

        num = v_freq.H*Rxx_inv_freq
        dem = v_freq.H*Rxx_inv_freq*v_freq
        r = num/(dem[0] + 1e-50)

        w_mvdr[:,i] = r
        
    for i in range(0,len(arr.Microphones)):
        arr.Microphones[i].signal.weight(w_mvdr[i])

    MVDR_Spectrum = np.sum([m.signal.weighted() for m in arr.Microphones],axis=0)

    return MVDR_Spectrum

#Invert the Rxx per frequency
#If diagonalLoad = 0 then uses a adaptative diagonalLoading
def invertRxx(Rxx, diagonalLoad):
    Rxx_inv = np.zeros((len(Rxx),len(Rxx),Rxx.shape[-1]), dtype=np.complex_)

    for i in range(0,Rxx.shape[-1]):
        Rxx_freq = np.matrix(Rxx[...,i], dtype=np.complex_)

        Rxx_real = np.real(Rxx_freq)/(np.real(Rxx_freq).max() + 1e-50)
        if diagonalLoad==0:
            dl = np.std(np.diag(Rxx_real))
        else:
            dl = diagonalLoad
        Rxx_real = Rxx_real + dl*np.eye(len(Rxx), dtype=np.complex_)

        Rxx_imag = np.imag(Rxx_freq)/(np.imag(Rxx_freq).max() + 1e-50)
        if diagonalLoad==0:
            dl = np.std(np.diag(Rxx_imag))
        else:
            dl = diagonalLoad
        Rxx_imag = Rxx_imag + dl*np.eye(len(Rxx), dtype=np.complex_)

        #Lemma: inv(A + U*C*V) = inv(A) - inv(A)*U*inv(inv(C) + V*inv(A)*U)*V*inv(A)
        A = Rxx_real
        U = np.matrix(np.eye(len(Rxx)), dtype=np.complex_)
        C = 1j * Rxx_imag
        V = np.matrix(np.eye(len(Rxx)), dtype=np.complex_)
        Rxx_freq_inv = inv(A) - inv(A)*U*inv(inv(C) + V*inv(A)*U)*V*inv(A)

        Rxx_inv[...,i] = Rxx_freq_inv
    
    return Rxx_inv

#Calculate the covariance by multiplying its spectrum
def DirectCovariance(arr, Sn):
    lFreq = Sn.shape[1]
    Rxx = np.zeros((len(Sn),len(Sn),lFreq), dtype=np.complex_)

    for i in range(0,lFreq):
        Rxx[...,i] = np.matrix(Sn[:,i])*np.matrix(Sn[:,i]).H

    return  Rxx

#Calculate the covariance by diving the 
def IterativeCovariance(arr, Sn, parts):
    #If no parts, calculate directly
    if parts==0:
        return DirectCovariance(arr, Sn)

    lFreq = Sn.shape[1]
    Rxx = np.zeros((len(Sn),len(Sn),lFreq))

    #Defines a min step to get better results
    min_step = math.trunc(((arr.Microphones[-1].dist_source - arr.Microphones[0].dist_source)/arr.c)*arr.fs*3.5)
    
    sigs_time = np.array([toTime(np.array(s)[0]) for s in Sn])
    lTime = len(sigs_time[0])

    if parts == 0:
        parts = math.trunc(lTime/min_step)

    step = math.trunc(lTime/parts)
    sigs_time = np.array([np.resize(s,len(s)+2*step) for s in sigs_time])
    sigs_time[:,0:2*step] = 0
    sigs_time = np.array([np.roll(s,-step) for s in sigs_time])

    if step<min_step:
        raise Exception("Min step not satisfied!")

    Rxx = np.zeros((len(Sn),len(Sn),lFreq), dtype=np.complex_)

    rStep = []
    for i in range(0,parts+1):
        pos_forward = (i+1)*step #Current position for forward stepping
        pos_backward = len(sigs_time[0]) - pos_forward #Current position for backward stepping

        tFor_pos = sigs_time[...,pos_forward-step:pos_forward+step]*np.hanning(2*step)
        fFor_pos = np.matrix([toFreq(t,lTime) for t in tFor_pos], dtype=np.complex_)
        RFor_pos = np.array([fFor_pos[...,j]*fFor_pos[...,j].H/(2*fFor_pos.shape[0]*(parts+1)) for j in range(0,fFor_pos.shape[1])], dtype=np.complex_)
        rStep.append(RFor_pos.T)

        tBack_pos = sigs_time[...,pos_backward-step:pos_backward+step]*np.hanning(2*step)
        fBack_pos = np.matrix([toFreq(t,lTime) for t in tBack_pos], dtype=np.complex_)
        RBack_pos = np.array([fBack_pos[...,j]*fBack_pos[...,j].H/(2*fBack_pos.shape[0]*(parts+1)) for j in range(0,fBack_pos.shape[1])], dtype=np.complex_)
        rStep.append(RBack_pos.T)
    
    Rxx[...,:] = np.sum(rStep,axis=0)
    return Rxx
