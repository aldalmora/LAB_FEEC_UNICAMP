import glob
import pickle as pk

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage.filters import gaussian_filter1d
from scipy.signal import savgol_filter

from .Beamforming import (DirectCovariance, IterativeCovariance,
                          applyDAS, applyMVDR, calcManifoldVector, invertRxx)
from .CustomPlots import Plot_Mics_Time, plotFreq, plotTime
from .Utils import DeslocSig

#Test the parameters in MVDR beamformer and stores in a pickle file
def Test(file_name, Original, arr, partsAC, diagonalLoading):
    DAS_Inc = applyDAS(arr, True)

    DAS_Refl = applyDAS(arr, False)
    Dls = DeslocSig(DAS_Refl,calcManifoldVector(arr, False))
    Sn = np.matrix(Dls)

    Rxx_Sn_Inc = IterativeCovariance(arr, Sn, partsAC)
    Conds_Rxx_Inc = np.array([np.linalg.cond(r) for r in np.array(Rxx_Sn_Inc[...,:].T)], dtype=np.complex_)

    Rxx_Sn_Inv_Inc = invertRxx(Rxx_Sn_Inc, diagonalLoading)
    Conds_Rxx_Inv_Inc = np.array([np.linalg.cond(r) for r in np.array(Rxx_Sn_Inv_Inc[...,:].T)], dtype=np.complex_)

    MVDR_Inc = applyMVDR(arr, True, Rxx_Sn_Inv_Inc)

    Dls = DeslocSig(MVDR_Inc,calcManifoldVector(arr, True))
    Sn = np.matrix(Dls)

    Rxx_Sn_Refl = IterativeCovariance(arr, Sn, partsAC)
    Conds_Rxx_Refl = np.array([np.linalg.cond(r) for r in np.array(Rxx_Sn_Refl[...,:].T)], dtype=np.complex_)

    Rxx_Sn_Inv_Refl = invertRxx(Rxx_Sn_Refl,diagonalLoading)
    Conds_Rxx_Inv_Refl = np.array([np.linalg.cond(r) for r in np.array(Rxx_Sn_Inv_Refl[...,:].T)], dtype=np.complex_)

    MVDR_Refl = applyMVDR(arr, False, Rxx_Sn_Inv_Refl)
       
    #'Rxx_Sn_Inv_Inc' 		: Rxx_Sn_Inv_Inc,
    #'Rxx_Sn_Refl' 			: Rxx_Sn_Refl,
    #'Rxx_Sn_Inv_Refl' 		: Rxx_Sn_Inv_Refl,
    data = {
        'fs'                    : arr.fs,
        'divisoesAC'            : partsAC,
        'diagonalLoading'       : diagonalLoading,
        'Original' 				: Original,
        'DAS_Inc' 				: DAS_Inc,
        'DAS_Refl' 				: DAS_Refl,   
        'Rxx_Sn_Inc' 			: Rxx_Sn_Inc,    
        'Conds_Rxx_Inc' 		: Conds_Rxx_Inc,
        'Conds_Rxx_Inv_Inc' 	: Conds_Rxx_Inv_Inc,
        'MVDR_Inc' 				: MVDR_Inc,
        'Conds_Rxx_Refl' 		: Conds_Rxx_Refl,
        'Conds_Rxx_Inv_Refl' 	: Conds_Rxx_Inv_Refl,
        'MVDR_Refl' 			: MVDR_Refl
        }
    with open('data/' + file_name + '.pickle', 'wb') as f:
        pk.dump(data, f, pk.HIGHEST_PROTOCOL)


#Reads a pickle file and presents its data
def Show_Data(file_name, xmin, xmax):
    files = ["data/"+file_name+'.pickle']
    data = []

    for a in files:
        with open(a, 'rb') as f:
            T = pk.load(f)
            data.append(T)
            f.close()
        
    dado = data[0]

    fs = dado['fs']
    Original = dado['Original']
    DAS_Inc = dado['DAS_Inc']
    DAS_Refl = dado['DAS_Refl']

    #Rxx_Sn_Inc = dado['Rxx_Sn_Inc']
    #Rxx_Sn_Inv_Inc = dado['Rxx_Sn_Inv_Inc']
    #Conds_Rxx_Inc = dado['Conds_Rxx_Inc']
    #Conds_Rxx_Inv_Inc = dado['Conds_Rxx_Inv_Inc']
    MVDR_Inc = dado['MVDR_Inc']

    #Rxx_Sn_Refl = dado['Rxx_Sn_Refl']
    #Rxx_Sn_Inv_Refl = dado['Rxx_Sn_Inv_Refl']
    #Conds_Rxx_Refl = dado['Conds_Rxx_Refl']
    #Conds_Rxx_Inv_Refl = dado['Conds_Rxx_Inv_Refl']
    MVDR_Refl = dado['MVDR_Refl']


    plt.figure(figsize=(9,5))
    plt.title('Original')
    plt.xlim(xmin,xmax)
    plotFreq(fs,Original,'b')

    #plt.figure()
    #plt.title('DAS Incidente')
    #plotFreq(fs,DAS_Inc,'b')

    #plt.figure()
    #plt.title('DAS Refletido')
    #plotFreq(fs,DAS_Refl,'b')

    plt.figure()
    plt.title('MVDR(Incidente)')
    plt.xlim(xmin,xmax)
    plotFreq(fs,MVDR_Inc,'b')

    #plt.figure()
    #plt.title('Condition Number Sn (MVDR Inc)')
    #nmr = [np.trace(Rxx_Sn_Inc[...,i])/np.linalg.norm(Rxx_Sn_Inc[...,i]) for i in range(0,len(Conds_Rxx_Inc)) ]
    #plotFreq(fs,Conds_Rxx_Inv_Refl,'b')

    #plt.figure()
    #plt.title('Condition Number Inverse Sn (MVDR Inc)')
    #plotFreq(fs,Conds_Rxx_Inv_Inc,'b')

    plt.figure()
    plt.title('MVDR(Refletido)')
    plt.xlim(xmin,xmax)
    plotFreq(fs,MVDR_Refl,'b')
        
    plt.figure(figsize=(9,5))
    plt.title('Coeficiente de Reflexão')
    rate = (np.absolute(MVDR_Refl)/np.absolute(MVDR_Inc))
    rate = pd.DataFrame(rate).rolling(300,min_periods=1,center=True).median()._values
    plt.xlabel('Frequência(Hz)')
    plt.ylabel('Amplitude')
    plt.xlim(xmin,xmax)
    plt.semilogx(np.linspace(0,fs/2-1,len(MVDR_Refl)),rate,'b')

    plt.figure()
    plt.title('MVDR(Incidente)')
    plotTime(fs,MVDR_Inc,'b')

    plt.figure()
    plt.title('MVDR(Refletido)')
    plotTime(fs,MVDR_Refl,'b')

    plt.show()