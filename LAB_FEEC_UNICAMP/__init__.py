from .Array import *
from .Utils  import *
from .Microphone import *
from .Beamforming   import *
from .CustomPlots  import *
from .Test_Show import *


#Specific code for the project.
#Load the measurements of the 4-microphone linear array stored in ITA-Toolbox format
arranjoGde = []
arranjoMed = []
arranjoPqn = []
def LoadITAs():
    #
    #Velocidade do Som
    c_ = 343.30
    #FS das medidas
    fs_ = 44100

    #Distâncias  GDE
    #amos - mic1: 3cm
    #mic1 - mic2: 9cm
    #mic2 - mic3: 27cm
    #mic3 - mic4: 18cm
    #mic4 - autf: nome do arquivo
    medicoesGde = [ ['ITAs/ArranjoGde/029cm.ita',  0.03,0.09,0.27,0.18,0.29],
                    ['ITAs/ArranjoGde/054cm.ita',  0.03,0.09,0.27,0.18,0.54],
                    ['ITAs/ArranjoGde/079cm.ita',  0.03,0.09,0.27,0.18,0.79],
                    ['ITAs/ArranjoGde/105cm.ita',  0.03,0.09,0.27,0.18,1.05],
                    ['ITAs/ArranjoGde/131cm.ita',  0.03,0.09,0.27,0.18,1.31],
                    ['ITAs/ArranjoGde/232,5cm.ita',0.03,0.09,0.27,0.18,2.325]]
    
    for m in medicoesGde:
        arranjoGde.append(Array.fromIta(c=c_,fs=fs_,arq=m[0],dis_mic_amos=m[1],dis_mic_autf=m[5]))

    #Distâncias MED
    #amos - mic1: 34cm
    #mic1 - mic2: 14cm
    #mic2 - mic3: 21cm
    #mic3 - mic4: 7cm
    #mic4 - autf: nome do arquivo
    medicoesMed = [ ['ITAs/ArranjoMed/047cm.ita',  0.34,0.14,0.21,0.07,0.47],
                    ['ITAs/ArranjoMed/063cm.ita',  0.34,0.14,0.21,0.07,0.63],
                    ['ITAs/ArranjoMed/083cm.ita',  0.34,0.14,0.21,0.07,0.83],
                    ['ITAs/ArranjoMed/104cm.ita',  0.34,0.14,0.21,0.07,1.04],
                    ['ITAs/ArranjoMed/122cm.ita',  0.34,0.14,0.21,0.07,1.22],
                    ['ITAs/ArranjoMed/145cm.ita',  0.34,0.14,0.21,0.07,1.45],
                    ['ITAs/ArranjoMed/161cm.ita',  0.34,0.14,0.21,0.07,1.61]]
    
    for m in medicoesMed:
        arranjoMed.append(Array.fromIta(c=c_,fs=fs_,arq=m[0],dis_mic_amos=m[1],dis_mic_autf=m[5]))

    #Distâncias Pqn
    #amos - mic1: 3cm
    #mic1 - mic2: 3cm
    #mic2 - mic3: 9cm
    #mic3 - mic4: 6cm
    #mic4 - autf: nome do arquivo
    medicoesPqn = [ ['ITAs/ArranjoPqn/065cm.ita',  0.03,0.03,0.09,0.06,0.65],
                    ['ITAs/ArranjoPqn/092,5cm.ita',0.03,0.03,0.09,0.06,0.925],
                    ['ITAs/ArranjoPqn/120cm.ita',  0.03,0.03,0.09,0.06,1.20],
                    ['ITAs/ArranjoPqn/149cm.ita',  0.03,0.03,0.09,0.06,1.49],
                    ['ITAs/ArranjoPqn/174,5cm.ita',0.03,0.03,0.09,0.06,1.745],
                    ['ITAs/ArranjoPqn/268cm.ita',  0.03,0.03,0.09,0.06,2.68]]
    
    for m in medicoesPqn:
        arranjoPqn.append(Array.fromIta(c=c_,fs=fs_,arq=m[0],dis_mic_amos=m[1],dis_mic_autf=m[5]))

