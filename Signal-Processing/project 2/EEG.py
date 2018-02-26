import numpy as np
from  scipy.signal import butter, buttord     # funkcje do projektowania filtrów  
from  scipy.signal import cheby1, cheb1ord    # funkcje do projektowania filtrów 
from  scipy.signal import cheby2, cheb2ord    # funkcje do projektowania filtrów 
from  scipy.signal import ellip, ellipord,filtfilt     # funkcje do projektowania filtrów eliptycznych
import pylab as py
from funkcje import sortandplot,sortandplot2,compare,compare2,parametry_atomu,TFRPlot,tfr_atomu,rekonstrukcja_atomu#,skalogram

file=open('inb14_fragment.bin', 'rb')
s = np.fromfile(file, dtype='float32')
Fs=128
nq=Fs/2.
wp=11/(Fs/2.)
ws=16/(Fs/2.)
gpass = 1
gstop = 25
rp = gpass
rs = gstop
analog=0
T = np.arange(0, 600, 1./Fs)
t=2
f = np.arange(0.01,Fs/2,0.01)
#########Filtrowanie
##filtruj sygnał filtrem pasmowo-przepustowym w paśmie 11-16 Hz
###BUTTER##
[n1,Wn1]=buttord(wp,ws, gpass, gstop)
[b1,a1]=butter(n1,[11./nq,16./nq],btype='bandpass')
filter_1_1 = filtfilt(b1,a1,s)
###Czebysze1##
[n2,Wn2]=cheb1ord(wp, ws, gpass, gstop, analog)
[b2,a2]=cheby1(n2, gpass, [11./nq,16./nq], btype='bandpass',  output='ba')
filter_1_2 = filtfilt(b2,a2,s)
###Czebysze2###
[n3,Wn3]=cheb2ord(wp, ws, gpass, gstop, analog=0);
[b3,a3]=cheby2(n3, gstop, [11./nq,16./nq], btype='bandpass',  output='ba')
#filter_1_2 = filtfilt(b3,a3,s)
filter_1_3 = filtfilt(b3,a3,s)
###eliptyczny##
[n4,Wn4]=ellipord(wp, ws, rp,rs);
[b4,a4]=ellip(n4, rp, rs, [11./nq,16./nq], btype='bandpass', analog=0, output='ba')
filter_1_4 = filtfilt(b4,a4,s)

py.figure(1)
compare(a1,b1,a2,b2,a3,b3,a4,b4,f,t,Fs,"Bandpass filter--")
#podniesione do kwaratu
filter_2_1,filter_2_2,filter_2_3,filter_2_4=filter_1_1**2,filter_1_2**2,filter_1_3**2,filter_1_4**2


#przefiltruj filtrem dolnoprzepustowym z częstości�
 odcięcia 10 Hz
wp1=9/(Fs/2.)#testowanie
ws1=10/(Fs/2.)#testowanie

[n1,Wn1]=buttord(wp1,ws, gpass, gstop,analog=0)
[b1,a1]=butter(n1,Wn1,btype='lowpass',output='ba',analog=0)
filter_3_1 = filtfilt(b1,a1,np.sqrt(filter_2_1))
###Czebysze1##
[n2,Wn2]=cheb1ord(wp1, ws1, gpass, gstop, analog)
[b2,a2]=cheby1(n2, gpass, Wn2, btype='lowpass',  output='ba')
filter_3_2 = filtfilt(b2,a2,np.sqrt(filter_2_2))
###Czebysze2###
[n3,Wn3]=cheb2ord(wp1, ws1, gpass, gstop, analog=0);
[b3,a3]=cheby2(n3, gstop, Wn3, btype='lowpass',  output='ba')
#filter_1_2 = filtfilt(b3,a3,s)
filter_3_3 = filtfilt(b3,a3,np.sqrt(filter_2_3))
###eliptyczny##
[n4,Wn4]=ellipord(wp1, ws1, rp,rs);
[b4,a4]=ellip(n4, rp, rs, Wn4, btype='lowpass', analog=0, output='ba')
filter_3_4 = filtfilt(b4,a4,np.sqrt(filter_2_4))

py.figure(2)
compare(a1,b1,a2,b2,a3,b3,a4,b4,f,t,Fs,"Lowpass filter--")

py.figure(3)
py.subplot(221)
py.ylabel("[micro V]")
filter_4_1= sortandplot(filter_3_1,5,64,T,s,Fs,"butter") #funkcja zaimportowana
print("Wrzecion dla butter jest: {}".format(len(filter_4_1)))

py.subplot(222)
filter_4_2= sortandplot(filter_3_2,5,64,T,s,Fs,"Czebyszew 1") #funkcja zaimportowana
print("Wrzecion dla czebyszewa 1 jest: {}".format(len(filter_4_2)))

py.subplot(223)
py.ylabel("[micro V]")
py.xlabel("czas [s]")
filter_4_3= sortandplot(filter_3_3,5,64,T,s,Fs,"Czebyszew 2") #funkcja zaimportowana
print("Wrzecion dla czebyszewa 2 jest: {}".format(len(filter_4_3)))

py.subplot(224)
py.xlabel("czas [s]")
filter_4_4= sortandplot(filter_3_4,5,64,T,s,Fs,"Eliptyczny") #funkcja zaimportowana
print("Wrzecion dla eliptyczny jest: {}".format(len(filter_4_4)))

py.figure(4)
py.subplot(221)
compare2(a1,b1,f,T,Fs,"--butter")
py.subplot(222)
compare2(a2,b2,f,T,Fs,"--Czebyszew 1")
py.subplot(223)
compare2(a3,b3,f,T,Fs,"--Czebyszew 2")
py.subplot(224)
compare2(a4,b4,f,T,Fs,"--Eliptyczny")


"""
Zatem widać, że filtr czebyszew 2 ma najbardziej gładk�
 odpowiedź impulsow�
 oraz najmniejsz�
 amplitudę
po filtrowaniu filtrem pasmowo-przepustowym (Figura 1), jak i 
najmniej chaotyczn�
 odpowiedź impulsow�
 po przefiltrowaniu filtrem
dolno-przepustowym(Figura 3). 
Natomiast filtr eliptyczny także w dwóch przypadkach okazał się dobry,
więc porównano dodatkowo opóźnienie fazowe, i w końcu wybrano jednak filtr eliptyczny ponieważ ma on
najmniejsze opóźnienie fazowe (strome zbocze) jak i nie wycina wielu potrzebnych fragmentów sygnału. 
Filtr eliptyczny daję mniejsz�
 liczbę wrzecion co sugeruję, że
w Czebyszewie 2 mogło dojść do zaburzę które zostały zarejestrowane jako wrzeciona.
"""
"""
##################
##     Drugi    ##
## klasyfikator ##
##################

import book_reader as br
import os

a=np.array(s,'float32')
f=open('inb14_fragment.dat','wb')
a.tofile(f)
f.close()


PlikSygnalu='inb14_fragment.dat' 

liczbaProbek_w_Epoce=Fs*20
liczbaKanalow=1
wybraneKanaly=1
wybraneEpoki="1-30"
maxIteracji=50
procentEnergii=95.
energyError=0.1
algorytm='SMP'
    
fo = open('inb14_fragment.set', "wt")
fo.write('# OBLIGATORY PARAMETERS\n')
fo.write('nameOfDataFile            '+PlikSygnalu+'\n')
fo.write('nameOfOutputDirectory     ./\n')
fo.write('writingMode               CREATE \n')  
fo.write('samplingFrequency         '+str(Fs)+'\n')
fo.write('numberOfChannels          '+str(liczbaKanalow)+'\n')
fo.write('selectedChannels          '+str(wybraneKanaly)+'\n')
fo.write('numberOfSamplesInEpoch    '+str(liczbaProbek_w_Epoce)+ '\n')
fo.write('selectedEpochs            '+str(wybraneEpoki)+'\n')
fo.write('typeOfDictionary          OCTAVE_FIXED\n')
fo.write('energyError               '+str(energyError)+' 100.0 \n')
fo.write('randomSeed                auto \n')
fo.write('reinitDictionary          NO_REINIT_AT_ALL \n')
fo.write('maximalNumberOfIterations '+str(maxIteracji)+'\n')
fo.write('energyPercent             '+str(procentEnergii)+'\n')
fo.write('MP                        '+algorytm+'\n')
fo.write('scaleToPeriodFactor       1.0 \n')
fo.write('pointsPerMicrovolt        1.0 \n')
fo.write('\n# ADDITIONAL PARAMETERS\n')
fo.write('normType                  L2 \n') 
fo.write('diracInDictionary         YES \n')
fo.write('gaussInDictionary         YES \n')
fo.write('sinCosInDictionary        YES \n')
fo.write('gaborInDictionary         YES \n')
fo.close()
    
os.system('empi-win64.exe inb14_fragment.set')
book = br.BookImporter('inb14_fragment_smp.b')
py.figure(5)
py.subplot(212)
py.plot(T,s,'k')
for a in range (1,31):
    for atom in book.atoms[a]:
        if (atom['params']['f']*book.fs/2>=11 
                and atom['params']['f']*book.fs/2<=16 
                and atom['params']['scale']/book.fs>=.5):
            py.plot(T[atom['params']['t']+(a-1)*20*Fs:atom['params']['t']+(a-1)*20*Fs+atom['params']['scale']],s[atom['params']['t']+(a-1)*20*Fs:atom['params']['t']+atom['params']['scale']+(a-1)*20*Fs],'r')
            py.title('Wrzeciona snu z drugiego klasyfikatora')

py.subplot(211)
filter_4_4= sortandplot(filter_3_4,5,64,T,s,Fs,"Eliptyczny")
py.show()
rekonstrukcja = np.zeros(len(s)) 
mapaEnergii = np.zeros((128.,len(s))) 

for a in range (1,31):
    for atom in book.atoms[a]:     
        f_Hz, A, faza, t0, skala  = parametry_atomu(book, atom)
        t, atom_czas = rekonstrukcja_atomu(book, atom)
        t,f, atom_tfr = tfr_atomu (book, atom, 128.)
        mapaEnergii[:,(a-1)*Fs*20:a*Fs*20] += atom_tfr
        rekonstrukcja[(a-1)*Fs*20:a*Fs*20] += atom_czas


py.figure(6)
"""
py.subplot(411)
#py.xlim([0,20])
py.ylim([1,3])
t1,f,p=skalogram(s,50,1.0,20,0.5,Fs)
TFRPlot(p,t1,f,Fs,'Skalogram')
"""

py.subplot(411)
py.plot(T,s)
py.title("sygnal orginalny")

py.subplot(412)
py.ylim([1,3])
py.xlim([0,600])
TFRPlot(mapaEnergii, t, f, Fs=128,title ='Mapa czas-czestosc z MP') 

py.subplot(413)
py.xlim([0,600])
py.ylim([0,2])
filter_4_4= sortandplot2(filter_3_4,5,64,T,s,Fs,"Eliptyczny")

py.subplot(414)
for a in range (1,31):
    for atom in book.atoms[a]:
        if (atom['params']['f']*book.fs/2>=11 
                and atom['params']['f']*book.fs/2<=16 
                and atom['params']['scale']/book.fs>=.5):
            py.plot(T[atom['params']['t']+(a-1)*20*Fs:atom['params']['t']+(a-1)*20*Fs+atom['params']['scale']],s[atom['params']['t']+(a-1)*20*Fs:atom['params']['t']+atom['params']['scale']+(a-1)*20*Fs],'r')
            py.ylim([0,2])            
            py.title('Wrzeciona snu z drugiego klasyfikatora')
"""
"""
Zgodność wyników obu klasyfikatorów zależy przede wszystkim od odpowiedniego doboru parametru energyError.
Czym mniejszy energyError tym lepsze oszacowani liczby wrzecion. Parametr ten jest kluczowy do znajdywania wrzecion
metod�
 MP. Jest głównym regulatorem kompromisu stosunku użytej pamięci obliczeniowej do jakości dekompozycji.


##########################################################
##Energy Error## ilość wrzecion## ~czas obliczenia [min]##
##########################################################
##    0.1     ##      11       ##          6            ##
##    0.01    ##      20       ##          25           ##
########################################################## 

Parametr energyError jest tak naprawdę kluczowy dla znajdowania kolejnych wrzecion metod�
 dekompozycji MP.
Reguluje on kompromis stosunku użytej pamięci obliczeniowej do jakości dekompozycji. Czym wartość się zbliża
do zero tym większa liczba funkcji słowniku, skutkuj�
c w lepsz�
 jakość dekompozycji MP kosztem zajmowania
ogromnej pamięci obliczeniowej.

Główn�
 wad�
 zastosowania klasyfikatora opartego o dekompozycję MP zatem jest czas obliczeniowy.
Aby uzyskać bardziej precyzyjny wynik warto zwiększać liczbę interacji, ponieważ wrzeciona snu mog�
 
mieć mniejsz�
 energię niż inne struktury sygnału, ale w tym samym czasie obliczenia by były jeszcze 
bardziej uci�
żliwe dla systemu(pamięć obliczeniowa) jak i by trwały jeszcze dłużej.

Pierwszy klasyfikator także ma swoje wady, a do ich określenia potrzebna jest analiza własności danego filtru.
Wybrany filtr eliptyczny ma wadę tak�
 że może gubić wrzeciona snu jeżeli dudnienia s�
 w paśmie zaporowym.

W pierwszym klasyfikatorze uzyskano :
#####################################
##   Filtr     ##  ilość wrzecion  ##
#####################################
## butter      ##       15         ##
##czebyszewa 1 ##        9         ##
##czebyszew 2  ##       15         ##
##eliptyczny   ##       11         ##
#####################################
Zatem ilość wrzecion znalezionych pierwszym klasyfikatorem poprzez wybrany filtr (eliptyczny),
zgadza się z ilości�
 wrzecion znalezionych drugim klasyfikatorem na poziomie wartości energyError = 0.1.
Natomiast przy energyError 0.01 wrzecion jest 20 co mogło by sugerować wybranie 
innego filtru(np. czebyszew II) jak pocz�
tkowo była przeprowadzona analiza.

Podsumowuj�
c, zgodność dwóch klasyfikatorów zależy od wybranych parametrów klasyfikatora. Patrzac na figure 6,
widać, że pomimo zgodności ilości to wrzeciona nie występuj�
 w tym samym miejscu w sygnale sugeruj�
c że,
wrzeciona bardzo zalęża od techniki pomimo faktu, że w teorii powinny znajdywać się w tym samym miejscu.

Typu filtru dl pierwszego, a wartości energyError dla drugiego. 
Natomiast dekompozycja MP byłaby najlepszym klasyfikatorem odnajdowania wrzecion snu, gdyby dysponowało się ogromn�
 pamięci�
 operacyjn�
.
"""