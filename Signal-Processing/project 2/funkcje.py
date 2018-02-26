from  scipy.signal import freqz  #funkcja obliczająca funkcję systemu
from  scipy.signal import lfilter #filtfilt # funkcje do aplikowania filtrów
import numpy as np
import pylab as py

def group(L):
    first = last = L[0]
    for n in L[1:]:
        if n - 1 == last: 
            last = n
        else: 
            yield first, last
            first = last = n
    yield first, last 
def sortandplot(x,min_val,length,T,s,Fs,title):
    py.plot(T,s,'k')
    momenty=[]    
    poz = np.where(x>min_val)
    q = list(group(poz[0]))
    w=poz[0]
    dt=1./Fs
    for i in q:
        if (np.diff(i)+1)>length:
            momenty.append(x[min(i):max(i)+1])  
    b = []
    tr=" dla filtru-- "
    for i in range (len(w)-1):
        if (w[i]+1==w[i+1]):
            b.append(w[i])
        else:
            if len(b)>=length:
                py.plot(np.array(b)*dt,s[b],'c')
                py.title('Wrzecion snu z pierszego klasyfikatora jest: '+str(len(momenty))+str(tr) +str(title))
            b = []
    return momenty

def sortandplot2(x,min_val,length,T,s,Fs,title):
    #py.plot(T,s,'k')
    momenty=[]    
    poz = np.where(x>min_val)
    q = list(group(poz[0]))
    w=poz[0]
    dt=1./Fs
    for i in q:
        if (np.diff(i)+1)>length:
            momenty.append(x[min(i):max(i)+1])  
    b = []
    tr=" dla filtru-- "
    for i in range (len(w)-1):
        if (w[i]+1==w[i+1]):
            b.append(w[i])
        else:
            if len(b)>=length:
                py.plot(np.array(b)*dt,s[b],'c')
                py.title('Wrzecion snu z pierszego klasyfikatora jest: '+str(len(momenty))+str(tr) +str(title))
            b = []
    return momenty


def compare(a1,b1,a2,b2,a3,b3,a4,b4,f,T,Fs,title):
    t = np.arange(-T, T, 1/Fs)
    # obliczamy odpowiedź impulsową
    x = np.zeros(len(t))
    x[len(t)//2] = 1 # impuls
    y1 = lfilter(b1,a1,x)
    y2 = lfilter(b2,a2,x)
    y3 = lfilter(b3,a3,x)
    y4 = lfilter(b4,a4,x)      
    
    py.subplot(411)    
    py.title(str(title)+'butter')    
    py.plot(t, x,t,y1) 
    py.xlim([-T/2,T])
    py.grid('on')
    py.ylim([-0.1,0.1])
    py.ylabel('Amplitude [dB]')
    
    py.subplot(412) 
    py.title('czebyszew 1') 
    py.plot(t, x,t,y2)    
    py.xlim([-T/2,T])
    py.grid('on')
    py.ylim([-0.1,0.1])
    py.ylabel('Amplitude [dB]')
    py.subplot(413)    
    py.title('czebyszew 2')    
    py.plot(t, x,t,y3)
    py.xlim([-T/2,T])
    py.grid('on')
    py.ylim([-0.1,0.1])    
    py.ylabel('Amplitude [dB]')
    py.subplot(414) 
    py.title('eliptyczny') 
    py.plot(t, x,t,y4) 
    py.xlim([-T/2,T])
    py.ylim([-0.1,0.1])
    py.grid('on')
    py.xlabel('Częstość [Hz]')
    py.ylabel('Amplitude [dB]')    
    
def compare2(a,b,f,T,Fs,title):
    # oś częstości przeliczamy na radiany
    w = 2*np.pi* f/Fs
    # obliczamy transmitancję
    w, h = freqz(b,a,w)     
    # obliczamy fazę i "rozwijamy" ją   
    faza = np.unwrap(np.angle(h))
    # obliczamy opóźnienie fazowe
    opoznienieFazowe = -faza/w
    py.title('opóźnienie fazowe'+str(title))
    py.plot(f, opoznienieFazowe)
    py.ylabel('próbki')
    py.xlabel('Częstość [Hz]')
    py.grid('on')

def parametry_atomu(book, atom):
   f_Hz  = atom['params']['f']*book.fs/2     
   A     = atom['params']['amplitude']       
   faza  = atom['params']['phase']           
   t0    = atom['params']['t']/book.fs       
   skala = atom['params']['scale']/book.fs   
   return f_Hz, A, faza, t0, skala  
   
def TFRPlot(TFR, t_mapy, f_mapy, Fs=128,title =''):
    df = f_mapy[1]-f_mapy[0]
    dt = t_mapy[1]-t_mapy[0]
    py.xlim((t_mapy.min(), t_mapy.max()))
    py.imshow(TFR,aspect='auto',origin='lower',interpolation='nearest', 
                  extent=(t_mapy.min()-dt/2,t_mapy.max()+dt/2,f_mapy.min()-df/2,f_mapy.max()+df/2))
    py.title(title)
    py.show()
    
def tfr_atomu (book, atom, N_czestosci):
   Fs=128.
   f_Hz, A, faza, t0, skala = parametry_atomu(book, atom)
   t = np.arange(0,book.epoch_s/book.fs,1/book.fs)
   f = np.linspace(0, Fs / 2, book.fs)
   rec_t = np.zeros((1,book.epoch_s))
   rec_f = np.zeros((N_czestosci,1))
   rec_t[0,:] = np.exp(-np.pi*((t-t0)/skala)**2)     
   rec_f[:,0] = np.exp(-np.pi*((f-f_Hz)*skala)**2)   
   tfr_atom = np.kron(rec_t,rec_f) 
   tfr_atom/= np.sum(np.sum(tfr_atom))  
   tfr_atom *= atom['params']['modulus']**2 
   return t, f, tfr_atom

def rekonstrukcja_atomu(book, atom):
   f_Hz, A, faza, t0, skala = parametry_atomu(book, atom)
   t = np.arange(0,book.epoch_s/book.fs,1/book.fs)     
   rekonstrukcja = A * np.exp(-np.pi*((t-t0)/skala)**2)*np.cos(2*np.pi*f_Hz*(t-t0)+faza) 
   return t, rekonstrukcja

def skalogram(x,  w=7.0, MinF = 1.0 ,MaxF = 64, df=1.0, Fs=128.):
    T= len(x)/Fs
    M = len(x)
    t = np.arange(0,T,1./Fs)
    freqs = np.arange(MinF,MaxF,df)
    P = np.zeros((len(freqs),M))
    X = np.fft.fft(x)
    for i,f in enumerate(freqs):
        s = T*f/(2*w)
        xx = np.linspace(-s * 2 * np.pi, s * 2 * np.pi, M)
        psi = np.exp(1j * w * xx) * np.exp(-0.5 * (xx**2)) * np.pi**(-0.25)   
        Psi = np.fft.fft(psi)
        Psi /= np.sqrt(np.sum(Psi*Psi.conj()))  
        tmp = np.fft.fftshift(np.fft.ifft(X*Psi))  
        P[i,:] = (tmp*tmp.conj()).real 
    return t, freqs, P