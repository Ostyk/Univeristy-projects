import numpy as np
import pylab as py
from numpy.fft import rfft, rfftfreq
from scipy.io.wavfile import read


def periodogram(sig, Fs, okno):
    okno = okno / np.linalg.norm(okno)  # normalizacja okna
    s = sig * okno  # okienkowanie sygnalu
    S = rfft(s)
    P = S * S.conj()
    P = P / Fs
    P = P.real
    P = 20 * np.log10(np.abs(P))  # zmiana na skale decybelowa
    if len(s) % 2 == 0:  # moc z ujemnej części widma
        P[1:-1] *= 2
    else:
        P[1:] *= 2
    F = rfftfreq(len(s), 1. / Fs)
    return F, P


def spectrogram(sig, Fs, okno, N_wsp):
    N_okna = len(okno)  # nieskrocone okna #liczba kawalkow --> ilosc periodogramow
    N_cut = N_okna - N_wsp  # skrocone okna
    n = int((len(sig) - N_okna) / (N_cut) + 1)  # liczba kawalkow
    spectra = np.zeros(
        (int(N_okna / 2 + 1), n))  # bierzemy prawdziwa czesc .+1 aby sie zgadzal wymiar(linijka wczesniej n =)
    for i in range(n):  # n=360
        f, spectra[:, i] = periodogram(sig[i * N_cut:i * N_cut + N_okna], Fs, okno)  # wpakowujemy widma w kolumny
    np.place(spectra, spectra < 0, [0])  # zerujemy wszystkie te mniejze od zera
    psd_time = np.sum(spectra, axis=0)  # sumujemy kolumny, to jest moc w czasie
    psd_freq = np.sum(spectra, axis=1)  # sumujemy wiersze, to jest moc w czestosci
    return spectra, psd_time, psd_freq, n


def plot(Fs, sig, N_wsp, okno, N_okna, t_min, t_max, f_min, f_max, title, n_fig):
    a, b, c, n = spectrogram(sig, Fs, okno, N_wsp)
    time_axis = np.linspace(t_min, t_max, n)
    freq_axis = np.linspace(f_min, f_max, (int(N_okna / 2 + 1)))
    py.figure(n_fig)
    py.subplot(222)
    py.imshow(a, aspect='auto', cmap='YlOrRd', clim=(0, 120), origin='lower', extent=(t_min, t_max, f_min, f_max))
    py.title(title)
    py.subplot(224)
    py.title('moc w czasie')
    py.plot(time_axis, b)
    py.axis([min(time_axis), max(time_axis), min(b), max(b)])
    py.subplot(221)
    py.title('moc w czestosci')
    py.axis([min(c), max(c), f_min, f_max])
    py.plot(c, freq_axis)
    py.savefig(title + ".png")
    py.show()


Fs1, sig1 = read("Candy_Dulfer_-_Lily_Was_Here.wav")
Fs2, sig2 = read("Rupert_Blaise_-_06_-_What_A_Wonderful_World.wav")
songs = ["Candy_Dulfer-Lily Was Here", "Rupert_Blaise-What a Wonderful World"]
x = [[Fs1, sig1, 1, songs[0]], [Fs2, sig2, 2, songs[1]]]
for Fs, sig, n_fig, title in x:
    okno = np.blackman(4410)
    N_okna = len(okno)
    N_wsp = 10
    t_min = 0
    t_max = len(sig) / Fs
    f_min = 0
    f_max = Fs / 2
    plot(Fs, sig, N_wsp, okno, N_okna, t_min, t_max, f_min, f_max, title, n_fig)
