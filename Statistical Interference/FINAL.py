# -*- coding: utf-8 -*-
import scipy.optimize as opt
import pylab as py
import numpy as np
import scipy.stats as st

xdata,ydata = np.loadtxt("dane.txt")
N=len(ydata)
def funkcja_do_fitowania(x,y0, A, x0, w, a):
	y =y0+ A*(x-x0)*np.sin(w*(x-x0)**a)
	return y
def funkcja_bledu(x, y, funkcja, params, err):
	y_fit = funkcja(x, *params) # aktualne wartosci y z dopasowania
	residuum = y-y_fit # residua wchodzą do sumy kwadratów z wagą odwrotnie proporcjonalną do standardowego odchylenia
	residuum_wazone = residuum/ err 
	return residuum_wazone	
def funkcja_bledu_dla_wielomianow(x, y, wspolczynniki1, err1):
	W = np.poly1d(wspolczynniki1)
	y_fit1 = W(x) # aktualne wartości y z dopasowania
	residuum1 = y-y_fit1 # residua wchodzą do sumy kwadratów z wagą odwrotnie proporcjonalną do standardowego odchylenia
	residuum_wazone1 = residuum1/ err1 
	return residuum_wazone1

params_final, covar = opt.curve_fit(funkcja_do_fitowania,xdata,ydata,[0,1,-0.5,1,2],0.05)
print("Dopasowane parametry funkcji: ",params_final)
print("Macierz kowariancji funkcji :\n",covar)
A,B,C,D,E = params_final # dopasowane parametry
residua = funkcja_bledu(xdata, ydata, funkcja_do_fitowania, params_final, 0.05)# tym razem residua już są podzielone przez standardowe odchylenie, każde przez swoje
W, p =st.shapiro(residua)
print('Test normalnosci residuow funkcji : p = %.3f'%(p))

params_final_wielomian, covar1=np.polyfit(xdata, ydata, deg=15,cov=True) 
print("Dopasowane wspolczynniki wielomianu : ",params_final_wielomian)
print("Macierz kowariancji :\n",covar)
residua_wielomian = funkcja_bledu_dla_wielomianow(xdata, ydata, params_final_wielomian, 0.05)# tym razem residua juz sa podzielone przez standardowe odchylenie, każde przez swoje
W, p =st.shapiro(residua_wielomian)
print('Test normalnosci residuow wielomianu: p = %.3f'%(p))

residua_teoretyczne_chi2=np.sum(residua**2)
chi2_t=np.zeros(26)+residua_teoretyczne_chi2
wielomiany_chi2=np.zeros(26)

for i in range(0,26):
    params_final_wielomian=np.polyfit(xdata, ydata,i)
    W = np.poly1d(params_final_wielomian)
    residua_w=ydata-W(xdata)
    wielomiany_chi2[i]=np.sum(residua_w**2)/0.05**2

p_wielomiany_chi2=np.zeros(26)
for i in range(0,26):
	if wielomiany_chi2[i] < N - 5:
		p_wielomiany_chi2[i]=st.chi2.sf(wielomiany_chi2[i], N-5)
	else:
		p_wielomiany_chi2[i] = 1 - st.chi2.sf(wielomiany_chi2[i], N-5)
p_chi_t=st.chi2.cdf(residua_teoretyczne_chi2, N-5)
p_chi2_t=np.zeros(26)+p_chi_t #do rysowania

py.subplot(2,2,1)# RYSUJEMY
W=np.poly1d(params_final_wielomian)
py.plot(xdata,ydata,".k",label="dane")
py.plot(xdata, W(xdata),"b-",label="wielomian W(x) rzedu 15")
py.plot(xdata, funkcja_do_fitowania(xdata,A,B,C,D,E),"g-",label="funkcja F(x)") 
py.legend(loc = 2)
py.xlabel('x')
py.ylabel('y, W(x), F(x)')
py.title(u'dane dopasowane do krzywej')

py.subplot(2,2,2)
py.plot(xdata, residua,"g-",label="funkcja F(x)") # residua funkcji
py.plot(xdata, residua_wielomian,"b-",label="wielomian W(x) rzedu 15") # residua wielomianu
py.xlabel('x')
py.ylabel('y-W(x), y-F(x)')
py.title(u'Wykres residuów')
py.legend(loc = 1)

py.subplot(2,2,3)
py.plot(np.arange(0,26),np.ones(26)+residua_teoretyczne_chi2,'g--', label="funkcja F(x)") #żeby wyrysowało prostą 
py.plot(np.arange(0,26),wielomiany_chi2 ,'b.-', label="wielomian W(x) rzedu n")
py.yscale('log')
py.title('statystyka chi2')
py.legend(loc = 1)
py.xlabel('n (stopien wielomianu)')
py.ylabel('Chi2_fit')
py.title('Chi2')

py.subplot(2,2,4)
p_wielomiany_chi2,=py.plot(np.arange(0,26),1-p_wielomiany_chi2, 'b.-', label="wielomian W(x) rzedu n")
chi2_t,=py.plot(np.arange(0,26),p_chi2_t, 'g--', label="funkcja F(x)")
py.legend(loc = 3)
py.xlabel('n (stopien wielomianu)')
py.ylabel('Chi2_fit')
py.title('Prawdopodobienstwo uzyskania statystyki chi2')
py.show()