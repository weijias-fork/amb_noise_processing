import numpy as np
from scipy.io import loadmat as lm
import matplotlib.pyplot as plt
#from scipy.ndimage.filters import gaussian_filter as gs
from scipy.signal import gaussian
import math
np.seterr(divide='ignore')


x = lm('/Users/tzompantli/Documents/GitHub/amb_noise_processing/cc.mat')
freq = x['freq'][0]
signal = x['smoothed'][0]
dist=.160
alpha=25

print (len(signal), len(np.real(np.fft.fftshift(np.fft.ifft(signal)))))
dkj


fqa =  np.arange(.1,25,.25) # demander  a quoi sert ce param
nzeropad = 20000

nfa = len(fqa)
nf = len(freq)
df = freq[1]-freq[0]

nt = 2*(nf-1)+nzeropad
f = np.zeros(nfa)
dt = 1/((2*(nf-1)+nzeropad)*df)
t = np.arange(0, dt*(2*(nf-1)+nzeropad-1), dt)

v00 = dist/t

print (v00, np.shape(v00))
anv = 4*np.sqrt(1/dist)
nt2 = np.min([math.floor(nt/2)+1, np.where(v00>.050)[0][-1]])
nt1 = np.where(v00<1.000)[0][0]

nkeep   = 3
tresh   = .01 #diference in velocity to determine jump in DC
npoints = 5#??? max number points in jump
perc    = 50 #   - minimal length of of output segment vs freq. range
t0R     = np.zeros((nfa,nkeep))




t0wR    = np.ones((nfa,nt2))
phiR	= np.ones((nfa,nt2))
dphiR = np.ones((nfa,nt2))


#anv : ancho de la ventana (ici independant de la distance)
nw      = np.round(anv/df)        # ancho de la ventana en termino de puntos
nw      = np.min([nw,nf])
if np.mod(nw,2)==0: #                % even only to make the center ok
     nw += 1

nw2     = (nw-1)/2
# hw      = blackmanharris(nw);   %construcion de la ventana mobile de base
#hw      =  gaussian(nw,8.5)
from scipy.signal import blackmanharris
hw      =  blackmanharris(int(nw))

# % ciclo en ventana en frecuencia
# % figure(10);hold on;
for j in np.arange(0, nfa, 1):
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #% construcion del filtro con troncatura en los bordes %
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #% se tendria que agregar una interpolacion de hw para mejor caer sobre
    #% la frecuencia deseada
#% wi = interp1(X.f, tap, Fqi, 'pchip');
#% [~, ij] = max(wi);
#% f(j) = Fqi(ij);

    #% look for the index of the frequency
    ic = np.argmin(np.abs(fqa[j] - freq))
#    %construct the windows
    i0w = ic - nw2
    ifw = ic + nw2
    if i0w<1:
        #%left truncated
        if i0w == 0:
            i0w = -1
        tmpw = hw[int(-i0w):int(nw)]
        if len(tmpw)>nf:
            tmpw=tmpw[:int(nf)]

        tap = np.concatenate((tmpw, np.zeros(nf-len(tmpw))))
    elif ifw > nf:
        #%right truncated
        tmpw = hw[:int(nw-(ifw-nf))]
        tap = np.concatenate((np.zeros(nf-len(tmpw)),tmpw))
    else:
        #%full
        tap  = np.concatenate((np.zeros(int(i0w-1)),hw, np.zeros(int(nf-(i0w-1)-nw))))

    #% frecuencia del maximo de la ventana mobil
    #% [~, ij] = max(tap);
    #% correction for not symetric taper
    ij = np.where(np.cumsum(tap) >= np.sum(tap) / 2)[0][0]

    #    % frequency of the  ventana
    f[j] = freq[ij] + df / 2

    #% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    #% selecion de la rebanada del espectro %
    #% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    sp1 = 1j * np.real(signal)*tap

    #% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    #% construction de la envolvente de la forma de onda %
    #% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    #% crepe + analitic
    if nzeropad == 0:
        sptotm = []
    else:
        sptotm = np.zeros((nzeropad))

    sptot0 = sp1[:int(nf)]
    sptotf = np.zeros(int(nf) - 2)
    sptot = 2 * np.concatenate((sptot0, sptotm, sptotf))

    isp = np.fft.ifft(sptot)
    sig = np.abs(isp)
    phi = np.unwrap(np.angle(isp))
    envolv = sig[:nt2]

    t0wR[j,:]= t0wR[j]* envolv #;%./max(envolv);%somme logarithmic campillo
    phiR[j,:]= phi[:nt2]#%./max(envolv);
    dphiR[j,:]=np.gradient(phiR[j])





fac0=1e-3


dnt=np.max(np.round(nt2/600))

tmp = t0wR[:int(nt2)]
tmp = tmp[::int(dnt/4)]
x = f
y = v00[:int(nt2)]


print(np.shape(tmp), np.shape(f), np.shape(v00[:int(nt2)]))

x, y  = np.meshgrid(x,y)

print()

z = np.log( (fac0*np.abs(tmp)) +1 )#[:,1:dnt:nt2]))
#plt.imshow(Z, aspect='auto')

print (np.shape(x), np.shape(y), np.shape(z))

print (x, y, z)
plt.xlim(0,24)
plt.ylim(0,1)

plt.pcolor(x,y,z.T,cmap=plt.jet())
plt.show()
