import numpy as np
from findmax import findmax

zeta = np.linspace(-0.025, -0.2, 26)
SigS = np.logspace(-3, -1, 26)

lab = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12",
       "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24",
       "25", "26"]
N = len(lab)
sigm = np.zeros((N,N))
taum = np.zeros((N,N))

zeta = np.linspace(-0.025, -0.2, N)
SigS = np.logspace(-3, -1, N)

zetam = zeta[0:N]
SigSm = SigS[0:N]
zm, sm = np.meshgrid(zetam, SigSm)
safe = []

for ii in range(N):
    for jj in range(N):
        dat = np.loadtxt("num/tot/spec_"+lab[ii]+"_"+lab[jj]+"_num.txt",
                         dtype=np.complex_)
        w = dat[:,0]
        sig = dat[:,1]
        sigmax, tau = findmax(w.real, sig.imag, 0)
        sigm[ii,jj] = sigmax
        taum[ii,jj] = tau
        safe.append([zm[ii,jj], sm[ii,jj], sigmax, tau])
    np.savetxt("tau_max_num.txt", safe)