import numpy as np
from findmax import findmax
import conductivity_sl as sl

N = 26
zeta = np.linspace(-0.025, -0.2, 26)
SigS = np.logspace(-3, -1, 26)

lab = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12",
       "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24",
       "25", "26"]

zetam = zeta[0:N]
SigSm = SigS[0:N]
zm, sm = np.meshgrid(zetam, SigSm)
safe = []

w = np.logspace(-4,12,1600)
sigm = np.zeros((N,N))
taum = np.zeros((N,N))

for ii in range(N):
    for jj in range(N):
        sig = sl.conductivity_sl(w, 1, 293, 80, 5e-8, 5e-9, 9e-5/2, 1e-5/2, 2e-6,
                                 2e-7, zm[ii,jj], sm[ii,jj])
        np.savetxt("spectra/sl/spec_"+lab[ii]+"_"+lab[jj]+"_sl.txt",
                   np.transpose([w, sig]))
        sigmax, tau = findmax(w, sig.imag, 0)
        sigmaxh, tauh = findmax(w, sig.imag, 1)
        if tau == tauh:
            sigm[ii,jj] = np.nan
            taum[ii,jj] = np.nan
            safe.append([zm[ii,jj], sm[ii,jj], np.nan, np.nan])
        else:
            sigm[ii,jj] = sigmax
            taum[ii,jj] = tau
            safe.append([zm[ii,jj], sm[ii,jj], sigmax, tau])
        print((N*ii+jj+1)/N**2)
    np.savetxt("tau_max_sl.txt", safe)