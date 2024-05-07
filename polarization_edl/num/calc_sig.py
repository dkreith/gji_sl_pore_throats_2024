import numpy as np
import matplotlib.pyplot as plt

lab = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12",
       "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24",
       "25", "26"]
N = len(lab)

for ii in range(N):
    for jj in range(N):
        SigS = lab[ii]
        zeta = lab[jj]        
        
        dl = np.loadtxt("DL/spec_"+SigS+"_"+zeta+"_d.txt", skiprows=5,
                        dtype=np.complex_)
        sl = np.loadtxt("SL/spec_"+SigS+"_"+zeta+"_s.txt", skiprows=8,
                        dtype=np.complex_)

        w = dl[:,0].real
        sig_d = dl[:,1]
        sig_s = sl[1:len(sl)]

        sig = sig_d+sig_s
        
        np.savetxt("tot/spec_"+SigS+"_"+zeta+"_num.txt",
                   np.transpose([w, sig]))

"""
fig, ax = plt.subplots(2,1)
ax[0].semilogx(w, sig.real)
ax[0].semilogx(w, sig_d.real)
ax[0].semilogx(w, sig_s.real)
ax[1].semilogx(w, sig.imag)
ax[1].semilogx(w, sig_d.imag)
ax[1].semilogx(w, sig_s.imag)
"""