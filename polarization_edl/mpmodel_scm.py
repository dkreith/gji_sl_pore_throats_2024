import numpy as np
from findmax import findmax
import conductivity_sl as sl

strlist = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11",
           "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22",
           "23", "24", "25", "26", "27", "28", "29", "30"]
c0_c = np.logspace(-3,-1,30)
c0_s = np.logspace(-4,-1,30)
Sig0 = np.logspace(-1,np.log10(0.15),30)
pH = np.linspace(5,9,30)

# Import clay data
zeta_c = np.zeros((30,30))
SigS_c = np.zeros((30,30))
test_c = np.zeros((30,30))

for ii in range(30):
    dat = np.loadtxt("pH_c0_Sig0/clay/Sig_"+strlist[ii]+".txt", skiprows=5)
    for jj in range(30):
        zeta_c[ii,jj] = dat[jj,1]
        SigS_c[ii,jj] = dat[jj,2]
        test_c[ii,jj] = dat[jj,3]
        
# Import quartz data
zeta_s = np.zeros((30,30))
SigS_s = np.zeros((30,30))
test_s = np.zeros((30,30))

for ii in range(30):
    dat = np.loadtxt("pH_c0_Sig0/quartz/pH_"+strlist[ii]+".txt", skiprows=5)
    for jj in range(30):
        zeta_s[ii,jj] = dat[jj,1]
        SigS_s[ii,jj] = dat[jj,2]
        test_s[ii,jj] = dat[jj,3]

w = np.logspace(-2,12,1400)

safe_c = []
sigm_c = np.zeros((30,30))
taum_c = np.zeros((30,30))

for ii in range(30):
    for jj in range(30):
        sig = sl.conductivity_sl(w, 1e3*c0_c[jj], 293, 80, 5e-8, 5e-9, 9e-5/2,
                                 1e-5/2, 2e-6, 2e-7, zeta_c[ii,jj],
                                 SigS_c[ii,jj])
        sigmax, tau = findmax(w, sig.imag, 0)
        sigmaxh, tauh = findmax(w, sig.imag, 1)
        if tau == tauh:
            sigm_c[ii,jj] = np.nan
            taum_c[ii,jj] = np.nan
            safe_c.append([c0_c[jj], Sig0[ii], zeta_c[ii,jj], SigS_c[ii,jj],
                           np.nan, np.nan])
        else:
            sigm_c[ii,jj] = sigmax
            taum_c[ii,jj] = tau
            safe_c.append([c0_c[jj], Sig0[ii], zeta_c[ii,jj], SigS_c[ii,jj],
                           sigmax, tau])
        print(["Clay", (30*ii+jj+1)/30**2])
    np.savetxt("tau_max_clay.txt", safe_c)

safe2_c = []
sigm2_c = np.zeros((30,30))
taum2_c = np.zeros((30,30))

n = 2.5

for ii in range(30):
    for jj in range(30):
        lambdaD = np.sqrt((8.85e-12*80*1.38e-23*293)/
                          (2e3*1.602e-19*96485.33*c0_c[jj]))
        sig = sl.conductivity_sl(w, 1e3*c0_c[jj], 293, 80, 5e-8, 5e-9,
                                 5e-6/2, 5e-7/2, n*2*lambdaD, 2*lambdaD,
                                 zeta_c[ii,jj], SigS_c[ii,jj])
        sigmax, tau = findmax(w, sig.imag, 0)
        sigmaxh, tauh = findmax(w, sig.imag, 1)
        if tau == tauh:
            sigm2_c[ii, jj] = np.nan
            taum2_c[ii, jj] = np.nan
            safe2_c.append([c0_c[jj], Sig0[ii], zeta_c[ii,jj],
                           SigS_c[ii,jj], np.nan, np.nan])
        else:
            sigm2_c[ii, jj] = sigmax
            taum2_c[ii, jj] = tau
            safe2_c.append([c0_c[jj], Sig0[ii], zeta_c[ii,jj],
                           SigS_c[ii,jj], sigmax, tau])
        print(["Clay (small)", (30*ii+jj+1)/30**2])
    np.savetxt("tau_max_clay_small.txt", safe2_c)

safe_s = []
sigm_s = np.zeros((30,30))
taum_s = np.zeros((30,30))

for ii in range(30):
    for jj in range(30):
        sig = sl.conductivity_sl(w, 1e3*c0_s[jj], 293, 80, 5e-8, 5e-9, 9e-5/2,
                                 1e-5/2, 2e-6, 2e-7, zeta_s[ii,jj],
                                 SigS_s[ii,jj])
        sigmax, tau = findmax(w, sig.imag, 0)
        sigmaxh, tauh = findmax(w, sig.imag, 1)
        if tau == tauh:
            sigm_s[ii,jj] = np.nan
            taum_s[ii,jj] = np.nan
            safe_s.append([c0_s[jj], pH[ii], zeta_s[ii,jj], SigS_s[ii,jj],
                           np.nan, np.nan])
        else:
            sigm_s[ii,jj] = sigmax
            taum_s[ii,jj] = tau
            safe_s.append([c0_s[jj], pH[ii], zeta_s[ii,jj], SigS_s[ii,jj],
                           sigmax, tau])
        print(["Quartz", (30*ii+jj+1)/30**2])
    np.savetxt("tau_max_quartz.txt", safe_s)