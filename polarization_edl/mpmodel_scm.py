import numpy as np
import matplotlib.pyplot as plt
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

plt.rcParams["figure.figsize"] = (20,8)
plt.rcParams.update({'font.size': 20})

fig1, ax1 = plt.subplots(1,2)
c10 = ax1[0].pcolor(c0_c, Sig0, zeta_c, vmin=-0.17, vmax=0)
c11 = ax1[1].pcolor(c0_s, pH, zeta_s, vmin=-0.17, vmax=0)
ax1[0].set_xscale("log")
ax1[0].set_yscale("log")
ax1[1].set_xscale("log")
ax1[0].set_xlabel("NaCl concentration [mol/L]")
ax1[0].set_ylabel("Surface charge dens. [C/m$^2$]")
ax1[1].set_xlabel("NaCl concentration [mol/L]")
ax1[1].set_ylabel("pH [-]")
cax10 = fig1.colorbar(c10, ax=ax1[0])
cax11 = fig1.colorbar(c11, ax=ax1[1])
cax10.set_label("Zeta potential [V]")
cax11.set_label("Zeta potential [V]")
fig1.tight_layout()

fig2, ax2 = plt.subplots(1,2)
c20 = ax2[0].pcolor(c0_c, Sig0, np.log10(SigS_c), vmin=-5.5, vmax=0)
c21 = ax2[1].pcolor(c0_s, pH, np.log10(SigS_s), vmin=-5.5, vmax=0)
ax2[0].set_xscale("log")
ax2[0].set_yscale("log")
ax2[1].set_xscale("log")
ax2[0].set_xlabel("NaCl concentration [mol/L]")
ax2[0].set_ylabel("Surface charge dens. [C/m$^2$]")
ax2[1].set_xlabel("NaCl concentration [mol/L]")
ax2[1].set_ylabel("pH [-]")
cax20 = fig2.colorbar(c20, ax=ax2[0])
cax21 = fig2.colorbar(c21, ax=ax2[1])
cax20.set_label("Log. SL surf. charge dens. [C/m$^2$]")
cax21.set_label("Log. SL surf. charge dens. [C/m$^2$]")
fig2.tight_layout()


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
            safe_c.append([c0_c[ii], Sig0[jj], zeta_c[ii,jj], SigS_c[ii,jj],
                           np.nan, np.nan])
        else:
            sigm_c[ii,jj] = sigmax
            taum_c[ii,jj] = tau
            safe_c.append([c0_c[jj], Sig0[ii], zeta_c[ii,jj], SigS_c[ii,jj],
                           sigmax, tau])
        print(["Clay", (30*ii+jj+1)/30**2])
    np.savetxt("tau_max_clay.txt", safe_c)
        
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

fig3, ax3 = plt.subplots(1,2)
c30 = ax3[0].pcolor(c0_c, Sig0, np.log10(taum_c))
c31 = ax3[1].pcolor(c0_s, pH, np.log10(taum_s))
ax3[0].set_xscale("log")
ax3[0].set_yscale("log")
ax3[1].set_xscale("log")
ax3[0].set_xlabel("NaCl concentration [mol/L]")
ax3[0].set_ylabel("Surface charge dens. [C/m$^2$]")
ax3[1].set_xlabel("NaCl concentration [mol/L]")
ax3[1].set_ylabel("pH [-]")
ax3[0].set_title("(a)", loc="left", fontsize=36)
ax3[1].set_title("(b)", loc="left", fontsize=36)
cax30 = fig3.colorbar(c30, ax=ax3[0])
cax31 = fig3.colorbar(c31, ax=ax3[1])
cax30.set_label("Log. Relaxation time [s]")
cax31.set_label("Log. Relaxation time [s]")
fig3.tight_layout()