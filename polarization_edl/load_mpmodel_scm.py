import numpy as np
import matplotlib.pyplot as plt
import conductivity_sl as sl
from findmax import findmax
# subplot_mosaic statt subfigures

conc_c = np.zeros((30, 30))
Sig0_c = np.zeros((30, 30))
sigm_c = np.zeros((30, 30))
taum_c = np.zeros((30, 30))
zeta_c = np.zeros((30, 30))
SigS_c = np.zeros((30, 30))

dat_c = np.loadtxt("tau_max_clay.txt")
for ii in range(30):
    for jj in range(30):
        conc_c[ii][jj] = dat_c[jj+ii*30][0]
        Sig0_c[ii][jj] = dat_c[jj+ii*30][1]
        zeta_c[ii][jj] = dat_c[jj+ii*30][2]
        SigS_c[ii][jj] = dat_c[jj+ii*30][3]
        sigm_c[ii][jj] = dat_c[jj+ii*30][4]
        taum_c[ii][jj] = dat_c[jj+ii*30][5]

conc_s = np.zeros((30, 30))
pH_s = np.zeros((30, 30))
sigm_s = np.zeros((30, 30))
taum_s = np.zeros((30, 30))
zeta_s = np.zeros((30, 30))
SigS_s = np.zeros((30, 30))

dat_s = np.loadtxt("tau_max_quartz.txt")
for ii in range(30):
    for jj in range(30):
        conc_s[ii][jj] = dat_s[jj+ii*30][0]
        pH_s[ii][jj] = dat_s[jj+ii*30][1]
        zeta_s[ii][jj] = dat_s[jj+ii*30][2]
        SigS_s[ii][jj] = dat_s[jj+ii*30][3]
        sigm_s[ii][jj] = dat_s[jj+ii*30][4]
        taum_s[ii][jj] = dat_s[jj+ii*30][5]

plt.rcParams["figure.figsize"] = (15,12)
plt.rcParams.update({"font.size": 18})
cmap = "viridis"

fig, ax = plt.subplots(2,2)
c00 = ax[0,0].pcolor(1e3*conc_c, Sig0_c, np.log10(sigm_c), cmap=cmap)
cp00z = ax[0,0].contour(1e3*conc_c, Sig0_c, zeta_c, 4, linestyles="solid", colors="k")
cp00s = ax[0,0].contour(1e3*conc_c, Sig0_c, SigS_c, 4, linestyles="dashed", colors="w")
c01 = ax[0,1].pcolor(1e3*conc_s, pH_s, np.log10(sigm_s), cmap=cmap)
cp01z = ax[0,1].contour(1e3*conc_s, pH_s, zeta_s, 4, linestyles="solid", colors="k")
#cp01s = ax[0,1].contour(conc_s, pH_s, SigS_s, [2e-6,2e-5,2e-4,2e-3,2e-2,2e-1], linestyles="dashed", colors="w")
cp01s = ax[0,1].contour(1e3*conc_s, pH_s, np.log10(SigS_s/2), [-6,-5,-4,-3,-2,-1], linestyles="dashed", colors="w")

c10 = ax[1,0].pcolor(1e3*conc_c, Sig0_c, np.log10(taum_c), cmap=cmap)
cp10z = ax[1,0].contour(1e3*conc_c, Sig0_c, zeta_c, 4, linestyles="solid", colors="k")
cp10s = ax[1,0].contour(1e3*conc_c, Sig0_c, SigS_c, 4, linestyles="dashed", colors="w")
c11 = ax[1,1].pcolor(1e3*conc_s, pH_s, np.log10(taum_s), cmap=cmap)
cp11z = ax[1,1].contour(1e3*conc_s, pH_s, zeta_s, 4, linestyles="solid", colors="k")
#cp11s = ax[1,1].contour(conc_s, pH_s, SigS_s, [2e-6,2e-5,2e-4,2e-3,2e-2,2e-1], linestyles="dashed", colors="w")
cp11s = ax[1,1].contour(1e3*conc_s, pH_s, np.log10(SigS_s/2), [-6,-5,-4,-3,-2,-1], linestyles="dashed", colors="w")

ax[0,0].set_xscale("log")
#ax[0,0].set_yscale("log")
ax[0,1].set_xscale("log")
ax[1,0].set_xscale("log")
#ax[1,0].set_yscale("log")
ax[1,1].set_xscale("log")

"""
ax[0,0].set_xlabel("NaCl-Concentration [mol/L]")
ax[0,0].set_ylabel("Abs. Surface Charge Dens. [C/m$^2$]")
ax[0,1].set_xlabel("NaCl-Concentration [mol/L]")
ax[0,1].set_ylabel("pH [-]")
ax[1,0].set_xlabel("NaCl-Concentration [mol/L]")
ax[1,0].set_ylabel("Abs. Surface Charge Dens. [C/m$^2$]")
ax[1,1].set_xlabel("NaCl-Concentration [mol/L]")
ax[1,1].set_ylabel("pH [-]")
"""

ax[0,0].set_xlabel("$c_\infty$ [mol/m$^3$]")
ax[0,0].set_ylabel("|$\u03A3_0$| [C/m$^2$]")
ax[0,1].set_xlabel("$c_\infty$ [mol/m$^3$]")
ax[0,1].set_ylabel("$pH$ [-]")
ax[1,0].set_xlabel("$c_\infty$ [mol/m$^3$]")
ax[1,0].set_ylabel("|$\u03A3_0$| [C/m$^2$]")
ax[1,1].set_xlabel("$c_\infty$ [mol/m$^3$]")
ax[1,1].set_ylabel("$pH$ [-]")

ax[0,0].set_title("(a)", loc="left", fontsize=22)
ax[0,1].set_title("(b)", loc="left", fontsize=22)
ax[1,0].set_title("(c)", loc="left", fontsize=22)
ax[1,1].set_title("(d)", loc="left", fontsize=22)

cax00 = fig.colorbar(c00, ax=ax[0,0])
cax01 = fig.colorbar(c01, ax=ax[0,1])
#cax00.set_label("Log. Max. Norm. Imag. Cond. [-]")
#cax01.set_label("Log. Max. Norm. Imag. Cond. [-]")
cax00.set_label("log($\u03C3''_{max}/\u03C3_0$) [-]")
cax01.set_label("log($\u03C3''_{max}/\u03C3_0$) [-]")
cax10 = fig.colorbar(c10, ax=ax[1,0])
cax11 = fig.colorbar(c11, ax=ax[1,1])
#cax10.set_label("Log. Relaxation time [s]")
#cax11.set_label("Log. Relaxation time [s]")
cax10.set_label("log($\u03C4$) [s]")
cax11.set_label("log($\u03C4$) [s]")

ax[0,0].clabel(cp00z, fontsize=15, inline=True)
ax[0,0].clabel(cp00s, fontsize=15, inline=True)
ax[0,1].clabel(cp01z, fontsize=15, inline=True)
ax[0,1].clabel(cp01s, fontsize=15, inline=True, fmt=r"$2 \times 10^{%1i}$")
ax[1,0].clabel(cp10z, fontsize=15, inline=True)
ax[1,0].clabel(cp10s, fontsize=15, inline=True)
ax[1,1].clabel(cp11z, fontsize=15, inline=True)
ax[1,1].clabel(cp11s, fontsize=15, inline=True, fmt=r"$2 \times 10^{%1i}$")

fig.tight_layout()

fig.savefig("Figure_09.png",dpi=300,bbox_inches="tight")

w = np.logspace(-2, 2)

col = plt.cm.viridis(np.linspace(0.1, 1, 7))
col_ind = 0

wmax = []

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (10, 8)

fig1, ax1 = plt.subplots()
#for ii in [11, 13, 15, 17, 19, 21, 23]:
for ii in [15, 17, 19, 21, 23, 25, 27]:
    #sig = sl.conductivity_sl(w, 1e3*conc_s[0, ii], 293, 80, 5e-8, 5e-9, 9e-5/2,
    #                         1e-5/2, 2e-6, 2e-7, zeta_s[0, ii],
    #                         SigS_s[0, ii])
    sig = sl.conductivity_sl(w, 1e3*conc_s[4, ii], 293, 80, 5e-8, 5e-9, 9e-5/2,
                             1e-5/2, 2e-6, 2e-7, zeta_s[4, ii],
                             SigS_s[4, ii])
    sigmax, tau = findmax(w, sig.imag, 0)
    wmax.append([1/tau, sigmax])
    ax1.loglog(w, sig.imag, color=col[col_ind], linewidth=2.5)
    col_ind = col_ind+1
wmax = np.array(wmax)
ax1.loglog(wmax[:, 0], wmax[:, 1], "ok", fillstyle="none", markeredgewidth=1.5,
           markersize=8)
#ax1.set_xlabel("Angular Frequency [rad/s]")
#ax1.set_ylabel("Norm. Imag. Conductivity [-]")
ax1.set_xlabel("$\u03C9$ [rad/s]")
ax1.set_ylabel("$\u03C3''/\u03C3_0$[-]")
ax1.legend(["$c_\infty = 3.56$ mol/m$^3$", "$c_\infty = 5.74$ mol/m$^3$",
            "$c_\infty = 9.24$ mol/m$^3$", "$c_\infty = 14.87$ mol/m$^3$",
            "$c_\infty = 23.95$ mol/m$^3$", "$c_\infty = 38.57$ mol/m$^3$",
            "$c_\infty = 62.10$ mol/m$^3$"])
fig1.tight_layout()

fig1.savefig("Figure_10.pdf",dpi=300,bbox_inches="tight")
