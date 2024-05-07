import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 13})

dat_num = np.loadtxt("tau_max_num.txt")
N_num = 26
zeta_num = np.zeros((N_num, N_num))
SigS_num = np.zeros((N_num, N_num))
taum_num = np.zeros((N_num, N_num))
sigm_num = np.zeros((N_num, N_num))
for ii in range(N_num):
    for jj in range(N_num):
        zeta_num[ii][jj] = dat_num[ii*N_num+jj][0]
        SigS_num[ii][jj] = dat_num[ii*N_num+jj][1]
        sigm_num[ii][jj] = dat_num[ii*N_num+jj][2]
        taum_num[ii][jj] = dat_num[ii*N_num+jj][3]
        
dat_sl = np.loadtxt("tau_max_sl.txt")
N_sl = 26
zeta_sl = np.zeros((N_sl, N_sl))
SigS_sl = np.zeros((N_sl, N_sl))
taum_sl = np.zeros((N_sl, N_sl))
sigm_sl = np.zeros((N_sl, N_sl))
for ii in range(N_sl):
    for jj in range(N_sl):
        zeta_sl[ii][jj] = dat_sl[ii*N_sl+jj][0]
        SigS_sl[ii][jj] = dat_sl[ii*N_sl+jj][1]
        sigm_sl[ii][jj] = dat_sl[ii*N_sl+jj][2]
        taum_sl[ii][jj] = dat_sl[ii*N_sl+jj][3]

dat_mm = np.loadtxt("tau_max_mm.txt")
N_mm = 26
zeta_mm = np.zeros((N_mm, N_mm))
SigS_mm = np.zeros((N_mm, N_mm))
taum_mm = np.zeros((N_mm, N_mm))
sigm_mm = np.zeros((N_mm, N_mm))
for ii in range(N_mm):
    for jj in range(N_mm):
        zeta_mm[ii][jj] = dat_mm[ii*N_mm+jj][0]
        SigS_mm[ii][jj] = dat_mm[ii*N_mm+jj][1]
        sigm_mm[ii][jj] = dat_mm[ii*N_mm+jj][2]
        taum_mm[ii][jj] = dat_mm[ii*N_mm+jj][3]

zeta = np.linspace(-0.025,-0.2,26)
SigS = np.logspace(-3,-1,26)
zeta_num = np.linspace(-0.025,-0.2,26)
zeta_sl = np.linspace(-0.025,-0.2,26)
zeta_mm = np.linspace(-0.025,-0.2,26)
SigS_num = np.logspace(-3,-1,26)
SigS_sl = np.logspace(-3,-1,26)
SigS_mm = np.logspace(-3,-1,26)

fig, [[ax1, ax2], [ax3, ax4], [ax5, ax6]]  = plt.subplots(3,2)
fig.set_figwidth(12)
fig.set_figheight(15)

cmap = "viridis"
shading = "nearest"

pl1 = ax1.pcolor(zeta_num, SigS_num, np.log10(sigm_num),
                     vmin=-4.5, vmax=-1.5,
                     shading=shading, cmap=cmap)
pl2 = ax2.pcolor(zeta_num, SigS_num, np.log10(taum_num),
                     vmin=-1.1, vmax=0.3,
                     shading=shading, cmap=cmap)
pl3 = ax3.pcolor(zeta_num, SigS_num, np.log10(sigm_sl),
                     vmin=-4.5, vmax=-1.5,
                     shading=shading, cmap=cmap)
pl4 = ax4.pcolor(zeta_num, SigS_num, np.log10(taum_sl),
                     vmin=-1.1, vmax=0.3,
                     shading=shading, cmap=cmap)
pl5 = ax5.pcolor(zeta_num, SigS_num, np.log10(sigm_mm),
                     vmin=-4.5, vmax=-1.5,
                     shading=shading, cmap=cmap)
pl6 = ax6.pcolor(zeta_num, SigS_num, np.log10(taum_mm),
                     vmin=-1.1, vmax=0.3,
                     shading=shading, cmap=cmap)

ax1.plot(zeta[0], SigS[0], "*", c="orangered", markersize=12)
ax1.plot(zeta[3], SigS[22], "o", c="white")
ax1.plot(zeta[20], SigS[22], "o", c="white")
ax1.plot(zeta[3], SigS[5], "o", c="white")
ax1.plot(zeta[20], SigS[5], "o", c="white")

ax2.plot(zeta[0], SigS[0], "*", c="orangered", markersize=12)
ax2.plot(zeta[3], SigS[22], "o", c="white")
ax2.plot(zeta[20], SigS[22], "o", c="white")
ax2.plot(zeta[3], SigS[5], "o", c="white")
ax2.plot(zeta[20], SigS[5], "o", c="white")

ax3.plot(zeta[0], SigS[0], "*", c="orangered", markersize=12)
ax3.plot(zeta[3], SigS[22], "o", c="white")
ax3.plot(zeta[20], SigS[22], "o", c="white")
ax3.plot(zeta[3], SigS[5], "o", c="white")
ax3.plot(zeta[20], SigS[5], "o", c="white")

ax4.plot(zeta[0], SigS[0], "*", c="orangered", markersize=12)
ax4.plot(zeta[3], SigS[22], "o", c="white")
ax4.plot(zeta[20], SigS[22], "o", c="white")
ax4.plot(zeta[3], SigS[5], "o", c="white")
ax4.plot(zeta[20], SigS[5], "o", c="white")

ax5.plot(zeta[0], SigS[0], "*", c="orangered", markersize=12)
ax5.plot(zeta[3], SigS[22], "o", c="white")
ax5.plot(zeta[20], SigS[22], "o", c="white")
ax5.plot(zeta[3], SigS[5], "o", c="white")
ax5.plot(zeta[20], SigS[5], "o", c="white")

ax6.plot(zeta[0], SigS[0], "*", c="orangered", markersize=12)
ax6.plot(zeta[3], SigS[22], "o", c="white")
ax6.plot(zeta[20], SigS[22], "o", c="white")
ax6.plot(zeta[3], SigS[5], "o", c="white")
ax6.plot(zeta[20], SigS[5], "o", c="white")

ax1.set_xlim(-0.0209, -0.2037)
ax2.set_xlim(-0.0209, -0.2037)
ax3.set_xlim(-0.0209, -0.2037)
ax4.set_xlim(-0.0209, -0.2037)
ax5.set_xlim(-0.0209, -0.2037)
ax6.set_xlim(-0.0209, -0.2037)

ax1.set_yscale("log")
ax2.set_yscale("log")
ax3.set_yscale("log")
ax4.set_yscale("log")
ax5.set_yscale("log")
ax6.set_yscale("log")

ax1.set_title("(a)", loc="left")
ax2.set_title("(b)", loc="left")
ax3.set_title("(c)", loc="left")
ax4.set_title("(d)", loc="left")
ax5.set_title("(e)", loc="left")
ax6.set_title("(f)", loc="left")

ax1.set_xticklabels([])
ax2.set_xticklabels([])
ax2.set_yticklabels([])
ax3.set_xticklabels([])
ax4.set_xticklabels([])
ax4.set_yticklabels([])
ax6.set_yticklabels([])

"""
ax1.set_ylabel("SL Surf.-Charge Dens. [C/m$^2$]")
ax3.set_ylabel("SL Surf.-Charge Dens. [C/m$^2$]")
ax5.set_ylabel("SL Surf.-Charge Dens. [C/m$^2$]")
ax5.set_xlabel("d-Plane Pot. [V]")
ax6.set_xlabel("d-Plane Pot. [V]")
"""
ax1.set_ylabel("$\u03A3_S$ [C/m$^2$]")
ax3.set_ylabel("$\u03A3_S$ [C/m$^2$]")
ax5.set_ylabel("$\u03A3_S$ [C/m$^2$]")
ax5.set_xlabel("$\u03C6_d$ [V]")
ax6.set_xlabel("$\u03C6_d$ [V]")

cbs = fig.colorbar(pl1, ax=[ax1, ax3, ax5], shrink=0.5)
cbt = fig.colorbar(pl2, ax=[ax2, ax4, ax6], shrink=0.5)
#cbs.set_label('Log. Max. Norm. Imag. Cond. [-]')
#cbt.set_label('Log. Relaxation Time [s]')
cbs.set_label("log($\u03C3''_{max}/\u03C3_0$) [-]")
cbt.set_label("log($\u03C4$) [s]")

fig.savefig("Figure_07.png",dpi=300,bbox_inches="tight")

fig2, [[ax7, ax8],[ax9, ax0]] = plt.subplots(2,2)
fig2.set_figwidth(10)
fig2.set_figheight(8)

pl7 = ax7.pcolor(zeta_num, SigS_num, (sigm_sl-sigm_num)/sigm_num,
                 cmap="seismic", vmin=-5, vmax=5, shading=shading)
pl8 = ax8.pcolor(zeta_num, SigS_num, np.log10(taum_sl)-np.log10(taum_num),
                 cmap="seismic", vmin=-0.8, vmax=0.8, shading=shading)
pl9 = ax9.pcolor(zeta_num, SigS_num, (sigm_mm-sigm_num)/sigm_num,
                 cmap="seismic", vmin=-5, vmax=5, shading=shading)
pl0 = ax0.pcolor(zeta_num, SigS_num, np.log10(taum_mm)-np.log10(taum_num),
                 cmap="seismic", vmin=-0.8, vmax=0.8, shading=shading)

ax7.set_xlim(-0.0209, -0.2037)
ax8.set_xlim(-0.0209, -0.2037)
ax9.set_xlim(-0.0209, -0.2037)
ax0.set_xlim(-0.0209, -0.2037)

ax7.set_yscale("log")
ax8.set_yscale("log")
ax9.set_yscale("log")
ax0.set_yscale("log")

ax7.set_title("(a)", loc="left")
ax8.set_title("(b)", loc="left")
ax9.set_title("(c)", loc="left")
ax0.set_title("(d)", loc="left")

"""
ax7.set_ylabel("SL Surf.-Charge Dens. [C/m$^2$]")
ax9.set_ylabel("SL Surf.-Charge Dens. [C/m$^2$]")
ax9.set_xlabel("d-Plane Pot. [V]")
ax0.set_xlabel("d-Plane Pot. [V]")
"""

ax7.set_ylabel("$\u03A3_S$ [C/m$^2$]")
ax9.set_ylabel("$\u03A3_S$ [C/m$^2$]")
ax9.set_xlabel("$\u03C6_d$ [V]")
ax0.set_xlabel("$\u03C6_d$ [V]")

cbs2 = fig2.colorbar(pl7, ax=[ax7, ax9], shrink=0.5)
cbt2 = fig2.colorbar(pl8, ax=[ax8, ax0], shrink=0.5)
#cbs2.set_label('Rel. Diff. of Max. Imag. Cond. [-]')
#cbt2.set_label('Log. Diff. of Rel. Time [s]')
cbs2.set_label("$\u0394_{rel}\u03C3''_{max}$ [-]")
cbt2.set_label("$\u0394_{log}\u03C4$ [s]")

ax7.set_xticklabels([])
ax8.set_xticklabels([])
ax8.set_yticklabels([])
ax0.set_yticklabels([])

fig2.savefig("Figure_08.png",dpi=300,bbox_inches="tight")