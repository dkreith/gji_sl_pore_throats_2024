import numpy as np
import matplotlib.pyplot as plt

def find_bound(x,y):

    xmin = np.min(x)
    xmax = np.max(x)

    ii_min = 0
    ii_max = len(x)
    for ii in range(len(y)):
        if y[ii] <= xmin and y[ii+1] > xmin:
            ii_min = ii
        elif y[ii-1] < xmax and y[ii] >= xmin:
            ii_max = ii
    
    return [ii_min, ii_max]

plt.rcParams.update({'font.size': 16})

sig0_n = np.loadtxt("spectra/num/spec_01_01_num.txt",
                    dtype=np.complex_)
sig0_s = np.loadtxt("spectra/sl/spec_01_01_sl.txt",
                    dtype=np.complex_)
sig0_m = np.loadtxt("spectra/mm/spec_01_01_mm.txt",
                    dtype=np.complex_)

plt.rcParams["figure.figsize"] = (10,10)

i0s, j0s = find_bound(sig0_n[:,0], sig0_s[:,0])
i0m, j0m = find_bound(sig0_n[:,0], sig0_m[:,0])

"""
fig1, ax1 = plt.subplots(2,1)
ax1[0].semilogx(sig0_n[:,0].real, sig0_n[:,1].real, "xk")
ax1[0].semilogx(sig0_s[i0s:j0s,0].real, sig0_s[i0s:j0s,1].real, "r-")
ax1[0].semilogx(sig0_m[i0m:j0m,0].real, sig0_m[i0m:j0m,1].real, "b--")
ax1[0].fill_between([2.5e1, 2e4], 0.09377, 0.10115, alpha=0.4, color="silver")
ax1[0].set_xlim([5e-3, 2e4])
ax1[0].set_ylim([0.09377, 0.10115])


ax1[1].loglog(sig0_n[:,0].real, sig0_n[:,1].imag, "xk")
ax1[1].loglog(sig0_s[i0s:j0s,0].real, sig0_s[i0s:j0s,1].imag, "r-")
ax1[1].loglog(sig0_m[i0m:j0m,0].real, sig0_m[i0m:j0m,1].imag, "b--")
ax1[1].plot([0.431138, 0.431138],[5e-7, 6.69586e-5], color="gold")
ax1[1].plot([5e-3, 0.431138],[6.69586e-5, 6.69586e-5], color="gold")
ax1[1].loglog(sig0_n[:,0].real, sig0_n[:,1].imag, "xk")
ax1[1].loglog(sig0_s[i0s:j0s,0].real, sig0_s[i0s:j0s,1].imag, "r-")
ax1[1].loglog(sig0_m[i0m:j0m,0].real, sig0_m[i0m:j0m,1].imag, "b--")
ax1[1].fill_between([2.5e1, 2e4], 5e-7, 4e-3, alpha=0.4, color="silver")
ax1[1].set_xlim([5e-3, 2e4])
ax1[1].set_ylim([5e-7,4e-3])

ax1[1].legend(["Num. model", "New anal. model", "Old anal. model"])

ax1[0].set_xlabel("Angular Frequency [rad/s]")
ax1[0].set_ylabel("Norm. Real Conductivity [-]")
ax1[1].set_xlabel("Angular Frequency [rad/s]")
ax1[1].set_ylabel("Norm. Imag. Conductivity [-]")
fig1.tight_layout()
"""

fig1, ax1 = plt.subplot_mosaic([["real", "norm"],["imag", "imag"]],
                               layout="tight")
ax1["real"].semilogx(sig0_s[i0s:j0s,0].real, sig0_s[i0s:j0s,1].real, "r-")
ax1["real"].semilogx(sig0_n[:,0].real, sig0_n[:,1].real, "xk")
ax1["real"].semilogx(sig0_m[i0m:j0m,0].real, sig0_m[i0m:j0m,1].real, "b--")
ax1["real"].fill_between([2.5e1, 2e4], 0.09277, 0.10215, alpha=0.4, color="silver")
ax1["real"].set_xlim([5e-3, 2e4])
ax1["real"].set_ylim([0.09277, 0.10215])

ax1["norm"].semilogx(sig0_s[i0s:j0s,0].real,
                     sig0_s[i0s:j0s,1].real/sig0_s[i0s,1].real, "r-")
ax1["norm"].semilogx(sig0_n[:,0].real, sig0_n[:,1].real/sig0_n[0,1].real, "xk")
ax1["norm"].semilogx(sig0_m[i0m:j0m,0].real,
                     sig0_m[i0m:j0m,1].real/sig0_m[i0m,1].real, "b--")
ax1["norm"].fill_between([2.5e1, 2e4], 0.9995, 1.0035, alpha=0.4, color="silver")
ax1["norm"].set_xlim([5e-3, 2e4])
ax1["norm"].set_ylim([0.9995, 1.0035])

ax1["imag"].loglog(sig0_s[i0s:j0s,0].real, sig0_s[i0s:j0s,1].imag, "r-")
ax1["imag"].loglog(sig0_n[:,0].real, sig0_n[:,1].imag, "xk")
ax1["imag"].loglog(sig0_m[i0m:j0m,0].real, sig0_m[i0m:j0m,1].imag, "b--")
ax1["imag"].plot([0.431138, 0.431138],[5e-7, 6.69586e-5], color="gold")
ax1["imag"].plot([5e-3, 0.431138],[6.69586e-5, 6.69586e-5], color="gold")
ax1["imag"].loglog(sig0_n[:,0].real, sig0_n[:,1].imag, "xk")
ax1["imag"].loglog(sig0_s[i0s:j0s,0].real, sig0_s[i0s:j0s,1].imag, "r-")
ax1["imag"].loglog(sig0_m[i0m:j0m,0].real, sig0_m[i0m:j0m,1].imag, "b--")
ax1["imag"].fill_between([2.5e1, 2e4], 5e-7, 4e-3, alpha=0.4, color="silver")
ax1["imag"].set_xlim([5e-3, 2e4])
ax1["imag"].set_ylim([5e-7,4e-3])

ax1["imag"].legend(["New semi-analytical model", "Numerical model",
                    "Bücker & Hördt (2013)"])

ax1["real"].set_title("(a)", loc="left")
ax1["norm"].set_title("(b)", loc="left")
ax1["imag"].set_title("(c)", loc="left")

ax1["real"].set_xlabel("$\u03C9$ [rad/s]")
ax1["real"].set_ylabel("$\u03C3'/\u03C3_0$ [-]")
ax1["norm"].set_xlabel("$\u03C9$ [rad/s]")
ax1["norm"].set_ylabel("$\u03C3'/\u03C3_{DC}$ [-]")
ax1["imag"].set_xlabel("$\u03C9$ [rad/s]")
ax1["imag"].set_ylabel("$\u03C3''/\u03C3_0$ [-]")
fig1.tight_layout()

fig1.savefig("Figure_05.png",dpi=300,bbox_inches="tight")

sig1_n = np.loadtxt("spectra/num/spec_23_04_num.txt",
                    dtype=np.complex_)
sig1_s = np.loadtxt("spectra/sl/spec_23_04_sl.txt",
                    dtype=np.complex_)
sig1_m = np.loadtxt("spectra/mm/spec_23_04_mm.txt",
                    dtype=np.complex_)

sig2_n = np.loadtxt("spectra/num/spec_23_21_num.txt",
                    dtype=np.complex_)
sig2_s = np.loadtxt("spectra/sl/spec_23_21_sl.txt",
                    dtype=np.complex_)
sig2_m = np.loadtxt("spectra/mm/spec_23_21_mm.txt",
                    dtype=np.complex_)

sig3_n = np.loadtxt("spectra/num/spec_06_04_num.txt",
                    dtype=np.complex_)
sig3_s = np.loadtxt("spectra/sl/spec_06_04_sl.txt",
                    dtype=np.complex_)
sig3_m = np.loadtxt("spectra/mm/spec_06_04_mm.txt",
                    dtype=np.complex_)

sig4_n = np.loadtxt("spectra/num/spec_06_21_num.txt",
                    dtype=np.complex_)
sig4_s = np.loadtxt("spectra/sl/spec_06_21_sl.txt",
                    dtype=np.complex_)
sig4_m = np.loadtxt("spectra/mm/spec_06_21_mm.txt",
                    dtype=np.complex_)

plt.rcParams["figure.figsize"] = (12,8)

i1s, j1s = find_bound(sig1_n[:,0], sig1_s[:,0])
i2s, j2s = find_bound(sig2_n[:,0], sig2_s[:,0])
i3s, j3s = find_bound(sig3_n[:,0], sig3_s[:,0])
i4s, j4s = find_bound(sig4_n[:,0], sig4_s[:,0])

i1m, j1m = find_bound(sig1_n[:,0], sig1_m[:,0])
i2m, j2m = find_bound(sig2_n[:,0], sig2_m[:,0])
i3m, j3m = find_bound(sig3_n[:,0], sig3_m[:,0])
i4m, j4m = find_bound(sig4_n[:,0], sig4_m[:,0])

fig2, ax2 = plt.subplots(2,2)

ax2[0,0].loglog(sig1_s[i1s:j1s,0].real, sig1_s[i1s:j1s,1].imag, "r-")
ax2[0,0].loglog(sig1_n[:,0].real, sig1_n[:,1].imag, "xk")
ax2[0,0].loglog(sig1_m[i1m:j1m,0].real, sig1_m[i1m:j1m,1].imag, "b--")

ax2[0,1].loglog(sig2_s[i2s:j2s,0].real, sig2_s[i2s:j2s,1].imag, "r-")
ax2[0,1].loglog(sig2_n[:,0].real, sig2_n[:,1].imag, "xk")
ax2[0,1].loglog(sig2_m[i2m:j2m,0].real, sig2_m[i2m:j2m,1].imag, "b--")

ax2[1,0].loglog(sig3_s[i3s:j3s,0].real, sig3_s[i3s:j3s,1].imag, "r-")
ax2[1,0].loglog(sig3_n[:,0].real, sig3_n[:,1].imag, "xk")
ax2[1,0].loglog(sig3_m[i3m:j3m,0].real, sig3_m[i3m:j3m,1].imag, "b--")

ax2[1,1].loglog(sig4_s[i4s:j4s,0].real, sig4_s[i4s:j4s,1].imag, "r-")
ax2[1,1].loglog(sig4_n[:,0].real, sig4_n[:,1].imag, "xk")
ax2[1,1].loglog(sig4_m[i4m:j4m,0].real, sig4_m[i4m:j4m,1].imag, "b--")

ax2[0,1].legend(["New semi-analytical model", "Numerical model",
                 "Bücker & Hördt (2013)"],
                loc=4)

ax2[0,0].set_title("(a) Spectrum 1", loc="left")
ax2[0,1].set_title("(b) Spectrum 2", loc="left")
ax2[1,0].set_title("(c) Spectrum 3", loc="left")
ax2[1,1].set_title("(d) Spectrum 4", loc="left")

ax2[0,0].set_ylabel("$\u03C3''/\u03C3_0$ [-]")
ax2[1,0].set_ylabel("$\u03C3''/\u03C3_0$ [-]")
ax2[1,0].set_xlabel("$\u03C9$ [rad/s]")
ax2[1,1].set_xlabel("$\u03C9$ [rad/s]")
fig2.tight_layout()

fig2.savefig("Figure_06.pdf",dpi=300,bbox_inches="tight")