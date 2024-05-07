import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 16})
plt.rcParams["figure.figsize"] = (7,10)
col = plt.cm.viridis(np.linspace(0,0.95,3))

mtm_100_c = np.loadtxt("phreeqc/montmorillonite_Sig0_-100_c.txt", skiprows=5)
mtm_010_c = np.loadtxt("phreeqc/montmorillonite_Sig0_-125_c.txt", skiprows=5)
mtm_001_c = np.loadtxt("phreeqc/montmorillonite_Sig0_-150_c.txt", skiprows=5)

mtm_100_p = np.loadtxt("phreeqc/montmorillonite_Sig0_-100_p.txt")
mtm_010_p = np.loadtxt("phreeqc/montmorillonite_Sig0_-125_p.txt")
mtm_001_p = np.loadtxt("phreeqc/montmorillonite_Sig0_-150_p.txt")

fig1, ax1 = plt.subplots(3,1)
cc1, = ax1[0].semilogx(mtm_001_c[:,0],mtm_001_c[:,1], color=col[0])
cc2, = ax1[0].semilogx(mtm_010_c[:,0],mtm_010_c[:,1], color=col[1])
cc3, = ax1[0].semilogx(mtm_100_c[:,0],mtm_100_c[:,1], color=col[2])
cp1, = ax1[0].semilogx(mtm_001_p[:,0],mtm_001_p[:,1], "o", color=col[0])
cp2, = ax1[0].semilogx(mtm_010_p[:,0],mtm_010_p[:,1], "o", color=col[1])
cp3, = ax1[0].semilogx(mtm_100_p[:,0],mtm_100_p[:,1], "o", color=col[2])

ax1[1].loglog(mtm_001_c[:,0],mtm_001_c[:,2], color=col[0])
ax1[1].loglog(mtm_010_c[:,0],mtm_010_c[:,2], color=col[1])
ax1[1].loglog(mtm_100_c[:,0],mtm_100_c[:,2], color=col[2])
ax1[1].loglog(mtm_001_p[:,0],mtm_001_p[:,2], "o", color=col[0])
ax1[1].loglog(mtm_010_p[:,0],mtm_010_p[:,2], "o", color=col[1])
ax1[1].loglog(mtm_100_p[:,0],mtm_100_p[:,2], "o", color=col[2])

ax1[2].semilogx(mtm_001_c[:,0],mtm_001_c[:,3], color=col[0])
ax1[2].semilogx(mtm_010_c[:,0],mtm_010_c[:,3], color=col[1])
ax1[2].semilogx(mtm_100_c[:,0],mtm_100_c[:,3], color=col[2])
ax1[2].semilogx(mtm_001_p[:,0],mtm_001_p[:,3], "o", color=col[0])
ax1[2].semilogx(mtm_010_p[:,0],mtm_010_p[:,3], "o", color=col[1])
ax1[2].semilogx(mtm_100_p[:,0],mtm_100_p[:,3], "o", color=col[2])

ax1[0].set_ylabel("$\u03C6_d$ [V]")
ax1[1].set_ylabel("$\u03A3_S$ [C/m$^2$]")
ax1[2].set_ylabel("$f_Q$ [-]")
ax1[2].set_xlabel("$c_0$ [mol/L]")

ax1[0].legend([(cc1, cp1), (cc2, cp2), (cc3, cp3)],
              ["$\u03A3_0$ = -100 mC/m$^2$", "$\u03A3_0$ = -125 mC/m$^2$",
               "$\u03A3_0$ = -150 mC/m$^2$"])
ax1[0].set_ylim([-0.24, 0])

fig1.tight_layout()

fig1.savefig("Figure_02.pdf",dpi=300,bbox_inches="tight")

qtz_100_c = np.loadtxt("phreeqc/quartz_c0_1e-1_c.txt", skiprows=5)
qtz_010_c = np.loadtxt("phreeqc/quartz_c0_1e-2_c.txt", skiprows=5)
qtz_001_c = np.loadtxt("phreeqc/quartz_c0_1e-3_c.txt", skiprows=5)

qtz_100_p = np.loadtxt("phreeqc/quartz_c0_1e-1_p.txt")
qtz_010_p = np.loadtxt("phreeqc/quartz_c0_1e-2_p.txt")
qtz_001_p = np.loadtxt("phreeqc/quartz_c0_1e-3_p.txt")


fig2, ax2 = plt.subplots(3,1)
qc1, = ax2[0].plot(qtz_001_c[:,0],qtz_001_c[:,1], color=col[0])
qc2, = ax2[0].plot(qtz_010_c[:,0],qtz_010_c[:,1], color=col[1])
qc3, = ax2[0].plot(qtz_100_c[:,0],qtz_100_c[:,1], color=col[2])
qp1, = ax2[0].plot(qtz_001_p[:,0],qtz_001_p[:,1], "o", color=col[0])
qp2, = ax2[0].plot(qtz_010_p[:,0],qtz_010_p[:,1], "o", color=col[1])
qp3, = ax2[0].plot(qtz_100_p[:,0],qtz_100_p[:,1], "o", color=col[2])

ax2[1].semilogy(qtz_001_c[:,0],qtz_001_c[:,2], color=col[0])
ax2[1].semilogy(qtz_010_c[:,0],qtz_010_c[:,2], color=col[1])
ax2[1].semilogy(qtz_100_c[:,0],qtz_100_c[:,2], color=col[2])
ax2[1].semilogy(qtz_001_p[:,0],qtz_001_p[:,2], "o", color=col[0])
ax2[1].semilogy(qtz_010_p[:,0],qtz_010_p[:,2], "o", color=col[1])
ax2[1].semilogy(qtz_100_p[:,0],qtz_100_p[:,2], "o", color=col[2])

ax2[2].plot(qtz_001_c[:,0],qtz_001_c[:,3], color=col[0])
ax2[2].plot(qtz_010_c[:,0],qtz_010_c[:,3], color=col[1])
ax2[2].plot(qtz_100_c[:,0],qtz_100_c[:,3], color=col[2])
ax2[2].plot(qtz_001_p[:,0],qtz_001_p[:,3], "o", color=col[0])
ax2[2].plot(qtz_010_p[:,0],qtz_010_p[:,3], "o", color=col[1])
ax2[2].plot(qtz_100_p[:,0],qtz_100_p[:,3], "o", color=col[2])

ax2[0].set_ylabel("$\u03C6_d$ [V]")
ax2[1].set_ylabel("$\u03A3_S$ [C/m$^2$]")
ax2[2].set_ylabel("$f_Q$ [-]")
ax2[2].set_xlabel("pH [-]")

ax2[0].legend([(qc1, qp1), (qc2, qp2), (qc3, qp3)],
              ["$c_0$ = 1e-3 mol/L", "$c_0$ = 1e-2 mol/L", "$c_0$ = 1e-1 mol/L"])
ax2[0].set_ylim([-0.19,0])

fig2.tight_layout()

fig2.savefig("Figure_03.pdf",dpi=300,bbox_inches="tight")