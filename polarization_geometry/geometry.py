import numpy as np
import conductivity_sl as sl
import conductivity_mm as mm
import findmax

def zeta2sigma(zeta, c0, epsr, T, mod):
    
    # Constants
    kB = 8.6173e-5  # Boltzmanm's constant [eV/K]
    F = 96485  # Faraday's constant [C/mol]
    eps0 = 8.85e-12  # vacuum permittivity [C/Vm]
    
    # Inverse Debye length
    k = np.sqrt(2*c0*F/epsr/eps0/kB/T)
    
    # Calculate diffuse layer surface charge density
    if mod == 'd':
        SigD = -4*F*c0*np.sinh(zeta/2/kB/T)/k
    elif mod == 'd+':
        SigD = (2*F*c0/k)*(np.exp(-zeta/(2*kB*T))-1)
    
    return SigD

# Parameters
c0 = 1
T = 293
mu = 5e-8
muS = 5e-9
epsr = 80
zeta = -0.075
SigS = zeta2sigma(zeta, c0, epsr, T, 'd') # SigS = SigD --> fQ = 0.5

w = np.logspace(-4,8)

step = 100

L1 = np.logspace(-7,-3,step)
L2 = np.logspace(-8,-4,step)
r1x = np.logspace(-8,-4,step)
r1y = np.logspace(-7.5,-4,step)
r2 = np.logspace(-9,-5,step)

L1_mesh = np.zeros((step,step))
L2_mesh = np.zeros((step,step))
r1x_mesh = np.zeros((step,step))
r1y_mesh = np.zeros((step,step))
r2_mesh = np.zeros((step,step))

sig_L1_L2_sl = np.zeros((step,step))
tau_L1_L2_sl = np.zeros((step,step))
sig_L1_L2_mm = np.zeros((step,step))
tau_L1_L2_mm = np.zeros((step,step))

sig_L1_r1_sl = np.zeros((step,step))
tau_L1_r1_sl = np.zeros((step,step))
sig_L1_r1_mm = np.zeros((step,step))
tau_L1_r1_mm = np.zeros((step,step))

sig_r1_r2_sl = np.zeros((step,step))
tau_r1_r2_sl = np.zeros((step,step))
sig_r1_r2_mm = np.zeros((step,step))
tau_r1_r2_mm = np.zeros((step,step))

for ii in range(step):
    for jj in range(step):
        sigma_sl = sl.conductivity_sl(w, c0, T, epsr, mu, muS, L1[ii], L2[jj],
                                    100e-9, 10e-9, zeta, SigS)
        sigma_mm = mm.conductivity_mm(w, T, c0, mu, epsr, L1[ii], L2[jj],
                                      100e-9, 10e-9, zeta, SigS)
        sig_l_sl, tau_l_sl = findmax.findmax(w,sigma_sl.imag,0)
        sig_h_sl, tau_h_sl = findmax.findmax(w,sigma_sl.imag,1)
        sig_mm, tau_mm = findmax.findmax(w,sigma_mm.imag,0)
        if sig_l_sl == sig_h_sl:
            sig_L1_L2_sl[ii][jj] = np.nan
            tau_L1_L2_sl[ii][jj] = np.nan
        else:
            sig_L1_L2_sl[ii][jj] = sig_l_sl
            tau_L1_L2_sl[ii][jj] = tau_l_sl
        sig_L1_L2_mm[ii][jj] = sig_mm
        tau_L1_L2_mm[ii][jj] = tau_mm
        L1_mesh[ii][jj] = L1[ii]
        L2_mesh[ii][jj] = L2[jj]
        print(["L1/L2", (100*ii+jj+1)/100**2])
    np.savetxt('geometry/L1_L2/SL_sig_L1_L2.txt', sig_L1_L2_sl) 
    np.savetxt('geometry/L1_L2/SL_tau_L1_L2.txt', tau_L1_L2_sl)
    np.savetxt('geometry/L1_L2/MM_sig_L1_L2.txt', sig_L1_L2_mm) 
    np.savetxt('geometry/L1_L2/MM_tau_L1_L2.txt', tau_L1_L2_mm)
    np.savetxt('geometry/L1_L2/L1_mesh.txt', L1_mesh)
    np.savetxt('geometry/L1_L2/L2_mesh.txt', L2_mesh)

for ii in range(step):
    for jj in range(step):
        sigma_sl = sl.conductivity_sl(w, c0, T, epsr, mu, muS, L1[ii], 50e-6,
                                    r1y[jj], 30e-9, zeta, SigS)
        sigma_mm = mm.conductivity_mm(w, T, c0, mu, epsr, L1[ii], 50e-6,
                                      r1y[jj], 30e-9, zeta, SigS)
        sig_l_sl, tau_l_sl = findmax.findmax(w,sigma_sl.imag,0)
        sig_h_sl, tau_h_sl = findmax.findmax(w,sigma_sl.imag,1)
        sig_mm, tau_mm = findmax.findmax(w,sigma_mm.imag,0)
        if sig_l_sl == sig_h_sl:
            sig_L1_r1_sl[ii][jj] = np.nan
            tau_L1_r1_sl[ii][jj] = np.nan
        else:
            sig_L1_r1_sl[ii][jj] = sig_l_sl
            tau_L1_r1_sl[ii][jj] = tau_l_sl
        sig_L1_r1_mm[ii][jj] = sig_mm
        tau_L1_r1_mm[ii][jj] = tau_mm
        r1y_mesh[ii][jj] = r1y[jj]
        print(["L1/r1", (100*ii+jj+1)/100**2])
    np.savetxt('geometry/L1_r1/SL_sig_L1_r1.txt', sig_L1_r1_sl) 
    np.savetxt('geometry/L1_r1/SL_tau_L1_r1.txt', tau_L1_r1_sl)
    np.savetxt('geometry/L1_r1/MM_sig_L1_r1.txt', sig_L1_r1_mm) 
    np.savetxt('geometry/L1_r1/MM_tau_L1_r1.txt', tau_L1_r1_mm)
    np.savetxt('geometry/L1_r1/L1_mesh.txt', L1_mesh)
    np.savetxt('geometry/L1_r1/r1y_mesh.txt', r1y_mesh)


for ii in range(step):
    for jj in range(step):
        if r1x[ii] < r2[jj]:
            sig_r1_r2_sl[ii][jj] = np.nan
            tau_r1_r2_sl[ii][jj] = np.nan
            sig_r1_r2_mm[ii][jj] = np.nan
            tau_r1_r2_mm[ii][jj] = np.nan
        else:
            sigma_sl = sl.conductivity_sl(w, c0, T, epsr, mu, muS, 1e-5, 1e-7,
                                       r1x[ii], r2[jj], zeta, SigS)
            sigma_mm = mm.conductivity_mm(w, T, c0, mu, epsr, 1e-5, 1e-7,
                                      r1x[ii], r2[jj], zeta, SigS)
            sig_l_sl, tau_l_sl = findmax.findmax(w,sigma_sl.imag,0)
            sig_h_sl, tau_h_sl = findmax.findmax(w,sigma_sl.imag,1)
            sig_mm, tau_mm = findmax.findmax(w,sigma_mm.imag,0)
            if sig_l_sl == sig_h_sl or sig_l_sl < 0:
                sig_r1_r2_sl[ii][jj] = np.nan
                tau_r1_r2_sl[ii][jj] = np.nan
            else:
                sig_r1_r2_sl[ii][jj] = sig_l_sl
                tau_r1_r2_sl[ii][jj] = tau_l_sl
            sig_r1_r2_mm[ii][jj] = sig_mm
            tau_r1_r2_mm[ii][jj] = tau_mm
        r1x_mesh[ii][jj] = r1x[ii]
        r2_mesh[ii][jj] = r2[jj]
        print(["r1/r2", (100*ii+jj+1)/100**2])
    np.savetxt('geometry/r1_r2/SL_sig_r1_r2.txt', sig_r1_r2_sl) 
    np.savetxt('geometry/r1_r2/SL_tau_r1_r2.txt', tau_r1_r2_sl)
    np.savetxt('geometry/r1_r2/MM_sig_r1_r2.txt', sig_r1_r2_mm) 
    np.savetxt('geometry/r1_r2/MM_tau_r1_r2.txt', tau_r1_r2_mm)
    np.savetxt('geometry/r1_r2/r1x_mesh.txt', r1x_mesh)
    np.savetxt('geometry/r1_r2/r2_mesh.txt', r2_mesh)       