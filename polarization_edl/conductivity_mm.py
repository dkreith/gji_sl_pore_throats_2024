import numpy as np
import cmath
import effmob

def conductivity_mm(w, T, c0, mu, epsr, L1, L2, r1, r2, zeta, SigS):
    
    # Redefine pore length according to MM model
    L1 = 2*L1
    L2 = 2*L2
    
    # Constants
    kB = 8.6173e-5  # Boltzmanm's constant [eV/K]
    F = 96485  # Faraday's constant [C/mol]
    
    # Calculate equilibrium conductivity
    sigma0 = 2*F*mu*c0
    
    fQ = partcoeff(SigS, zeta, c0, T, epsr, 'f')
    
    # Diffuse layer mobilities
    mup1, mun1, mup2, mun2 = effmob.effmob(T, c0, mu, epsr, r1, r2, zeta, fQ)
    
    # Normalize mobilities in narrow pore to cross section ratio
    R = r2**2/r1**2
    mup2 = R*mup2
    mun2 = R*mun2
    
    # Derived values
    Dp1 = kB*T*mup1
    Dp2 = kB*T*mup2
        
    tp1 = mup1/(mup1+mun1)
    tp2 = mup2/(mup2+mun2)
    tn1 = mun1/(mup1+mun1)
    tn2 = mun2/(mup2+mun2)
    
    A = L1/L2
    B = Dp1/Dp2
    S1 = tn1/tp1
    S2 = tn2/tp2
        
    sigma = []
    
    for ii in range(len(w)):
        X1 = cmath.sqrt(1j*w[ii]/2/tn1/Dp1)/2*L1
        X2 = cmath.sqrt(1j*w[ii]/2/tn2/Dp2)/2*L2
        
        # Calculate MM impedance
        Z = L1/mup1/c0/F*(tp1+B/A*tp2+(S2-S1)**2/
                          (X1*S1/tp2**2/tp1/cmath.tanh(X1)+
                           A*X2*S2/B/tp1**2/tp2/cmath.tanh(X2)))
    
        sigma.append((L1+L2)/Z)
        
    # Normalize to equilibrium conductivity
    sig = np.array(sigma)/sigma0
    
    return sig

def partcoeff(SigS, zeta, c0, T, epsr, mod):
    
    if mod == 'f':
        SigD = zeta2sigma(zeta, c0, epsr, T, 'd+')
    elif mod == 'p':
        SigD = zeta2sigma(zeta, c0, epsr, T, 'd')
        
    fQ = SigS/(SigS+SigD)
    
    return fQ

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