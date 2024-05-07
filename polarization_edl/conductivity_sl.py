import numpy as np
import cmath
import effmob

def conductivity_sl(w, c0, T, epsr, mu, muS, L1, L2, r1, r2, zeta, SigS):

    # Constants
    kB = 8.6173e-5  # Boltzmanm's constant [eV/K]
    F = 96485  # Faraday's constant [C/mol]
    eps0 = 8.85e-12  # vacuum permittivity [C/Vm]
    
    # Calculate equilibrium conductivity
    sigma0 = 2*F*mu*c0

    # Diffuse layer mobilities
    mup1, mun1, mup2, mun2 = effmob.effmob(T, c0, mu, epsr, r1, r2, zeta, 0)

    # Calculate effective Stern-layer mobilities and concentrations
    muS1 = muS
    muS2 = muS
    cS1 = 2*SigS/r1/F
    cS2 = 2*SigS/r2/F

    # Calculate diffusivities using Nernst-Einstein equation
    Dp1 = kB*T*mup1
    Dp2 = kB*T*mup2
    Dn1 = kB*T*mun1
    Dn2 = kB*T*mun2
    DS1 = kB*T*muS1
    DS2 = kB*T*muS2

    # Define parameters
    k = np.sqrt(2*c0*F/kB/T/eps0/epsr)  # Inverse Debye length [1/m]
    a = r2**2/r1**2
    M = 1 + (k*SigS/2/F/c0/np.cosh(zeta/2/kB/T))
    xi1 = cS1/c0
    xi2 = cS2/c0

    sigma = []

    # Iterate over frequencies
    for ii in range(len(w)):
        # Create matrix of differential equation
        A1 = np.matrix([[1j*w[ii]/Dp1+k**2/2, -k**2/2, k**2/2],
                        [-k**2/2, 1j*w[ii]/Dn1+k**2/2, -k**2/2],
                        [xi1*k**2/M/2, -xi1*k**2/M/2,
                         1j*w[ii]/M/DS1+xi1*k**2/M/2]])
        A2 = np.matrix([[1j*w[ii]/Dp2+k**2/2, -k**2/2, k**2/2],
                        [-k**2/2, 1j*w[ii]/Dn2+k**2/2, -k**2/2],
                        [xi2*k**2/M/2, -xi2*k**2/M/2,
                         1j*w[ii]/M/DS2+xi2*k**2/M/2]])

        # Solve for eigenvalues and calculate lambda
        eigval1, eigvec1 = np.linalg.eig(A1)
        eigval2, eigvec2 = np.linalg.eig(A2)
        lamb11, lamb21, lamb31 = np.sqrt(eigval1)
        lamb12, lamb22, lamb32 = np.sqrt(eigval2)

        # Define and solve equation systems for coefficient ratio
        v = F*c0/eps0/epsr
        CRRHS1 = np.array([-mup1*v, mun1*v])
        CRRHS2 = np.array([-mup2*v, mun2*v])

        CREQS11 = np.matrix([[1j*w[ii]-Dp1*lamb11**2+mup1*v,
                              -mup1*v],
                             [-mun1*v,
                              1j*w[ii]-Dn1*lamb11**2+mun1*v]])
        CREQS21 = np.matrix([[1j*w[ii]-Dp1*lamb21**2+mup1*v,
                              -mup1*v],
                             [-mun1*v,
                              1j*w[ii]-Dn1*lamb21**2+mun1*v]])
        CREQS31 = np.matrix([[1j*w[ii]-Dp1*lamb31**2+mup1*v,
                              -mup1*v],
                             [-mun1*v,
                              1j*w[ii]-Dn1*lamb31**2+mun1*v]])
        
        CREQS12 = np.matrix([[1j*w[ii]-Dp2*lamb12**2+mup2*v,
                              -mup2*v],
                             [-mun2*v,
                              1j*w[ii]-Dn2*lamb12**2+mun2*v]])
        CREQS22 = np.matrix([[1j*w[ii]-Dp2*lamb22**2+mup2*v,
                              -mup2*v],
                             [-mun2*v,
                              1j*w[ii]-Dn2*lamb22**2+mun2*v]])
        CREQS32 = np.matrix([[1j*w[ii]-Dp2*lamb32**2+mup2*v,
                              -mup2*v],
                             [-mun2*v,
                              1j*w[ii]-Dn2*lamb32**2+mun2*v]])

        alph11, beta11 = np.linalg.solve(CREQS11, CRRHS1)
        alph21, beta21 = np.linalg.solve(CREQS21, CRRHS1)
        alph31, beta31 = np.linalg.solve(CREQS31, CRRHS1)

        alph12, beta12 = np.linalg.solve(CREQS12, CRRHS2)
        alph22, beta22 = np.linalg.solve(CREQS22, CRRHS2)
        alph32, beta32 = np.linalg.solve(CREQS32, CRRHS2)

        # Define abbreviations
        gamm11 = L1/lamb11*(1-cmath.tanh(lamb11*L1)/L1/lamb11)
        gamm21 = L1/lamb21*(1-cmath.tanh(lamb21*L1)/L1/lamb21)
        gamm31 = L1/lamb31*(1-cmath.tanh(lamb31*L1)/L1/lamb31)

        gamm12 = L2/lamb12*(1-cmath.tanh(lamb12*L2)/L2/lamb12)
        gamm22 = L2/lamb22*(1-cmath.tanh(lamb22*L2)/L2/lamb22)
        gamm32 = L2/lamb32*(1-cmath.tanh(lamb32*L2)/L2/lamb32)

        deltp11 = 2*Dp1*(L1+L2)*lamb11/k**2/(a*Dp2-Dp1)
        deltp21 = 2*Dp1*(L1+L2)*lamb21/k**2/(a*Dp2-Dp1)
        deltp31 = 2*Dp1*(L1+L2)*lamb31/k**2/(a*Dp2-Dp1)

        deltp12 = 2*a*Dp2*(L1+L2)*lamb12/k**2/(a*Dp2-Dp1)
        deltp22 = 2*a*Dp2*(L1+L2)*lamb22/k**2/(a*Dp2-Dp1)
        deltp32 = 2*a*Dp2*(L1+L2)*lamb32/k**2/(a*Dp2-Dp1)

        deltn11 = 2*Dn1*(L1+L2)*lamb11/k**2/(a*Dn2-Dn1)
        deltn21 = 2*Dn1*(L1+L2)*lamb21/k**2/(a*Dn2-Dn1)
        deltn31 = 2*Dn1*(L1+L2)*lamb31/k**2/(a*Dn2-Dn1)

        deltn12 = 2*a*Dn2*(L1+L2)*lamb12/k**2/(a*Dn2-Dn1)
        deltn22 = 2*a*Dn2*(L1+L2)*lamb22/k**2/(a*Dn2-Dn1)
        deltn32 = 2*a*Dn2*(L1+L2)*lamb32/k**2/(a*Dn2-Dn1)

        deltS11 = 2*M*DS1*(L1+L2)*lamb11/k**2/(a*xi2*DS2-xi1*DS1)
        deltS21 = 2*M*DS1*(L1+L2)*lamb21/k**2/(a*xi2*DS2-xi1*DS1)
        deltS31 = 2*M*DS1*(L1+L2)*lamb31/k**2/(a*xi2*DS2-xi1*DS1)

        deltS12 = 2*a*M*DS1*(L1+L2)*lamb12/k**2/(a*xi2*DS2-xi1*DS1)
        deltS22 = 2*a*M*DS1*(L1+L2)*lamb22/k**2/(a*xi2*DS2-xi1*DS1)
        deltS32 = 2*a*M*DS1*(L1+L2)*lamb32/k**2/(a*xi2*DS2-xi1*DS1)

        # Define matrix elements
        Kp11 = ((alph11-beta11+1)*gamm11+alph11*deltp11)/cmath.tanh(lamb11*L1)
        Kp21 = ((alph21-beta21+1)*gamm21+alph21*deltp21)/cmath.tanh(lamb21*L1)
        Kp31 = ((alph31-beta31+1)*gamm31+alph31*deltp31)/cmath.tanh(lamb31*L1)

        Kp12 = ((alph12-beta12+1)*gamm12-alph12*deltp12)/cmath.tanh(lamb12*L2)
        Kp22 = ((alph22-beta22+1)*gamm22-alph22*deltp22)/cmath.tanh(lamb22*L2)
        Kp32 = ((alph32-beta32+1)*gamm32-alph32*deltp32)/cmath.tanh(lamb32*L2)

        Kn11 = ((alph11-beta11+1)*gamm11-beta11*deltn11)/cmath.tanh(lamb11*L1)
        Kn21 = ((alph21-beta21+1)*gamm21-beta21*deltn21)/cmath.tanh(lamb21*L1)
        Kn31 = ((alph31-beta31+1)*gamm31-beta31*deltn31)/cmath.tanh(lamb31*L1)

        Kn12 = ((alph12-beta12+1)*gamm12+beta12*deltn12)/cmath.tanh(lamb12*L2)
        Kn22 = ((alph22-beta22+1)*gamm22+beta22*deltn22)/cmath.tanh(lamb22*L2)
        Kn32 = ((alph32-beta32+1)*gamm32+beta32*deltn32)/cmath.tanh(lamb32*L2)

        KS11 = ((alph11-beta11+1)*gamm11+deltS11)/cmath.tanh(lamb11*L1)
        KS21 = ((alph21-beta21+1)*gamm21+deltS21)/cmath.tanh(lamb21*L1)
        KS31 = ((alph31-beta31+1)*gamm31+deltS31)/cmath.tanh(lamb31*L1)

        KS12 = ((alph12-beta12+1)*gamm12-deltS12)/cmath.tanh(lamb12*L2)
        KS22 = ((alph22-beta22+1)*gamm22-deltS22)/cmath.tanh(lamb22*L2)
        KS32 = ((alph32-beta32+1)*gamm32-deltS32)/cmath.tanh(lamb32*L2)

        # Create and solve matrix for calculation of coefficients
        alphmax = max(np.abs([alph11, alph21, alph31, alph12, alph22, alph32]))
        betamax = max(np.abs([beta11, beta21, beta31, beta12, beta22, beta32]))
        Kpmax = max(np.abs([Kp11, Kp21, Kp31, Kp12, Kp22, Kp32]))
        Knmax = max(np.abs([Kn11, Kn21, Kn31, Kn12, Kn22, Kn32]))
        KSmax = max(np.abs([KS11, KS21, KS31, KS12, KS22, KS32]))

        C = np.matrix([np.array([alph11, alph21, alph31,
                                 alph12, alph22, alph32])/alphmax,
                       np.array([beta11, beta21, beta31,
                                 beta12, beta22, beta32])/betamax,
                       np.array([1, 1, 1, np.sqrt(a), np.sqrt(a), np.sqrt(a)]),
                       np.array([Kp11, Kp21, Kp31, Kp12, Kp22, Kp32])/Kpmax,
                       np.array([Kn11, Kn21, Kn31, Kn12, Kn22, Kn32])/Knmax,
                       np.array([KS11, KS21, KS31, KS12, KS22, KS32])/KSmax])

        c = np.array([0, 0, 0, -1/2/Kpmax, -1/2/Knmax, -1/2/KSmax])

        C11, C21, C31, C12, C22, C32 = np.linalg.solve(C, c)

        # Calculate effective complex conductivity
        Z = ((alph11-beta11+1)*gamm11/cmath.tanh(lamb11*L1)/L1*C11 +
             (alph21-beta21+1)*gamm21/cmath.tanh(lamb21*L1)/L1*C21 +
             (alph31-beta31+1)*gamm31/cmath.tanh(lamb31*L1)/L1*C31 -
             (alph12-beta12+1)*gamm12/cmath.tanh(lamb12*L2)/L2*C12 -
             (alph22-beta22+1)*gamm22/cmath.tanh(lamb22*L2)/L2*C22 -
             (alph32-beta32+1)*gamm32/cmath.tanh(lamb32*L2)/L2*C32)

        sigma.append(2*eps0*epsr*(k**2/4/(L1+L2) *
                                  (L1*(Dp1+Dn1+xi1*DS1) +
                                   a*L2*(Dp2+Dn2+xi2*DS2) +
                                   2*L1*L2*Z*(a*Dp2-Dp1+a*Dn2-Dn1 +
                                              a*xi2*DS2-xi1*DS1)) +
                                  (a*Dp2-Dp1)*(alph11*C11 +
                                               alph21*C21+alph31*C31) -
                                  (a*Dn2-Dn1)*(beta11*C11 +
                                               beta21*C21+beta31*C31) +
                                  M*(a*np.sqrt(a)*DS2-DS1)*(C11+C21+C31)))

    # Normalize to equilibrium conductivity
    sig = np.array(sigma)/sigma0
    
    return sig
