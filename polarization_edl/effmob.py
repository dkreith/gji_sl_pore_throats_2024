import numpy as np
from scipy import special
from scipy import integrate

def effmob(T, c0, mu, epsr, r1, r2, zeta, fQ):
    c1 = meanc(r1, c0, epsr, T, zeta, fQ)
    c2 = meanc(r2, c0, epsr, T, zeta, fQ)
        
    mp1 = mu*c1[0]/c0
    mn1 = mu*c1[1]/c0
    mp2 = mu*c2[0]/c0
    mn2 = mu*c2[1]/c0

    m = [mp1, mn1, mp2, mn2]
    return m


def meanc(r0, c0, epsr, T, zeta, fQ):
    # Constants
    kB = np.float64(8.6173e-5)  # Boltzmanm's constant [eV/K]
    F = np.float64(96485)  # Faraday's constant [C/mol]
    eps0 = np.float64(8.85e-12)  # vacuum permittivity [C/Vm]

    # Calculate Debye screening length [m]
    Ld = 1/np.sqrt(2*c0*F/(eps0*epsr*kB*T))

    # Average ion concentrations
    if r0 < 100*Ld:
        # For small pore radii, use potential in cylindrical pore
        def integrandp(r):
            return (r*np.exp(-zeta*special.jv(0, 1j*r/Ld) /
                             (special.jv(0, 1j*r0/Ld)*kB*T)))

        def integrandn(r):
            return (r*np.exp(zeta*special.jv(0, 1j*r/Ld) /
                             (special.jv(0, 1j*r0/Ld)*kB*T)))
    else:
        # For large pore radii, use approximation of plane surface
        def integrandp(r):
            return r*np.exp(-zeta*np.exp((r-r0)/Ld)/(kB*T))

        def integrandn(r):
            return r*np.exp(zeta*np.exp((r-r0)/Ld)/(kB*T))
    
    steps = 1024
    
    def integration_c(r0, steps):
        meshr = (r0 - np.logspace(np.log(1e-20),
                                  np.log(r0), num=steps, base=np.e))
    
        arr_intp = []
        arr_intn = []
        for ii in range(len(meshr)):
            arr_intp.append(integrandp(meshr[ii]))
            arr_intn.append(integrandn(meshr[ii]))
        arr_intp = np.array(arr_intp)
        arr_intn = np.array(arr_intn)

        # Numerical integration of real part
        intp_re = integrate.simps(arr_intp.real, meshr)
        # Numerical integration of imaginary part
        intp_im = integrate.simps(arr_intp.imag, meshr)
        intp = intp_re + 1j * intp_im
        # Numerical integration of real part
        intn_re = integrate.simps(arr_intn.real, meshr)
        # Numerical integration of imaginary part
        intn_im = integrate.simps(arr_intn.imag, meshr)
        intn = intn_re + 1j * intn_im
        # Average cation concentration [mol/m^3]
        cp_mean = np.abs(2*c0/(r0**2)*intp)
        cp_mean = 1/(1-fQ)*(cp_mean-c0*fQ)
        # Average anion concentration [mol/m^3]
        cn_mean = np.abs(2*c0/(r0**2)*intn)
        
        if np.isnan(cp_mean) or np.isnan(cn_mean):
            steps = int(steps/2)
            cp_mean, cn_mean = integration_c(r0,steps)
        
        return [cp_mean, cn_mean]
    
    cp_mean, cn_mean = integration_c(r0,steps)
    return [cp_mean, cn_mean]