# constants describing air

cpd  = 1005.7 # isobaric specific heat of dry air at constant pressure [J/(kg K)]
cpv = 1864. # isobaric specific heat of water vapour at constant pressure and 300K [J/(kg K)]
rhod  = 1.2 # specific mass of dry air [kg/m**3] for standard atmosphere
p0 = 101325.0 # reference pressure (Pa)
Rd = 287. # specific gas constant for dry air [J kg^-1 K^-1]
Rv = 462. # specific gas constant for water vapour [J kg^-1 K^-1]

g = 9.81 # gravitational accelleration

def rhov(T,p,sh=0.):
    '''purpose: calculate the density from pressure temperature and specific humidity'''
    R = Rd*(1.-sh) + Rv*sh
    return p/R/T

def t_from_theta(ptemp,pressure):
    return(ptemp*(pressure/p0)**(Rd/cpd))

def theta_from_t(ptemp,pressure):
    return(ptemp*(pressure/p0)**(Rd/cpd))

def t_from_thetam(PT,p,w = 0.,qv = None):
    if qv != None:
        w = qv/(1.-qv)
    alpha = Rd*(1.-0.28*w)/cpd
    return PT/((p/p0)**alpha)

#hypsometric equation for moist with constant moist potential temperature
#  input:
#  w: mixing ratio (optional)
def z_const_theta(z_0, p_0, p_1, pt, w = 0.0):
    alpha = Rd*(1-0.28*w)/cpd
    #result = z_0 +  R_spec * ttemp(pt,p_0)/g * log(p_0/p_1)
    return z_0 - (p_1**(alpha) - p_0**(alpha))/(alpha) * Rd *pt / g / (p0**alpha)
