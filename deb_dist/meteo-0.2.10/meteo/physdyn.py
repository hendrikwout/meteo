from numpy import log, nan, sign, pi, sqrt, arctan

g = 9.81 # gravitational aCCelleration ms^-2
krm = 0.40 # von Karman constant


# constants for the flux-profile relationships
gammaM = 16.0
gammaH = 16.0
alphaM = 0.25
alphaH=0.5
AA = 6.1
BB = 2.5
CC = 5.3
DD = 1.1
p2 = pi/2.

def Lmo(ustar,thetastar,theta0=293.15):
    return ustar**2 * theta0/(krm*g*thetastar)

# change to make it applicable for zref = z0
def funm(z,zref,zeta):
    funm = log(z/zref) - psim(zeta/zref*z) + psim(zeta)
    print ('funm: ',funm )
    return funm

def funh(z,zref,zeta):
    funh = log(z/zref) - psih(zeta/zref*z) + psih(zeta)
    print ('funh: ',funh)
    return funh

def ustr(uref,cm):
    # # if we do not follow convention, but more logical...
    # return sqrt(cm)*abs(uref)
    
    # following the convention
    return sqrt(cm)*abs(uref)

def thetastr(thetaref,thetas, uref,ch,ustar):
    return ch*(thetaref - thetas)*uref/ustar

# In order to calculate the 10m wind and 2m temperature, your only need
# zref, uref, thetaref, cm, ch
def uprofile(z,ustar,zref,zeta,uref):
    print(ustar,uref)
    # calculate wind profile from meteorological surface parameters
    # 
    # # if we do not follow convention, but more logical... (see also 'ustr()' )
    # return uref + uref*ustar/krm*funm(z,zref,zeta)

    # following the convention:
    return uref + sign(uref)*ustar/krm*funm(z,zref,zeta)
def thetaprofile(z,thetastar,zref,zeta,thetaref):
    # calculate pot. temperature profile from meteorological surface parameters
    return thetaref + thetastar/krm*funh(z,zref,zeta)

def CM(z,z0,ustar,thetastar):
    Lm = Lmo(ustar,thetastar)
    return krm**2/funm**2
def CH(z,z0,z0h,Lm=None,ustar=nan,thetastar=nan):
    if Lm != None:
        Lmdef = Lm
    else:
        Lmdef = Lmo(ustar,thetastar)
    return krm**2/(funm*funh)

# integrated flux-profile relationship for momentum
def psim(zeta):
    if (zeta < 0):
        # unstable case 
        #   Paulson (1970) integrated flux-profile relationships derived from 
        #   Dyer (1967) /Businger (1966)
        XX=(1.-gammaM*zeta)**alphaM
        return 2.*log((1.+XX)/2.)+log((1.+XX**2.)/2.)-2.*arctan(XX)+p2
    else :
        # stable case
        # Cheng and Brutsaert 2007
        return -AA*log(zeta+(1.+zeta**BB)**(1./BB))
# integrated flux-profile relationship for heat
def psih(zeta):
    if (zeta < 0):
        # unstable case 
        #   Paulson (1970) integrated flux-profile relationships derived from 
        #   Dyer (1967) /Businger (1966)
        XX=(1.-gammaH*zeta)**alphaH
        return 2.*log((1.+XX)/2.)
    else :
        # stable case
        # Cheng and Brutsaert 2005
        return -CC*log(zeta+(1.+zeta**DD)**(1./DD))

# todo: check sign with regard to cm and ch
def zetaC(z,cm,ch,u,theta,thetas,theta0=None):
    if theta0 == None:
        theta0 = thetas
    # zeta = z k g thetast / (ust**2 theta0)
    # zeta = z k g thetast ust / (ust**3 theta0)
    # zeta = z k g -ch ua (theta - thetas) / (ust**3 theta0)
    zeta = z* krm* g* ch*abs(u) *(theta - thetas) / ((cm**(0.5) *abs(u))**3 * theta0)
    return zeta

# no RSL
def zetaCRiB(cm,ch,RiB):
    # zeta = RiB * funm**2/funh
    # zeta = RiB * funm**3/(funh*funm)
    # zeta = RiB * k**3/cm**1.5/(k**2/ch)
    # zeta = RiB * k/cm**(1.5)/(1/ch)
    zeta = RiB*krm*ch/cm**1.5
    return zeta
