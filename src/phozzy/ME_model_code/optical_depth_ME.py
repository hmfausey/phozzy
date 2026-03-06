
#Optical depth package for GRB210905A


##############################################################################
#################################IMPORTS######################################
##############################################################################

import numpy as np

##############################################################################
############################CONSTANTS(in cgs)#################################
##############################################################################

#Basic constants
c = 2.99792458 * 10**(10) #cm/s
e = 4.8032068 * 10**(-10) #esu or g^(1/2)cm^(3/2)s^(-1)
m_e = 9.1093897 * 10**(-28) #g
m_p = 1.6726231*10**(-24)
G = 6.67259 * 10**(-8) #dyn cm^2 g^(-2)
pi = np.pi
lambda_alpha = 1215.67*10**(-8) #in cm
nu_alpha = c/lambda_alpha #1/s

#From Totani et al. paper
f_alpha = 0.4162
gugl = 3
Lambda_cla =(8 * pi**2 * e**2)/(3*m_e*c * lambda_alpha**2)
#Lambda_cla = 1.503*10**(9)
#print(Lambda_cla) ##Sanity Check

Lambda_alpha = 3 * (gugl)**(-1) * f_alpha * Lambda_cla
#Lambda_alpha = 6.25 * 10**(8) #1/s


#Cosmological parameters (Using Plank 2018 results, using specifically TT, TE, EE+lowE+lensing results, as are used as default in the paper)
H0 = 67.4 #km/s/Mpc
H0_cgs = H0 * (1000/1)*(100/1)*(1/(3.086*10**(24)))
Omega_M = 0.315
h = H0/100
Omega_B = 0.0224/(h**2)
Omega_k = -0.011
Omega_lam = 0.6847
Yp = 0.2454


### Add XH1 as constant here###
x_H1 =0

##############################################################################
################################FUNCTIONS#####################################
##############################################################################

def sigma_alpha(nu):
    ##From Totani et. al. 2014 Eq. 1/ Miralda Escude 1998 Eq. 8
    constant = (3*(lambda_alpha**2)*f_alpha*Lambda_cla)/(8*pi)
    numerator = Lambda_alpha * (nu/nu_alpha)**(4)
    denominator = 4*(pi**2)*(nu - nu_alpha)**(2) + ((Lambda_alpha**2)*(nu/nu_alpha)**(6))/4
    
    return constant * (numerator/denominator)


#most important for our purposes
def T_DLA(logNH, lam_obs, z_DLA): 
    #From Totani pg 4 col 1
    nu_obs = c/lam_obs
    return 10**(logNH) * sigma_alpha(nu_obs*(1+z_DLA)) #
    
def I(x):
    #From Miralda Escude Eq 12/ Totani Eq 3
    result = ((x**(9/2))/(1-x)) + (9/7)*x**(7/2) + 9/5*x**(5/2) + 3*x**(3/2) + 9*x**(1/2) - (9/2)*np.log((1 + (x**(1/2)))/(1 - (x**(1/2))))
    return result

def T_GP(z):
    #Gunn peterson optical depth from Miralda Escude appendix and Totani Eq 4

    p_crit = (3*(H0_cgs**2))/(8*pi*G)
    full = ((3*f_alpha*Lambda_cla*(lambda_alpha**3)*p_crit*Omega_B*(1-Yp))/(8*pi*m_p*H0_cgs*(Omega_M**(1/2))))*((1+z)**(3/2))

    return full

#most important for our purposes
def T_IGM(x_H1, z_host, z_IGMu, z_IGMl, lam_observed): 
    #Combines all relevant optical depths and other contributing factors into the IGM optical depth
    T_0 = T_GP(z_host)
    R_alpha = (Lambda_alpha*lambda_alpha)/(4*pi*c)
    #R_alpha = 2.02 * 10**(-8)
    z_obs_eq = lam_observed/lambda_alpha
    term1 = ((x_H1 * R_alpha * T_0)/(pi)) * ((z_obs_eq)/(1 + z_host))**(3/2)
    term2 = I((1 + z_IGMu)/z_obs_eq) - I((1 + z_IGMl)/z_obs_eq)
    return term1*term2

def H(z):
    #From Miralda Escude 1998 pg 16 col 2
    Hz = H0 * np.sqrt(Omega_M*(1+z)**(3) + Omega_lam + Omega_k*(1+z)**(2))
    return Hz

def T_IGM_McQuinn(x_H1, Rb, z_IGMu, lam_obs):
    
    term1 = 900 * x_H1 * ((1+z_IGMu)/8)**(3/2) #in km/s
    
    nu_z = (c/lam_obs)*(1+z_IGMu)
    term2 = (H(z_IGMu)* Rb)/(1+z_IGMu) - c*10**(-5)* ((nu_z - nu_alpha)/nu_alpha)
    
    return term1/term2


def Teff_total(lam,logNH, x_H1, z_l, z_host,z_DLA, z_u):
    T_IGM1 = T_IGM(x_H1, z_host, z_u, z_l, lam)
    
    T_DLA1 = T_DLA(logNH, lam, z_DLA)
    
    T_eff = T_IGM1 + T_DLA1
    
    return T_eff



#T_IGM(0.2, 6.318, 6, 6.318, 9000 * 10**(-8))
