import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import quad
from scipy.special import gamma
import matplotlib.pyplot as plt 
import plots

def kappa(p, lam, d):
    """kappa_p(lambda) in manuscript"""
    return lam*((p - 1)*lam + d - p)

def rhs(theta, y, p, lam, d):
    """Right hand side of the first order system y = (f, f')."""
    f, fp = y
    k = kappa(p, lam, d)
    D = fp*fp + lam*lam*f*f         
    num = ((d - 2)*fp/np.tan(theta) + k*f)*D + (p - 2)*lam*lam*f*fp*fp
    den = (p - 1)*fp*fp + lam*lam*f*f     
    return [fp, -num/den]

def shoot(lam, p, d, eps=1e-6, rtol=1e-12, atol=1e-14, dense=False):
    """Integrate ODE for spherical part, given a choice of lambda
    Starts at theta = eps using the Taylor expansion near theta=0
    """
    k = kappa(p, lam, d)       
    y0 = [1.0 - 1.0/d - k*eps*eps/(2.0*d), -k*eps/d]
    sol = solve_ivp(rhs, [eps, np.pi/2], y0, args=(p, lam, d), method='DOP853', rtol=rtol, atol=atol, dense_output=dense)

    if not sol.success:
        raise RuntimeError(f"IVP failed: p={p}, lam={lam}, d={d}: {sol.message}")
    return sol

def find_lam(p, d, eps=1e-6):
    """Use shooting method to find lambda""" 

    #Initialization
    a = 4/3  #Aronsson solution exponent
    b = 2.0  #For p=2
    lam = (a+b)/2

    #binary search
    for i in range(50):
        sol = shoot(lam, p, d)
        fp = sol.y[1,-1]
        if fp > 0:
            b = lam
        else:
            a = lam
        lam = (a+b)/2

    return lam

def sigma(d):
    """sigma_d from manuscript"""
    return gamma(d/2)/(np.sqrt(np.pi)*gamma((d-1)/2))

def spherical_mean(sol, d, eps=1e-6):
    """Computes the mean of axially symmetric function over sphere"""
    #Contribution from [\eps,pi/2]
    integral = quad(lambda t: sol.sol(t)[0]*np.sin(t)**(d-2), eps, np.pi/2,limit=200)[0]
    #Contribution from [0, eps]
    integral  += sol.sol(eps)[0] * eps ** (d - 1) / (d - 1)
    return 2*sigma(d)*integral

def C(p, d, lam, sol, eps=1e-6):
    mean = spherical_mean(sol, d, eps)
    return 0.5*(p-2)*(1-1/d + sol.sol(np.pi/2)[0]) + d*(d+2)/(lam+d)*mean

def H(m,d,lam):
    return gamma(d/2)*gamma((m+lam)/2)/(gamma(m/2)*gamma((d+lam)/2))

def Cmd(p, m, d, lam, sol, eps=1e-6):
    mean = spherical_mean(sol, m, eps)
    return 0.5*(p-2)*(1-1/m + sol.sol(np.pi/2)[0]) + d*(d+2)/(lam+d)*H(m,d,lam)*mean

print("Plotting several profiles...")
d = 5
h = 0.01
t = np.arange(0,np.pi/2,h)
plt.figure()
for p in [3,5,10,50,1000]:
    lam = find_lam(p,d)
    sol = shoot(lam, p, d, dense=True)
    f = sol.sol(t)[0]
    plt.plot(t,f,label='$p=%d, \\lambda=%.2f$'%(p,lam))

plt.legend()
plt.xlabel('$\\theta$')
plt.ylabel('$f_p(\\theta)$')
plots.savefig('profiles.pdf',grid=True,axis=True)

print("Computing C(p) and lambda(p)...")
p_vals = np.arange(2,50.5,0.5)
Cp_fig, Cp_ax = plt.subplots()
lam_fig, lam_ax = plt.subplots()
for d in [3,4,5]:
    print('    d=%d...'%d)
    Cp_vals = []
    lam_vals = []
    for p in p_vals:
        lam = find_lam(p,d)
        sol = shoot(lam, p, d, dense=True)
        Cp = C(p,d,lam,sol)
        lam_vals += [lam]
        Cp_vals += [Cp]
    Cp_vals = np.array(Cp_vals)
    lam_vals = np.array(lam_vals)

    Cp_ax.plot(p_vals,Cp_vals,label='$d=%d$'%d)
    lam_ax.plot(p_vals,lam_vals,label='$d=%d$'%d)

plt.sca(Cp_ax)
plt.legend()
plt.xlabel('$p$')
plt.ylabel('$\\mathcal C(p)$')
plots.savefig('Cp.pdf',grid=True,axis=True)

plt.sca(lam_ax)
plt.legend()
plt.xlabel('$p$')
plt.ylabel('$\\lambda(p)$')
plots.savefig('lam.pdf',grid=True,axis=True)

print("Computing crossing points of C(p)...")
for d in range(3,11):
    a = 3
    b = 5
    Cp = -1
    while Cp < 0:
        b = 2*b
        lam = find_lam(b,d)
        sol = shoot(lam, b, d, dense=True)
        Cp = C(b,d,lam,sol)

    p = (a+b)/2
    for i in range(50):
        lam = find_lam(p,d)
        sol = shoot(lam, p, d, dense=True)
        Cp = C(p,d,lam,sol)
        if Cp > 0:
            b = p
        else:
            a = p
        p = (a+b)/2

    print('d=%d, p=%.5f'%(d,p)) 

print("Computing C_{3 -> d}(p)...")
p_vals = np.arange(2,61)
plt.figure()
m = 3
for d in [4,5,6]:
    print('    d=%d...'%d)
    Cp_vals = []
    lam_vals = []
    for p in p_vals:
        lam = find_lam(p,m)
        sol = shoot(lam, p, m, dense=True)
        Cp = Cmd(p,m,d,lam,sol)
        lam_vals += [lam]
        Cp_vals += [Cp]
    Cp_vals = np.array(Cp_vals)
    lam_vals = np.array(lam_vals)

    plt.plot(p_vals,Cp_vals,label='$d=%d$'%d)

plt.legend()
plt.xlabel('$p$')
plt.ylabel('$\\mathcal C_{3\\mapsto d}(p)$')
plots.savefig('Cmdp.pdf',grid=True,axis=True)

print("Computing crossing points of C_{3 -> d}(p)...")
for d in range(4,11):
    a = 3
    b = 5
    Cp = -1
    while Cp < 0:
        b = 2*b
        lam = find_lam(b,m)
        sol = shoot(lam, b, m, dense=True)
        Cp = Cmd(b,m,d,lam,sol)

    p = (a+b)/2
    for i in range(50):
        lam = find_lam(p,m)
        sol = shoot(lam, p, m, dense=True)
        Cp = Cmd(p,m,d,lam,sol)
        if Cp > 0:
            b = p
        else:
            a = p
        p = (a+b)/2

    print('d=%d, p=%.5f'%(d,p)) 



