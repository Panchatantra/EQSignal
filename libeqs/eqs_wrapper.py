from ctypes import POINTER, c_int, c_double, byref
import numpy as np

libeqs = np.ctypeslib.load_library("libeqs", ".")

# subroutine acc2vd(a,v,d,n,dt,v0,d0)
libeqs.acc2vd.argtypes = [
    POINTER(c_double), POINTER(c_double), POINTER(c_double),
    POINTER(c_int), POINTER(c_double), POINTER(c_double), POINTER(c_double)
]

# subroutine ratacc2vd(a,v,d,n,dt,v0,d0)
libeqs.ratacc2vd.argtypes = [
    POINTER(c_double), POINTER(c_double), POINTER(c_double),
    POINTER(c_int), POINTER(c_double), POINTER(c_double), POINTER(c_double)
]

# subroutine spectrum(acc,n,dt,zeta,P,np,SPA,SPI,SM)
libeqs.spectrum.argtypes = [
    POINTER(c_double), POINTER(c_int), POINTER(c_double), POINTER(c_double),
    POINTER(c_double), POINTER(c_int), POINTER(c_double), POINTER(c_int), POINTER(c_int)
]

# subroutine spectrumavd(acc,n,dt,zeta,P,np,SPA,SPI,SPV,SPD,SPEV,SM)
libeqs.spectrumavd.argtypes = [
    POINTER(c_double), POINTER(c_int), POINTER(c_double), POINTER(c_double),
    POINTER(c_double), POINTER(c_int), POINTER(c_double), POINTER(c_int),
    POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_int)
]

# subroutine spmu(acc,n,dt,zeta,P,nP,SPA,SPI,SPV,SPD,SPE,mu,model,tol,maxiter,uy,rk)
libeqs.spmu.argtypes = [
    POINTER(c_double), POINTER(c_int), POINTER(c_double), POINTER(c_double),
    POINTER(c_double), POINTER(c_int), POINTER(c_double), POINTER(c_int),
    POINTER(c_double), POINTER(c_double), POINTER(c_double),
    POINTER(c_double), POINTER(c_int), POINTER(c_double), POINTER(c_int),
    POINTER(c_double), POINTER(c_double)
]

# subroutine r(acc,n,dt,zeta,P,ra,rv,rd,SM)
libeqs.r.argtypes = [
    POINTER(c_double), POINTER(c_int), POINTER(c_double), POINTER(c_double),
    POINTER(c_double), POINTER(c_double), POINTER(c_double),
    POINTER(c_double), POINTER(c_int)
]

# subroutine rnl(acc,n,dt,zeta,P,ra,rv,rd,SM,mu)
libeqs.rnl.argtypes = [
    POINTER(c_double), POINTER(c_int), POINTER(c_double), POINTER(c_double),
    POINTER(c_double), POINTER(c_double), POINTER(c_double),
    POINTER(c_double), POINTER(c_int), POINTER(c_double)
]