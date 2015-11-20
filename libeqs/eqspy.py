from enum import Enum
from ctypes import POINTER, c_int, c_double, byref
from eqs_wrapper import libeqs
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')


def plotTH(t, acc, vel, dsp, title="Time History Curves"):

    plt.figure(title,(8,6))
    plt.subplot(3,1,1)
    plt.plot(t,acc)
    plt.grid()
    plt.subplot(3,1,2)
    plt.plot(t,vel)
    plt.grid()
    plt.subplot(3,1,3)
    plt.plot(t,dsp)
    plt.grid()
    plt.show()

class EQSignal(object):
    """docstring for EQSignal"""
    v0 = 0.0
    d0 = 0.0
    def __init__(self, acc=np.zeros(1024), dt=0.02):

        super(EQSignal, self).__init__()
        self.acc = np.asarray(acc)
        self.dt = dt
        self.n = len(self.acc)
        self.t = np.linspace(0.0, self.dt*self.n-self.dt, self.n)
        self.vel = np.zeros(self.n)
        self.dsp = np.zeros(self.n)

        self.generate_ptr()


    @staticmethod
    def from_file(fname, dt=0.02):
        acc = np.loadtxt(fname)
        return EQSignal(acc,dt)


    def allocate(self):
        self.t = np.linspace(0.0, self.dt*self.n-self.dt, self.n)
        self.vel = np.zeros(self.n)
        self.dsp = np.zeros(self.n)


    def generate_ptr(self):

        self.acc_ptr = self.acc.ctypes.data_as(POINTER(c_double))
        self.vel_ptr = self.vel.ctypes.data_as(POINTER(c_double))
        self.dsp_ptr = self.dsp.ctypes.data_as(POINTER(c_double))

        self.n_ptr = byref(c_int(self.n))
        self.dt_ptr = byref(c_double(self.dt))
        self.v0_ptr = byref(c_double(self.v0))
        self.d0_ptr = byref(c_double(self.d0))


    def acc2vd(self):

        libeqs.acc2vd(self.acc_ptr, self.vel_ptr, self.dsp_ptr,
                        self.n_ptr, self.dt_ptr, self.v0_ptr, self.d0_ptr)


    def response(self, zeta=0.05, P=2.0, SM=0, mu=None):

        zeta_ptr = byref(c_double(zeta))
        P_ptr = byref(c_double(P))
        SM_ptr = byref(c_int(SM))

        ra = np.zeros(self.n)
        rv = np.zeros(self.n)
        rd = np.zeros(self.n)

        bl = np.zeros(self.n)

        ra_ptr = ra.ctypes.data_as(POINTER(c_double))
        rv_ptr = rv.ctypes.data_as(POINTER(c_double))
        rd_ptr = rd.ctypes.data_as(POINTER(c_double))

        if mu == None:
            libeqs.r(self.acc_ptr, self.n_ptr, self.dt_ptr, zeta_ptr, P_ptr,
                        ra_ptr, rv_ptr, rd_ptr, SM_ptr)
        else:
            mu_ptr = byref(c_double(mu))
            libeqs.rnl(self.acc_ptr, self.n_ptr, self.dt_ptr, zeta_ptr, P_ptr,
                        ra_ptr, rv_ptr, rd_ptr, SM_ptr, mu_ptr)

        w = 2.0*np.pi/P
        k = w*w
        c = 2.0*w*zeta
        
        rf = -ra - c*rv
        
        return self.t, ra, rv, rd, rf


    def spectrum(self):
        pass


    def plotTH(self):
        plotTH(self.t, self.acc, self.vel, self.dsp)

class Response(object):
    """docstring for Response"""
    def __init__(self, eqs, zeta=0.05, P=2.0):
        super(Response, self).__init__()
        self.acc = eqs.acc
        self.t = eqs.t
        self.dt = eqs.dt
        self.n = eqs.n
        self.zeta = zeta
        self.P = P

        self.ra = np.zeros(self.n)
        self.rv = np.zeros(self.n)
        self.rd = np.zeros(self.n)

        self.generate_ptr()
    
    def generate_ptr(self):

        self.acc_ptr = self.acc.ctypes.data_as(POINTER(c_double))
        self.t_ptr = self.t.ctypes.data_as(POINTER(c_double))
        self.ra_ptr = self.ra.ctypes.data_as(POINTER(c_double))
        self.rv_ptr = self.rv.ctypes.data_as(POINTER(c_double))
        self.rd_ptr = self.rd.ctypes.data_as(POINTER(c_double))

        self.dt_ptr = byref(c_double(self.dt))
        self.n_ptr = byref(c_int(self.n))
        self.zeta_ptr = byref(c_double(self.zeta))
        self.P_ptr = byref(c_double(self.P))


class EQSpectra(object):
    """docstring for EQSpectra"""
    def __init__(self, eqs=EQSignal(), zeta=0.05, P=np.logspace(-1.3,1,30),
                    RM=3, pseudo=0):
        super(EQSpectra, self).__init__()
        self.acc = eqs.acc
        self.n = eqs.n
        self.dt = eqs.dt
        self.zeta = zeta
        self.P = P
        self.nP = len(P)
        self.RM = RM + pseudo*3

        self.SPA = np.zeros(self.nP)
        self.SPV = np.zeros(self.nP)
        self.SPD = np.zeros(self.nP)
        self.SPE = np.zeros(self.nP)
        self.SPI = np.zeros(self.nP, dtype=np.int32)

        self.generate_ptr()

    def generate_ptr(self):

        self.acc_ptr = self.acc.ctypes.data_as(POINTER(c_double))
        self.P_ptr = self.P.ctypes.data_as(POINTER(c_double))
        self.SPA_ptr = self.SPA.ctypes.data_as(POINTER(c_double))
        self.SPV_ptr = self.SPV.ctypes.data_as(POINTER(c_double))
        self.SPD_ptr = self.SPD.ctypes.data_as(POINTER(c_double))
        self.SPE_ptr = self.SPE.ctypes.data_as(POINTER(c_double))
        self.SPI_ptr = self.SPI.ctypes.data_as(POINTER(c_int))

        self.n_ptr = byref(c_int(self.n))
        self.nP_ptr = byref(c_int(self.nP))
        self.dt_ptr = byref(c_double(self.dt))
        self.zeta_ptr = byref(c_double(self.zeta))
        self.RM_ptr = byref(c_int(self.RM))

    def calc(self):
        libeqs.spectrumavd(self.acc_ptr,self.n_ptr,self.dt_ptr,
                        self.zeta_ptr,self.P_ptr,self.nP_ptr,self.SPA_ptr,
                        self.SPI_ptr,self.SPV_ptr,self.SPD_ptr,
                        self.SPE_ptr, self.RM_ptr)

    def calc_nl(self,mu=2.0,model=0,rk=0.0,tol=0.001,maxiter=100):

        mu_ptr = byref(c_double(mu))
        model_ptr = byref(c_int(model))
        rk_ptr = byref(c_double(rk))
        tol_ptr = byref(c_double(tol))
        maxiter_ptr = byref(c_int(maxiter))
        
        uy = np.zeros(self.nP)
        uy_ptr = uy.ctypes.data_as(POINTER(c_double))
        
        libeqs.spmu(self.acc_ptr,self.n_ptr,self.dt_ptr,
                        self.zeta_ptr,self.P_ptr,self.nP_ptr,self.SPA_ptr,
                        self.SPI_ptr,self.SPV_ptr,self.SPD_ptr,
                        self.SPE_ptr,mu_ptr,model_ptr,tol_ptr,
                        maxiter_ptr,uy_ptr,rk_ptr)


    def plotSPA(self,xs="log",show=True):

        plt.figure("SPA",(8,6))
        plt.plot(self.P,self.SPA,label="$ \zeta $ = %.2f"%self.zeta)
        plt.xscale(xs)
        plt.grid()
        if show:
            plt.legend(loc=0)
            plt.show()

    def plotSP(self,xs="log",show=True):

        plt.figure("EQSignal --- Spectra",(8,6))
        plt.subplot(2,2,1)
        plt.plot(self.P,self.SPA,label="$ \zeta $ = %.2f"%self.zeta)
        plt.xscale(xs)
        plt.grid()
        if show: plt.legend(loc=0)
        plt.subplot(2,2,2)
        plt.plot(self.P,self.SPV,label="$ \zeta $ = %.2f"%self.zeta)
        plt.xscale(xs)
        plt.grid()
        if show: plt.legend(loc=0)
        plt.subplot(2,2,3)
        plt.plot(self.P,self.SPD,label="$ \zeta $ = %.2f"%self.zeta)
        plt.xscale(xs)
        plt.grid()
        if show: plt.legend(loc=0)
        plt.subplot(2,2,4)
        plt.plot(self.P,self.SPE,label="$ \zeta $ = %.2f"%self.zeta)
        plt.xscale(xs)
        plt.grid()
        if show: plt.legend(loc=0)
        if show:
            plt.show()
        


if __name__ == '__main__':

    eqs = EQSignal.from_file("../data/AWX09-1.txt")
    # eqs.acc2vd()
    # eqs.plotTH()
    # t,ra,rv,rd,rf = eqs.response(mu=2.0)
    # plotTH(t,ra,rv,rd)

    # plt.figure("Hysteretic",(8,6))
    # plt.plot(rd,rf)
    # plt.show()

    # sp = EQSpectra(eqs, 0.05)
    # sp.calc()
    # sp.plotSP(show=False)

    # sp = EQSpectra(eqs, 0.02)
    # sp.calc()
    # sp.plotSP(show=False)

    sp = EQSpectra(eqs, zeta=0.01)
    sp.calc()
    # sp.calc_nl(mu=2.0)
    sp.plotSP()