#ifndef EQSIGNAL_H
#define EQSIGNAL_H

#include <QVector>
#include <QString>
#include "Spectra.h"
#include <complex>

struct Response
{
    double zeta, P;
    double *ra, *rv, *rd, *ku;
    double fy, uy;
};

class EQSignal
{

public:

    EQSignal();
    EQSignal(int N, double DT, double V0=0.0, double D0=0.0);
    EQSignal(double *a, int N, double DT, double V0=0.0, double D0=0.0);
    ~EQSignal();

    void readtxt(const char *filename, double DT=0.02, bool NORM=false);
    void readtxt(QString filename, double DT=0.02, bool NORM=false);
    void readnga(const char *filename, bool NORM=false);
    void readnga(QString filename, bool NORM=false);

    void resample(int r);
    void interpolate();
    void recover();
    void trim(int ind1, int ind2);
    void confirm();
    void norm();
    void calcAriasIntensity();

    void copyAccFrom(EQSignal *eqs);

    int *autoTrimEdges(int method=0, double thd1=0.02, double thd2=0.98, bool EZ=true);

    void savetxt(const char *filename);
    void savetxt(QString filename);
    void savecsv(const char *filename);
    void savecsv(QString filename);

    void savetxtsp(const char *filename);
    void savetxtsp(QString filename);
    void savecsvsp(const char *filename);
    void savecsvsp(QString filename);

    void savetxtsp(const char *filename,int indsp);
    void savetxtsp(QString filename,int indsp);
    void savecsvsp(const char *filename,int indsp);
    void savecsvsp(QString filename,int indsp);

	void setAcc(double *a);
    void modifyAcc(double *a, int ind1, int ind2);

    void a2vd(bool raw=false, bool rat=false);
    void detrend(int method, int oh=2, int ol=0, bool raw=true);
    void align(int method, int ntp=4, int oh=1, int ol=0, bool raw=true, bool EWZ=false);
    void filt(int ftype, int order, double f1, double f2, bool raw=true);

    void setupSP(int NSP, int np, double P1, double P2, int dm, int sm, bool pseudo, double *Zetas);
    void setupSP(int NSP, int np, double *p, int sm, bool pseudo, double *Zetas);

    void calcSP(bool allSP=false);

    void calcSP(int i, bool allSP=false);

    void calcFFT();
    void calcPSD(double olr=0.5, bool win=true);
    void setSPT(double Tg, double amax, double scale);

    void fitSP(int i, double tol, int mit, int fm, double peak0);

    double getDR();
    double getAR();

    int getN() {return n;}
    int getNsp() {return nsp;}
    double getDt() {return dt;}
    double getV0() {return v0;}
    double getD0() {return d0;}

    double getPeakAcc() {return peak(acc,n);}
    double getPeakVel() {return peak(vel,n);}
    double getPeakDisp() {return peak(disp,n);}
    
    double getRMSAcc() {return rms(acc,n);}
    double getRMSVel() {return rms(vel,n);}
    double getRMSDisp() {return rms(disp,n);}

    double *getT()    {return t;}
    double *getAcc0() {return acc0;}
    double *getAcc()  {return acc;}
    double *getVel()  {return vel;}
    double *getDisp() {return disp;}
    double *getTd()   {return td;}
    double *getTv()   {return tv;}
    double *getTa()   {return ta;}
    double *getRd()   {return rd;}
    double *getRv()   {return rv;}
    double *getRa()   {return ra;}
    double *getZeta()   {return zeta;}
    double *getIa()   {return Ia;}
    double *getIv()   {return Iv;}
    double *getId()   {return Id;}
    double *getFreqs()   {return freqs;}
    double *getAmpf()   {return ampf;}
    double *getAngf()   {return angf;}
    double *getDAngf()   {return dangf;}
    double *getPsd()   {return psd;}
    double *getFpsd()   {return fpsd;}
    double *getRF();
    double **getTHData();
    double **getSPData(int i);
    double ***getSPData();

    double **getEnergy();

    Spectra *getSP() {return sp;}
    Spectra *getSP(int i) {return sp+i;}

    void response(double zeta, double P, int method=2);
    void responseNL(double zeta, double P, int method, double *cp);

    QVector<double> qGetT()    {return A2QV(t,n);}
    QVector<double> qGetAcc0() {return A2QV(acc0,n);}
    QVector<double> qGetAcc()  {return A2QV(acc,n);}
    QVector<double> qGetVel()  {return A2QV(vel,n);}
    QVector<double> qGetDisp() {return A2QV(disp,n);}
    QVector<double> qGetTa()   {return A2QV(ta,n);}
    QVector<double> qGetTv()   {return A2QV(tv,n);}
    QVector<double> qGetTd()   {return A2QV(td,n);}
	QVector<double> qGetRa()   {return A2QV(ra,n);}
    QVector<double> qGetRv()   {return A2QV(rv,n);}
    QVector<double> qGetRd()   {return A2QV(rd,n);}
    QVector<double> qGetIa()   {return A2QV(Ia,n);}
    QVector<double> qGetIv()   {return A2QV(Iv,n);}
    QVector<double> qGetId()   {return A2QV(Id,n);}
    QVector<double> qGetFreqs()  { return A2QV(freqs, nfft / 2 + 1); }
	QVector<double> qGetAmpf()   { return A2QV(ampf, nfft / 2 + 1); }
    QVector<double> qGetAngf()   { return A2QV(angf, nfft / 2 + 1); }
    QVector<double> qGetDAngf()   { return A2QV(dangf, nfft / 2 + 1); }
    QVector<double> qGetPsd()   { return A2QV(psd, npsd/2 + 1); }
    QVector<double> qGetFpsd()   { return A2QV(fpsd, npsd/2 + 1); }

    QVector<double> qGetRF()   {return A2QV(getRF(),n);}

    Response getRes() {return res;}

    void endAlign(int ntp=8, bool raw=true, int IZC=1);
private:
    int n, nsp, nfft, npsd;
    double dt, v0, d0;
    double *t, *acc, *vel, *disp;
    double *acc0, *Ia, *Iv, *Id;
    double *td, *tv, *ta;
    double *rd, *rv, *ra;
	Spectra *sp;
    double *zeta;
    std::complex<double> *af;
    double *freqs, *ampf, *angf, *dangf, *psd, *fpsd;
    Response res;
    void reallocateTH();

};

#endif // EQSIGNAL_H
