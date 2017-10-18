#ifndef SPECTRA_H
#define SPECTRA_H

#include "eqs.h"
#include <QtCore/QVector>

class Spectra
{
public:
	Spectra();
    Spectra(int np, double P1, double P2, int dm=0, int sm=1, bool pseudo=false, double Zeta=0.05);
    Spectra(int np, double *p, int sm=1, bool pseudo=false, double Zeta=0.05);
    //	~Spectra();

	int getNP()       { return NP; }
	double getZeta()  { return zeta; }
	double *getP()    { return P; }
	double *getSPA()  { return SPA; }
	double *getSPV()  { return SPV; }
	double *getSPD()  { return SPD; }
	double *getSPE()  { return SPE; }
    double *getSPT() { return SPT; }
    double **getSPData();

    void setP(double *p, int np);
    void setNP(int np);
    void setSM(int sm = 1) {SM = sm;}
    void setDM(int dm = 0) {DM = dm;}
    void setZeta(double z = 0.05) {zeta = z;}

    void setSPT(double *spt);
    void setSPT(double *p, double *spt, int np);
    void setSPT(double Tg, double PAF=2.25, double scale=1.0);
    
    void setAcc(double *Acc, int N, double DT);
	void calc(bool all = false);

	void calcNL(double mu, int model, double rk, double alpha);

    double *fitSP(double tol, int mit, int fm, double peak0, int kpb=1);

    void fitError(double &Emax, double &Emean, double &CV);

    QVector<double> qGetP()    {return A2QV(P,NP);}
    QVector<double> qGetSPA()  {return A2QV(SPA,NP);}
    QVector<double> qGetSPV()  {return A2QV(SPV,NP);}
    QVector<double> qGetSPD()  {return A2QV(SPD,NP);}
    QVector<double> qGetSPE()  {return A2QV(SPE,NP);}
    QVector<double> qGetSPT()  {return A2QV(SPT,NP);}

private:
    double *P, *SPA, *SPV, *SPD, *SPE, *SPT;
    int *SPI;
	int NP, SM, DM;
    double zeta, p1, p2;
    int n;
    double *acc;
    double dt;

    void reallocate();
};

#endif // SPECTRA_H
