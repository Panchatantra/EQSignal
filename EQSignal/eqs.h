#ifndef NUMEXT_H
#define NUMEXT_H

#include <stddef.h>
#include <math.h>
#include <complex>
#include <QVector>
#include <QList>

#define PI  3.141592653589793
#define PI2 9.869604401089358

using namespace std;

extern "C"
{
    extern void acc2vd(double *a, double *v, double *d, int *n, double *dt, double *v0, double *d0);
    extern void ariasintensity(double *a,double *Ia,int *n);
    extern void energy(double *acc, int *n, double *dt, double *zeta, double *P, double *ra, double *rv, double *rd, double *Ek, double *Ez, double *Es, double *Eh, double *Ein, double *ku);
    extern void fft(double *in, std::complex<double> *out, int *n);
    extern void fftfreqs(int *Nfft,double *fs,double *freqs);
    extern void fftresample(double *a,int *n,int *r,double *ar,int *nr);
    extern void fitspectrum(double *acc,int *n,double *dt,double *zeta,double *P,int *nP,double *SPT,double *a,double *tol,int *mit,int *fm, int *kpb);
    extern void ifft(std::complex<double> *in, double *out, int *n);
    extern void polyblc(double *a, int *n, int *oh, int *ol, double *dt, double *v0, double *d0);
    extern void polydetrend(double *a, int *n, int *oh, int *ol);
	extern void r(double *acc, int *n, double *dt, double *zeta, double *P, double *ra, double *rv, double *rd, double *ku, int *SM);
	extern void rfreq(double *acc, int *n, double *dt, double *zeta, double *P, double *ra, double *rv, double *rd);
	extern void rmixed(double *acc, int *n, double *dt, double *zeta, double *P, double *ra, double *rv, double *rd);
	extern void rnmk(double *acc, int *n, double *dt, double *zeta, double *P, double *ra, double *rv, double *rd);
    extern void ratacc2vd(double *a, double *v, double *d, int *n, double *dt, double *v0, double *d0);
	extern void rnl(double *acc, int *n, double *dt, double *zeta, double *P, double *ra, double *rv, double *rd, double *ku, int *SM, double *cp);
    extern void spectrum(double *acc,int *n,double *dt,double *zeta,double *P,int *np,double *SPA,int *SPI,int *SM);
    extern void spectrumavd(double *acc, int *n, double *dt, double *zeta, double *P, int *np, double *SPA, int *SPI, double *SPV, double *SPD, double *SPE, int *SM);
    // ÑÓÐÔ·´Ó¦Æ×
	extern void spmu(double *acc,int *n,double *dt,double *zeta,double *P,int *nP,double *SPA,int *SPI,double *SPV,double *SPD,double *SPE,double *mu,int *model,
                     double *uy,double *rk,double *alpha);
    extern void targetdc(double *a, double *td, int *n, int *tp, int *ntp, int *ph, int *pl, double *dt, double *v0, double *d0);
    extern void targetdvac(double *a, double *td, double *tv, double *ta, int *n, int *tp, int *ntp, int *ph, int *pl, double *dt, double *v0, double *d0);
    extern void targetdvc(double *a, double *td, double *tv, int *n, int *tp, int *ntp, int *ph, int *pl, double *dt, double *v0, double *d0);
    extern void welch(double *a,int *n,int *m,double *olr,double *psd,int *win);

	extern void bwfilt(double *acc, int n, double dt, int ftype, int order, double f1, double f2);

    extern void spectrum_endur(double *acc, int *n, double *dt, double *zeta, double *P, int *np, int *DI, int *nd, double *SPA, int *SPI, int *SM);
	extern void initartwave(double *a, int *n, double *dt, double *zeta, double *P, int *nP, double *SPT);
	extern void whitenoise(double *a, int *n, double *dt);
    extern void adjustspectra_endur(double *acc, int *n, double *dt, double *zeta, double *P, int *nP, int *DI, int *nD, double *SPAT1, double *a, double *tol, int *mit, int *kpb);
	extern void adjustspectra_md(double *acc, int *n, double *dt, double *zeta, double *P, int *nP, int *nD, double *SPAT, int *nT, double *a, double *tol, int *mit, int *kpb);
}

int static inline nextpow2(int n)
{
	int m = 1;

	while (m<n)
	{
		m *= 2;
	}

	return m;
}

int static inline nextpow(int n, int base=2)
{
	int m = 1;

	while (m<n)
	{
		m *= base;
	}

	return m;
}

bool static inline iseven(int n)
{
	return (n%2 == 0);
}

template<typename T>
T static inline mean(T *a, int n)
{
	double m = 0.0;
	for (int i = 0; i < n; i++) m += a[i];
	return m / (T)n;
}

template<typename T>
T static inline rms(T *a, int n)
{
	T m = 0.0;
	for (int i = 1; i < n; i++) m += (a[i] * a[i] + a[i-1] * a[i-1]);
	m *= 0.5;
	return sqrt(m / (T)(n - 1));
}

template<typename T=double>
T static inline peak(T *a, int n)
{
    T m = 0.0;
	for (int i = 0; i < n; i++)
	{
		if (fabs(a[i]) > fabs(m)) m = a[i];
	}
	return m;
}

template<typename T=double>
int static inline peakLoc(T *a, int n)
{
    int m = 0;
    for (int i = 0; i < n; i++)
    {
		if (fabs(a[i]) > fabs(a[m])) m = i;
    }
    return m;
}

template<typename T = double>
T static inline max(T *a, int n)
{
	T m = 0.0;
	for (int i = 0; i < n; i++)
	{
		if ( a[i] > m ) m = a[i];
	}
	return m;
}

template<typename T = double>
T static inline min(T *a, int n)
{
	T m = 0.0;
	for (int i = 0; i < n; i++)
	{
		if (a[i] < m) m = a[i];
	}
	return m;
}

template<typename T = double>
void static inline fill(T *x, int n, T val)
{
	for (int i = 0; i < n; i++)
	{
		x[i] = val;
	}
}

void static inline autoScale(double *x, int n, double &rl, double &ru, double &dx)
{
	double xmin = min(x, n);
	double xmax = max(x, n);

    double peak = fmax(xmax, -xmin);
	double base = pow(10.0, floor(log10(peak)))*0.25;

    rl = -base;
    ru =  base;

    while ( rl > xmin ) rl -= base;
    while ( ru < xmax ) ru += base;

	dx = base;
	if ( (int)ceil(peak / dx) < 2 ) { dx *= 0.4; base *= 0.4; }
	while ( (int)round(peak / dx) > 3 ) dx += base;
	 
}

void static inline normalize(double *a, int n)
{
	double m = peak(a, n);
	for (int i = 0; i < n; i++) a[i] = a[i] / fabs(m);
}

void static inline scale(double *a, int n, double factor)
{
    for (int i = 0; i < n; i++) a[i] = a[i] * factor;
}

double static inline *scaled(double *a, int n, double factor)
{
    double *b = new double[n];
    for (int i = 0; i < n; i++) b[i] = a[i] * factor;
    return b;
}

template <typename T = double>
T static inline *zeros(int n)
{
    T *res = new T[n];
	for (int i = 0; i < n; i++)
	{
		res[i] = 0.0;
	}
	return res;
}

template <typename T=double>
T static inline *ones(int n)
{
    T *res = new T[n];
	for (int i = 0; i < n; i++)
	{
		res[i] = 1.0;
	}
	return res;
}

double static inline *linspace(double s1, double s2, int n)
{
	double *res = new double[n];
	double ds = (s2 - s1) / (double)(n - 1);
	for (int i = 0; i < n; i++)
	{
		res[i] = s1 + i*ds;
	}
	return res;
}

double static inline *logspace(double s1, double s2, int n)
{
	double *res = new double[n];
	double ls1 = log10(s1);
	double ls2 = log10(s2);
	double lds = (ls2 - ls1) / (double)(n - 1);
	double tmp;

	for (int i = 0; i < n; i++)
	{
		tmp = ls1 + i*lds;
		res[i] = pow(10.0, tmp);
	}
	return res;
}

template <typename T=double>
T static inline *arraySlice(T *a, int ind1, int ind2)
{
    int n = ind2-ind1+1;
    T *res = new T[n];

#pragma omp parallel for
	for (int i = ind1; i <= ind2; ++i)
    {
		res[i - ind1] = a[i];
    }
    return res;
}

template <typename T=double>
QVector<T> static inline arraySlice(QVector<T> a, int ind1, int ind2)
{
    int n = ind2-ind1+1;
    QVector<T> res = QVector<T>(n);

#pragma omp parallel for
    for (int i = ind1; i <= ind2; ++i)
    {
        res[i - ind1] = a[i];
    }
    return res;
}

template <typename T = double>
void static inline arrayCopy(T *a, T *b, int n)
{
#pragma omp parallel for
	for (int i = 0; i <= n; ++i)
	{
		b[i] = a[i];
	}
}


template <typename T=double>
QVector<T> static inline A2QV(T *a, int n)
{
    QVector<T> qv(n);
    qCopy(a,a+n,qv.begin());

    return qv;
}

template <typename T = double>
void static inline peakScale(T *a, int n, T peak0)
{
    T p = fabs(peak(a,n));

#pragma omp parallel for
	for (int i = 0; i < n; i++)
	{
		a[i] = a[i]*peak0/p;
	}

}

template <typename T = double>
void static inline hist(QVector<T> a, T rl, T ru,
                        int N, QVector<T> &x, QVector<T> &y)
{
    T dx = (ru - rl)/(double)N;

    for (int I=0; I<N; I++) x[I] = rl + dx*0.5 + dx*I;

    int indx;

    for (int i=0; i<a.count(); i++)
    {
        indx = (int)floor((a[i]-rl)/dx);
        y[indx] += 1;
    }
}

template<typename T>
void static inline bilinearDetrend(T *a, int n)
{
    T slope1, slope2, errMin;
    int I;
    T base;
    T *err = new T[n];

    err[0] = 0.0;
    slope2 = (a[n-1] - a[0])/(n-1);

    for (int j=0; j<n; j++)
    {
        base = a[0] + slope2*j;
        err[0] += (base-a[j])*(base-a[j]);
    }

    err[0] = sqrt(err[0]/n);

#pragma omp parallel for private(slope1,slope2,base) shared(a,n,err)
    for (int i=1; i<n; i++)
    {

        slope1 = (a[i] - a[0])/(i);
        slope2 = (a[n-1] - a[i])/(n-1-i);

        err[i] = 0.0;
        for (int j=0; j<i; j++)
        {
            base = a[0] + slope1*(j);
            err[i] += (base-a[j])*(base-a[j]);
        }
        for (int j=i; j<n; j++)
        {
            base = a[i] + slope2*(j-i);
            err[i] += (base-a[j])*(base-a[j]);
        }

        err[i] = sqrt(err[i]/n);
    }

    errMin = err[0];
    I = 0;
    for (int i=1; i<n; i++)
    {
        if (err[i]<errMin)
        {
            errMin = err[i];
            I = i;
        }
    }

    double a0 = a[0];
    double ai = a[I];
    slope1 = (ai - a0)/(I);
    slope2 = (a[n-1] - ai)/(n-1-I);

    for (int j=0; j<I; j++)
    {
        a[j] -= a0 + slope1*(j);
    }
    for (int j=I; j<n; j++)
    {
        a[j] -= ai + slope2*(j-I);
    }

}

template<typename T>
void static inline endLinearDetrend(T *a, int n, int IZC=1)
{
    int I0 = n-1;
    int ZC = IZC;
    for (int i=n-1; i>=0; i--)
    {
        if (a[i]*a[i-1]<0.0)
        {
            I0 = i;
            ZC--;
        }
        if ( ZC==0 ) break;
    }

    if (I0 == n-1) {
        a[I0] = 0.0;
        return;
    }

    double ts = n-1-I0;
    double slp = (a[n-1])/ts;

    for (int i = I0; i < n; i++)
    {
        a[i] -= slp*(i-I0);
    }
}

template<typename T>
void static inline beginLinearDetrend(T *a, int n, int IZC=1)
{
    int I0 = 0;
    int ZC = IZC;
    for (int i=0; i<n; i++)
    {
        if (a[i]*a[i+1]<0.0)
        {
            I0 = i;
            ZC--;
        }
        if ( ZC==0 ) break;
    }

    if (I0 == 0) {
        a[0] = 0.0;
        return;
    }

    double ts = I0;
    double slp = (a[0])/ts;

    for (int i = 0; i <= I0; i++)
    {
        a[i] -= slp*(I0);
    }
}

#endif // NUMEXT_H

