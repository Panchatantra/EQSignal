#ifndef NUMEXT_H
#define NUMEXT_H

#include <stddef.h>
#include <math.h>
#include <complex>
#include <QtCore/QVector>

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
    extern void fitspectrum(double *acc,int *n,double *dt,double *zeta,double *P,int *nP,double *SPT,double *a,double *tol,int *mit,int *fm);
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
    extern void spmu(double *acc,int *n,double *dt,double *zeta,double *P,int *nP,double *SPA,int *SPI,double *SPV,double *SPD,double *SPE,double *mu,int *model,
                     double *rtol,int *maxiter,double *uy,double *rk,double *alpha);
    extern void targetdc(double *a, double *td, int *n, int *tp, int *ntp, int *ph, int *pl, double *dt, double *v0, double *d0);
    extern void targetdvac(double *a, double *td, double *tv, double *ta, int *n, int *tp, int *ntp, int *ph, int *pl, double *dt, double *v0, double *d0);
    extern void targetdvc(double *a, double *td, double *tv, int *n, int *tp, int *ntp, int *ph, int *pl, double *dt, double *v0, double *d0);
    extern void welch(double *a,int *n,int *m,double *olr,double *psd,int *win);

	extern void bwfilt(double *acc, int n, double dt, int ftype, int order, double f1, double f2);
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

template <typename T=double>
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

template <typename T = double>
void static inline *arrayCopy(T *a, T *b, int n)
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
void static inline *peakScale(T *a, int n, T peak0)
{
    T p = fabs(peak(a,n));

#pragma omp parallel for
	for (int i = 0; i < n; i++)
	{
		a[i] = a[i]*peak0/p;
	}

}

#endif // NUMEXT_H

