#include "DspFilters/Dsp.h"

extern "C" __declspec(dllexport) void bwfilt(double *acc, int n, double dt, int ftype, int order, double f1, double f2)
{

	double *Data[1];
	Data[0] = acc;

	switch (ftype)
	{
	case 0:
	{
		Dsp::SimpleFilter < Dsp::Butterworth::HighPass<8>, 1 > f;
		f.setup(order, 1.0 / dt, f1);
		f.process(n, Data);
		break;
	}
	case 1:
	{
		Dsp::SimpleFilter < Dsp::Butterworth::LowPass<8>, 1 > f;
		f.setup(order, 1.0 / dt, f2);
		f.process(n, Data);
		break;
	}
	case 2:
	{
		Dsp::SimpleFilter < Dsp::Butterworth::BandPass<8>, 1 > f;
		f.setup(order, 1.0 / dt, (f1 + f2) / 2.0, (f2 - f1) / 2.0);
		f.process(n, Data);
		break;
	}
	case 3:
	{
		Dsp::SimpleFilter < Dsp::Butterworth::BandStop<8>, 1 > f;
		f.setup(order, 1.0 / dt, (f1 + f2) / 2.0, (f2 - f1) / 2.0);
		f.process(n, Data);
		break;
	}

	default:
		break;
	}

	if (order % 2 == 1) for (int i = 0; i < n; i++) acc[i] = -acc[i];
}