#pragma once

#include <ctime>
#include "model.hpp"


class Integrator {
protected:
	long double eps = 1e-8l;
public:
	Integrator(long double eps) : eps(eps) {};
	virtual void run(model_t& system)=0;
};

class DormandPrinceIntegrator : public Integrator {
protected:
	const Vector<long double> c = Vector<long double>({
			0.0l,	0.2l, 0.3l, 0.8l, 8.0l/9.0l, 1.0l, 1.0l
		});

	const Matrix<long double> a = Matrix<long double>(7, 7, {
			0,				0,					0,				0,				0,					0,			0,
			1.0l/5.0l,		0,					0,				0,				0,					0,			0,
			3.0l/40.0l,		9.0l/40.0l,			0,				0,				0,					0,			0,
			44.0l/45.0l,		-56.0l/15.0l,			32.0l/9.0l,		0,				0,					0,			0,
			19372.0l/6561.0l, -25360.0l/2187.0l,	64448.0l/6561.0l, -212.0l/729.0l,	0,					0,			0,
			9017.0l/3168.0l,	-355.0l/33.0l,		46732.0l/5247.0l,	49.0l/176.0l,		-5103.0l/18656.0l,	0,			0,
			35.0l/384.0l,		0,					500.0l/1113.0l,	125.0l/192.0l,	-2187.0l/6784.0l,		11.0l/84.0l,	0
		});

	const Vector<long double> b = Vector<long double>({
			35.0l / 384.0l,	0,	500.0l / 1113.0l,	125.0l / 192.0l,	-2187.0l / 6784.0l,	11.0l / 84.0l,	0
		});

	const Vector<long double> b1 = Vector<long double>({
			5179.0l/57600.0l,	0,	7571.0l/16695.0l,	393.0l/640.0l,	-92097.0l/339200.0l,	187.0l/2100.0l,	1.0l/40.0l
		});

	Vector<Vector<long double>> k{ 7 };
public:
	DormandPrinceIntegrator(long double eps) : Integrator(eps) {};
	
	virtual void run(model_t& system) override;
};