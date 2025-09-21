#pragma once
#include <math.h>
#include <fstream>
#include <iomanip>
#include "funcm.hpp"
#include "quartenion.hpp"


class model_t {
   
protected:
	Matrix<long double> res;
	long double sample_inc, t0, t1;
	Vector<long double> x0;
public:
	model_t(const Vector< long double>& vec,long double t0,long double t1,long double inc);

	void load_res2file(const char* filename);
	double get_t0() const noexcept { return t0; };
	double get_t1() const noexcept { return t1; };
	double get_step() const noexcept { return sample_inc; };
	Vector<long double> get_init() const noexcept { return x0; };
	Matrix<long double> get_result() const noexcept { return res; };

	virtual void add_result(const Vector<long double>& X, double t);
	virtual Vector<long double> get_right(const Vector<long double>& X, long double t) const;
};

class earth_move_model : public model_t {
protected:
	const double mu_s = 132712.43994e15;
public:
	earth_move_model(const Vector<long double>& vec, long double t0, long double t1, long double inc);

	Vector<long double> get_right(const Vector<long double>& X, long double t) const override;
};

// l = 1m
class sundial_model : public earth_move_model {
protected:
	const double Re = 6371300;
	const double Ω = 7.292115e-5;
	double φ, λ, date;
	double s_0;
public:
	sundial_model(double φ_, double λ_, double date_);
	
	double get_siderial_time(double Y, double M, double D, double h, double m, double s) const noexcept;

	void add_result(const Vector<long double>& X, double t) override;
};

class blag_time_model : public earth_move_model {
protected:
	enum day_state {
		sunrise,
		sunset,
	};

	const double Re = 6371300;
	const double Ω = 7.292115e-5;
	const double φ = 0.97302106;
	const double λ = 0.65624380;
	const double s_0 = 1.75659;
	const int UTC_n = 3;

	double time_v = 0.0;
	double time_z = 0.0;
	day_state state = day_state::sunset;
public:
	blag_time_model();

	void add_result(const Vector<long double>& X, double t) override;
};