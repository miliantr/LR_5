#include "model.hpp"

model_t::model_t(const Vector<long double>& vec,long double t0,long double t1,long double inc) : x0(vec), sample_inc(inc), t0(t0), t1(t1) {};

void model_t::add_result(const Vector<long double>& X, double t) {
	res.push_row(X);
}

void model_t::load_res2file(const char* filename) {
	std::ofstream f(filename, std::ios::trunc);

	if (!f.is_open())
		throw std::logic_error("load res2file");
	else
		std::cout << "Create file" << '\n';

	for (uint64_t count = 0u; count < res.rows(); ++count)
		f << res.get_rows(count) << '\n';

	f.close();
}

Vector<long double> model_t::get_right(const Vector<long double>& X, long double t) const {
	Vector<long double> dX(X.dimension());

	long double mu = 0.012277471l;
	long double mu_ = 1 - mu;
	long double D1 = pow((X.at(0) + mu) * (X.at(0) + mu) + X.at(2) * X.at(2), 3.0l / 2.0l);
	long double D2 = pow((X.at(0) - mu_) * (X.at(0) - mu_) + X.at(2) * X.at(2), 3.0l / 2.0l);

	dX.at(0) = X.at(1);
	dX.at(1) = X.at(0) + 2 * X.at(3) - mu_ * (X.at(0) + mu) / D1 - mu * (X.at(0) - mu_) / D2;
	dX.at(2) = X.at(3);
	dX.at(3) = X.at(2) - 2 * X.at(1) - mu_ * X.at(2) / D1 - mu * X.at(2) / D2;

	return dX;
};

earth_move_model::earth_move_model(const Vector<long double>& vec, long double t0, long double t1, long double inc) : model_t(vec, t0, t1, inc) {};

Vector<long double> earth_move_model::get_right(const Vector<long double>& X, long double t) const {
	Vector<long double> dX(X.dimension());

	long double modul = sqrt(pow(X.at(0), 2) + pow(X.at(1), 2) + pow(X.at(2), 2));

	dX.at(0) = X.at(3);
	dX.at(1) = X.at(4);
	dX.at(2) = X.at(5);
	dX.at(3) = -mu_s * X.at(0) / pow(modul, 3);
	dX.at(4) = -mu_s * X.at(1) / pow(modul, 3);
	dX.at(5) = -mu_s * X.at(2) / pow(modul, 3);

	return dX;
};

sundial_model::sundial_model(double φ_, double λ_, double date_) : φ(φ_), λ(λ_), date(date_),
earth_move_model(Vector<long double>({ -2.6005047996994e10, 1.32621705709054e11, 5.7523888683657e10, -2.9832953e4, -4.715287e3, -2.043123e3 }), 2460310.50 * 86400.0, (date_+ 1.0) * 86400.0, 60.0)
{
	s_0 = get_siderial_time(2024, 1, 1, 0, 0, 0);
	s_0 = wrap_angle(2 * math_const::π * s_0 / 86400.0); // угол ориентации гринвичского меридиана 
};

double sundial_model::get_siderial_time(double Y, double M, double D, double h, double m, double s) const noexcept { // время звёздное

	double JD = get_JDN(Y, M, D, h, m, s);

	double t_star = int(JD - 2451544.5) / 36525.0;//t*

	double sg_0 = 24110.54841 + 8640184.812866 * t_star + 0.093104 * pow(t_star, 2) - 6.2e-6 * pow(t_star, 3);//гринвическое начальное звёздное время места

	return sg_0;
};

void sundial_model::add_result(const Vector<long double>& X, double t) {
	using namespace math_const;

	//test day
	if (t <= date * 86400.0)
		return;

	//part 0
	t -= get_t0();

	//part 1
	Vector<long double> coord{ 3 };

	double s = s_0 + Ω * t + λ;

	s = wrap_angle(s);

	coord.at(0) = Re * cos(φ) * cos(s);
	coord.at(1) = Re * cos(φ) * sin(s);
	coord.at(2) = Re * sin(φ);

	auto ort_r = coord / coord.length();//нормированный вектор

	//part 2
	Vector<long double> earth_r({ X.at(0), X.at(1), X.at(2) });
	auto ort_earth_r = earth_r / earth_r.length();

	double angle = acos(ort_earth_r * ort_r);

	if (angle <= π / 2)
		return;
	else {
		angle -= π;
	}

	//part 3
	auto earth_r_star = (-1 / (ort_earth_r * ort_r)) * ort_earth_r;

	//part 4
	auto vec_shadow = ort_r + earth_r_star;

	Matrix<long double> M{ 3,3,{-sin(φ) * cos(s), -sin(φ) * sin(s), cos(φ), -cos(φ) * cos(s), -cos(φ) * sin(s), -sin(φ), -sin(s), cos(s), 0} };

	auto finally = M * vec_shadow;

	auto vec_r = finally.concat(Vector<long double>({ -angle, t - (date - 2460310.50) * 86400.0 }));

	//add result
	res.push_row(vec_r);
};

blag_time_model::blag_time_model() :
	earth_move_model(Vector<long double>({ -2.6005047996994e10, 1.32621705709054e11, 5.7523888683657e10, -2.9832953e4, -4.715287e3, -2.043123e3 }), 2460310.50 * 86400.0, (2460310.50 + 365.0) * 86400.0, 60.0) {};

void blag_time_model::add_result(const Vector<long double>& X, double t) {
	using namespace math_const;
	

	//part 0
	t -= get_t0();

	int day = t / 86400;

	//part 1
	Vector<long double> coord{ 3 };

	double s = s_0 + Ω * t + λ;

	s = wrap_angle(s);

	coord.at(0) = Re * cos(φ) * cos(s);
	coord.at(1) = Re * cos(φ) * sin(s);
	coord.at(2) = Re * sin(φ);

	auto ort_r = coord / coord.length();

	//part 2
	
	Vector<long double> earth_r({ X.at(0), X.at(1), X.at(2) });
	auto ort_earth_r = earth_r / earth_r.length();

	double angle = acos(ort_earth_r * ort_r);

	double time{};

	if (t - day * 86400.0 + UTC_n * 3600.0 > 86400.0)
		time = t - day * 86400.0 + UTC_n * 3600.0 - 86400.0;
	else
		time = t - day * 86400.0 + UTC_n * 3600.0;

	if (angle < π / 2) {
		if (state == day_state::sunrise) {
			state = day_state::sunset;

			time_z = time;

			//add result
			res.push_row(Vector<long double>({ time_v, time_z}));
		}
		return;
	}
	else {
		if (state == day_state::sunset) {
			state = day_state::sunrise;

			time_v = time;

		}
		angle -= π;
	}
}
