#include "integrator.hpp"

void DormandPrinceIntegrator::run(model_t& system) {

	//clock_t start_time = clock();

	long double h = 1e-5l;
	long double h_new = h;
	long double t0 = 0;
	long double t1 = 10;
	long double step = 0.00001;
	long double t = t0 + step;
	Vector<long double> x0 = system.get_init();

	long double v{ 1 };
	long double u;
	while (1 + v > 1) {
		u = v;
		v /= 2;
	}

	while (t0 < t1) {
		h = h_new;

		k.at(0) = system.get_right(x0, t0);
		k.at(1) = system.get_right(x0 + h * a(1, 0) * k.at(0), t0 + c.at(1) * h);
		k.at(2) = system.get_right(x0 + h * (a(2, 0) * k(0) + a(2, 1) * k(1)), t0 + c.at(2) * h);
		k.at(3) = system.get_right(x0 + h * (a(3, 0) * k(0) + a(3, 1) * k(1) + a(3, 2) * k(2)), t0 + c.at(3) * h);
		k.at(4) = system.get_right(x0 + h * (a(4, 0) * k(0) + a(4, 1) * k(1) + a(4, 2) * k(2) + a(4, 3) * k(3)), t0 + c.at(4) * h);
		k.at(5) = system.get_right(x0 + h * (a(5, 0) * k(0) + a(5, 1) * k(1) + a(5, 2) * k(2) + a(5, 3) * k(3) + a(5, 4) * k(4)), t0 + c.at(5) * h);
		k.at(6) = system.get_right(x0 + h * (a(6, 0) * k(0) + a(6, 1) * k(1) + a(6, 2) * k(2) + a(6, 3) * k(3) + a(6, 4) * k(4) + a(6, 5) * k(5)), t0 + c.at(6) * h);

		Vector<long double> x1 = x0 + h * (b(0) * k(0) + b(1) * k(1) + b(2) * k(2) + b(3) * k(3) + b(4) * k(4) + b(5) * k(5));
		Vector<long double> x1_ = x0 + h * (b1(0) * k(0) + b1(1) * k(1) + b1(2) * k(2) + b1(3) * k(3) + b1(4) * k(4) + b1(5) * k(5) + b1(6) * k(6));

		long double new_eps{ 0 };

		for (uint64_t count = 0u; count < x0.dimension(); ++count) {
			double max = std::max({ 1e-5l, abs(x0.at(count)), abs(x1.at(count)), 2.0l * u / eps });
			new_eps += pow(h * (x1.at(count) - x1_.at(count)) / max, 2.0l);
		}
		new_eps /= x0.dimension();
		new_eps = sqrt(new_eps);

		h_new = h / std::max(0.1l, std::min(5.0l, pow(new_eps / eps, (1.0l / 5.0l)) / 0.9l));

		if (new_eps > eps)
			continue;

		while ((t < t0 + h) && (t <= t1 + step)) {
			Vector<long double> d(6);
			long double theta = (t - t0) / h;

			d.at(0) = theta * (1 + theta * (-1337.0l / 480.0l + theta * (1039.0l / 360.0l + theta * (-1163.0l / 1152.0l))));
			d.at(1) = 0;
			d.at(2) = 100.0l * pow(theta, 2.0l) * (1054.0l / 9275.0l + theta * (-4682.0l / 27825.0l + theta * (379.0l / 5565.0l))) / 3.0l;
			d.at(3) = -5.0l * pow(theta, 2.0l) * (27.0l / 40.0l + theta * (-9.0l / 5.0l + theta * (83.0l / 96.0l))) / 2.0l;
			d.at(4) = 18225.0l * pow(theta, 2.0l) * (-3.0l / 250.0l + theta * (22.0l / 375.0l + theta * (-37.0l / 600.0l))) / 848.0l;
			d.at(5) = -22.0l * pow(theta, 2.0l) * (-3.0l / 10.0l + theta * (29.0l / 30.0l + theta * (-17.0l / 24.0l))) / 7.0l;

			Vector<long double> res = x0 + h * (d(0) * k(0) + d(1) * k(1) + d(2) * k(2) + d(3) * k(3) + d(4) * k(4) + d(5) * k(5));
			system.add_result(res, t);
			t += step;
		}

		t0 += h;
		x0 = x1_;
	}

	//clock_t end_time = clock();

	//std::cout << "Elapsed time: " << -(start_time - end_time) / 1000000.0 << " s" << '\n';
}