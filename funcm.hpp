#pragma once
#include "constants.hpp"
#include <iostream>
#include <math.h>

inline double rad(double angle) {
	return angle * math_const::π / 180.0;
}

template<typename T>
T wrap_angle(T theta)
{
	using namespace math_const;

	const T modded = fmod(theta, (T)2.0 * (T)π);
	return (modded > (T)π) ?
		(modded - (T)2.0 * (T)π) :
		modded;
}

inline double get_JDN(int Y, int M, int D, int h, int m, int s) {
	int a = (14 - M) / 12;
	int year = Y + 4800 - a;
	int month = M + 12 * a - 3;
	int JDN = D + ((153 * month + 2) / 5) + 365 * year + year / 4 - year / 100 + year / 400 - 32045;
	double JD = JDN + (h - 12) / 24.0 + m / 1440.0 + s / 86400.0;
	
	return JD;
}


inline long double Legendre(long double arg, long double n,long double m) {
	long double answer{};

	if (n == m && n != 0) {

		long double delta;
		if (m - 1 == 0)
			delta = 0.5;
		else
			delta = 1;

		answer = Legendre(arg, n - 1, m - 1) * cos(arg) * sqrt((2 * n + 1) / (2 * n * delta));
	}
	else if (n==m && n == 0) {
		answer = 1.0;
	}

	if (n > m) {
		answer = Legendre(arg, n - 1, m) * sin(arg) * sqrt((4.0 * n * n - 1.0) / (n * n - m*m)) - Legendre(arg, n - 2, m) * sqrt((((n - 1) * (n - 1) - m*m) * (2 * n + 1)) / ((n * n - m*m) * (2 * n - 3)));
	}
	else if (n < m) {
		answer = 0.0;
	}
	
	return answer;
}

