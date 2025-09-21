#pragma once
#include "vector.hpp"

class Quartenion 
{
private:
	double q;
	Vector<double>* Q;
public:
	Quartenion();
	Quartenion(double l0, double l1, double l2, double l3);
	Quartenion(double phi,const Vector<double>& e);
	Quartenion(const Quartenion& quar);
	Quartenion& operator=(const Quartenion& quar);
	~Quartenion();

	double scal() const;
	Vector<double> vec() const;
	Matrix<double> toRotateMatrix() const;
	Quartenion& normalization();
	Quartenion conj() const;

	Quartenion operator-() const;
	Quartenion operator+(const Quartenion& quar) const;
	Quartenion operator-(const Quartenion& quar) const;
	Quartenion operator*(const Quartenion& quar) const;
	Quartenion operator*(const Vector<double>& vec) const;
	Quartenion operator*(double s) const;
	Quartenion operator!() const;
	friend Quartenion operator*(double s, const Quartenion& quar);
	friend std::ostream& operator<<(std::ostream& out, const Quartenion& quar);
};
