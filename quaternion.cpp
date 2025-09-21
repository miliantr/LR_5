#include "quaternion.hpp"

Quartenion::Quartenion() :q(0) {Q = new Vector<double>({ 0, 0, 0 });}

Quartenion::Quartenion(double l0, double l1, double l2, double l3): q(l0){Q = new Vector<double>({ l1,l2,l3 });}

Quartenion::Quartenion(double phi,const Vector<double>& e) : q(cos(phi / 2))
{
	Vector<double> temp{ e };
	Q= new Vector<double>(temp.normalization() * sin(phi / 2));
}

Quartenion::Quartenion(const Quartenion& quar): q(quar.q)
{
	Q = new Vector<double>({0,0,0});
	(*Q) = quar.vec();
}

Quartenion& Quartenion::operator=(const Quartenion& quar)
{
	if (this != &quar) {
		q = quar.q;
		(*Q) = quar.vec();
	}

	return *this;
}

Quartenion::~Quartenion(){delete Q;}

Quartenion Quartenion::operator-() const
{
	Vector<double> temp = -vec();
	return Quartenion(-q, temp.at(0), temp.at(1), temp.at(2));
};

Quartenion Quartenion::operator+(const Quartenion& quar) const 
{
	Vector<double> temp = quar.vec() + vec();
	return Quartenion(q + quar.q, temp.at(0), temp.at(1), temp.at(2));
};

Quartenion Quartenion::operator-(const Quartenion& quar) const {return (*this) + (-quar);};

Quartenion Quartenion::operator*(const Quartenion& quar) const 
{
	double scal = this->scal() * quar.scal() - vec() * quar.vec();
	Vector<double> part = this->scal() * quar.vec() + quar.scal() * this->vec();
	Vector<double> last = part + (this->vec() ^ quar.vec());
	return { scal, last.at(0), last.at(1), last.at(2) };
};

Quartenion Quartenion::operator*(double s) const
{
	Vector<double> temp = this->vec() * s;
	return Quartenion(q*s,temp.at(0), temp.at(1), temp.at(2));
};

std::ostream& operator<<(std::ostream& out, const Quartenion& quar) 
{
	out << '(' << quar.q << ' ' << quar.Q->at(0) << ' ' << quar.Q->at(1) << ' ' << quar.Q->at(2) << ')';

	return out;
};

Quartenion operator*(double s, const Quartenion& quar){return quar*s;};

Quartenion Quartenion::operator!() const
{
	Quartenion temp = conj();
	double norm = q * q + vec() * vec();

	temp.q /= norm;
	(*temp.Q) = temp.vec() / norm;

	return temp;
};

double Quartenion::scal() const{return q;};

Vector<double> Quartenion::vec() const {return *Q;}

Quartenion& Quartenion::normalization()
{
	double norm = q * q + vec()*vec();

	q /= sqrt(norm);
	*Q = vec() / sqrt(norm);

	return *this;
}

Quartenion Quartenion::conj() const
{
	Vector<double> temp = vec();
	temp = -temp;
	return Quartenion(q,temp.at(0), temp.at(1), temp.at(2));
}

Matrix<double> Quartenion::toRotateMatrix() const
{
	Matrix<double> output(3, 3);
	
	double norm = sqrt(q * q + vec() * vec());

	double q0 = q / norm;
	double q1 = Q->at(0)/norm;
	double q2 = Q->at(1)/norm;
	double q3 = Q->at(2)/norm;

	output.at(0, 0) = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
	output.at(0, 1) = 2 * (q1 * q2 - q0 * q3);
	output.at(0, 2) = 2 * (q1 * q3 + q0 * q2);
	output.at(1, 0) = 2 * (q2 * q1 + q0 * q3);
	output.at(1, 1) = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
	output.at(1, 2) = 2 * (q2 * q3 - q0 * q1);
	output.at(2, 0) = 2 * (q3 * q1 - q0 * q2);
	output.at(2, 1) = 2 * (q3 * q2 + q0 * q1);
	output.at(2, 2) = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
	
	return output;
};


Quartenion Quartenion::operator*(const Vector<double>& vec) const 
{
	Quartenion temp(0, vec.at(0), vec.at(1), vec.at(2));
	return (*this) * temp;
};