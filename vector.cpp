#include "vector.hpp"
#include "quartenion.hpp"

template<typename T>
Quartenion Vector<T>::operator*(const Quartenion& quar) const {
	Quartenion q(0, _data.at(0), _data.at(1), _data.at(2));
	return q * quar;
}

template<typename T>
Vector<T> Vector<T>::rotate(double phi, const Vector<T>& axis) const {
	Quartenion quat(phi, axis);
	Quartenion output = quat * (*this) * (quat.conj());
	return output.vec();
}

template<typename T>
Vector<T> Vector<T>::rotateByQuartenion(const Quartenion& L) const{
	Quartenion temp(L);
	temp.normalization();
	return (temp*(*this)*(!temp)).vec();
};

template class Vector<double>;