#pragma once 

#include <vector>
#include "matrix.hpp"
#include "quartenion.hpp"

class Quartenion;

template<typename T>
class Vector {
protected:
	std::vector<T> _data;
public:
	Vector();
	Vector(uint64_t size);
	Vector(const std::vector<T>& vec);
	Vector(const Vector<T>& vec);
	Vector<T>& operator=(const Vector<T>& rval);
	~Vector();
	
	T& operator[](int64_t index);
	T operator[](int64_t index) const;

	T& operator()(int64_t index);
	T operator()(int64_t index) const;
	
	T& at(const int index);
	T at(const int index) const;

	int dimension() const noexcept;
	void resize(uint64_t size);
	long double cross(const Vector<T>& vec) const;
	long double length() const noexcept;
	Vector<T> vec_cross(const Vector<T>& vec) const;
	Vector<T> rotateByRodrigFormula(double phi, Vector<T> axis) const;
	Vector<T> rotate(double phi, const Vector<T>& axis) const;
	Vector<T> rotateByQuartenion(const Quartenion& L) const;
	Vector<T>& normalization();
	Vector<T>& push_back(const T& element) noexcept;
	Vector<T>& concat(const Vector<T>& vec) noexcept;
	Vector<T>& add(const Vector<T>& vec);
	template<typename S> Vector<T>& scale_mult(const S& s) noexcept;
	template<typename S> Vector<T>& mat_mult(const Matrix<S>& mat) noexcept;

	long double operator*(const Vector<T>& vec) const;
	Vector<T> operator-() const;
	Vector<T> operator+(const Vector<T>& vec) const;
	Vector<T> operator-(const Vector<T>& vec) const;
	Vector<T> operator*(double s) const;
	Vector<T> operator*(const Matrix<T>& mat) const;
	Vector<T> operator^(const Vector<T>& vec) const;
	Quartenion operator*(const Quartenion& quar) const;
	template<typename S> Vector<T> operator/(const S& arg);

	template<typename S> friend bool operator>(const Vector<S>& lval, const Vector<S>& rval);
	template<typename S> friend bool operator<(const Vector<S>& lval, const Vector<S>& rval);
	template<typename S> friend bool operator==(const Vector<S>& lval, const Vector<S>& rval);
	template<typename S> friend Vector<S> operator*(double s, const Vector<S>& vec);
	template<typename S> friend std::ostream& operator<<(std::ostream& out, const Vector<S>& vec);
	template<typename S> friend std::ofstream& operator<<(std::ofstream& out, const Vector<S>& vec);
};

template<typename T>
Vector<T>::Vector(const std::vector<T>& vec) : _data(vec) {};

template<typename T>
Vector<T>::Vector() {
	_data.reserve(3);
};

template<typename T>
Vector<T>::Vector(uint64_t size) {
	_data = std::vector<T>(size, 0);
};

template<typename T>
Vector<T>::Vector(const Vector<T>& vec) {
	_data = vec._data;
};

template<typename T>
Vector<T>& Vector<T>::operator=(const Vector<T>& rval) {
	if (this != &rval)
		this->_data = rval._data;
	return *this;
};

template<typename T>
void Vector<T>::resize(uint64_t size) {
	_data.resize(size);
};

template<typename T>
T& Vector<T>::operator[](int64_t index) {
	return this->at(index);
};

template<typename T>
T Vector<T>::operator[](int64_t index) const {
	return this->at(index);
};

template<typename T>
T& Vector<T>::operator()(int64_t index) {
	return this->at(index);
};

template<typename T>
T Vector<T>::operator()(int64_t index) const {
	return this->at(index);
};

template<typename T>
Vector<T>::~Vector() {
	_data.clear();
};

template<typename T>
Vector<T>& Vector<T>::normalization() {
	for (auto& el : _data)
		el /= length();
	return *this;
};

template<typename T>
Vector<T>& Vector<T>::push_back(const T& element) noexcept {
	_data.push_back(element);

	return *this;
};

template<typename T>
Vector<T>& Vector<T>::concat(const Vector<T>& vec) noexcept {
	for (const auto& el : vec._data)
		_data.push_back(el);

	return *this;
};

template<typename T>
Vector<T>& Vector<T>::add(const Vector<T>& vec) {
	if (this->dimension() != vec.dimension())
		throw std::logic_error("Dimension vectors are other!");

	for (uint64_t count = 0u; count < this->dimension(); ++count)
		_data.at(count) += vec._data.at(count);

	return *this;
};

template<typename T>
long double Vector<T>::cross(const Vector<T>& vec) const {
	if (this->dimension() != vec.dimension())
		throw std::logic_error("Dimension vectors are other!");

	long double sum{};

	for (uint64_t count = 0u; count < _data.size(); ++count)
		sum += _data.at(count) * vec._data.at(count);

	return sum;
};

template<typename T>
int Vector<T>::dimension() const noexcept {
	return _data.size();
};

template<typename T>
T& Vector<T>::at(const int index) {
	if (index >= _data.size())
		throw std::out_of_range("at vec");

	return _data.at(index);
};

template<typename T>
T Vector<T>::at(const int index) const {
	if (index >= _data.size())
		throw std::out_of_range("at vec");

	return _data.at(index);
};

template<typename T>
Vector<T> Vector<T>::vec_cross(const Vector<T>& vec) const {
	if (this->dimension() != 3 || vec.dimension() != 3)
		throw std::logic_error("Vec cross vector dimension must be 3!");

	T c_1 = at(1) * vec.at(2) - at(2) * vec.at(1);
	T c_2 = at(2) * vec.at(0) - at(0) * vec.at(2);
	T c_3 = at(0) * vec.at(1) - at(1) * vec.at(0);

	return Vector<T>{ {c_1, c_2, c_3}};
};

template<typename T>
template<typename S>
Vector<T>& Vector<T>::scale_mult(const S& s) noexcept {
	for (uint64_t count = 0u; count < this->dimension(); ++count)
		this->at(count) *= s;

	return *this;
};

template<typename T>
template<typename S>
Vector<T>& Vector<T>::mat_mult(const Matrix<S>& mat) noexcept {
	Vector<T> temp{ *this };

	for (uint64_t count = 0u; count < this->dimension(); ++count) {
		this->at(count) = 0u;
		for (uint64_t col = 0u; col < this->dimension(); ++col) {
			this->at(count) += temp.at(col) * mat(col, count);
		}
	}

	return *this;
};

template<typename T>
long double Vector<T>::length() const noexcept {
	long double len{};

	for (const auto& el : _data)
		len += el*el;

	return sqrt(len);
};

template<typename T>
Vector<T> Vector<T>::operator-() const {
	Vector<T> temp{ *this };
	for (uint64_t count = 0u; count < temp.dimension(); ++count)
		temp.at(count) = -temp.at(count);
	return Vector<T>(temp);
};

template<typename T>
Vector<T> Vector<T>::operator+(const Vector<T>& vec) const {
	if (this->dimension() != vec.dimension())
		throw std::logic_error("n+m");

	Vector<T> temp{ *this };
	for (uint64_t count = 0u; count < temp.dimension(); ++count)
		temp.at(count) += vec.at(count);

	return Vector<T>(temp);
};

template<typename T>
Vector<T> Vector<T>::operator-(const Vector<T>& vec) const {
	return this->operator+(-vec);
};

template<typename T>
long double Vector<T>::operator*(const Vector<T>& vec) const {
	return this->cross(vec);
};

template<typename T>
Vector<T> Vector<T>::operator*(double s) const {
	Vector<T> temp{ *this };

	for (uint64_t count = 0u; count < temp.dimension(); ++count)
		temp.at(count) *= s;

	return Vector<T>(temp);
};

template<typename T>
Vector<T> Vector<T>::operator*(const Matrix<T>& mat) const {
	Vector<T> temp{ *this };

	for (uint64_t count = 0u; count < this->dimension(); ++count) {
		temp.at(count) = 0u;
		for (uint64_t col = 0u; col < this->dimension(); ++col) {
			temp.at(count) += this->at(col) * mat(col, count);
		}
	}

	return Vector<T>(temp);
};

template<typename T>
inline Vector<T> Vector<T>::operator^(const Vector<T>& vec) const {
	return this->vec_cross(vec);
}

template<typename T>
Vector<T> Vector<T>::rotateByRodrigFormula(double phi, Vector<T> axis) const {
	axis.normalization();
	Vector<T> output{ *this };
	output = output * cos(phi) + (axis ^ output) * sin(phi) + axis * (axis * output) * (1 - cos(phi));
	return output;
}

template<typename T>
template<typename S>
Vector<T> Vector<T>::operator/(const S& arg) {
	Vector<T> output{ *this };
	for (auto& el : output._data)
		el = el / arg;
	return output;
};


template<typename S>
inline std::ostream& operator<<(std::ostream& out, const Vector<S>& vec) {
	for (uint64_t count = 0u; count < vec.dimension(); ++count) {
		if (count == vec.dimension() - 1u)
			out << vec.at(count);
		else
			out << vec.at(count) << ' ';

	}

	return out;
};

template<typename S>
inline std::ofstream& operator<<(std::ofstream& out, const Vector<S>& vec) {
	for (uint64_t count = 0u; count < vec.dimension(); ++count) {
		if (count == vec.dimension() - 1u)
			out << vec.at(count);
		else
			out << vec.at(count) << ' ';

	}

	return out;
}

template<typename S>
inline bool operator>(const Vector<S>& lval, const Vector<S>& rval) {
	return lval.length() > rval.length();
};

template<typename S>
inline bool operator<(const Vector<S>& lval, const Vector<S>& rval) {
	return lval.length() < rval.length();
};

template<typename S>
inline bool operator==(const Vector<S>& lval, const Vector<S>& rval) {
	if (rval.dimension() != lval.dimension())
		return false;

	for (uint64_t count = 0u; count < lval.dimension(); ++count)
		if (lval.at(count) != rval.at(count))
			return false;

	return true;
}

template<typename S>
inline Vector<S> operator*(double s, const Vector<S>& vec) {
	return vec * s;
}
