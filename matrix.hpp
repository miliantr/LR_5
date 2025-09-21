#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include "vector.hpp"

template<typename T>
class Vector;

template<typename T>
inline size_t count_inverse(const std::vector<T>& vec) {
	size_t count{};
	for (uint64_t i = 0u; i < vec.size(); ++i)
		for (uint64_t j = i + 1u; j < vec.size(); ++j)
			if (vec.at(i) > vec.at(j))
				++count;
	return count;
}

template<typename T>
inline void generate_permutation(std::vector<T>& vec, const size_t size) {
	vec.resize(size);
	std::iota(vec.begin(), vec.end(), 0);
}

template<typename T>
class Matrix {
private:
	uint64_t _rows;
	uint64_t _cols;
	std::vector<T> _data;
public:
	Matrix() noexcept;
	Matrix(uint64_t rows, uint64_t cols) noexcept;
	Matrix(uint64_t rows, uint64_t cols, const Vector<T>& vec) noexcept;
	Matrix(uint64_t rows, uint64_t cols, const std::vector<T>& vec) noexcept;
	Matrix(const std::vector<T>& vec);
	Matrix(const Matrix<T>& mat);
	Matrix<T>& operator=(const Matrix<T>& mat);
	~Matrix();

	T& at(int64_t row, int64_t col);
	T at(int64_t row, int64_t col) const;
	T& operator()(int64_t row, int64_t col);
	T operator()(int64_t row, int64_t col) const;

	Matrix<T> operator-() const;
	Matrix<T> operator+(const Matrix<T>& mat) const;
	Matrix<T> operator*(const Matrix<T>& mat) const;
	Vector<T> operator*(const Vector<T>& vec) const;
	Matrix<T> operator!() const;
	template<typename S> Matrix<T> operator*(const S& s) const;
	template<typename S> friend std::ostream& operator<<(std::ostream& out, const Matrix<S>& mat);
	template<typename S, typename U> friend Matrix<U> operator*(const S& val, const Matrix<U>& mat);

	void resize(uint64_t rows, uint64_t cols);
	void set_rows(int64_t row, const Vector<T>& vec);
	uint64_t cols() const;
	uint64_t rows() const;
	double determinate() const;
	Vector<T> get_rows(uint64_t row) const;
	Vector<T> get_cols(uint64_t col) const;
	Matrix<T>& transpose() noexcept;
	Matrix<T>& multiply(const Matrix<T>& mat) noexcept;
	Matrix<T>& swap_rows(int64_t i, uint64_t j);
	Matrix<T>& push_row(const Vector<T>& vec);
	Matrix<T>& push_col(const Vector<T>& vec);
	static Matrix<double> E(uint64_t size);
};

template<typename T>
Matrix<T>::Matrix() noexcept : _rows(0), _cols(0) {};

template<typename T>
Matrix<T>::Matrix(uint64_t rows, uint64_t cols) noexcept : _rows(rows), _cols(cols) {
	_data.resize(_rows * _cols);
}

template<typename T>
Matrix<T>::Matrix(uint64_t rows, uint64_t cols, const Vector<T>& vec) noexcept : Matrix(rows, cols) {
	for (uint64_t row = 0u; row < _rows; ++row)
		for (uint64_t col = 0u; col < _cols; ++col) {
			this->at(row, col) = vec.at(row * _cols + col);
		}
}

template<typename T>
Matrix<T>::Matrix(uint64_t rows, uint64_t cols, const std::vector<T>& vec) noexcept : _rows(rows), _cols(cols) {
	_data = vec;
}

template<typename T>
Matrix<T>::Matrix(const std::vector<T>& vec) :_rows(1), _cols(vec.size()) {
	_data = vec;
}

template<typename T>
Matrix<T>::Matrix(const Matrix<T>& mat) :_rows(mat._rows), _cols(mat._cols) {
	_data = mat._data;
}

template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& mat) {
	if (this != &mat) {
		this->_data = mat._data;
		this->_rows = mat._rows;
		this->_cols = mat._cols;
	}
	return *this;
}

template<typename T>
Matrix<T>::~Matrix() {
	_data.clear();
}

template<typename T>
uint64_t Matrix<T>::cols() const {
	return _cols;
}

template<typename T>
uint64_t Matrix<T>::rows() const {
	return _rows;
}

template<typename T>
void Matrix<T>::resize(uint64_t rows, uint64_t cols) {
	_rows = rows;
	_cols = cols;
	_data.resize(rows * cols);
}

template<typename S>
inline std::ostream& operator<<(std::ostream& out, const Matrix<S>& mat) {
	for (uint64_t row = 0u; row < mat._rows; ++row) {
		for (uint64_t col = 0u; col < mat._cols; ++col) {
			out << mat._data.at(row * mat._cols + col) << " ";
		}
		out << std::endl;
	}

	return out;
}

template<typename T>
Matrix<T>& Matrix<T>::transpose() noexcept {
	Matrix<T> temp(_cols, _rows);

	for (uint64_t row = 0u; row < temp._rows; ++row)
		for (uint64_t col = 0u; col < temp._cols; ++col)
			temp(row, col) = this->at(col, row);

	_cols = temp._cols;
	_rows = temp._rows;

	*this = temp;

	return *this;
}

template<typename T>
T& Matrix<T>::at(int64_t row, int64_t col) {
	if (row >= _rows || col >= _cols)
		throw std::out_of_range("at");

	return _data.at(row * _cols + col);
}

template<typename T>
T Matrix<T>::at(int64_t row, int64_t col) const {
	if (row >= _rows || col >= _cols)
		throw std::out_of_range("at");

	return _data.at(row * _cols + col);
}

template<typename T>
T& Matrix<T>::operator()(int64_t row, int64_t col) {
	return this->at(row, col);
}

template<typename T>
T Matrix<T>::operator()(int64_t row, int64_t col) const {
	return this->at(row, col);
}

template<typename T>
Matrix<T>& Matrix<T>::multiply(const Matrix<T>& mat) noexcept {
	Matrix<T> output{ _rows, mat._cols, };

	for (uint64_t row = 0u; row < _rows; ++row) {
		for (uint64_t col = 0u; col < mat._cols; ++col) {

			for (uint64_t count = 0u; count < _cols; ++count) {
				output.at(row, col) += this->at(row, count) * mat.at(count, col);
			}
		}
	}

	*this = output;

	return *this;
}

template<typename T>
Matrix<T> Matrix<T>::operator-() const {
	Matrix<T> temp{ *this };
	for (auto& el : temp._data)
		el = -el;
	return temp;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& mat) const {
	if (this->_cols != mat._cols && this->_rows != mat._rows)
		throw std::logic_error("Dimensional matrixs must be also");

	Matrix<T> temp{ *this };
	for (uint64_t row = 0u; row < mat._rows; ++row)
		for (uint64_t col = 0u; col < mat._cols; ++col)
			temp.at(row, col) += mat.at(row, col);

	return temp;
}

template<typename T>
template<typename S>
Matrix<T> Matrix<T>::operator*(const S& s) const {
	Matrix<T> temp{ *this };
	for (uint64_t row = 0u; row < _rows; ++row)
		for (uint64_t col = 0u; col < _cols; ++col)
			temp.at(row, col) *= s;

	return temp;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& mat) const {
	Matrix<T> temp{ *this };
	return temp.multiply(mat);
}

template<typename T>
Vector<T> Matrix<T>::operator*(const Vector<T>& vec) const {
	Matrix<T> temp{ *this };
	temp.transpose();
	return vec * temp;
}

template<typename S, typename U>
inline Matrix<U> operator*(const S& val, const Matrix<U>& mat) {
	return mat * val;
}

template<typename T>
double Matrix<T>::determinate() const {
	if (_rows != _cols)
		throw std::logic_error("determinate");

	double answer{};

	std::vector<int> perm;
	generate_permutation(perm, _rows);

	do {
		double temp{ 1 };

		for (uint64_t i = 0u; i < _rows; ++i)
			temp *= at(i, perm.at(i));
		answer += pow(-1, count_inverse(perm)) * temp;

	} while (std::next_permutation(perm.begin(), perm.end()));

	return answer;
}

template<typename T>
Matrix<double> Matrix<T>::E(uint64_t size) {
	Matrix<double> E(size, size);
	for (uint64_t count = 0u; count < size; ++count)
		E.at(count, count) = 1;
	return E;
}

template<typename T>
Vector<T> Matrix<T>::get_rows(uint64_t row) const {
	if (row >= _rows)
		throw std::logic_error("get_rows");

	Vector<T> output;
	for (uint64_t col = 0u; col < _cols; ++col)
		output.push_back(this->at(row, col));
	return output;
}

template<typename T>
Vector<T> Matrix<T>::get_cols(uint64_t col) const {
	if (col >= _cols)
		throw std::logic_error("get_cols");
	Vector<T> output;
	for (uint64_t row = 0u; row < _rows; ++row)
		output.push_back(this->at(row, col));
	return output;
}

template<typename T>
void Matrix<T>::set_rows(int64_t row, const Vector<T>& vec) {
	if (row >= _rows && vec.dimension() != _cols)
		throw std::logic_error("set_rows");

	for (int64_t col = 0u; col < _cols; ++col)
		this->at(row, col) = vec.at(col);
}

template<typename T>
Matrix<T>& Matrix<T>::swap_rows(int64_t i, uint64_t j) {
	Vector<T> temp = this->get_rows(i);
	this->set_rows(i, this->get_rows(j));
	this->set_rows(j, temp);

	return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::push_row(const Vector<T>& vec) {
	if (_cols == 0 && _rows == 0) {
		_cols = vec.dimension();
		_rows = 1;
		_data.resize(_cols * _rows);

		for (uint64_t count = 0u; count < _cols; ++count)
			_data.at(count) = vec.at(count);

		return *this;
	}

	if (vec.dimension() != _cols)
		throw std::logic_error("push row");

	for (uint64_t count = 0u; count < _cols; ++count)
		_data.push_back(vec.at(count));

	++_rows;

	return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::push_col(const Vector<T>& vec) {
	if (vec.dimension() != _rows)
		throw std::logic_error("push col");

	(*this).transpose();

	push_row(vec);

	return (*this).transpose();

}

template<typename T>
Matrix<T> Matrix<T>::operator!() const {
	if (_cols != _rows)
		throw std::logic_error("Inverse matrix");

	Matrix<T> output(*this);

	for (uint64_t count = 0u; count < _rows; ++count) {
		std::vector<T> col(_rows, 0);
		col.at(count) = 1;
		output.push_col(col);
	}

	for (uint64_t count = 0u; count < _rows; ++count) {

		uint64_t row_with_biggest_element = count;

		for (uint64_t row = count; row < _rows; ++row)
			if (output.at(row, count) > abs(output.at(row_with_biggest_element, count)))
				row_with_biggest_element = row;

		if (!output.at(row_with_biggest_element, count))
			throw std::logic_error("Vyroshdena");

		output.swap_rows(count, row_with_biggest_element);

		Vector<T> temp = output.get_rows(count);
		temp = temp / temp(count);
		output.set_rows(count, temp);

		for (uint64_t row = count + 1u; row < _rows; ++row) {
			T coeff = output.at(row, count);
			Vector<T> temp_temp = temp * (-coeff);
			Vector<T> other_row = output.get_rows(row);
			other_row = other_row + temp_temp;
			output.set_rows(row, other_row);
		}
	}

	for (int64_t count = _rows - 1; count > 0; --count) {

		Vector<T> temp = output.get_rows(count);
		temp = temp / temp(count);
		output.set_rows(count, temp);

		for (int64_t row = count - 1; row >= 0; --row) {
			T coeff = output.at(row, count);
			Vector<T> temp_temp = temp * (-coeff);
			Vector<T> other_row = output.get_rows(row);
			other_row = other_row + temp_temp;
			output.set_rows(row, other_row);
		}
	}

	Matrix<T> inv(_rows, _rows);

	for (uint64_t row = 0u; row < _rows; ++row)
		for (uint64_t col = 0u; col < _cols; ++col) {
			inv.at(row, col) = output.at(row, col + _rows);
		}

	return inv;
}