#define __Rational__
template <typename T>
class Rational {
	void normalize() {
		T g = gcd(numerator, denominator);
		numerator /= g;
		denominator /= g;
		if (denominator < 0) {
			numerator = -numerator;
			denominator = -denominator;
		}
	}
	static T _gcd(T a, T b) {
		return b ? _gcd(b, a % b) : a;
	}
public:
	static T gcd(T a, T b) { return std::abs(_gcd(a, b)); }
	T numerator;
	T denominator;
	Rational(T n = (T)0, T d = (T)1) : numerator(n), denominator(d) {};
	operator bool() const { return numerator; }
	Rational& operator+=(const Rational r) {
		if (this == &r) {
			if (denominator & 1) numerator *= 2;
			else denominator /= 2;
			return *this;
		}
		this->numerator *= r.denominator;
		this->numerator += r.numerator * this->denominator;
		this->denominator *= r.denominator;
		normalize();
		return *this;
	}
	Rational& operator-=(const Rational r) {
		if (this == &r) {
			numerator = 0;
			denominator = 1;
			return *this;
		}
		this->numerator *= r.denominator;
		this->numerator -= r.numerator * this->denominator;
		this->denominator *= r.denominator;
		normalize();
		return *this;
	}
	Rational& operator*=(const Rational r) {
		if (this == &r) {
			numerator *= numerator;
			denominator *= denominator;
			return *this;
		}
		this->numerator *= r.denominator;
		this->denominator *= r.denominator;
		normalize();
		return *this;
	}
};

#include <ostream>
template <typename T>
std::ostream& operator<<(std::ostream& o, const Rational<T>& r) {
	return std::cout << "{" << r.numerator << "/" << r.denominator << "}";
}

template <typename T>
Rational<T> operator+(const Rational<T>& l, const Rational<T>& r) {
	T n = l.numerator * r.denominator + r.numerator * l.denominator;
	T d = l.denominator * r.denominator;
	T g = Rational<T>::gcd(n, d);
	return Rational<T>(n/g, d/g);
}

template <typename T>
Rational<T> operator-(const Rational<T>& l, const Rational<T>& r) {
	T n = l.numerator * r.denominator - r.numerator * l.denominator;
	T d = l.denominator * r.denominator;
	T g = Rational<T>::gcd(n, d);
	return Rational<T>(n / g, d / g);
}

template <typename T>
Rational<T> operator*(const Rational<T>& l, const Rational<T>& r) {
	T n = l.numerator * r.numerator;
	T d = l.denominator * r.denominator;
	T g = Rational<T>::gcd(n, d);
	return Rational<T>(n / g, d / g);
}