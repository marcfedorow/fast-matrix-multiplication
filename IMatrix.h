#ifndef __IMatrix__
#define __IMatrix__
#include <stdexcept>
template <typename __matrix_elem_type>
class IMatrix {
protected:
	size_t row, col;
	virtual void clear() = 0;
	virtual void re_alloc(size_t r, size_t c) = 0;
public:
	IMatrix(size_t r = 1000, size_t c = 1000) : row(r), col(c) {};
	IMatrix(const IMatrix<__matrix_elem_type>& r) : row(r.row), col(r.col) {};
	virtual IMatrix<__matrix_elem_type>& operator=(const IMatrix<__matrix_elem_type>& r) = 0;
	virtual ~IMatrix() {};

	virtual void set(size_t r, size_t c, __matrix_elem_type value, void(*)(void) = [](void) {throw std::out_of_range("set failure: out of range"); }) = 0;
	virtual __matrix_elem_type get(size_t r, size_t c, __matrix_elem_type(*)(void) = [](void)->__matrix_elem_type {throw std::out_of_range("get failure: out of range");}) const = 0;

	virtual IMatrix<__matrix_elem_type>& operator+=(const IMatrix<__matrix_elem_type>& r);
	virtual IMatrix<__matrix_elem_type>& operator-=(const IMatrix<__matrix_elem_type>& r);
	virtual IMatrix<__matrix_elem_type>& operator*=(const IMatrix<__matrix_elem_type>& r);
	virtual bool operator==(const IMatrix<__matrix_elem_type>& r) const;

	virtual void fill_n(size_t n = 1000, __matrix_elem_type(*)(void) = [](void)->__matrix_elem_type {return (rand()) - RAND_MAX / 2; }) final;

	//for operators = and <<
	virtual size_t getrow() const final { return row; }
	virtual size_t getcol() const final { return col; }
};

#include "IMatrix.ipp"
#endif // !__IMatrix__
