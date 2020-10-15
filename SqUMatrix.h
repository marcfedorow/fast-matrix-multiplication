#ifndef __SqUMatrix__
#define __SqUMatrix__
#include "IMatrix.h"
template <typename __matrix_elem_type>
class SqUMatrix : public IMatrix<__matrix_elem_type> {
	static size_t extend_size(size_t, size_t);
	void assign(const IMatrix<__matrix_elem_type>& r);
public:
	__matrix_elem_type* matrix;
protected:
	virtual void clear() override;
	virtual void re_alloc(size_t r, size_t c) override;
public:
	static SqUMatrix<__matrix_elem_type> strassen(const SqUMatrix<__matrix_elem_type>& l, const SqUMatrix<__matrix_elem_type>& r, int&, bool w = false);
	static SqUMatrix<__matrix_elem_type> mul(const SqUMatrix<__matrix_elem_type>& l, const SqUMatrix<__matrix_elem_type>& r, int& th, bool t = true);
	static SqUMatrix<__matrix_elem_type> GPUmul(const SqUMatrix<__matrix_elem_type>& l, const SqUMatrix<__matrix_elem_type>& r, int& th, bool t = false);
	explicit SqUMatrix(size_t s);
	SqUMatrix(SqUMatrix<__matrix_elem_type>&& r) noexcept;
	SqUMatrix(const SqUMatrix<__matrix_elem_type>& r);
	SqUMatrix(const IMatrix<__matrix_elem_type>& r);
	virtual ~SqUMatrix() override;
	virtual IMatrix<__matrix_elem_type>& operator=(SqUMatrix<__matrix_elem_type>&& r) noexcept;
	virtual IMatrix<__matrix_elem_type>& operator=(const SqUMatrix<__matrix_elem_type>& r);
	virtual IMatrix<__matrix_elem_type>& operator=(const IMatrix<__matrix_elem_type>& r);
	virtual void set(size_t r, size_t c, __matrix_elem_type value, void(*)(void) = [](void) {throw std::out_of_range("set failure: out of range"); }) override;
	virtual __matrix_elem_type get(size_t r, size_t c, __matrix_elem_type(*)(void) = [](void)->__matrix_elem_type {throw std::out_of_range("get failure: out of range"); }) const override;
	virtual SqUMatrix<__matrix_elem_type> operator*(const IMatrix<__matrix_elem_type>& r) const;
};
#include "SqUMatrix.ipp"
#endif // !__SqUMatrix__