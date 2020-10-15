#ifndef __UMatrix__
#define __UMatrix__
#include "IMatrix.h"
template <typename __matrix_elem_type>
class UMatrix : public IMatrix<__matrix_elem_type> {
	void assign(const IMatrix<__matrix_elem_type>& r);
	__matrix_elem_type* matrix;
protected:
	virtual void clear() override;
	virtual void re_alloc(size_t r, size_t c) override;
public:
	UMatrix(size_t r, size_t c);
	UMatrix(UMatrix<__matrix_elem_type>&& r) noexcept;
	UMatrix(const UMatrix<__matrix_elem_type>& r);
	UMatrix(const IMatrix<__matrix_elem_type>& r);
	virtual ~UMatrix() override;
	virtual IMatrix<__matrix_elem_type>& operator=(UMatrix<__matrix_elem_type>&& r) noexcept;
	virtual IMatrix<__matrix_elem_type>& operator=(const UMatrix<__matrix_elem_type>& r);
	virtual IMatrix<__matrix_elem_type>& operator=(const IMatrix<__matrix_elem_type>& r);

	virtual void set(size_t r, size_t c, __matrix_elem_type value, void(*)(void) = [](void) {throw std::out_of_range("set failure: out of range"); }) override;
	virtual __matrix_elem_type get(size_t r, size_t c, __matrix_elem_type(*)(void) = [](void)->__matrix_elem_type {throw std::out_of_range("get failure: out of range"); }) const override;

	virtual UMatrix<__matrix_elem_type> operator+(const IMatrix<__matrix_elem_type>& r) const;
	virtual UMatrix<__matrix_elem_type> operator*(const IMatrix<__matrix_elem_type>& r) const;
};
#include "UMatrix.ipp"
#endif // !__UMatrix__