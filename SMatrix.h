#ifndef __SMatrix__
#define __SMatrix__
#include "IMatrix.h"
template <typename __matrix_elem_type>
class SMatrix : public IMatrix<__matrix_elem_type> {
	void assign(const IMatrix<__matrix_elem_type>& r);
	struct node {
		size_t r, c;
		__matrix_elem_type value;
		node* next;
	};
	node* head;
protected:
	void clear() override;
	void re_alloc(size_t r, size_t c) override;
public:
	SMatrix(size_t r, size_t c);
	SMatrix(SMatrix<__matrix_elem_type>&& r) noexcept;
	SMatrix(const SMatrix<__matrix_elem_type>& r);
	SMatrix(const IMatrix<__matrix_elem_type>& r);
	virtual ~SMatrix() override;
	virtual IMatrix<__matrix_elem_type>& operator=(SMatrix<__matrix_elem_type>&& r) noexcept;
	virtual IMatrix<__matrix_elem_type>& operator=(const SMatrix<__matrix_elem_type>& r);
	virtual IMatrix<__matrix_elem_type>& operator=(const IMatrix<__matrix_elem_type>& r);

	virtual void set(size_t r, size_t c, __matrix_elem_type value, void(*)(void) = [](void) {throw std::out_of_range("set failure: out of range"); }) override;
	virtual __matrix_elem_type get(size_t r, size_t c, __matrix_elem_type(*)(void) = [](void)->__matrix_elem_type {throw std::out_of_range("get failure: out of range"); }) const override;

	virtual SMatrix<__matrix_elem_type> operator+(const IMatrix<__matrix_elem_type>& r) const;
	virtual SMatrix<__matrix_elem_type> operator*(const IMatrix<__matrix_elem_type>& r) const;
};
#include "SMatrix.ipp"
#endif // !__SMatrix__