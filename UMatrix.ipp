#ifdef __UMatrix__

template <typename __matrix_elem_type>
void UMatrix<__matrix_elem_type>::assign(const IMatrix<__matrix_elem_type>& r) {
	matrix = new __matrix_elem_type[this->row * this->col]{};
	for (size_t i = 0; i < this->row; ++i) {
		for (size_t j = 0; j < this->col; ++j) {
			set(i, j, r.get(i, j));
		}
	}
}

template <typename __matrix_elem_type>
void UMatrix<__matrix_elem_type>::clear() {
	delete[] matrix;
	matrix = nullptr;
}

template <typename __matrix_elem_type>
void UMatrix<__matrix_elem_type>::re_alloc(size_t r, size_t c) {
	clear();
	this->row = r, this->col = c;
	matrix = new __matrix_elem_type[r * c]{};
}

template <typename __matrix_elem_type>
UMatrix<__matrix_elem_type>::UMatrix(size_t r, size_t c) {
	this->row = r, this->col = c;
	matrix = new __matrix_elem_type[this->row * this->col]{};
}

template <typename __matrix_elem_type>
UMatrix<__matrix_elem_type>::UMatrix(UMatrix<__matrix_elem_type>&& r) noexcept : IMatrix<__matrix_elem_type>(r) {
	matrix = r.matrix;
	r.matrix = nullptr;
}

template <typename __matrix_elem_type>
UMatrix<__matrix_elem_type>::UMatrix(const UMatrix<__matrix_elem_type>& r) : IMatrix<__matrix_elem_type>(r) {
	assign(r);
}

template <typename __matrix_elem_type>
UMatrix<__matrix_elem_type>::UMatrix(const IMatrix<__matrix_elem_type>& r) : IMatrix<__matrix_elem_type>(r) {
	assign(r);
}

template <typename __matrix_elem_type>
UMatrix<__matrix_elem_type>::~UMatrix() {
	clear();
}

template <typename __matrix_elem_type>
IMatrix<__matrix_elem_type>& UMatrix<__matrix_elem_type>::operator=(UMatrix<__matrix_elem_type>&& r) noexcept {
	if (this == &r) return *this;
	matrix = r.matrix;
	r.matrix = nullptr;
	return *this;
}

template <typename __matrix_elem_type>
IMatrix<__matrix_elem_type>& UMatrix<__matrix_elem_type>::operator=(const IMatrix<__matrix_elem_type>& r) {
	if (this == &r) return *this;
	clear();
	assign(r);
	return *this;
}

template <typename __matrix_elem_type>
IMatrix<__matrix_elem_type>& UMatrix<__matrix_elem_type>::operator=(const UMatrix<__matrix_elem_type>& r) {
	if (this == &r) return *this;
	clear();
	assign(r);
	return *this;
}

template <typename __matrix_elem_type>
void UMatrix<__matrix_elem_type>::set(size_t r, size_t c, __matrix_elem_type value, void(*f)(void)) {
	if (r > this->row || c > this->col) f();
	matrix[r * this->col + c] = value;
}

template <typename __matrix_elem_type>
__matrix_elem_type UMatrix<__matrix_elem_type>::get(size_t r, size_t c, __matrix_elem_type(*def)(void)) const {
	if (r > this->row || c > this->col) return def();
	return matrix[r * this->col + c];
}

template <typename __matrix_elem_type>
UMatrix<__matrix_elem_type> UMatrix<__matrix_elem_type>::operator+(const IMatrix<__matrix_elem_type>& r) const {
	UMatrix res(*this);
	res += r;
	return res;
}

template <typename __matrix_elem_type>
UMatrix<__matrix_elem_type> UMatrix<__matrix_elem_type>::operator*(const IMatrix<__matrix_elem_type>& r) const {
	UMatrix res(*this);
	res *= r;
	return res;
}

#endif // __UMatrix__