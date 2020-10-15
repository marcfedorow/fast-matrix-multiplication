#ifdef __IMatrix__

template <typename __matrix_elem_type>
IMatrix<__matrix_elem_type>& IMatrix<__matrix_elem_type>::operator+=(const IMatrix<__matrix_elem_type>& r) {
	if (row != r.row || col != r.col) throw std::invalid_argument("Matrix addition error: sizes do not match");
	for (size_t i = 0; i < row; ++i) {
		for (size_t j = 0; j < col; ++j) {
			set(i, j, get(i, j) + r.get(i, j));
		}
	}
	return *this;
}

template <typename __matrix_elem_type>
IMatrix<__matrix_elem_type>& IMatrix<__matrix_elem_type>::operator-=(const IMatrix<__matrix_elem_type>& r) {
	if (row != r.row || col != r.col) throw std::invalid_argument("Matrix addition error: sizes do not match");
	for (size_t i = 0; i < row; ++i) {
		for (size_t j = 0; j < col; ++j) {
			set(i, j, get(i, j) - r.get(i, j));
		}
	}
	return *this;
}

template <typename __matrix_elem_type>
IMatrix<__matrix_elem_type>& IMatrix<__matrix_elem_type>::operator*=(const IMatrix<__matrix_elem_type>& r) {
	if (col != r.row) throw std::invalid_argument("Matrix multiplication error: sizes do not match");
	size_t size = row * r.col;
	if (size == 0 || size / row != r.col || size / r.col != row) throw std::bad_array_new_length();
	__matrix_elem_type* m = new __matrix_elem_type[size]{};
	for (size_t i = 0; i < row; ++i) {
		for (size_t j = 0; j < r.col; ++j) {
			for (size_t k = 0; k < col; ++k) {
				m[i * col + j] += get(i, k) * r.get(k, j);
			}
		}
	}
	//clear();
	re_alloc(row, r.col);
	for (size_t i = 0; i < row; ++i) {
		for (size_t j = 0; j < col; ++j) {
			set(i, j, m[i * col + j]);
		}
	}
	delete[] m;
	return *this;
}

template <typename __matrix_elem_type>
bool IMatrix<__matrix_elem_type>::operator==(const IMatrix<__matrix_elem_type>& r) const {
	if (this == &r) return true;
	if (col != r.col || row != r.row) return false;
	for (size_t i = 0; i < row; ++i) {
		for (size_t j = 0; j < r.col; ++j) {
			if (get(i, j) != r.get(i, j)) return false;
		}
	}
	return true;
}

#include <set>
#include <utility>
template <typename __matrix_elem_type>
void IMatrix<__matrix_elem_type>::fill_n(size_t n, __matrix_elem_type(*f)(void)) {
	if (n > row * col) {
		n = row * col;
	}
	std::set<std::pair<int, int>> s;
	while (s.size() < n) {
		size_t r = rand() % row, c = rand() % col;
		s.insert(std::make_pair(r, c));
	}
	for (auto &p : s) set(p.first, p.second, f());
}

#include <ostream>
template <typename __matrix_elem_type>
std::ostream& operator<<(std::ostream& o, const IMatrix<__matrix_elem_type>& m) {
	size_t r = m.getrow(), c = m.getcol();
	for (size_t i = 0; i < r; ++i) {
		std::cout << "[";
		for (size_t j = 0; j < c; ++j) {
			std::cout << m.get(i, j) << ' ';
		}
		std::cout << "\b]\n";
	}
	return o;
}
#endif // __IMatrix__