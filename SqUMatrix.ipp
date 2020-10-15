#include "SqUMatrix.h"
#ifdef __SqUMatrix__

template <typename __matrix_elem_type>
inline size_t SqUMatrix< __matrix_elem_type>::extend_size(size_t r, size_t c) {
	size_t s = (r > c) ? r : c;
	if (!(s & (s - 1))) return s;
	size_t l = 1;
	while (s >>= 1) ++l;
	return size_t(1) << l;
}

template <typename __matrix_elem_type>
void SqUMatrix<__matrix_elem_type>::assign(const IMatrix<__matrix_elem_type>& r) {
	size_t s = extend_size(r.getrow(), r.getcol());
	matrix = new __matrix_elem_type[s * s]();
	for (size_t i = 0; i < this->row; ++i) {
		for (size_t j = 0; j < this->col; ++j) {
			set(i, j, r.get(i, j));
		}
	}
}

template <typename __matrix_elem_type>
void SqUMatrix<__matrix_elem_type>::clear() {
	delete[] matrix;
	matrix = nullptr;
}

template <typename __matrix_elem_type>
void SqUMatrix<__matrix_elem_type>::re_alloc(size_t r, size_t c) {
	clear();
	this->row = r, this->col = c;
	size_t s = extend_size(r, c);
	matrix = new __matrix_elem_type[s * s]{};
}

template <typename __matrix_elem_type>
SqUMatrix<__matrix_elem_type>::SqUMatrix(size_t size) : IMatrix<__matrix_elem_type>(size, size) {
	size_t s = extend_size(size, size);
	matrix = new __matrix_elem_type[s * s]();
}

template <typename __matrix_elem_type>
SqUMatrix<__matrix_elem_type>::SqUMatrix(SqUMatrix<__matrix_elem_type>&& r) noexcept : IMatrix<__matrix_elem_type>(r) {
	matrix = r.matrix;
	r.matrix = nullptr;
}

template <typename __matrix_elem_type>
SqUMatrix<__matrix_elem_type>::SqUMatrix(const SqUMatrix<__matrix_elem_type>& r) : IMatrix<__matrix_elem_type>(r) {
	assign(r);
}

template <typename __matrix_elem_type>
SqUMatrix<__matrix_elem_type>::SqUMatrix(const IMatrix<__matrix_elem_type>& r) : IMatrix<__matrix_elem_type>(r) {
	assign(r);
}

template <typename __matrix_elem_type>
SqUMatrix<__matrix_elem_type>::~SqUMatrix() {
	clear();
}

template <typename __matrix_elem_type>
IMatrix<__matrix_elem_type>& SqUMatrix<__matrix_elem_type>::operator=(SqUMatrix<__matrix_elem_type>&& r) noexcept {
	if (this == &r) return *this;
	matrix = r.matrix;
	r.matrix = nullptr;
	return *this;
}

template <typename __matrix_elem_type>
IMatrix<__matrix_elem_type>& SqUMatrix<__matrix_elem_type>::operator=(const IMatrix<__matrix_elem_type>& r) {
	if (this == &r) return *this;
	clear();
	assign(r);
	return *this;
}

template <typename __matrix_elem_type>
IMatrix<__matrix_elem_type>& SqUMatrix<__matrix_elem_type>::operator=(const SqUMatrix<__matrix_elem_type>& r) {
	if (this == &r) return *this;
	clear();
	assign(r);
	return *this;
}

template <typename __matrix_elem_type>
void SqUMatrix<__matrix_elem_type>::set(size_t r, size_t c, __matrix_elem_type value, void(*f)(void)) {
	if (r > this->row || c > this->col) f();
	matrix[r * extend_size(this->col, this->col) + c] = value;
}

template <typename __matrix_elem_type>
__matrix_elem_type SqUMatrix<__matrix_elem_type>::get(size_t r, size_t c, __matrix_elem_type(*def)(void)) const {
	if (r > this->row || c > this->col) return def();
	return matrix[r * extend_size(this->col, this->col) + c];
}

template <typename __matrix_elem_type>
SqUMatrix<__matrix_elem_type> SqUMatrix<__matrix_elem_type>::operator*(const IMatrix<__matrix_elem_type>& r) const {
	if ((*this).col != r.getrow()) throw std::invalid_argument("operator* wrong sizes\n");
	SqUMatrix res(*this);
	res *= r; //here should be some smart thing that decides whitch of { strassen, mul, GPUmul } we should call
	return res;
}

#include "Strassen.h"

template<typename __matrix_elem_type>
inline SqUMatrix<__matrix_elem_type> SqUMatrix<__matrix_elem_type>::strassen(const SqUMatrix<__matrix_elem_type>& l, const SqUMatrix<__matrix_elem_type>& r, int& th, bool w) {
	auto f = w ? ::winograd<__matrix_elem_type> : ::strassen<__matrix_elem_type>;
	auto res = SqUMatrix<__matrix_elem_type>(0);
	delete[] res.matrix;
	res.row = l.getrow();
	res.col = r.getcol();
	size_t size = extend_size(res.row > res.col ? res.row : res.col, l.getcol());
	auto** lm = new __matrix_elem_type * [size]();
	for (size_t i = 0; i < size; ++i) { lm[i] = l.matrix + i * size; }
	auto** rm = new __matrix_elem_type * [size]();
	for (size_t i = 0; i < size; ++i) { rm[i] = r.matrix + i * size; }
	auto m = f(lm, rm, th, size);
	res.matrix = *m;
	delete[] lm;
	delete[] rm;
	delete[] m;
	return res;
}

template<typename __matrix_elem_type>
inline SqUMatrix<__matrix_elem_type> SqUMatrix<__matrix_elem_type>::mul(const SqUMatrix<__matrix_elem_type>& l, const SqUMatrix<__matrix_elem_type>& r, int& th, bool w) {
	if (th <= 0) th = 1;
	const int thr = th;
	const size_t size_tr = extend_size(r.col, r.row);

	auto** rt = new __matrix_elem_type * [size_tr]();
	rt[0] = new __matrix_elem_type[size_tr * size_tr]();
	for (size_t i = 1; i < size_tr; ++i) { rt[i] = rt[0] + i * size_tr; }

	auto* f = new std::future<void>[thr];
	const size_t r_row_dsc = r.row / thr;
	auto fnc = [&](size_t r_row_dsc, int t, const int thr, __matrix_elem_type** rt, const SqUMatrix<__matrix_elem_type>& r) {
		for (size_t i = t * r_row_dsc; i < (((t + 1) == thr) ? r.row : ((t + 1) * r_row_dsc)); ++i) {
			auto* ri = &r.matrix[i * size_tr];
			for (size_t j = 0; j < r.col; ++j) {
				rt[j][i] = ri[j];
			}
		}
	};

	for (int t = 0; t < thr; ++t) {
		f[t] = std::async(std::launch::async, fnc, r_row_dsc, t, thr, rt, std::ref(r));
	}
	for (int t = 0; t < thr; ++t) {
		f[t].get();
	}
	delete[] f;

	const size_t size_res = extend_size(l.col, r.row);
	auto** res = new __matrix_elem_type * [size_res]();
	res[0] = new __matrix_elem_type[size_res * size_res]();
	for (size_t i = 1; i < size_res; ++i) { res[i] = res[0] + i * size_res; }
	const size_t size_l = extend_size(l.col, l.row);

	f = new std::future<void>[thr];
	auto fnc2 = [](size_t r_row_dsc, int t, const int thr, __matrix_elem_type** rt, const SqUMatrix<__matrix_elem_type>& l, const SqUMatrix<__matrix_elem_type>& r, __matrix_elem_type** res, size_t size_l) {
		for (size_t i = t * r_row_dsc; i < (((t + 1) == thr) ? r.row : ((t + 1) * r_row_dsc)); ++i) {
			for (size_t j = 0; j < r.col; ++j) {
				register auto& cij = res[i][j];
				register auto* const ai = l.matrix + i * size_l;
				register auto* const btj = rt[j];
				for (size_t k = 0; k < l.col; ++k) {
					cij += ai[k] * btj[k];
				}
			}
		}
	};

	for (int t = 0; t < thr; ++t) {
		f[t] = std::async(std::launch::async, fnc2, r_row_dsc, t, thr, rt, std::ref(l), std::ref(r), res, size_l);
	}
	for (int t = 0; t < thr; ++t) {
		f[t].get();
	}
	delete[] f;
	delete[] rt[0];
	delete[] rt;
	SqUMatrix<__matrix_elem_type> result(0);
	result.clear();
	result.row = l.row;
	result.col = r.col;
	result.matrix = *res;
	delete[] res;
	return result;
}

#include <amp.h>
#include <amp_math.h>

template <typename __matrix_elem_type>
void MatrixMultiply(__matrix_elem_type* A, __matrix_elem_type* B, __matrix_elem_type* C, unsigned long size)
{
	concurrency::array_view <const __matrix_elem_type, 2> avA(size, size, A);
	concurrency::array_view <const __matrix_elem_type, 2> avB(size, size, B);
	concurrency::array_view <__matrix_elem_type, 2> avC(size, size, C);
	avC.discard_data();

	concurrency::parallel_for_each(avC.extent, [=](concurrency::index <2> idx) restrict(amp)
		{
			__matrix_elem_type sum = 0;

			for (unsigned long k = 0; k < size; ++k)
			{
				sum += avA(idx[0], k) * avB(k, idx[1]);
			}

			avC[idx] = sum;
		});
}

template<typename __matrix_elem_type>
inline SqUMatrix<__matrix_elem_type> SqUMatrix<__matrix_elem_type>::GPUmul(const SqUMatrix<__matrix_elem_type>& l, const SqUMatrix<__matrix_elem_type>& r, int& th, bool T)
{
	//it supposed that if (T) we'll transpond matrix & pass it to
	//sum += avA(idx[0], k) * avB(idx[1], k);
	//function, but it obviously takes too much time with no gain
	SqUMatrix<__matrix_elem_type> result(extend_size(l.row, r.col));
	MatrixMultiply(l.matrix, r.matrix, result.matrix, extend_size(l.row, r.col));
	return result;
}

#endif //__SqUMatrix__