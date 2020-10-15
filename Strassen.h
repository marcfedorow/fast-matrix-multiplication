#ifndef _Strassen_
#define _Strassen_

template <typename __matrix_elem_type>
__matrix_elem_type** cubic_mul_t(__matrix_elem_type** A, __matrix_elem_type** BT, size_t size) {
	__matrix_elem_type** C = new __matrix_elem_type * [size];
	C[0] = new __matrix_elem_type[size * size]();
	for (size_t i = 1; i < size; ++i) { C[i] = C[0] + i * size; }
	for (size_t i = 0; i < size; ++i) {
		for (size_t j = 0; j < size; ++j) {
			register auto cij = __matrix_elem_type();
			register auto* const ai = A[i];
			register auto* const btj = BT[j];
			for (size_t k = 0; k < size; ++k) {
				cij += ai[k] * btj[k];
			}
			C[i][j] = cij;
		}
	}
	return C;
}

#define c2da(n)										\
matrices[n] = new __matrix_elem_type * [size]();	\
matrices[n][0] = new __matrix_elem_type[size_sqr]();\
{													\
	register auto * const m_n = matrices[n];		\
	register auto const zero = m_n[0];				\
	register size_t sh  = size;						\
	for (size_t i = 1; i < size; ++i) {				\
		m_n[i] = zero + sh;							\
		sh += size;									\
	}												\
}

#define a1 l[i]
#define a2 l[i] + size
#define a3 l[i + size]
#define a4 l[i + size] + size
#define b1 r[i]
#define b2 r[i] + size
#define b3 r[i + size]
#define b4 r[i + size] + size

#define create(n, l, sign, r)						\
c2da(n);											\
for (size_t i = 0; i < size; ++i) {					\
	register auto * const m = matrices[n][i];		\
	register auto const * const la = l;				\
	register auto const * const ra = r;				\
	for (size_t j = 0; j < size; ++j) {				\
		m[j] = la[j] sign ra[j];					\
	}												\
}
#define create_t(n, l, sign, r)						\
c2da(n);											\
for (size_t i = 0; i < size; ++i) {					\
	register auto * m	= &matrices[n][0][i];		\
	register auto const * const la = l;				\
	register auto const * const ra = r;				\
	for (size_t j = 0; j < size; ++j) {				\
		*m = la[j] sign ra[j]; m += size;			\
	}												\
}
#define copy_shrink(n, f)							\
c2da(n);											\
for (size_t i = 0; i < size; ++i) {					\
	register auto * const m = matrices[n][i];		\
	register auto const * const ma = f;				\
	for (size_t j = 0; j < size; ++j) {				\
		m[j] = ma[j];								\
	}												\
}
#define copy_shrink_t(n, f)							\
c2da(n);											\
for (size_t i = 0; i < size; ++i) {					\
	register auto * m = &matrices[n][0][i];			\
	register auto const * const ma = f;				\
	for (size_t j = 0; j < size; ++j) {				\
		*m = ma[j]; m += size;						\
	}												\
}


size_t const lim = 1 << 9;
const bool last_step_transpond = true;
//set false and change line 94 when cubic_mul_GPU done well & is more exxicient
#include <future>
template <typename __matrix_elem_type>
__matrix_elem_type** strassen(__matrix_elem_type** l, __matrix_elem_type** r, int& th, size_t size)
noexcept { //if something got wrong, nothing could be done
	if (size < lim) return last_step_transpond
		? cubic_mul_t(l, r, size) //partial GPU multiplication not realized
		: (throw std::exception("noexcept will stop the program immediately"), nullptr); //cubic_mul_GPU(l, r, size);
	const bool transpond = (size == lim) && last_step_transpond;
	size /= 2;
	size_t size_sqr = size * size;

	__matrix_elem_type** matrices[14];

	create(0, a1, +, a4);
	create(1, a3, +, a4);
	//copy_shrink(2, a1);
	matrices[2] = l; //a1
	//copy_shrink(3, a4);
	matrices[3] = new __matrix_elem_type * [size]; //a4
	for (size_t i = 0; i < size; ++i) { matrices[3][i] = l[size + i] + size; }
	create(4, a1, +, a2);
	create(5, a3, -, a1);
	create(6, a2, -, a4);

	if (size * 2 == lim && transpond) {
		create_t(7, b1, +, b4);
		copy_shrink_t(8, b1);
		create_t(9, b2, -, b4);
		create_t(10, b3, -, b1);
		copy_shrink_t(11, b4);
		create_t(12, b1, +, b2);
		create_t(13, b3, +, b4);
	}
	else {
		create(7, b1, +, b4);
		//copy_shrink(8, b1);
		matrices[8] = r; //b1
		create(9, b2, -, b4);
		create(10, b3, -, b1);
		//copy_shrink(11, b4);
		matrices[11] = new __matrix_elem_type * [size]; //b4
		for (size_t i = 0; i < size; ++i) { matrices[11][i] = r[size + i] + size; }
		create(12, b1, +, b2);
		create(13, b3, +, b4);
	}

	__matrix_elem_type** p[7];
	if (th > 0) {
		std::future<__matrix_elem_type**> f[7];
		th -= 7;
		for (size_t i = 0; i < 7; ++i) {
			f[i] = std::async(&strassen<__matrix_elem_type>, matrices[i], matrices[i + 7], std::ref(th), size);
		}
		for (size_t i = 0; i < 7; ++i) {
			p[i] = f[i].get();
		}
		th += 7;
	}
	else {
		for (size_t i = 0; i < 7; ++i) {
			p[i] = strassen(matrices[i], matrices[i + 7], th, size);
		}
	}

	matrices[3][0] = nullptr;
	if (!(size * 2 == lim && transpond))
		matrices[11][0] = nullptr;
	for (size_t i = 0; i < 14; ++i) {
		if (i == 2 || i == 8 && !(size * 2 == lim && transpond)) continue;
		delete[] matrices[i][0];
		delete[] matrices[i];
	}

	size *= 2;
	__matrix_elem_type** res = new __matrix_elem_type * [size];
	res[0] = new __matrix_elem_type[size * size]();
	for (size_t i = 1; i < size; ++i) { res[i] = *res + i * size; }
	size /= 2;
	for (size_t i = 0; i < size; ++i) {
		register auto* const ri = res[i];
		register auto const* const p0i = p[0][i];
		register auto const* const p3i = p[3][i];
		register auto const* const p4i = p[4][i];
		register auto const* const p6i = p[6][i];
		for (size_t j = 0; j < size; ++j) {
			ri[j] = p0i[j] + p3i[j] - p4i[j] + p6i[j];
		}
	}
	for (size_t i = 0; i < size; ++i) {
		register auto* const ri = res[i] + size;
		register auto const* const p2i = p[2][i];
		register auto const* const p4i = p[4][i];
		for (size_t j = 0; j < size; ++j) {
			ri[j] = p2i[j] + p4i[j];
		}
	}
	for (size_t i = 0; i < size; ++i) {
		register auto* const ri = res[i + size];
		register auto const* const p1i = p[1][i];
		register auto const* const p3i = p[3][i];
		for (size_t j = 0; j < size; ++j) {
			ri[j] = p1i[j] + p3i[j];
		}
	}
	for (size_t i = 0; i < size; ++i) {
		register auto* const ri = res[i + size] + size;
		register auto const* const p0i = p[0][i];
		register auto const* const p1i = p[1][i];
		register auto const* const p2i = p[2][i];
		register auto const* const p5i = p[5][i];
		for (size_t j = 0; j < size; ++j) {
			ri[j] = p0i[j] - p1i[j] + p2i[j] + p5i[j];
		}
	}
	for (size_t i = 0; i < 7; ++i) {
		delete[] p[i][0];
		delete[] p[i];
	}
	return res;
}

#define s1 S[0][i]
#define s2 S[1][i]
#define s5 S[4][i]
#define s6 S[5][i]

template <typename __matrix_elem_type>
__matrix_elem_type** winograd(__matrix_elem_type** l, __matrix_elem_type** r, int& th,
	size_t size) noexcept { //if something got wrong, nothing could be done
	const bool transpond = (size == lim) && 1;
	size /= 2;
	size_t size_sqr = size * size;

	__matrix_elem_type** S[8];
#define matrices S
	if (th > 0) {
		th -= 3;
		std::future<void> f[3];
		f[0] = std::async(std::launch::async,
			[&]() {
				create(0, a3, +, a4);
				create(1, s1, -, a1);
				create(3, a2, -, s2);
			});
		f[1] = std::async(std::launch::async,
			[&]() {
				create(4, b2, -, b1);
				create(5, b4, -, s5);
				if (transpond) {
					create_t(7, s6, -, b3);
				} else {
					create(7, s6, -, b3);
				}
			});
		f[2] = std::async(std::launch::async,
			[&]() {
				create(2, a1, -, a3);
				if (transpond) {
					create_t(6, b4, -, b2);
				} else {
					create(6, b4, -, b2);
				}
				//may be splitted but what for?
			});
		for (size_t i = 0; i < 3; ++i) {
			f[i].get();
		}
		th += 3;
	} else {
		create(0, a3, +, a4);
		create(1, s1, -, a1);
		create(2, a1, -, a3);
		create(3, a2, -, s2);

		create(4, b2, -, b1);
		create(5, b4, -, s5);
		if (transpond) {
			create_t(6, b4, -, b2);
			create_t(7, s6, -, b3);
		} else {
			create(6, b4, -, b2);
			create(7, s6, -, b3);
		}
	}

#undef matrices

	__matrix_elem_type** a12 = new __matrix_elem_type * [size];
	for (size_t i = 0; i < size; ++i) { a12[i] = l[i] + size; }
	__matrix_elem_type** a22 = new __matrix_elem_type * [size];
	for (size_t i = 0; i < size; ++i) { a22[i] = l[size + i] + size; }

	__matrix_elem_type** b21 = nullptr;
	__matrix_elem_type** b22 = nullptr;
	if (!transpond) {
		b21 = new __matrix_elem_type * [size];
		for (size_t i = 0; i < size; ++i) { b21[i] = r[size + i]; }
		b22 = new __matrix_elem_type * [size];
		for (size_t i = 0; i < size; ++i) { b22[i] = r[size + i] + size; }
	}
	
	__matrix_elem_type** matrices[14];
	matrices[0] = S[1];
	matrices[1] = l; //a11
	matrices[2] = a12;
	matrices[3] = S[2];
	matrices[4] = S[0];
	matrices[5] = S[3];
	matrices[6] = a22;

	matrices[10] = S[6];
	matrices[13] = S[7];

	if (transpond) {
		copy_shrink_t(7, S[5][i]);
		copy_shrink_t(8, b1);
		copy_shrink_t(9, b3);
		copy_shrink_t(11, S[4][i])
		copy_shrink_t(12, b4);
	}
	else {
		matrices[7] = S[5];
		matrices[8] = r;
		matrices[9] = b21;
		matrices[11] = S[4];
		matrices[12] = b22;
	}
	

	__matrix_elem_type** p[7];
	if (th > 0) {
		std::future<__matrix_elem_type**> f[7];
		th -= 7;
		for (size_t i = 0; i < 7; ++i) {
			if (transpond) f[i] = std::async(&cubic_mul_t<__matrix_elem_type>, matrices[i], matrices[i + 7], size);
			else f[i] = std::async(&winograd<__matrix_elem_type>, matrices[i], matrices[i + 7], std::ref(th), size);
		}
		for (size_t i = 0; i < 7; ++i) {
			p[i] = f[i].get();
		}
		th += 7;
	}
	else {
		for (size_t i = 0; i < 7; ++i) {
			p[i] = transpond
				? cubic_mul_t(matrices[i], matrices[i + 7], size)
				: winograd(matrices[i], matrices[i + 7], th, size);
		}
	}
	
	if (transpond) {
		delete[] matrices[7][0];
		delete[] matrices[7];
		delete[] matrices[8][0];
		delete[] matrices[8];
		delete[] matrices[9][0];
		delete[] matrices[9];
		delete[] matrices[11][0];
		delete[] matrices[11];
		delete[] matrices[12][0];
		delete[] matrices[12];
	}

	delete[] a12;
	delete[] a22;
	delete[] b21;
	delete[] b22;

	for (size_t i = 0; i < 8; ++i) {
		delete[] S[i][0];
		delete[] S[i];
	}

	__matrix_elem_type** T[2];
#define matrices T
	create(0, p[0][i], +, p[1][i]);
	create(1, p[3][i], +, T[0][i]);
#undef matrices

	size *= 2;
	__matrix_elem_type** res = new __matrix_elem_type * [size];
	res[0] = new __matrix_elem_type[size * size]();
	for (size_t i = 1; i < size; ++i) { res[i] = *res + i * size; }
	size /= 2;
	if (th > 0) {
		th -= 4;
		std::future<void> f[4];
		f[0] = std::async(std::launch::async,
			[&]() {
				for (size_t i = 0; i < size; ++i) {
					register auto* const ri = res[i];
					register auto const* const p2i = p[1][i];
					register auto const* const p3i = p[2][i];
					for (size_t j = 0; j < size; ++j) {
						ri[j] = p2i[j] + p3i[j];
					}
				}
			});
		f[1] = std::async(std::launch::async,
			[&]() {
				for (size_t i = 0; i < size; ++i) {
					register auto* const ri = res[i] + size;
					register auto const* const t1i = T[0][i];
					register auto const* const p5i = p[4][i];
					register auto const* const p6i = p[5][i];
					for (size_t j = 0; j < size; ++j) {
						ri[j] = t1i[j] + p5i[j] + p6i[j];
					}
				}
			});
		f[2] = std::async(std::launch::async,
			[&]() {
				for (size_t i = 0; i < size; ++i) {
					register auto* const ri = res[i + size];
					register auto const* const t2i = T[1][i];
					register auto const* const p7i = p[6][i];
					for (size_t j = 0; j < size; ++j) {
						ri[j] = t2i[j] - p7i[j];
					}
				}
			});
		f[3] = std::async(std::launch::async,
			[&]() {
				for (size_t i = 0; i < size; ++i) {
					register auto* const ri = res[i + size] + size;
					register auto const* const t2i = T[1][i];
					register auto const* const p5i = p[4][i];
					for (size_t j = 0; j < size; ++j) {
						ri[j] = t2i[j] + p5i[j];
					}
				}
			});
		for (size_t i = 0; i < 4; ++i) f[i].get();
		th += 4;
	} else {
		for (size_t i = 0; i < size; ++i) {
			register auto* const ri = res[i];
			register auto const* const p2i = p[1][i];
			register auto const* const p3i = p[2][i];
			for (size_t j = 0; j < size; ++j) {
				ri[j] = p2i[j] + p3i[j];
			}
		}
		for (size_t i = 0; i < size; ++i) {
			register auto* const ri = res[i] + size;
			register auto const* const t1i = T[0][i];
			register auto const* const p5i = p[4][i];
			register auto const* const p6i = p[5][i];
			for (size_t j = 0; j < size; ++j) {
				ri[j] = t1i[j] + p5i[j] + p6i[j];
			}
		}
		for (size_t i = 0; i < size; ++i) {
			register auto* const ri = res[i + size];
			register auto const* const t2i = T[1][i];
			register auto const* const p7i = p[6][i];
			for (size_t j = 0; j < size; ++j) {
				ri[j] = t2i[j] - p7i[j];
			}
		}
		for (size_t i = 0; i < size; ++i) {
			register auto* const ri = res[i + size] + size;
			register auto const* const t2i = T[1][i];
			register auto const* const p5i = p[4][i];
			for (size_t j = 0; j < size; ++j) {
				ri[j] = t2i[j] + p5i[j];
			}
		}
	}
	for (size_t i = 0; i < 7; ++i) {
		delete[] p[i][0];
		delete[] p[i];
	}
	for (size_t i = 0; i < 2; ++i) {
		delete[] T[i][0];
		delete[] T[i];
	}
	return res;
}

#undef s1
#undef s2
#undef s5
#undef s6
#undef a1
#undef a2
#undef a3
#undef a4
#undef b1
#undef b2
#undef b3
#undef b4
#undef c2da
#undef create_t
#undef create

#endif //!_Strassen_