#define __MEMORY_LEAKS_CHECKER__ 0
#if __MEMORY_LEAKS_CHECKER__
	#define __CRTDBG_MAP_ALLOC
	#include <crtdbg.h>
	#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
	#define new DEBUG_NEW
#endif

#include <iostream>
//#include <cstdlib> //srand
//#include <memory> //smart pointers
//#include <thread>
#include <ctime>
#include "SMatrix.h"
#include "UMatrix.h"
#include "Rational.h"
#include "SqUMatrix.h"
using namespace std;

int __main() {
	using T = int;// Rational<int>;
	size_t size = 100;
	srand(1);
	SMatrix<T> s(size, size);
	s.fill_n(size);
	srand(1);
	UMatrix<T> u(size, size);
	u.fill_n(size);
	int th = std::thread::hardware_concurrency();
	//u.strassen(u, th);
	SMatrix<T> s1 = u; //copy constructor S(I)
	UMatrix<T> u1 = s; //copy constructor U(I)
	SMatrix<T> ss = s; //copy constructor S(S)
	UMatrix<T> uu = u; //copy constructor U(U)
	cout << "copy constructors: "
		<< (s == u) << (s1 == u) << (u1 == s) << (u1 == s1) << endl;
	s1 = u; u1 = s; //assignment operator S = I, U = I;
	cout << "assignment operators: " << (s1 == u) << (u1 == s);
	s1 = s; u1 = u; //assignment operator S = S, U = U;
	cout << (u1 == s) << (u1 == s1) << endl;
	
	srand(2);
	SMatrix<T> s2(size, size);
	s2.fill_n(size);

	srand(2);
	UMatrix<T> u2(size, size);
	u2.fill_n(size);

	unique_ptr<IMatrix<T>> sum[8] = {
		unique_ptr<IMatrix<T>>(new auto(s + s2)),
		unique_ptr<IMatrix<T>>(new auto(u + s2)),
		unique_ptr<IMatrix<T>>(new auto(s + u2)),
		unique_ptr<IMatrix<T>>(new auto(u + u2)),
		unique_ptr<IMatrix<T>>(new auto(s2 + s)),
		unique_ptr<IMatrix<T>>(new auto(s2 + u)),
		unique_ptr<IMatrix<T>>(new auto(u2 + s)),
		unique_ptr<IMatrix<T>>(new auto(u2 + u)),
	};

	for (int i = ((cout << "sum: "), 1); i < 8; ++i) {
		cout << (*sum[i - 1] == *sum[i]);
	} cout << endl;

	return 0;
}

int _main(int size) {
	using __matrix_elem_type = int;// Rational<long long>;
	cout << size << endl;
	size_t s = size;
	SqUMatrix<__matrix_elem_type> a(s);
	SqUMatrix<__matrix_elem_type> b(s);
	srand(static_cast<unsigned int>(time(0)));
	a.fill_n(100, []() {return __matrix_elem_type(rand() % 1024); });
	b.fill_n(100, []() {return __matrix_elem_type(rand() % 1024); });
	for (size_t i = 0; i < s; ++i) {
		a.set(i, i, 2);
		b.set(i, i, 3);
	}
	int th = std::thread::hardware_concurrency();
	auto cl = clock();
	auto c = SqUMatrix<__matrix_elem_type>::strassen(a, b, th);// SqUMatrix<__matrix_elem_type>::strassen(a, b, th);
	cout << clock() - cl << endl;
	cl = clock();
	auto d = SqUMatrix<__matrix_elem_type>::GPUmul(a, b, th);
	cout << clock() - cl << endl;
	cout << "equal S & WS: " << (c == d) << endl;
#define deep_check 0
#if deep_check
	auto naive = SqUMatrix<__matrix_elem_type>::mul(a, b, th);
	cout
		<< "naive-Sh: " << (naive == c) << endl
		<< "naive-WS: " << (naive == d) << endl
		;
	printf("   i |    j | naive | c val | d val\n");
	for (size_t i = 0; i < s; ++i) {
		for (size_t j = 0; j < s; ++j) {
			if (naive.get(i, j) == c.get(i, j) && c.get(i, j) == d.get(i, j)) continue;
			cout << i << '|' << j << '|' << naive.get(i, j) << c.get(i, j) << d.get(i, j);
		}
	}
#endif //deep_check
	return 0;
}

int main() {
	int s = (1 << 11) + 0;
	do try { _main(s); }
	catch (std::exception& e) {
		cout << e.what();
		return -1;
	} catch (...) {
		cout << "something was thrown";
		return -1;
	} while (s += 100, 0);
#if __MEMORY_LEAKS_CHECKER__
	_CrtDumpMemoryLeaks();
#endif
	return 0;
}