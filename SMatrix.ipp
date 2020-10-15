#ifdef __SMatrix__

template <typename __matrix_elem_type>
void SMatrix<__matrix_elem_type>::assign(const IMatrix<__matrix_elem_type>& r) {
	for (size_t i = this->row - 1; i != (size_t)-1; --i) {
		for (size_t j = this->col - 1; j != (size_t)-1; --j) {
			if (auto v = r.get(i, j)) set(i, j, v);
			//set(i, j, r.get(i, j));
		}
	}
}

template <typename __matrix_elem_type>
void SMatrix<__matrix_elem_type>::clear() {
	node* n = head;
	while (n) {
		node* tmp = n->next;
		delete n;
		n = tmp;
	}
	head = nullptr;
}

template <typename __matrix_elem_type>
void SMatrix<__matrix_elem_type>::re_alloc(size_t r, size_t c) {
	clear();
	this->row = r, this->col = c;
}

template <typename __matrix_elem_type>
SMatrix<__matrix_elem_type>::SMatrix(size_t r, size_t c) {
	this->row = r, this->col = c;
	head = nullptr;
}

template <typename __matrix_elem_type>
SMatrix<__matrix_elem_type>::SMatrix(SMatrix<__matrix_elem_type>&& r) noexcept : IMatrix<__matrix_elem_type>(r) {
	head = r.head;
	r.head = nullptr;
}

template <typename __matrix_elem_type>
SMatrix<__matrix_elem_type>::SMatrix(const SMatrix<__matrix_elem_type>& r) : IMatrix<__matrix_elem_type>(r) {
	assign(r);
}

template <typename __matrix_elem_type>
SMatrix<__matrix_elem_type>::SMatrix(const IMatrix<__matrix_elem_type>& r) : IMatrix<__matrix_elem_type>(r) {
	assign(r);
}

template <typename __matrix_elem_type>
IMatrix<__matrix_elem_type>& SMatrix<__matrix_elem_type>::operator=(SMatrix<__matrix_elem_type>&& r) noexcept {
	if (this == &r) return *this;
	head = r.head;
	r.head = nullptr;
	return *this;
}

template <typename __matrix_elem_type>
IMatrix<__matrix_elem_type>& SMatrix<__matrix_elem_type>::operator=(const SMatrix<__matrix_elem_type>& r) {
	if (this == &r) return *this;
	clear();
	assign(r);
	return *this;
}

template <typename __matrix_elem_type>
IMatrix<__matrix_elem_type>& SMatrix<__matrix_elem_type>::operator=(const IMatrix<__matrix_elem_type>& r) {
	if (this == &r) return *this;
	clear();
	assign(r);
	return *this;
}

template <typename __matrix_elem_type>
SMatrix<__matrix_elem_type>::~SMatrix() {
	clear();
}

template <typename __matrix_elem_type>
void SMatrix<__matrix_elem_type>::set(size_t r, size_t c, __matrix_elem_type value, void(*f)(void)) {
	if (r > this->row || c > this->col) f();
	if (!head) {
		if (!value) return;
		head = new node;
		head->r = r; head->c = c; head->value = value; head->next = nullptr;
		return;
	}
	if (head->r == r && head->c == c) {
		if (!value) {
			node* del = head;
			head = head->next;
			delete del;
			return;
		}
		head->value = value;
		return;
	}
	if (head->r > r || head->r == r && head->c > c) {
		if (!value) return;
		node* tmp = new node;
		tmp->r = r; tmp->c = c; tmp->value = value; tmp->next = head;
		//*tmp = { r, c, value, head };
		head = tmp;
		return;
	}
	node* cur;
	for (cur = head; cur->next; cur = cur->next) {
		if (cur->next->r == r && cur->next->c == c) {
			if (!value) {
				node* del = cur->next;
				cur->next = cur->next->next;
				delete del;
				return;
			}
			cur->next->value = value;
			return;
		}
		if (!cur->next || cur->next->r > r || cur->next->r == r && cur->next->c > c) {
			if (!value) return;
			node* tmp = new node;
			tmp->r = r; tmp->c = c; tmp->value = value;
			tmp->next = cur->next; cur->next = tmp;
			return;
		}
	}
	cur->next = new node;
	cur->next->r = r; cur->next->c = c; cur->next->value = value; cur->next->next = nullptr;
}

template <typename __matrix_elem_type>
__matrix_elem_type SMatrix<__matrix_elem_type>::get(size_t r, size_t c, __matrix_elem_type(*def)(void)) const {
	if (r > this->row || c > this->col) return def();
	for (node* cur = head; cur; cur = cur->next) {
		if (cur->r == r && cur->c == c) return cur->value;
		if (cur->r > r || cur->r == r && cur->c > c) return __matrix_elem_type();
	}
	return __matrix_elem_type();
}

template <typename __matrix_elem_type>
SMatrix<__matrix_elem_type> SMatrix<__matrix_elem_type>::operator+(const IMatrix<__matrix_elem_type>& r) const {
	SMatrix res(*this);
	res += r;
	return res;
}

template <typename __matrix_elem_type>
SMatrix<__matrix_elem_type> SMatrix<__matrix_elem_type>::operator*(const IMatrix<__matrix_elem_type>& r) const {
	SMatrix res(*this);
	res *= r;
	return res;
}

#endif // __SMatrix__