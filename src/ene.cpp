#include "ene.h"

pset_family ENECLASS::cpp_job() {

	chunk_num();
	std::vector<std::thread> vec;
	int t = THREADNUM > maxpset ? maxpset : THREADNUM;
	for (int i = 0; i < t - 1; i++) {
		vec.emplace_back(&ENECLASS::for_symmetry1, this);
	}
	for_symmetry1();
	for_each(vec.begin(), vec.end(), mem_fn(&std::thread::join));
	allflip();
	if (print_sym) outputSym();

	return sym;
}
pset_family ENECLASS::tbb_job() {
	chunk_num();
	int k = maxpset;
	tbb::parallel_for(0, k, 1, [=](unsigned long long i) {
		for_symmetry2(i);
		});
	allflip();
	if (print_sym) outputSym();
	return sym;
}
pset_family ENECLASS::openmp_job() {

	chunk_num();
#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < maxpset; i++) {
			for_symmetry2(i);
		}
	}
	allflip();
	if (print_sym) outputSym();

	return sym;
}

void ENECLASS::for_symmetry1() {
	pset fp, rp, k = set_new(NUMINPUTS * 2);
	int t = 0;

	unsigned long long pairid, chunk;
	int fend, rend;

	while ((pairid = aui.fetch_add(1)) < maxpset) {
		chunk = chunk_table[pairid];
		//chunk = pairid;
		fp = F->data + F->wsize * (chunk / col) * chunk_size;
		if (chunk / col == row - 1) fend = F->count % chunk_size;
		else fend = chunk_size;
		if (chunk % col == col - 1) rend = R->count % chunk_size;
		else rend = chunk_size;

		for (int i = 0; i < fend && !is_not_sym; i++, fp += F->wsize) {
			rp = R->data + R->wsize * (chunk % col) * chunk_size;
			for (int j = 0; j < rend && !is_not_sym; j++, rp += R->wsize) {
				if (OutputIntersection(fp, rp)) {
					int dist_on[2] = { -1, -1 };
					int dist = cdist123(fp, rp, dist_on);
					if (dist > 2) continue;
					else if (is_in_set(terminal_cond, dist_on[0])) {
						//skip.fetch_add(1);
						continue;
					}
					else if (dist == 2) {
						int position = dist_on[1] * 2;
						if (getvar(fp, dist_on[0]) == getvar(fp, dist_on[1])) {
							position += 1;
						}
						ddd[dist_on[0]]->to_remove(position);
					}
					else if (dist == 1) {
						k = set_fill(k, NUMINPUTS * 2);
						if (getvar(fp, dist_on[0]) == '1') {
							k = remove_symmetry(k, fp, rp);
						}
						else {
							k = remove_symmetry(k, rp, fp);
						}
						ddd[dist_on[0]]->to_and(k);
					}
					t++;
					if (t == Bat) {
						t = 0;
						flip();
						/*如果全部沒有對稱了就break*/
						is_not_sym = check_symmetry();
						if (is_not_sym) {
							break;
						}
					}
				}
			}
		}
	}
	FREE(k);
}
void ENECLASS::for_symmetry2(unsigned long long pairid) {
	pset fp, rp, k = set_new(NUMINPUTS * 2);
	int t = 0;

	unsigned long long chunk;
	int fend, rend;

	//cout << pairid << endl;
	//chunk = pairid;
	chunk = chunk_table[pairid];
	fp = F->data + F->wsize * (chunk / col) * chunk_size;
	if (chunk / col == row - 1) fend = F->count % chunk_size;
	else fend = chunk_size;
	if (chunk % col == col - 1) rend = R->count % chunk_size;
	else rend = chunk_size;

	for (int i = 0; i < fend && !is_not_sym; i++, fp += F->wsize) {
		rp = R->data + R->wsize * (chunk % col) * chunk_size;
		for (int j = 0; j < rend && !is_not_sym; j++, rp += R->wsize) {
			if (OutputIntersection(fp, rp)) {
				int dist_on[2] = { -1, -1 };
				int dist = cdist123(fp, rp, dist_on);
				if (dist > 2) continue;
				else if (is_in_set(terminal_cond, dist_on[0])) {
					//skip.fetch_add(1);
					continue;
				}
				else if (dist == 2) {
					int position = dist_on[1] * 2;
					if (getvar(fp, dist_on[0]) == getvar(fp, dist_on[1])) {
						position += 1;
					}
					ddd[dist_on[0]]->to_remove(position);
				}
				else if (dist == 1) {
					k = set_fill(k, NUMINPUTS * 2);
					if (getvar(fp, dist_on[0]) == '1') {
						k = remove_symmetry(k, fp, rp);
					}
					else {
						k = remove_symmetry(k, rp, fp);
					}
					ddd[dist_on[0]]->to_and(k);
				}
				t++;
				if (t == Bat) {
					t = 0;
					flip();
					/*如果全部沒有對稱了就break*/
					is_not_sym = check_symmetry();
					if (is_not_sym) {
						break;
					}
				}
			}
		}

	}
	//finish_flag = true;
	FREE(k);
}
void ENECLASS::flip() {
	pset p, last, q;

	for (int i = 0; i < NUMINPUTS; i++) {
		for (int j = i + 1; j < NUMINPUTS; j++) {
			p = ddd[j]->get_pset();
			{
				if (!is_in_set(p, i * 2)) {
					ddd[i]->to_remove(j * 2);
				}
				if (!is_in_set(p, i * 2 + 1)) {
					ddd[i]->to_remove(j * 2 + 1);
				}
			}
		}
	}
}
void ENECLASS::allflip() {
	pset p, last, q;

	for (int i = 0; i < NUMINPUTS; i++) {
		q = sym->data + sym->wsize * i;
		for (int j = i + 1; j < NUMINPUTS; j++) {
			p = sym->data + sym->wsize * j;
			{
				if (!is_in_set(p, i * 2)) {
					set_remove(q, j * 2);
				}
				if (!is_in_set(p, i * 2 + 1)) {
					set_remove(q, j * 2 + 1);
				}
				if (!is_in_set(q, j * 2)) {
					set_remove(p, i * 2);
				}
				if (!is_in_set(q, j * 2 + 1)) {
					set_remove(p, i * 2 + 1);
				}
			}
		}
	}
}
bool ENECLASS::check_symmetry() {
	bool terminal = true;
	std::vector<bool> check(NUMINPUTS, false);

	for (int i = 0; i < NUMINPUTS; i++) {
		if (!is_in_set(terminal_cond, i)) {
			for (int j = i + 1; j < NUMINPUTS; j++) {
				if (getvar(ddd[i]->get_pset(), j) != '?') {
					check[i] = true;
					check[j] = true;
					terminal = false;
				}
			}
			if (!check[i]) {
				set_insert(terminal_cond, i);
			}
		}
	}
	return terminal;
}
void ENECLASS::outputSym() {
	std::vector<int> times(NUMINPUTS + 1, 0);
	std::vector<int> isSym(NUMINPUTS, 1);
	int cnt = 0, count;
	for (int i = 0; i < NUMINPUTS; i++) {
		count = 0;
		pset p = ddd[i]->get_pset();
		for (int j = i; j < NUMINPUTS; j++) {
			if (isSym[j]) {
				if (is_in_set(p, j * 2)) {
					isSym[j] = 0;
					count++;
				}
				else if (is_in_set(p, j * 2 + 1)) {
					isSym[j] = 0;
					count++;
				}
			}
		}
		if (count != 0)
			times[count]++;
	}

	for (int i = 0; i <= NUMINPUTS; i++) {
		if (times[i] != 0) {
			printf("%d(%d)", times[i], i);
		}
	}
	printf("\n");
}
void ENECLASS::chunk_num() {
	std::vector<std::vector<unsigned long long>> v1(row, std::vector<unsigned long long>(col, -1));
	int cnt = 0;
	int c = 0, r = 0, i = 0, j = 0;

	while (cnt < maxpset) { //斜線
		v1[r][c] = cnt;
		cnt++;
		if (r + 1 >= row || c - 1 < 0) {
			if (i < col - 1) {
				i++;
				r = j;
				c = i;
			}
			else {
				j++;
				r = j;
				c = col - 1;
			}
		}
		else {
			r++;
			c--;
		}
	}
	cnt = 0;
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			chunk_table[v1[i][j]] = cnt;
			cnt++;
		}
	}
}
void ENECLASS::printSym() {
	for (int i = 0; i < NUMINPUTS; i++) {
		printf("[%5d]", i);
		ddd[i]->print_set();
	}
}

void ene(pset_family ENE, pset_family F, pset_family R, int option) {
	ENECLASS sym(F, R);
	if (option == 1) {
		sf_copy(ENE, sym.tbb_job());
	}
	else if (option == 2) {
		sf_copy(ENE, sym.openmp_job());
	}
	else {
		sf_copy(ENE, sym.cpp_job());
	}
}

pset set_dc(pset r, pset a) {
	unsigned int x;
	int tmp;

	if ((tmp = cube.inword) != -1) {
		x = (a[tmp] >> 1) & a[tmp];
		x = x & cube.inmask;
		r[tmp] = x | (x << 1);
	}
	for (int i = 1; i < tmp; i++) {
		x = (a[i] >> 1) & a[i];
		x = x & DISJOINT;
		r[i] = x | (x << 1);
	}
	return r;
}
int cdist123(pset a, pset b, int dist_on[2]) {
	unsigned int x;
	int tmp, ones; // ones用於存有多少1
	int counts = 0; // counts用於存到第幾個dist_on
	int dist = 0, n = 0x1;

	tmp = cube.inword;

	for (int j = 1; j < tmp; j++) {
		x = a[j] ^ b[j];
		x = (x >> 1) & x & DISJOINT;
		ones = count_ones(x);
		if (dist + ones > 2) return 3;
		else if (ones > 0) {
			dist += ones;
			for (int i = 0; i < BPI && ones; i += 2) {
				if (x & (n << i)) {
					dist_on[counts] = (j - 1) * varnum + i / 2;
					ones--;
					counts++;
				}
			}
		}
	}

	if (tmp != -1) {
		x = a[tmp] ^ b[tmp];
		x = (x >> 1) & x & cube.inmask;
		ones = count_ones(x);
		if (dist + ones > 2) return 3;
		else if (ones > 0) {
			dist += ones;
			for (int i = 0; i < BPI && ones; i += 2) {
				if (x & (n << i)) {
					dist_on[counts] = (tmp - 1) * varnum + i / 2;
					ones--;
					counts++;
				}
			}
		}
	}

	return dist;
}
pset remove_symmetry(pset r, pset one, pset two) {
	unsigned int x1, x2, x3, notdc;
	int tmp;
	pset dcpart = set_dc(set_new(cube.size), one);

	if ((tmp = cube.inword) != -1) {
		x1 = one[tmp] ^ two[tmp];
		x3 = one[tmp] | two[tmp];
		x3 = (x3 >> 1) & x3 & cube.inmask;
		x3 = x3 | (x3 << 1);
		if (dcpart[tmp]) {
			x2 = x1 & dcpart[tmp];
			x1 = x1 & (~dcpart[tmp]);
			r[tmp] = ((x2 & cube.inmask) << 1) | ((x2 >> 1) & cube.inmask) | x1 | (~x3);
		}
		else r[tmp] = x1 | (~x3);
	}

	for (int i = 1; i < tmp; i++) {
		x1 = one[i] ^ two[i];
		x3 = one[i] | two[i];
		x3 = (x3 >> 1) & x3 & DISJOINT;
		x3 = x3 | (x3 << 1);
		if (dcpart[i]) {
			x2 = x1 & dcpart[i];
			x1 = x1 & (~dcpart[i]);
			r[i] = ((x2 & DISJOINT) << 1) | ((x2 >> 1) & DISJOINT) | x1 | (~x3);
		}
		else r[i] = x1 | (~x3);
	}
	return r;
}
pset and_triangle(pset a, pset b) { //只and上三角形 有必要嗎
	return a;
}
void binary_out(unsigned int a) { //用於查看pset中的word為何
	int n = 0x00000001;
	for (int i = 0; i < BPI; i++) {
		if (a & n) std::cout << "1";
		else std::cout << "0";
		n <<= 1;
	}
}
void find_Symmetry(pset_family A) {
	std::vector<int> times(NUMINPUTS + 1, 0);
	std::vector<int> isSym(NUMINPUTS, 1);
	int cnt = 0, count;
	for (int i = 0; i < NUMINPUTS; i++) {
		count = 0;
		pset p = A->data + A->wsize * i;
		for (int j = i; j < NUMINPUTS; j++) {
			if (isSym[j]) {
				if (is_in_set(p, j * 2)) {
					isSym[j] = 0;
					count++;
				}
				else if (is_in_set(p, j * 2 + 1)) {
					isSym[j] = 0;
					count++;
				}
			}
		}
		if (count != 0)
			times[count]++;
	}

	for (int i = 0; i <= NUMINPUTS; i++) {
		if (times[i] != 0) {
			printf("%d(%d)", times[i], i);
		}
	}
}
void printSym(pset_family A) {
	for (int i = 0; i < NUMINPUTS; i++) {
		printf("[%5d]", i);
		pset p = A->data + A->wsize * i;
		std::cout << pbv1(p, NUMINPUTS * 2) << std::endl;
	}
}
void flip_Sym(pset_family A) {
	pset p, last, q;

	int cnt = 0;
	foreach_set(A, last, p) {
		q = A->data + cnt * A->wsize;
		for (int j = cnt + 1; j < NUMINPUTS; j++) {
			q += A->wsize;
			if (!is_in_set(q, cnt * 2)) {
				set_remove(p, j * 2);
			}
			if (!is_in_set(q, cnt * 2 + 1)) {
				set_remove(p, j * 2 + 1);
			}
		}
		cnt++;
	}
}
bool check_symmetry(pset tc, pset_family A) {
	bool terminal = true;
	std::vector<bool> check(NUMINPUTS, false);
	pset p;

	for (int i = 0; i < NUMINPUTS; i++) {
		if (!is_in_set(tc, i)) {
			p = A->data + A->wsize * i;
			for (int j = i + 1; j < NUMINPUTS; j++) {
				if (getvar(p, j) != '?') {
					check[i] = true;
					check[j] = true;
					terminal = false;
				}
			}
			if (!check[i])
				set_insert(tc, i);
		}
	}
	return terminal;
}
void printftriangle(pset_family A) {
	pset p;
	int i;

	foreachi_set(A, i, p) {
		printf("[%5d]", i);
		for (int j = 0; j < NUMINPUTS; j++) {
			if (j < i) {
				std::cout << "  ";
			}
			else {
				if (is_in_set(p, j * 2)) std::cout << "1";
				else std::cout << "0";
				if (is_in_set(p, j * 2 + 1)) std::cout << "1";
				else std::cout << "0";
			}
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}



