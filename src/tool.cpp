

#include "tool.h"

char getvar(pset a, int position) {
	char c = "?012"[GETINPUT(a, position)];
	return c;
}

int cdist1(pset a, pset b) {
	unsigned int x;
	int tmp, ones; // ones用於存有多少1
	int dist = 0, n = 0x1, var = -1;
	bool dist_is_one = true;

	tmp = cube.inword;

	for (int j = 1; j < tmp; j++) {
		x = a[j] ^ b[j];
		x = (x >> 1) & x & DISJOINT;
		ones = count_ones(x);

		if (dist + ones > 1) return -1;
		else {
			dist += ones;
			for (int i = 0; i < BPI && ones; i += 2) {
				if (x & (n << i)) {
					var = (j - 1) * BPI / 2 + i / 2;
					ones--;
				}
			}
		}
	}

	if (tmp != -1) {
		x = a[tmp] ^ b[tmp];
		x = (x >> 1) & x & cube.inmask;
		ones = count_ones(x);
		if (dist + ones > 1) return -1;
		else if (ones > 0) {
			dist += ones;
			for (int i = 0; i < BPI && ones; i += 2) {
				if (x & (n << i)) {
					var = (tmp - 1) * BPI / 2 + i / 2;
					ones--;
				}
			}
		}
	}

	return var;
}


/* check if pset A and pset B output is intersection*/
bool OutputIntersection(pset A, pset B) {
	//cost_temp = clock();
	unsigned int w = cube.first_word[cube.num_binary_vars]; /*word index*/
	unsigned int temp;
	temp = A[w] & B[w] & cube.var_mask[cube.num_binary_vars][w];
	if (temp > 0)
		return true;
	for (w += 1; w < LOOP(A) + 1; w++) {
		temp = A[w] & B[w];
		if (temp > 0) {
			return true;
		}
	}
	return false;
}

/* inverse a set family*/
pset_family sf_inverse(pset_family sf) {
	pset p, last;
	int w_size = ceil(sf->sf_size / 32.0) + 1;
	foreach_set(sf, last, p) {
		for (int w = 1; w < w_size; w++) {
			p[w] = ~p[w];
		}
	}
	return sf;
}


pset set_not(pset r, pset p) {
	r[0] = p[0];
	for (int w = 1; w <= LOOP(p); w++) {
		r[w] = ~p[w];
	}
	return r;
}

static void intcpy(unsigned int* d, unsigned int* s, long int n)
{
	int i;
	for (i = 0; i < n; i++) {
		*d++ = *s++;
	}
}

pset_family my_sf_copy(pset_family R, pset_family A)
{
	R->sf_size = A->sf_size;
	R->wsize = A->wsize;
	R->capacity = A->count;
	R->data = REALLOC(unsigned int, R->data, (long)R->capacity * R->wsize);
	R->count = A->count;
	R->active_count = A->active_count;
	intcpy(R->data, A->data, (long)A->wsize * A->count);
	return R;
}

/* sf_copy -- copy a set family */ /* 複製但不刪除原本的空間，但複製的是原資料的address*/
pset_family sf_copy2(pset_family R, pset_family A)
{
	R->sf_size = A->sf_size;
	R->wsize = A->wsize;
	R->capacity = A->count;
	R->data = REALLOC(unsigned int, R->data, (long)R->capacity * R->wsize);
	R->count = A->count;
	R->active_count = A->active_count;
	intcpy(R->data, A->data, (long)A->wsize * A->count);
	return R;
}

void intcpy2(unsigned int* dest, unsigned int* src, long int n)
{
	memcpy(dest, src, n);
}
/* 複製但不刪除原本的空間，連資料都是複製的*/
pset_family sf_copy3(pset_family R, pset_family A)
{
	R->sf_size = A->sf_size;
	R->wsize = A->wsize;
	R->capacity = A->count;
	R->data = REALLOC(unsigned int, R->data, (long)R->capacity * R->wsize);
	R->count = A->count;
	R->active_count = A->active_count;
	//intcpy2(R->data, A->data, (long)A->wsize * A->count);
	pset p1, p2;
	int index;
	foreachi_set(A, index, p1) {
		p2 = GETSET(R, index);
		set_copy(p2, p1);
	}
	return R;
}

/* ????*/
pset_family my_sf_append(pset_family A, pset_family B)
{
	long asize = A->count * A->wsize;
	long bsize = B->count * B->wsize;

	if (A->sf_size != B->sf_size) fatal("sf_append: sf_size mismatch");
	A->capacity = A->count + B->count;
	A->data = REALLOC(unsigned int, A->data, (long)A->capacity * A->wsize);
	intcpy(A->data + asize, B->data, bsize);
	A->count += B->count;
	A->active_count += B->active_count;
	return A;
}

/* sf_or -- form the "or" of all sets in a set family */
pset sf_or2(pset or , pset_family A)
{
	pset last, p;
	set_clear(or , A->sf_size);
	foreach_set(A, last, p)
		INLINEset_or(or , or , p);
	return or ;
}

pset or_exact_MvVar(pset r, pset p, int p_w_size, unsigned int first_word, unsigned int last_word, unsigned int first_bit, unsigned int last_bit) {
	unsigned int val = 0;
	for (unsigned int M = first_word, P = 1; M <= last_word; M++, P++) {
		if (P <= p_w_size)
			val |= (p[P] << first_bit);
		r[M] |= val;
		if (BPI - first_bit != 32)
			val = ((unsigned int)(p[P]) >> (BPI - first_bit)); // p[P]中剩餘的bits
		else
			val = 0;
	}
	return r;
}


pset get_mvcube_var(pset r, pset p, int p_w_size, unsigned int first_word, unsigned int last_word, unsigned int first_bit, unsigned int last_bit) {
	int P = first_word, R = 1, val;

	for (; P <= last_word && R <= p_w_size; P++, R++) {
		val = p[P] >> first_bit; /* var 的第一個bit對齊至bit 0*/
		if (P != last_word && first_bit != 0) /* 若var還有則將var下一個word的bits 補到前一個word中*/
			val |= (p[P + 1] << (BPI - first_bit));
		if (R == LOOP(r))  /* 提取的最後一個var需要做mask的處理*/
			val &= ~(0xFFFFFFFF << (((BPI - first_bit) + (last_bit + 1)) % 32));
		r[R] = val;
	}
	return r;
}

pset rshilft_or_disjoint(pset r, pset p) {
	r[0] = p[0];
	for (int i = 1; i <= LOOP(p); i++) {
		r[i] = (p[i] | (p[i] >> 1)) & DISJOINT;
	}
	return r;
}

pset set_inverse_disjoint(pset r, pset p) {
	r[0] = p[0];
	for (int i = 1; i <= LOOP(p); i++) {
		r[i] = (~p[i]) & DISJOINT;
	}
	return r;
}

int set_count_ones(pset p) {
	int cnt = 0;
	for (int i = 1; i <= LOOP(p); i++) {
		cnt += count_ones(p[i]);
	}
	return cnt;
}


// 適用ENE, KSIG_output 的pset分群
void pset_grouping(pset p, int var_size, int var_num, bool display) {
	std::vector<int> Group(var_num, -1); //(group_id)
	std::vector<bool> is_grouped(var_num, false);
	int group_id;
	
	pset_family sf_var = sf_new(var_num, var_size);
	pset pp = set_new(var_size);
	for (int i = 0; i < var_num; i++) {
		set_clear(pp, var_size);
		for (int j = 0; j < var_size; j++) {
			if (is_in_set(p, i * var_size + j))
				set_insert(pp, j);
		}
		sf_addset(sf_var, pp);
	}

	if (display) display_sf(sf_var);

	group_id = 0;
	for (int i = 0; i < var_num; i++) {
		if (is_grouped[i])
			continue;
		Group[i] = group_id; is_grouped[i] = true; 
		group_id++;
		for (int j = i+1; j < var_num; j++) {
			if (is_grouped[j])
				continue;
			if (setp_equal(GETSET(sf_var, i), GETSET(sf_var, j))) {
				Group[j] = Group[i]; is_grouped[j] = true;
			}
		}
	}
	if (display) {
		for (int& group_id : Group) {
			std::cout << group_id << " ";
		} std::cout << std::endl;
	}
	std::sort(Group.begin(), Group.end());
	std::vector<int> cnt_group(var_num, 0);

	int cnt_id = 0;

	for (int i = 0; i < Group.size(); i++) {
		if (i == 0) {
			cnt_id = 0;
		}
		else if (Group[i] == Group[i-1]) {
			cnt_id++;
		}
		else {
			cnt_group[cnt_id]++;
			cnt_id = 0;
		}
		if (i == Group.size() - 1) {
			cnt_group[cnt_id]++;
		}
	}
	if (display) {
		for (int& group_id : cnt_group) {
			std::cout << group_id << " ";
		} std::cout << std::endl;
	}
	int first = false;
	printf("output group: ");
	for (int i = cnt_group.size() - 1; i >= 0; i--) {
		if (cnt_group[i] != 0) {
			if (first) printf(" "); first = true;
			printf("%d(%d)", cnt_group[i], i+1);
		}
	} printf("\n");
	
	sf_free(sf_var);
}


pset string2pset(std::string str) {
	pset p = set_new(str.size());
	for (int i = 0; i < str.size(); i++) {
		if (str[i] == '1')
			set_insert(p, i);
	}
	return p;
}


/* sf_new -- allocate "num" sets of "size" elements each */
pset_family sf_new2(int num, int size)
{
	pset_family A;

	A = ALLOC(set_family_t, 1);
	A->sf_size = size;
	A->wsize = SET_SIZE(size);
	A->capacity = num;
	A->data = ALLOC(unsigned int, (long)A->capacity * A->wsize);
	A->count = 0;
	A->active_count = 0;
	return A;
}

/* sf_free -- free the storage allocated for a set family */
void sf_free2(pset_family A)
{
	FREE(A->data);
}


pset set_fill_disjoint(pset r, int size) {
	int i = LOOPINIT(size);
	*r = i;
	r[i] = ~DISJOINT;
	r[i] >>= i * BPI - size;
	while (--i > 0)
		r[i] = ~DISJOINT;
	return r;
}

/*
	vvvvvvvvvvvvvvvvvvvv修改至espresso中的cubestr.cvvvvvvvvvvvvvvvvvvvv
*/
/*
	cube_setup -- assume that the fields "num_vars", "num_binary_vars", and
	part_size[num_binary_vars .. num_vars-1] are setup, and initialize the
	rest of cube and cdata.

	If a part_size is < 0, then the field size is abs(part_size) and the
	field read from the input is symbolic.
*/
void cube_setup(cube_struct &cube)
{
	int i, var;
	pcube p;

	if (cube.num_binary_vars < 0 || cube.num_vars < cube.num_binary_vars)
		fatal("cube size is silly, error in .i/.o or .mv");

	cube.num_mv_vars = cube.num_vars - cube.num_binary_vars;
	cube.output = cube.num_mv_vars > 0 ? cube.num_vars - 1 : -1;

	cube.size = 0;
	cube.first_part = ALLOC(int, cube.num_vars);
	cube.last_part = ALLOC(int, cube.num_vars);
	cube.first_word = ALLOC(int, cube.num_vars);
	cube.last_word = ALLOC(int, cube.num_vars);
	for (var = 0; var < cube.num_vars; var++) {
		if (var < cube.num_binary_vars)
			cube.part_size[var] = 2;
		cube.first_part[var] = cube.size;
		cube.first_word[var] = WHICH_WORD(cube.size);
		cube.size += ABS(cube.part_size[var]);
		cube.last_part[var] = cube.size - 1;
		cube.last_word[var] = WHICH_WORD(cube.size - 1);
	}

	cube.var_mask = ALLOC(pset, cube.num_vars);
	cube.sparse = ALLOC(int, cube.num_vars);
	cube.binary_mask = new_cube();
	cube.mv_mask = new_cube();
	for (var = 0; var < cube.num_vars; var++) {
		p = cube.var_mask[var] = new_cube();
		for (i = cube.first_part[var]; i <= cube.last_part[var]; i++)
			set_insert(p, i);
		if (var < cube.num_binary_vars) {
			INLINEset_or(cube.binary_mask, cube.binary_mask, p);
			cube.sparse[var] = 0;
		}
		else {
			INLINEset_or(cube.mv_mask, cube.mv_mask, p);
			cube.sparse[var] = 1;
		}
	}
	if (cube.num_binary_vars == 0)
		cube.inword = -1;
	else {
		cube.inword = cube.last_word[cube.num_binary_vars - 1];
		cube.inmask = cube.binary_mask[cube.inword] & DISJOINT;
	}

	cube.temp = ALLOC(pset, CUBE_TEMP);
	for (i = 0; i < CUBE_TEMP; i++)
		cube.temp[i] = new_cube();
	cube.fullset = set_fill(new_cube(), cube.size);
	cube.emptyset = new_cube();

	cdata.part_zeros = ALLOC(int, cube.size);
	cdata.var_zeros = ALLOC(int, cube.num_vars);
	cdata.parts_active = ALLOC(int, cube.num_vars);
	cdata.is_unate = ALLOC(int, cube.num_vars);
}
/*
	setdown_cube -- free memory allocated for the cube/cdata structs
	(free's all but the part_size array)

	(I wanted to call this cube_setdown, but that violates the 8-character
	external routine limit on the IBM !)
*/
void setdown_cube(cube_struct& cube)
{
	int i, var;

	FREE(cube.first_part);
	FREE(cube.last_part);
	FREE(cube.first_word);
	FREE(cube.last_word);
	FREE(cube.sparse);

	free_cube(cube.binary_mask);
	free_cube(cube.mv_mask);
	free_cube(cube.fullset);
	free_cube(cube.emptyset);
	for (var = 0; var < cube.num_vars; var++)
		free_cube(cube.var_mask[var]);
	FREE(cube.var_mask);

	for (i = 0; i < CUBE_TEMP; i++)
		free_cube(cube.temp[i]);
	FREE(cube.temp);

	FREE(cdata.part_zeros);
	FREE(cdata.var_zeros);
	FREE(cdata.parts_active);
	FREE(cdata.is_unate);

	cube.first_part = cube.last_part = (int*)NULL;
	cube.first_word = cube.last_word = (int*)NULL;
	cube.sparse = (int*)NULL;
	cube.binary_mask = cube.mv_mask = (pcube)NULL;
	cube.fullset = cube.emptyset = (pcube)NULL;
	cube.var_mask = cube.temp = (pcube*)NULL;

	cdata.part_zeros = cdata.var_zeros = cdata.parts_active = (int*)NULL;
	cdata.is_unate = (bool*)NULL;
}
/*
	^^^^^^^^^^^^^^^^^^^^修改至espresso中的cubestr.c^^^^^^^^^^^^^^^^^^^^
*/



/*
	Display method
*/
void display_pset(pset p, int p_size) {
	for (int i = 0; i < p_size; i++) {
		printf("%d", is_in_set(p, i) ? 1 : 0);
	}printf("\n");
}

void display_pset_eachvar(pset p, int var_size, int var_num) {
	for (int i = 0; i < var_num; i++) {
		for (int j = 0; j < var_size; j++) {
			printf("%d", is_in_set(p, i * var_size + j) ? 1 : 0);
		}printf("\n");
	}
}

void display_binary(int p) {
	for (int i = 0; i < 32; i++) {
		printf("%d", p & (1 << WHICH_BIT(i)) ? 1 : 0);
	}printf("\n");
}

void display_pset_01X(pset p, int input_num, int output_num) {
	for (int i = 0; i < input_num; i++) {
		printf("%c", "?01-"[GETINPUT(p, i)]);
	}printf(" ");
	for (int i = 0; i < output_num; i++) {
		printf("%d", GETOUTPUT(p, i));
	}printf("\n");
}

void display_sf(pset_family sf) {
	pset p;
	int i;
	foreachi_set(sf, i, p) {
		printf("%d: ", i);  display_pset(p, sf->sf_size);
	}
}
void display_sf_eachvar(pset_family sf, int _size, int _num) {
	pset p;
	int i;
	foreachi_set(sf, i, p) {
		printf("%d:\n", i);  display_pset_eachvar(p, _size, _num);
	}
}
