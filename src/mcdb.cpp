#include "MCDB.h"

/* 找出literal-0 與literal-1 並且計算DC的個數*/
void MCDB::non_DC(pset lit01, pset lit10, pset _p) {
	int w = 1, mask;
	for (; w < cube.last_word[cube.num_binary_vars-1]; w++) {
		mask = _p[w] >> 1 & _p[w] & DISJOINT;
		mask |= mask << 1; /* mask 是DC的mask*/
		lit01[w] = _p[w] & ~mask;
		lit10[w] = (~_p[w]) & ~mask;
	}
	mask = _p[w] >> 1 & _p[w] & DISJOINT;
	mask |= mask << 1;
	lit01[w] = _p[w] & ~mask & cube.binary_mask[w];
	lit10[w] = (~_p[w]) & ~mask & cube.binary_mask[w];
	PUTSIZE(lit01, SIZE(_p));
	PUTSIZE(lit10, SIZE(_p));
}
/* count DC num*/
int MCDB::set_count_inputs_ones(pset _p) {
	int cnt = 0, i = cube.last_word[cube.num_binary_vars-1];
	int w = 1;
	unsigned int a;
	for (; w < cube.last_word[cube.num_binary_vars]; w++) {
		a = ((_p[w] >> 1 & _p[w]) & DISJOINT);
		cnt += count_ones(a);
	}
	a = ((_p[w] >> 1 & _p[w]) & DISJOINT) & cube.binary_mask[w];
	cnt += count_ones(a);
	
	return cnt;
}

/* 計算non DC 中lit 0 的個數*/
int MCDB::cnt_non_DC_0(pset _p) {
	int w = 1, cnt = 0;
	for (; w < cube.last_word[cube.num_binary_vars]; w++) {
		cnt += count_ones(_p[w] & Literal_0);
	}
	cnt += count_ones(_p[w] & cube.binary_mask[w] & Literal_0);
	return cnt;
}

/* 計算non DC 中lit 1 的個數*/
int MCDB::cnt_non_DC_1(pset _p) {
	int w = 1, cnt = 0;
	for (; w < cube.last_word[cube.num_binary_vars]; w++) {
		cnt += count_ones(_p[w] & Literal_1);
	}
	cnt += count_ones(_p[w] & cube.binary_mask[w] & Literal_1);
	return cnt;
}
/* 複製儲存literal position 的pset*/
/* a[i] = 0 表示位置已複製完 */
void MCDB::lit_pos_copy(pset r, pset a) {
	int i = 0;
	do r[i] = a[i]; while ((i == 1 || a[i] != 0) && ++i <= LOOP(a));
}

/* 將sf 中的 pset a 擴展成size倍數的pset r  ex: a=123, 3 => r=123123123*/
void MCDB::expand_sf_Q(unsigned int* first_word, unsigned int* last_word, unsigned int* first_bit, unsigned int* last_bit) {
	foreachi_set(sf_Q_non_DC, index1, p) {
		p2 = set_new(ex_Q_non_DC->sf_size);
		for (int var = 0; var < num_mv_var; var++) {
			or_exact_MvVar(p2, p, last_word[0], first_word[var], last_word[var], first_bit[var], last_bit[var]);
		}
		sf_addset(ex_Q_non_DC, p2);
	}
	foreachi_set(sf_q_non_DC, index1, p) {
		p2 = set_new(ex_q_non_DC->sf_size);
		for (int var = 0; var < num_mv_var; var++) {
			or_exact_MvVar(p2, p, last_word[0], first_word[var], last_word[var], first_bit[var], last_bit[var]);
		}
		sf_addset(ex_q_non_DC, p2);
	}
}

#define RULE_TEST(key, cost, cnt) {\
		sf_copy3(TABLE_PQ, R0_TABLE_PQ);\
		cnt_R1 = 0, cnt_R2 = 0, cnt_R3 = 0, cnt_R4 = 0;\
		EXECost_W(check_Rule(key), cost);\
		cnt = cal_TABLE_PQ();\
	}// cnt = (cnt_R1 + cnt_R2 + cnt_R3 + cnt_R4); 
void MCDB::run_RULE_TEST() {
	pset_family R0_TABLE_PQ = sf_new(TABLE_PQ->count, TABLE_PQ->sf_size);
	sf_copy3(R0_TABLE_PQ, TABLE_PQ);
	clock_t A = 0, B = A, C = A, D = A, AB = A, AC = A, AD = A, BC = A, BD = A, CD = A, ABC = A, ABD = A, ACD = A, BCD = A, ABCD = A;
	int cntA = 0, cntB = 0, cntC = 0, cntD = 0, cntAB = 0, cntAC = 0, cntAD = 0, cntBC = 0, cntBD = 0, cntCD = 0, cntABC = 0, cntABD = 0, cntACD = 0, cntBCD = 0, cntABCD = 0;

	RULE_TEST(0b0001, A, cntA);
	RULE_TEST(0b0010, B, cntB);
	RULE_TEST(0b0100, C, cntC);
	RULE_TEST(0b1000, D, cntD);
	RULE_TEST(0b0011, AB, cntAB);
	//RULE_TEST(0b0101, AC, cntAC);
	//RULE_TEST(0b1001, AD, cntAD);
	//RULE_TEST(0b0110, BC, cntBC);
	//RULE_TEST(0b1010, BD, cntBD);
	RULE_TEST(0b1100, CD, cntCD);
	RULE_TEST(0b0111, ABC, cntABC);
	//RULE_TEST(0b1011, ABD, cntABD);
	//RULE_TEST(0b1101, ACD, cntACD);
	//RULE_TEST(0b1110, BCD, cntBCD);
	RULE_TEST(0b1111, ABCD, cntABCD);
	
	printf("RULE TEST : A, B, C, D, AB, AC, AD, BC, BD, CD, ABC, ABD, ACD, BCD, ABCD\n");
	printf("RULE TEST cnt : %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", cntA, cntB, cntC, cntD, cntAB, cntAC, cntAD, cntBC, cntBD, cntCD, cntABC, cntABD, cntACD, cntBCD, cntABCD);
	printf("RULE TEST COST: %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf\n",
		(double)A / (double)CLOCKS_PER_SEC, (double)B / (double)CLOCKS_PER_SEC, (double)C / (double)CLOCKS_PER_SEC, (double)D / (double)CLOCKS_PER_SEC,
		(double)AB / (double)CLOCKS_PER_SEC, (double)AC / (double)CLOCKS_PER_SEC, (double)AD / (double)CLOCKS_PER_SEC, (double)BC / (double)CLOCKS_PER_SEC,
		(double)BD / (double)CLOCKS_PER_SEC, (double)CD / (double)CLOCKS_PER_SEC, (double)ABC / (double)CLOCKS_PER_SEC, (double)ABD / (double)CLOCKS_PER_SEC,
		(double)ACD / (double)CLOCKS_PER_SEC, (double)BCD / (double)CLOCKS_PER_SEC, (double)ABCD / (double)CLOCKS_PER_SEC);
	sf_free(R0_TABLE_PQ);
}