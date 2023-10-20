/*
	取名叫MvCube DataBase (MCDB) 但這裡實作的是MCQE，MCDB的可能並非完整

	MCQE MvCube Quick Eliminator 有4種Rule
	在計算MCQE前會先將PLA做一定的前處理 將PLA格式轉為類似MvCube的格式 但更為精簡
	之後才計算MCQE，故即使不需使用MCQE仍須使用此library

	使用Rule為一個int 的參數 第0個bit為1 表示使用Rule1 第1個bit為1 表示使用Rule2 以此類推 
*/

#ifndef MCDB_H
#define MCDB_H

#include <iostream>
#include <vector>
#include <utility> // for std::pair
#include "tool.h"

extern "C" {
	#include "espresso.h"
}

#define SORT_1(sf) {sf = sf_contain(sf);}
#define Literal_0 0x55555555
#define Literal_1 0xAAAAAAAA

#define EXECost_W(fct, cost) {cost = clock();fct; cost = clock()-cost;}

class MCDB {
private:
	bool phase_assignment; /* boolean matching 是否考慮input phaseassignment*/
	int num_mv_var, index1, index2, i, j, k, m, n; /*num_mv_var: input 個數, index1 index2: 暫存用變數*/
	int cnt_R0, cnt_R1, cnt_R2, cnt_R3, cnt_R4, cnt_R5;
	pset p, p2, last, last2; /*暫存用變數*/
	pset_family sf_P_non_DC, sf_p_non_DC, sf_Q_non_DC, sf_q_non_DC; /*只保留non don't care部分*/
	clock_t T_lit01, T_SORT, T_Out_IS, T_P_rela, T_Q_rela, T_R1, T_R2, T_R3, T_R4; /*記錄個步驟的執行時間*/
	std::vector<std::pair<int, int>> P_relation, Q_relation, P_re_relation, Q_re_relation; /* 紀錄U, V各自之間的包含關係 (i, j) 表示i 包含j*/

	//void non_DC(pset lit01, pset lit10, pset _p);
	
	void expand_sf_Q(unsigned int* first_word, unsigned int* last_word, unsigned int* first_bit, unsigned int* last_bit); /* 將sf 中的 pset a 擴展成size倍數的pset r  ex: a=123, 3 => r=123123123*/
	int cnt_non_DC_0(pset P);
	int cnt_non_DC_1(pset P);
	int set_count_inputs_ones(pset P);
public:
	void non_DC(pset lit01, pset lit10, pset _p);  //為了給paper.cpp使用


	pset_family P, Q;
	pset_family P_lit0_pos, P_lit1_pos; /* 每個uint初始全部為0，存一個位置，用SIZE來記錄有幾個位置*/
	pset_family ex_Q_non_DC; /* 指的是expand後的sf_Q_non_DC*/
	pset_family ex_q_non_DC; /* 指的是expand後的sf_Q_non_DC 的inverse (0=>1, 1=> 0)*/
	pset_family TABLE_PQ; //TABLE_PQ: P1, Q2
	
	MCDB() {};

	MCDB(pset_family f, pset_family g, int _num_mv_var, bool _phase_assignment) {
		setup(f, g, _num_mv_var, _phase_assignment);
	}

	void setup(pset_family f, pset_family g, int _num_mv_var, bool _phase_assignment) {
		phase_assignment = _phase_assignment;
		P = f;
		Q = g;
		num_mv_var = _num_mv_var;

		/* N個input var 2個bits = 2N*/
		sf_P_non_DC = sf_new(P->count, num_mv_var*2);
		sf_p_non_DC = sf_new(P->count, num_mv_var*2);
		sf_Q_non_DC = sf_new(Q->count, num_mv_var*2);
		sf_q_non_DC = sf_new(Q->count, num_mv_var*2);
		/* N個input var 2N個bits = 2N*N*/
		ex_Q_non_DC = sf_new(Q->count, num_mv_var*2 * num_mv_var);
		ex_q_non_DC = sf_new(Q->count, num_mv_var*2 * num_mv_var);
		/* onset 和offset配對的表 |on|*|off|*/
		TABLE_PQ = sf_new(P->count, Q->count);
		P_lit0_pos = sf_new(P->count, P->sf_size * 32);
		P_lit1_pos = sf_new(P->count, P->sf_size * 32);
		
		for (int i = 0; i < P->count; i++) {
			sf_addset(TABLE_PQ, set_new(Q->count));
			sf_addset(P_lit0_pos, set_new(P->sf_size * 32));
			sf_addset(P_lit1_pos, set_new(P->sf_size * 32));
		}
		int p_rel_max_size, q_rel_max_size;
		if (P->count > 1000)
			p_rel_max_size = 1000000;
		else
			p_rel_max_size = P->count * P->count;
		if (Q->count > 1000)
			q_rel_max_size = 1000000;
		else
			q_rel_max_size = Q->count * Q->count;

		P_relation.resize(p_rel_max_size);
		Q_relation.resize(q_rel_max_size);
		P_re_relation.resize(p_rel_max_size);
		Q_re_relation.resize(q_rel_max_size);
		cnt_R0 = 0, cnt_R1 = 0, cnt_R2 = 0, cnt_R3 = 0, cnt_R4 = 0, cnt_R5 = 0;
	}
	~MCDB() {
		sf_free(sf_P_non_DC);
		sf_free(sf_p_non_DC);
		sf_free(sf_Q_non_DC);
		sf_free(sf_q_non_DC);
		sf_free(ex_Q_non_DC);
		sf_free(ex_q_non_DC);
		sf_free(TABLE_PQ);
		sf_free(P_lit0_pos);
		sf_free(P_lit1_pos);
	};

	pset getPset(int set_index) {
		return GETSET(P, set_index);
	}

	pset getQset(int set_index) {
		return GETSET(Q, set_index);
	}

	void pre_run(unsigned int* first_word, unsigned int* last_word, unsigned int* first_bit, unsigned int* last_bit) {
		/* DC個數多到少排序
		   可以顯示P, Q各自的MIN,AVG,MAX,SD */
		EXECost_W(sort_P_Q(), T_SORT);
		T_lit01 = clock();
		/* 計算literal 0, literal 1 和DC個數*/
		cal_lit01_numDC();
		cal_lit0_lit1_pos();
		/* 將123 擴展成 123 123 123 123...*/
		expand_sf_Q(first_word, last_word, first_bit, last_bit);
		T_lit01 = clock() - T_lit01;
		EXECost_W(check_output_intersection(), T_Out_IS);
	}

	void check_Rule(int key) {
		fprintf(stderr, "MCQE key %d%d%d%d\n", (key & 1) > 0, (key & 2) > 0, (key & 4) > 0, (key & 8) > 0);
		if (key <= 0)
			return;
		EXECost_W(check_implies_relation(sf_P_non_DC, sf_p_non_DC, P_relation, P_re_relation), T_P_rela);
		EXECost_W(check_implies_relation(sf_Q_non_DC, sf_q_non_DC, Q_relation, Q_re_relation), T_Q_rela);

		if (key & 1)
			EXECost_W(check_Rule1(), T_R1);
		if (key & 2)
			EXECost_W(check_Rule2(), T_R2);
		if (key & 4)
			EXECost_W(check_Rule3(), T_R3);
		if (key & 8)
			EXECost_W(check_Rule4(), T_R4);
	}

	void pre_paper_run() {
		sort_P_Q();
	}

	/* DC個數多到少排序*/
	void sort_P_Q() {
		int P_sum_of_N_2 = 0, Q_sum_of_N_2 = 0, P_sum_of_N = 0, Q_sum_of_N = 0;
		int P_max = 0, P_min = 9999, Q_max = 0, Q_min = 9999;

		/* 計算P Q 的DC數量*/
		foreach_set(P, last, p) {
			PUTSIZE(p, set_count_inputs_ones(p));
		}
		foreach_set(Q, last, p) {
			PUTSIZE(p, set_count_inputs_ones(p));
		}
		/* 根據DC數量多到少排序*/
		SORT_1(P);
		SORT_1(Q);
	}
	/* 計算literal 0, literal 1 和DC個數*/
	void cal_lit01_numDC() {
		pset lit01 = set_new(num_mv_var * 2); // literal-0 對應到的mvcube
		pset lit10 = set_new(num_mv_var * 2); // literal-1 對應到的mvcube
		
		foreach_set(P, last, p) {
			/* N個input var 2個bits = 2N*/
			non_DC(lit01, lit10, p);
			sf_addset(sf_P_non_DC, lit01);
			sf_addset(sf_p_non_DC, lit10);
		}
		
		foreach_set(Q, last, p) {
			/* N個input var 2個bits = 2N*/
			non_DC(lit01, lit10, p);
			sf_addset(sf_Q_non_DC, lit01);
			sf_addset(sf_q_non_DC, lit10);
		}
		set_free(lit01);
		set_free(lit10);
	}
	/* 檢查P1 Q2 output是否相交*/
	void check_output_intersection() {
		pset p_talbe;
		int cnt_lit1;
		foreachi_set(P, index1, p) {
			p_talbe = GETSET(TABLE_PQ, index1);
			cnt_lit1 = cnt_non_DC_1(p);
			foreachi_set(Q, index2, p2) {
				if (OutputIntersection(p, p2)) {// && cnt_lit1 + cnt_non_DC_0(p2) < NUMINPUTS) {
					set_insert(p_talbe, index2);
					cnt_R0++;
				}
			}
		}
	}
	/*
		sf 之間的包含關係
		DC 多到少排序
		2 種包含關係 1. A 包含 B  2. A 包含 B'
	*/
	void check_implies_relation(pset_family sf, pset_family re_sf, std::vector<std::pair<int, int>>& relation, std::vector<std::pair<int, int>>& re_relation) {
		int cnt_relation = 0, cnt_re_relation = 0;
		pset p3;
		sf_active(sf);
		foreachi_set(sf, index1, p) {
			if ((index2 = index1 + 1) < sf->count) {
				for (p2 = p + sf->wsize, p3 = GETSET(re_sf, index2); index2 < sf->count; p2 += sf->wsize, p3 += re_sf->wsize, index2++) {
					// index1 包含 index2
					if (setp_implies(p2, p)) {
						relation[cnt_relation++] = std::make_pair(index1, index2);
					}
					// index2 包含 index1
					else if (setp_implies(p, p2)) {
						relation[cnt_relation++] = std::make_pair(index2, index1);
					}
					// index1 包含 re_index2
					else if (setp_implies(p3, p)) {
						re_relation[cnt_re_relation++] = std::make_pair(index1, index2);
					}
					// re_index2 包含 index1
					else if (setp_implies(p, p3)) {
						re_relation[cnt_re_relation++] = std::make_pair(index2, index1);
					}
				}
			}
		}
		relation.resize(cnt_relation);
		re_relation.resize(cnt_re_relation);
		relation.shrink_to_fit();
		re_relation.shrink_to_fit();
	}
	/* Rule1 同P不同Q (i, j) 包含(i, k)，故移除(i, j)*/
	void check_Rule1() {
		foreach_set(TABLE_PQ, last, p) {
			for (std::pair<int, int>& pair_q : Q_relation) {
				j = pair_q.first;
				k = pair_q.second;

				/*if (j > k)
					continue;*/

				if (!is_in_set(p, j) || !is_in_set(p, k))
					continue;

				set_remove(p, j);
				cnt_R1++;

				/* 如果P Q 未按照DC數量做多到少排序 則需額外處理*/
			}
		}
	}
	/* Rule2 不同P同Q(i, k) 包含(j, k)，故移除(i, k)*/
	void check_Rule2() {
		for (std::pair<int, int>& pair_p : P_relation) {
			p = GETSET(TABLE_PQ, pair_p.first);
			p2 = GETSET(TABLE_PQ, pair_p.second);
			for (k = 0; k < TABLE_PQ->sf_size; k++) {
				if (!is_in_set(p, k) || !is_in_set(p2, k))
					continue;
				set_remove(p, k);
				cnt_R2++;

				/* 如果P Q 未按照DC數量做多到少排序 則需額外處理*/
			}
		}
	}
	/* Rule3 不同U不同V(i, m) 包含(j, n)，故移除(i, m)*/
	void check_Rule3() {
		for (std::pair<int, int>& pair_p : P_relation) {
			p = GETSET(TABLE_PQ, pair_p.first);
			p2 = GETSET(TABLE_PQ, pair_p.second);
			for (std::pair<int, int>& pair_q : Q_relation) {
				m = pair_q.first;
				n = pair_q.second;
				if (!is_in_set(p, m) || !is_in_set(p2, n))
					continue;
				set_remove(p, m);
				cnt_R3++;

				/* 如果P Q 未按照DC數量做多到少排序 則需額外處理*/
			}
		}
	}
	/* 針對RE_RELATION
	   Rule4 不同U不同V(i, m) 包含(j, n)，故移除(i, m)*/
	void check_Rule4() {
		for (std::pair<int, int>& pair_p : P_re_relation) {
			//i = pair_p.first;
			//j = pair_p.second;
			p = GETSET(TABLE_PQ, pair_p.first);
			p2 = GETSET(TABLE_PQ, pair_p.second);
			for (std::pair<int, int>& pair_q : Q_re_relation) {
				m = pair_q.first;
				n = pair_q.second;
				if (!is_in_set(p, m) || !is_in_set(p2, n))
					continue;
				set_remove(p, m);
				cnt_R4++;

				/* 如果P Q 未按照DC數量做多到少排序 則需額外處理*/
			}
		}
	}

	void cal_lit0_lit1_pos() {
		pset p_lit0, p_lit1;
		int var, value, cnt_lit0, cnt_lit1;
		foreachi_set(P, index1, p) {
			p_lit0 = GETSET(P_lit0_pos, index1);
			p_lit1 = GETSET(P_lit1_pos, index1);
			cnt_lit0 = 1, cnt_lit1 = 1;
			for (var = 0; var < NUMINPUTS; var++) {
				value = GETINPUT(p, var);
				if (value == ZERO) {
					p_lit0[cnt_lit0++] = var;
				}
				else if (value == ONE) {
					p_lit1[cnt_lit1++] = var;
				}
			}
			PUTSIZE(p_lit0, cnt_lit0);
			PUTSIZE(p_lit1, cnt_lit1);
		}
	}

	void display() {
		printf("MCDB vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n");
		printf("P\n");
		foreachi_set(P, index1, p) {
			printf("%d:\n", index1);
			display_pset(p, cube.size);
			display_pset_01X(p, NUMINPUTS, NUMOUTPUTS);
			display_pset_01X(GETSET(sf_P_non_DC, index1), NUMINPUTS, NUMOUTPUTS);
		}
		printf("Q\n");
		foreachi_set(Q, index1, p) {
			printf("%d:\n", index1);
			display_pset(p, cube.size);
			display_pset_01X(p, NUMINPUTS, NUMOUTPUTS);
			display_pset_01X(GETSET(sf_Q_non_DC, index1), NUMINPUTS, NUMOUTPUTS);
			display_pset_01X(GETSET(sf_q_non_DC, index1), NUMINPUTS, NUMOUTPUTS);
		}
		printf("TABLE PQ\n");
		foreach_set(TABLE_PQ, last, p) {
			display_pset(p, TABLE_PQ->sf_size);
		}
		printf("MCDB ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
	}

	int cal_TABLE_PQ() {
		int cnt_ones = 0;
		foreach_set(TABLE_PQ, last, p) {
			cnt_ones += set_count_ones(p);
			//display_pset(p, TABLE_PQ->sf_size);
		}
		return cnt_ones;
	}

	void lit_pos_copy(pset r, pset a);

	void display_cost() {
		// P 與Q 各自之間的關係個數
		printf("P Q relation %d %d %d %d\n", P_relation.size(), Q_relation.size(), P_re_relation.size(), Q_re_relation.size());
		printf("MCQE Rule %d %d %d %d %d %d %d\n", TABLE_PQ->count * TABLE_PQ->sf_size, cnt_R0, cnt_R1, cnt_R2, cnt_R3, cnt_R4, cnt_R5);
		printf("NCQE Rule cost %.6lf  %.6lf  %.6lf  %.6lf  %.6lf  %.6lf  %.6lf  %.6lf  %.6lf\n",
			(double)T_SORT / (double)CLOCKS_PER_SEC, (double)T_lit01 / (double)CLOCKS_PER_SEC, 
			(double)T_Out_IS / (double)CLOCKS_PER_SEC, (double)T_P_rela / (double)CLOCKS_PER_SEC, (double)T_Q_rela / (double)CLOCKS_PER_SEC,
			(double)T_R1 / (double)CLOCKS_PER_SEC, (double)T_R2 / (double)CLOCKS_PER_SEC, (double)T_R3 / (double)CLOCKS_PER_SEC, (double)T_R4 / (double)CLOCKS_PER_SEC);
	}


	/* 以下是供原本paper作法使用*/
	bool check_output_intersection_pair(int index_P, int index_Q) {
		
		p = GETSET(P, index_P);
		p2 = GETSET(Q, index_Q);
		return OutputIntersection(p, p2);
	}
	/* 計算literal 0, literal 1 和DC個數*/
	/*
		pset lit01 = set_new(num_mv_var * 2); // literal-0 
		pset lit10 = set_new(num_mv_var * 2); // literal-1
	*/
	void cal_lit01_numDC_pair(int index_Q, pset lit01, pset lit10) {
		p = GETSET(Q, index_Q);
		
		non_DC(lit01, lit10, p);
	}
	void cal_lit0_lit1_pos_pair(int index_P, pset p_lit0, pset p_lit1) {
		p = GETSET(P, index_P);
		int var, value, cnt_lit0, cnt_lit1;

		cnt_lit0 = 1, cnt_lit1 = 1;
		for (var = 0; var < NUMINPUTS; var++) {
			value = GETINPUT(p, var);
			if (value == ZERO) {
				p_lit0[cnt_lit0++] = var;
			}
			else if (value == ONE) {
				p_lit1[cnt_lit1++] = var;
			}
		}
		PUTSIZE(p_lit0, cnt_lit0);
		PUTSIZE(p_lit1, cnt_lit1);
		
	}
	void run_RULE_TEST();
};

class paper_MCDB {
	int mv_num, mv_size;
	std::vector<pset_family> DB;

public:
	int cnt_skip;

	paper_MCDB(int _mv_num, int _mv_size) : mv_num(_mv_num), mv_size(_mv_size) {
		DB.resize(mv_size);

		for (int i = 0; i < mv_size; i++) {
			DB[i] = sf_new(1, mv_size);
		}
		cnt_skip = 0;
	}

	~paper_MCDB() {
		// 釋放動態配置的記憶體
		for (int i = 0; i < mv_num; i++) {
			sf_free(DB[i]);
		}
	}

	/* is Pset absorbed*/
	bool isPAbsorbed(pset _p) {
		pset p, last;
		pset_family temp_DB;

		//printf("%d %d\n", mv_num - 1, SIZE(_p));
		for (int i = mv_num-1; i >= (int)SIZE(_p); i--) {
			temp_DB = DB[i];
			foreach_set(temp_DB, last, p) {
				if (setp_implies(p, _p)) {
					sf_inactive(temp_DB);
					return true;
				}
				else if (setp_implies(_p, p)) {
					RESET(p, ACTIVE);
				}
			}
			sf_inactive(temp_DB);
		}
		temp_DB = DB[SIZE(_p)];
		SET(_p, ACTIVE);
		sf_addset(temp_DB, _p);//add pset into mcdb
		return false;
	}
	/* is Pset family absorbed*/
	bool isSFAbsorbed(pset_family sf) {
		pset p, last;
		bool isALLAbsorbed = true;


		foreach_set(sf, last, p) {
			if (!isPAbsorbed(p)) {
				isALLAbsorbed = false;
				//break;
			}
		}

		return isALLAbsorbed;
	}

};
#endif //MCDB_H



