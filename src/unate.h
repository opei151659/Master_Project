/*
	unate 特徵值的計算與比對，部分比對會在mainCplusplus中執行
*/

#ifndef UNATE_H
#define UNATE_H
#include <iostream>
#include <vector>
#include <time.h>
#include "tool.h"
extern "C" {
#include "espresso.h"   /* ESPRESSO.lib*/
}

#define output_non_support 0xFFFFFFFF
#define output_positive_unate 0x55555555
#define output_negative_unate 0xAAAAAAAA
#define output_binate 0

class Unate {
private:
	int num_inputs;
	int num_outputs;
	int num_mv_var;
	int mv_var_size;
	pset_family unate_table_A;
	pset_family unate_table_a; //A' inverse
	pset_family unate_table_B;
	pset_family unate_table_b; //B' inverse
public:
	Unate(int _num_in, int _num_out) {
		num_inputs = _num_in;
		num_outputs = _num_out;
		num_mv_var = num_inputs;
		mv_var_size = num_outputs * 2;
		unate_table_A = sf_new(num_mv_var, mv_var_size);
		unate_table_a = sf_new(num_mv_var, mv_var_size);
		unate_table_B = sf_new(num_mv_var, mv_var_size);
		unate_table_b = sf_new(num_mv_var, mv_var_size);
		for (int i = 0; i < num_inputs; i++) {
			sf_addset(unate_table_A, set_full(mv_var_size));
			sf_addset(unate_table_a, set_full(mv_var_size));
			sf_addset(unate_table_B, set_full(mv_var_size));
			sf_addset(unate_table_b, set_full(mv_var_size));
		}
	}
	~Unate() {
		sf_free(unate_table_A);
		sf_free(unate_table_a);
		sf_free(unate_table_B);
		sf_free(unate_table_b);
	}
	void check(pPLA PLA1, pPLA PLA2) {
		cal_unate(PLA1, unate_table_A, unate_table_a);
		cal_unate(PLA2, unate_table_B, unate_table_b);
	}
	/*
		檢查unate
		unate_table 是一個2維表 var_num * output_num   (每個output都用兩個bit表示)
		unate_table 初始化為全部為11
		ex:       Out1 Out2 Out3             00  binate
		    var0   11   11   11              01  negative unate
			var1   11   11   11              10  positive unate
			var2   11   11   11              11  non support

		使用onset(f^on) 和 offset(g^off)
		f^on_i(p_i), g^off_j(q_j)
		step1: 檢查p_i和q_j 的output是否相交，若否則跳過
		step2: 檢查p_i和q_j 的input距離是否等於1，若否則則跳過
		step3: 檢查p_i 在input距離為1位置的var是1還是0，1表示是Positive Unate (移除 NU bit)，0表示是Negative Unate (移除 PU bit)
		step4: 根據step3 的結果，將結果紀錄在input距離為1位置的var對應中有相交到的output上 (有兩張表PU和NU各一張)

	*/
	void cal_unate(pPLA PLA, pset_family unate_table_A, pset_family unate_table_a) {
		/* 每個var 對應到2 * NUMOUTPUTS 個bits  每兩個bits表示該var在對應output的unate關係*/
		/* 00: binate, 01: positive_unate, 10: negative_unage, 11: non-support */
		pset p1, p2, last1, last2, unate_A, unate_a;
		temp_psets t_ps(6);
		t_ps.psets[0] = set_new(NUMOUTPUTS * 2); //non-support 11
		t_ps.psets[1] = set_new(NUMOUTPUTS * 2); //binate 00
		t_ps.psets[2] = set_new(NUMOUTPUTS * 2); //positive unate 01
		t_ps.psets[3] = set_new(NUMOUTPUTS * 2); //negative unate 10
		t_ps.psets[4] = set_new(NUMOUTPUTS * 2); //output mask
		t_ps.psets[5] = set_new(NUMOUTPUTS * 2); //output temp mask
		pset NS = t_ps.psets[0], BN = t_ps.psets[1], PU = t_ps.psets[2], NU = t_ps.psets[3];
		pset output_mask = t_ps.psets[4];
		pset temp_output_mask = t_ps.psets[5];
		for (int word = 1; word < ceil((NUMOUTPUTS * 2.0) / 32.0); word++) {
			NS[word] = (unsigned int)output_non_support;
			PU[word] = (unsigned int)output_positive_unate;
			NU[word] = (unsigned int)output_negative_unate;
			BN[word] = (unsigned int)output_binate;
		}
		bool isPU, isNU, isOne;
		
		int var;
		foreach_set(PLA->F, last1, p1) {
			foreach_set(PLA->R, last2, p2) {
				// output 是否相交
				if (!OutputIntersection(p1, p2))
					continue;
				var = cdist1(p1, p2);
				if (var == -1) continue;
				if (GETINPUT(p1, var) == ONE)
					isOne = true; // 1 => 0
				else
					isOne = false; // 0 => 1

				set_fill(output_mask, NUMOUTPUTS * 2);
				// 確定哪幾個output相交
				for (int i = 0; i < NUMOUTPUTS; i++) {
					if (GETOUTPUT(p1, i) == 1 && GETOUTPUT(p2, i) == 1) {
						set_remove(output_mask, i * 2);
						set_remove(output_mask, i * 2 + 1);
					}
				}
				unate_A = GETSET(unate_table_A, var);
				unate_a = GETSET(unate_table_a, var);
				if (isOne) {  // PU
					set_and(unate_A, unate_A, set_or(temp_output_mask, PU, output_mask)); //remove NU bit
					set_and(unate_a, unate_a, set_or(temp_output_mask, NU, output_mask)); //remove PU bit
				}
				else { // NU
					set_and(unate_A, unate_A, set_or(temp_output_mask, NU, output_mask)); //remove PU bit
					set_and(unate_a, unate_a, set_or(temp_output_mask, PU, output_mask)); //remove NU bit
				}
			}
		}
	}
	/*
		使用unate
		定義: 當兩個var是相同的，則其var對應到各個output的PUNU關係必須相同
		      當var = var'時，則其var對應到各個output的PUNU關係必須相反
	*/
	/* notation1 適用表示法abcda'b'c'd'*/
	void var_match_notation1(pset p, int _var_size, bool INPA) {
		for (int varA = 0; varA < num_inputs; varA++) {
			for (int varB = 0; varB < num_inputs; varB++) {
				// 移除不符合的
				if (!setp_equal(GETSET(unate_table_A, varA), GETSET(unate_table_B, varB))) {
					set_remove(p, varA * _var_size + varB);
				}
				// 移除不符合的
				if (INPA && !setp_equal(GETSET(unate_table_A, varA), GETSET(unate_table_b, varB))) {
					set_remove(p, varA * _var_size + varB + num_inputs);
				}
			}
		}
	}
	/* notation1 適用表示法aa'bb'cc'dd'*/
	void var_match_notation2(pset p, int _var_size, bool INPA) {
		for (int varA = 0; varA < num_inputs; varA++) {
			for (int varB = 0; varB < num_inputs; varB++) {
				// 移除不符合的
				if (!setp_equal(GETSET(unate_table_A, varA), GETSET(unate_table_B, varB))) {
					set_remove(p, varA * _var_size + (varB*2+1));
				}
				// 移除不符合的
				if (INPA && !setp_equal(GETSET(unate_table_A, varA), GETSET(unate_table_b, varB))) {
					set_remove(p, varA * _var_size + (varB*2));
				}
			}
		}
	}

	void cnt_output_unateness(pset_family sf, std::vector<std::vector<int>>& out) {
		pset p, last;
		int out_var;
		foreach_set(sf, last, p) {
			for (out_var = 0; out_var < num_outputs; out_var++) {
				out[out_var][GETINPUT(p, out_var)]++;
			}
		}
	}

	void output_match(pset p_out, bool INPA) {
		// 4 個位置對應 00 01 10 11 在各個output上的個數
		std::vector<std::vector<int>> outA(NUMOUTPUTS, std::vector<int>(4, 0));
		std::vector<std::vector<int>> outB(NUMOUTPUTS, std::vector<int>(4, 0));
		pset p1, p2, last1, last2;
		int out_var;
		
		cnt_output_unateness(unate_table_A, outA);
		cnt_output_unateness(unate_table_B, outB);

		/*for (std::vector<int>& out_var : outA) {
			for (int& unate : out_var) {
				printf("%d ", unate);
			}printf("\n");
		}*/
		pset p_in_y = set_new(NUMOUTPUTS * NUMOUTPUTS);
		for (int i = 0; i < NUMOUTPUTS; i++) {
			for (int j = 0; j < NUMOUTPUTS; j++) {
				// ns = ns & binate = binate
				if (outA[i][0] == outB[j][0] && outA[i][3] == outB[j][3]) {
					if (outA[i][1] == outB[j][1] && outA[i][2] == outB[j][2]) {
						set_insert(p_in_y, i * NUMOUTPUTS + j);
					}
					else if (INPA && outA[i][1] + outA[i][2] == outB[j][1] + outB[j][2]) {
						set_insert(p_in_y, i * NUMOUTPUTS + j);
					}
				}
			}
		}

		set_and(p_out, p_out, p_in_y);
		set_free(p_in_y);
	}

	void display() {
		printf("unate table A: \n");
		display_sf(unate_table_A);
		printf("unate table a: \n");
		display_sf(unate_table_a);
		printf("unate table B: \n");
		display_sf(unate_table_B);
		printf("unate table b: \n");
		display_sf(unate_table_b);
	}
};


#endif //UNATE_H