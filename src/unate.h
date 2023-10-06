/*
	unate �S�x�Ȫ��p��P���A�������|�bmainCplusplus������
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
		�ˬdunate
		unate_table �O�@��2���� var_num * output_num   (�C��output���Ψ��bit���)
		unate_table ��l�Ƭ�������11
		ex:       Out1 Out2 Out3             00  binate
		    var0   11   11   11              01  negative unate
			var1   11   11   11              10  positive unate
			var2   11   11   11              11  non support

		�ϥ�onset(f^on) �M offset(g^off)
		f^on_i(p_i), g^off_j(q_j)
		step1: �ˬdp_i�Mq_j ��output�O�_�ۥ�A�Y�_�h���L
		step2: �ˬdp_i�Mq_j ��input�Z���O�_����1�A�Y�_�h�h���L
		step3: �ˬdp_i �binput�Z����1��m��var�O1�٬O0�A1��ܬOPositive Unate (���� NU bit)�A0��ܬONegative Unate (���� PU bit)
		step4: �ھ�step3 �����G�A�N���G�����binput�Z����1��m��var���������ۥ�쪺output�W (����i��PU�MNU�U�@�i)

	*/
	void cal_unate(pPLA PLA, pset_family unate_table_A, pset_family unate_table_a) {
		/* �C��var ������2 * NUMOUTPUTS ��bits  �C���bits��ܸ�var�b����output��unate���Y*/
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
				// output �O�_�ۥ�
				if (!OutputIntersection(p1, p2))
					continue;
				var = cdist1(p1, p2);
				if (var == -1) continue;
				if (GETINPUT(p1, var) == ONE)
					isOne = true; // 1 => 0
				else
					isOne = false; // 0 => 1

				set_fill(output_mask, NUMOUTPUTS * 2);
				// �T�w���X��output�ۥ�
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
		�ϥ�unate
		�w�q: ����var�O�ۦP���A�h��var������U��output��PUNU���Y�����ۦP
		      ��var = var'�ɡA�h��var������U��output��PUNU���Y�����ۤ�
	*/
	/* notation1 �A�Ϊ�ܪkabcda'b'c'd'*/
	void var_match_notation1(pset p, int _var_size, bool INPA) {
		for (int varA = 0; varA < num_inputs; varA++) {
			for (int varB = 0; varB < num_inputs; varB++) {
				// �������ŦX��
				if (!setp_equal(GETSET(unate_table_A, varA), GETSET(unate_table_B, varB))) {
					set_remove(p, varA * _var_size + varB);
				}
				// �������ŦX��
				if (INPA && !setp_equal(GETSET(unate_table_A, varA), GETSET(unate_table_b, varB))) {
					set_remove(p, varA * _var_size + varB + num_inputs);
				}
			}
		}
	}
	/* notation1 �A�Ϊ�ܪkaa'bb'cc'dd'*/
	void var_match_notation2(pset p, int _var_size, bool INPA) {
		for (int varA = 0; varA < num_inputs; varA++) {
			for (int varB = 0; varB < num_inputs; varB++) {
				// �������ŦX��
				if (!setp_equal(GETSET(unate_table_A, varA), GETSET(unate_table_B, varB))) {
					set_remove(p, varA * _var_size + (varB*2+1));
				}
				// �������ŦX��
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
		// 4 �Ӧ�m���� 00 01 10 11 �b�U��output�W���Ӽ�
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