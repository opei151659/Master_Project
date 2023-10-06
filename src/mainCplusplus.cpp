/*
	主要header的CPP檔
	包含輸入輸出特徵值的計算功能與比對過程
	包含MvCube的合法性檢查 分為notation1 與notation2 
		notation1 是abcda'b'c'd' 有實作合法性檢查 但未實作PLA轉乘notation1 的過程 (不可直接使用)
		notation2 是aa'bb'cc'dd' (實驗時使用這個) 理論上比notation1快，但要輸入個數很大效果才明顯

*/

#include "mainCplusplus.h"


void check_unate(pset p, pPLA PLA1, pPLA PLA2) {
	Unate char_unate(NUMINPUTS, NUMOUTPUTS);
	char_unate.check(PLA1, PLA2);
	//char_unate.display();
	if (notation == 1) char_unate.var_match_notation1(p, mvcube.mv_var_size, mvcube.is_input_phase_assigment);
	else if (notation == 2) char_unate.var_match_notation2(p, mvcube.mv_var_size, mvcube.is_input_phase_assigment);
	if (print_debug) { printf("unate:\n"); display_pset_eachvar(p, mvcube.mv_var_size, mvcube.num_mv_var); }
}

/* option 1: tbb, 2: openmp, 3: cpp*/
void check_ENE(pset p, pPLA PLA1, pPLA PLA2, int option) {
	pset_family ENE1 = sf_new(NUMINPUTS, NUMINPUTS * 2);
	pset_family ENE2 = sf_new(NUMINPUTS, NUMINPUTS * 2);
	ene(ENE1, PLA1->F, PLA1->R, option);
	//display_sf(ENE1);
	ene(ENE2, PLA2->F, PLA2->R, option);
	//display_sf(ENE2);
	vector<group_vars> group_E_1, group_NE_1, group_ENE_1, group_E_2, group_NE_2, group_ENE_2;
	vector<int> group_X_1, group_X_2;
	grouping_ENE(ENE1, group_E_1, group_NE_1, group_ENE_1, group_X_1);
	grouping_ENE(ENE2, group_E_2, group_NE_2, group_ENE_2, group_X_2);
	matching_ENE(group_E_1, group_NE_1, group_ENE_1, group_X_1, group_E_2, group_NE_2, group_ENE_2, group_X_2, p);
	if (print_debug) { printf("ENE:\n"); display_pset_eachvar(p, mvcube.mv_var_size, mvcube.num_mv_var); }
}

/* display KSO mode0
	for (int i = 0; i < KSIG_output1.size(); i++) {
		for (int j = 0; j < KSIG_output1[i].distx.size(); j++) {
			cout << KSIG_output1[i].distx[j];
			printf(" ");
		} printf("\n");
	}
*/
/* display KSI mode1
	for (int i = 0; i < KSI1.size(); i++) { //[output][input].dist[dist]
		for (int j = 0; j < KSI1[i].size(); j++) {
			for (int k = 0; k < KSI1[i][j].distx.size(); k++) {
				cout << KSI1[i][j].distx[k];
				printf(" ");
			} printf("\n");
		}
		printf("\n");
	}
*/
/* defualt: distN = N, option = 2 */
void check_KSIG(pset p, pPLA PLA1, pPLA PLA2, int distN, int mode) {
	if (mode == 0) { //KSO
		TABLE KSO1(NUMOUTPUTS, ksig_value(NUMINPUTS));
		TABLE KSO2(NUMOUTPUTS, ksig_value(NUMINPUTS));

		vector<TABLE> KSI_empty;

		cal_ksignature(PLA1, KSO1, KSI_empty, 0, 0); /* Mode 0: KSO, 1: KSI, 2: KSO+KSI*/
		cal_ksignature(PLA2, KSO2, KSI_empty, 0, 0); /* Mode 0: KSO, 1: KSI, 2: KSO+KSI*/

		matching_KSIG_output(p, KSO1, KSO2);
		if (print_debug) { printf("KSIG output:\n"); display_pset_eachvar(p, mvcube.num_output, mvcube.num_output); }
	}
	else if (mode == 1) { //KSI
		//KSI[output][input].dist[dist]
		vector<TABLE>  KSI1(NUMOUTPUTS, TABLE(NUMINPUTS, ksig_value(distN)));
		vector<TABLE>  KSI2(NUMOUTPUTS, TABLE(NUMINPUTS, ksig_value(distN)));

		TABLE KSO_empty;

		cal_ksignature(PLA1, KSO_empty, KSI1, distN, 1); /* Mode 0: KSO, 1: KSI, 2: KSO+KSI*/
		cal_ksignature(PLA2, KSO_empty, KSI2, distN, 1); /* Mode 0: KSO, 1: KSI, 2: KSO+KSI*/
		
		matching_KSIG(p, KSI1, KSI2);
		if (print_debug) { printf("KSIG:\n"); display_pset_eachvar(p, mvcube.mv_var_size, mvcube.num_mv_var); }
	}
	else if (mode == 2) {//KSI in output
		vector<TABLE>  KSI1(NUMOUTPUTS, TABLE(NUMINPUTS, ksig_value(distN)));
		vector<TABLE>  KSI2(NUMOUTPUTS, TABLE(NUMINPUTS, ksig_value(distN)));

		TABLE KSO_empty;

		cal_ksignature(PLA1, KSO_empty, KSI1, distN, 1); /* Mode 0: KSO, 1: KSI, 2: KSO+KSI*/
		cal_ksignature(PLA2, KSO_empty, KSI2, distN, 1); /* Mode 0: KSO, 1: KSI, 2: KSO+KSI*/

		matching_KSIG_output_with_KSIG_input(p, KSI1, KSI2);
		if (print_debug) { printf("KSIG:\n"); display_pset_eachvar(p, mvcube.num_output, mvcube.num_output); }
	}
	else if (mode == 3) {//KSO + KSI
		TABLE KSO1(NUMOUTPUTS, ksig_value(NUMINPUTS));
		TABLE KSO2(NUMOUTPUTS, ksig_value(NUMINPUTS));
		vector<TABLE>  KSI1(NUMOUTPUTS, TABLE(NUMINPUTS, ksig_value(distN)));
		vector<TABLE>  KSI2(NUMOUTPUTS, TABLE(NUMINPUTS, ksig_value(distN)));


		cal_ksignature(PLA1, KSO1, KSI1, distN, 2); /* Mode 0: KSO, 1: KSI, 2: KSO+KSI*/
		cal_ksignature(PLA2, KSO2, KSI2, distN, 2); /* Mode 0: KSO, 1: KSI, 2: KSO+KSI*/
		
		matching_KSIG_output(p, KSO1, KSO2);
		matching_KSIG_output_with_KSIG_input(p, KSI1, KSI2);

		if (print_debug) { printf("KSIG:\n"); display_pset_eachvar(p, mvcube.num_output, mvcube.num_output); }
	}
}



void check_cofactor(pset p, pPLA PLA1, pPLA PLA2) {
	pset cofactor_match = set_new(mvcube.mv_size);
	Cofactor Cof(PLA1->F);
	Cofactor Cof2(PLA2->F);
	Cof.cal_cofactor();
	//Cof.display();
	Cof2.cal_cofactor();
	//Cof.display();
	Cof.matching_cofactor(cofactor_match, Cof2.var_cof, mvcube.is_input_phase_assigment, notation);
	set_and(p, p, cofactor_match);
	if (print_debug) { printf("Cofactor:\n"); display_pset_eachvar(p, mvcube.mv_var_size, mvcube.num_mv_var); }
	set_free(cofactor_match);
}





void check_Unate_output(pset p_out, pPLA PLA1, pPLA PLA2) {
	Unate char_unate(NUMINPUTS, NUMOUTPUTS);
	char_unate.check(PLA1, PLA2);
	//char_unate.display();

	char_unate.output_match(p_out, mvcube.is_input_phase_assigment);
}

void check_character_value(int key, pset totality, pPLA PLA1, pPLA PLA2) {
	fprintf(stderr, "CV key %d%d%d%d\n", (key & 1) > 0, (key & 2) > 0, (key & 4) > 0, (key & 8) > 0);
	if (key & 1) {
		EXECost(check_unate(totality, PLA1, PLA2), "PUNU check");
		printf("cnt_unate %d\n", set_count_ones(totality));
	}
	if (key & 2) {
		EXECost(check_ENE(totality, PLA1, PLA2, 1), "ENE check"); /* option 1: tbb, 2: openmp, 3: cpp*/
		printf("cnt_ENE %d\n", set_count_ones(totality));
	}
	if (key & 4) {

		EXECost(check_KSIG(totality, PLA1, PLA2, NUMINPUTS, 1), "KSIG check"); /*default: distN = N, mode = 1 KSI*/
		printf("cnt_KSI_n %d\n", set_count_ones(totality));
	}
		
	if (key & 8)
		EXECost(check_cofactor(totality, PLA1, PLA2), "Cofactor check");
}

#define CV_TEST(key, cost, cnt) {\
		set_copy(temp_totality, totality);\
		EXECost_W(check_character_value(key, temp_totality, PLA1, PLA2), cost);\
		cnt =  set_count_ones(temp_totality);\
	}
void run_InputSIG_TEST(pset totality, pPLA PLA1, pPLA PLA2) {
	pset temp_totality = set_new(mvcube.mv_size);
	int init_cnt = set_count_ones(totality);
	clock_t A = 0, B = A, C = A, D = A, AB = A, AC = A, AD = A, BC = A, BD = A, CD = A, ABC = A, ABD = A, ACD = A, BCD = A, ABCD = A;
	int cntA = 0, cntB = 0, cntC = 0, cntD = 0, cntAB = 0, cntAC = 0, cntAD = 0, cntBC = 0, cntBD = 0, cntCD = 0, cntABC = 0, cntABD = 0, cntACD = 0, cntBCD = 0, cntABCD = 0;

	CV_TEST(0b0001, A, cntA);
	CV_TEST(0b0010, B, cntB);
	CV_TEST(0b0100, C, cntC);
	CV_TEST(0b0111, ABC, cntABC);
	//CV_TEST(0b1000, D, cntD);  //fprintf(stderr, "finish D %.6lf %.6lf\n", cntD, (double)D / (double)CLOCKS_PER_SEC);
	/*CV_TEST(0b1000, D, cntD);  fprintf(stderr, "finish D %.6lf %.6lf\n", cntD, (double)D / (double)CLOCKS_PER_SEC);
	CV_TEST(0b0011, AB, cntAB);
	CV_TEST(0b0101, AC, cntAC);
	CV_TEST(0b1001, AD, cntAD); fprintf(stderr, "finish AD %.6lf %.6lf\n", cntAD, (double)AD / (double)CLOCKS_PER_SEC);
	CV_TEST(0b0110, BC, cntBC);
	CV_TEST(0b1010, BD, cntBD); fprintf(stderr, "finish BD %.6lf %.6lf\n", cntBD, (double)BD / (double)CLOCKS_PER_SEC);
	CV_TEST(0b1100, CD, cntCD); fprintf(stderr, "finish CD %.6lf %.6lf\n", cntCD, (double)CD / (double)CLOCKS_PER_SEC);
	CV_TEST(0b0111, ABC, cntABC);
	CV_TEST(0b1011, ABD, cntABD); fprintf(stderr, "finish ABD %.6lf %.6lf\n", cntABD, (double)ABD / (double)CLOCKS_PER_SEC);
	CV_TEST(0b1101, ACD, cntACD); fprintf(stderr, "finish ACD %.6lf %.6lf\n", cntD, (double)ACD / (double)CLOCKS_PER_SEC);
	CV_TEST(0b1110, BCD, cntBCD); fprintf(stderr, "finish BCD %.6lf %.6lf\n", cntD, (double)BCD / (double)CLOCKS_PER_SEC);*/
	//CV_TEST(0b1111, ABCD, cntABCD); //fprintf(stderr, "finish ABCD %.6lf %.6lf\n", cntD, (double)ABCD / (double)CLOCKS_PER_SEC);

	printf("%s finish CV TEST\n", PLA1->filename);
	printf("CV TEST : A, B, C, D, AB, AC, AD, BC, BD, CD, ABC, ABD, ACD, BCD, ABCD\n");
	printf("CV TEST cnt : %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", cntA, cntB, cntC, cntD, cntAB, cntAC, cntAD, cntBC, cntBD, cntCD, cntABC, cntABD, cntACD, cntBCD, cntABCD);
	printf("CV TEST COST: %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf\n",
		(double)A / (double)CLOCKS_PER_SEC, (double)B / (double)CLOCKS_PER_SEC, (double)C / (double)CLOCKS_PER_SEC, (double)D / (double)CLOCKS_PER_SEC,
		(double)AB / (double)CLOCKS_PER_SEC, (double)AC / (double)CLOCKS_PER_SEC, (double)AD / (double)CLOCKS_PER_SEC, (double)BC / (double)CLOCKS_PER_SEC,
		(double)BD / (double)CLOCKS_PER_SEC, (double)CD / (double)CLOCKS_PER_SEC, (double)ABC / (double)CLOCKS_PER_SEC, (double)ABD / (double)CLOCKS_PER_SEC,
		(double)ACD / (double)CLOCKS_PER_SEC, (double)BCD / (double)CLOCKS_PER_SEC, (double)ABCD / (double)CLOCKS_PER_SEC);


	set_free(temp_totality);
}

/*
	ENE 分群
	根據建立好的ENE關係表將各個var之間的E / NE /ENE /X 關係做分群
	1. E  分群: 2個一群，最多不超過C N 取N-1 群
	2. NE 分群: 1~N個一群，最多不超過N 群
	3. ENE分群: 1~N個一群，最多不超過N 群
	4. X  分群: 1~N個一群，最多1 群
	step 1: 建立所有var之間的ENE關係表 var*var (2個bit分別代表E /NE 關係)  (ene())
	step 2: 根據ENE關係表做分群(為每個var儲存其對應到的E/ NE/ ENE群的id) (vars_group.grouping)
	step 3: 根據step2 每個var對應的E/ NE/ ENE群id 進行整理
	step 4: 將NE群最大化(合併NE-E-NE(NE間若是有E的關係則可以合併))
*/
void grouping_ENE(pset_family sf, vector<group_vars>& group_E, vector<group_vars>& group_NE, vector<group_vars>& group_ENE, vector<int>& group_X) {
	// E 的分群 最多有C N 取N-1 種
	group_E.resize(NUMINPUTS * (NUMINPUTS - 1) / 2);
	for (int i = 0; i < NUMINPUTS * (NUMINPUTS - 1) / 2; i++) {
		group_E[i].vars.resize(2); // 每一個分群最多只會有兩個var
	}
	// NE 的分群 最多NUMINPUTS 種
	group_NE.resize(NUMINPUTS);
	E_NE_group vars_group;
	int i, var;
	pset p;

	//對variable分群 E/ NE/ E+NE/ X
	foreachi_set(sf, i, p) {
		for (int j = i + 1; j < NUMINPUTS; j++) {
			vars_group.grouping(i, j, GETINPUT(p, j), group_E);
		}
	}

	//vars_group.display();
	//整理分群的結果
	//printf("(%d, %d)\n", group_E.size(), group_NE.size());
	vars_group.classify(group_E, group_NE, group_ENE, group_X);

	//printf("(%d, %d)\n", group_E.size(), group_NE.size());
	//將NEgroup最大化(合併NE-E-NE)
	vars_group.NE_E_NE_grouping(sf, group_E, group_NE);


	sort(group_E.begin(), group_E.end(), greater());
	sort(group_NE.begin(), group_NE.end(), greater());
	if (print_debug) {
		for (int i = 0; i < group_E.size(); i++) {
			printf(" E(%3d): ", i);
			group_E[i].display();
		}
		for (int i = 0; i < group_NE.size(); i++) {
			printf("NE(%3d): ", i);
			group_NE[i].display();
		}
		for (int i = 0; i < group_ENE.size(); i++) {
			printf("ENE(%3d): ", i);
			group_ENE[i].display();
		}
		printf("X: (");
		for (int i = 0; i < group_X.size(); i++) {
			printf(" %d", group_X[i]);
		}
		printf(")\n");
	}
}

/*
	ENE 比對
	輸入E /NE /ENE /X 分群的結果 (有小到大sort過)
	分兩部分處理
	1. 不考慮input phase assigment
		a.  E <=> E  (這裡的E皆沒有NE的關係)
		    (1). E 沒有 equivalence relation(不能做額外的化簡)
		    (2). 每個PLA1的E group 中的var 需要對應到所有PLA2的E group 中的var
		b.  NE <=> NE (NE需經過E的最大化)
		    (1). 每個PLA1的NE 會對應到所有group大小相同的PLA2的NE
			(2). NE 具有equivalence relation，故每個PLA1的NE group 中的var 僅需要對應到PLA2的NE group 中的一個var
		c.  ENE <=> ENE
		    (1). ENE 具有equivalence relation，故每個PLA1的ENE group 中的var 僅需要對應到PLA2的ENE group 中的一個var
		d.  X <=> X
			(1). 每個PLA1的X group 中的var 需要對應到所有PLA2的X group 中的var
	2. 考慮input phase assigment
	    a. NE <=> E
			(1). 若PLA1的NE group的大小為2，則每個PLA1的NE group 中的var 需要對應到PLA2的E group 中的一個var (因E可以透過phase assigment變成NE))
		b. E <=> NE
		    (1). 同NE <=> E

*/
#define ENE_match_notation1(p, pos, num_mv_var, phase_assigment) {set_insert(p, pos); if (phase_assigment) set_insert(p, pos + num_mv_var);}
#define ENE_match_notation2(p, pos, phase_assigment) {set_insert(p, pos*2+1); if (phase_assigment) set_insert(p, pos*2);}
void matching_ENE(vector<group_vars>& group_E_1, vector<group_vars>& group_NE_1, vector<group_vars>& group_ENE_1, vector<int>& group_X_1
	, vector<group_vars>& group_E_2, vector<group_vars>& group_NE_2, vector<group_vars>& group_ENE_2, vector<int>& group_X_2, pset totality) {
	// 不考慮invers
	temp_psets t_ps(1), t_ps2(NUMINPUTS);
	t_ps.psets[0] = set_new(mvcube.mv_size);
	pset ENE_mapping = t_ps.psets[0]; // ENE的var對應結果 初始全為0
	for (int i = 0; i < NUMINPUTS; i++)
		t_ps2.psets[i] = set_new(mvcube.mv_var_size);

	int pre_NE_size = 0, cur_g2 = 0, var2;
	// ***************** E <=> E no equivalence relation *****************
	/*for (int var = 0; var < 2; var++) {
		set_clear(t_ps2.psets[var], mvcube.mv_var_size);
	}*/
	set_clear(t_ps2.psets[0], mvcube.mv_var_size);
	//PLA2 中所有具有E關係的var
	for (int g2 = 0; g2 < group_E_2.size(); g2++) {
		for (int var = 0; var < group_E_2[g2].vars.size(); var++) {
			if (notation == 1) { ENE_match_notation1(t_ps2.psets[0], group_E_2[g2].vars[var], NUMINPUTS, mvcube.is_input_phase_assigment); }
			else if (notation == 2) { ENE_match_notation2(t_ps2.psets[0], group_E_2[g2].vars[var], mvcube.is_input_phase_assigment); }
		}
	}
	//printf("E<=>E: "); display_pset(t_ps2.psets[0], mvcube.mv_var_size);
	//PLA1 中所有具有E關係的都有可能對應到PLA2 所有具有E關係的可能
	for (int g = 0; g < group_E_1.size(); g++) {
		for (int var = 0; var < group_E_1[g].vars.size(); var++) {
			var2 = group_E_1[g].vars[var];
			or_exact_MvVar(ENE_mapping, t_ps2.psets[0], mvcube.mv_var_word_size, mvcube.first_word[var2], mvcube.last_word[var2], mvcube.first_bit[var2], mvcube.last_bit[var2]);
		}
	}

	// ***************** NE <=> NE *****************
	// 計算同個group size會對應到的1to1可能
	pre_NE_size = 0; cur_g2 = 0;
	for (int g = 0; g < group_NE_1.size(); g++) {
		// 同數量的group_NE 會對應到相同的var可能
		if (pre_NE_size != group_NE_1[g].vars.size()) {
			pre_NE_size = group_NE_1[g].vars.size();
			for (int var = 0; var < pre_NE_size; var++) {
				set_clear(t_ps2.psets[var], mvcube.mv_var_size);
			}
			for (int g2 = cur_g2; g2 < group_NE_2.size(); g2++) {
				if (group_NE_1[g].vars.size() < group_NE_2[g2].vars.size())
					continue;
				else if (group_NE_1[g].vars.size() == group_NE_2[g2].vars.size()) {
					// 1 to 1
					for (int var = 0; var < group_NE_1[g].vars.size(); var++) {
						if (notation == 1) { ENE_match_notation1(t_ps2.psets[var], group_NE_2[g2].vars[var], NUMINPUTS, mvcube.is_input_phase_assigment); }
						else if (notation == 2) { ENE_match_notation2(t_ps2.psets[var], group_NE_2[g2].vars[var], mvcube.is_input_phase_assigment); }
					}
					cur_g2 = g2 + 1;
				}
				else
					break;
			}
			// ***************** 考慮invers  NE <=> E (size = 2) *****************
			if (mvcube.is_input_phase_assigment && group_NE_1[g].vars.size() == 2) {
				for (int g2 = 0; g2 < group_E_2.size(); g2++) {
					for (int var = 0; var < group_E_2[g2].vars.size(); var++) {
						if (notation == 1) { ENE_match_notation1(t_ps2.psets[var], group_E_2[g2].vars[var], NUMINPUTS, mvcube.is_input_phase_assigment); }
						else if (notation == 2) { ENE_match_notation2(t_ps2.psets[var], group_E_2[g2].vars[var], mvcube.is_input_phase_assigment); }
					}
				}
			}
		}
		
		for (int var = 0; var < group_NE_1[g].vars.size(); var++) {
			var2 = group_NE_1[g].vars[var];
			or_exact_MvVar(ENE_mapping, t_ps2.psets[var], mvcube.mv_var_word_size, mvcube.first_word[var2], mvcube.last_word[var2], mvcube.first_bit[var2], mvcube.last_bit[var2]);
		}
	}


	// ***************** 考慮invers E <=> NE (size = 2) *****************
	if (mvcube.is_input_phase_assigment) {
		for (int var = 0; var < 2; var++) {
			set_clear(t_ps2.psets[var], mvcube.mv_var_size);
		}
		for (int g2 = 0; g2 < group_NE_2.size() && group_NE_2[g2].vars.size() >= 2; g2++) {
			if (group_NE_2[g2].vars.size() == 2) {
				for (int var = 0; var < group_NE_2[g2].vars.size(); var++) {
					if (notation == 1) { ENE_match_notation1(t_ps2.psets[var], group_NE_2[g2].vars[var], NUMINPUTS, mvcube.is_input_phase_assigment); }
					else if (notation == 2) { ENE_match_notation2(t_ps2.psets[var], group_NE_2[g2].vars[var], mvcube.is_input_phase_assigment); }
				}
			}
		}
		for (int g = 0; g < group_E_1.size(); g++) {
			for (int var = 0; var < group_E_1[g].vars.size(); var++) {
				var2 = group_E_1[g].vars[var];
				or_exact_MvVar(ENE_mapping, t_ps2.psets[var], mvcube.mv_var_word_size, mvcube.first_word[var2], mvcube.last_word[var2], mvcube.first_bit[var2], mvcube.last_bit[var2]);
			}
		}
	}

	// ***************** ENE <=> ENE *****************
	int pre_ENE_size = 0; cur_g2 = 0;
	for (int g = 0; g < group_ENE_1.size(); g++) {
		// 同數量的group_ENE 會對應到相同的var可能
		if (pre_ENE_size != group_ENE_1[g].vars.size()) {
			pre_ENE_size = group_ENE_1[g].vars.size();
			for (int var = 0; var < pre_ENE_size; var++) {
				set_clear(t_ps2.psets[var], mvcube.mv_var_size);
			}
			for (int g2 = cur_g2; g2 < group_ENE_2.size(); g2++) {
				if (group_ENE_1[g].vars.size() < group_ENE_2[g2].vars.size())
					continue;
				else if (group_ENE_1[g].vars.size() == group_ENE_2[g2].vars.size()) {
					// 1 to 1
					for (int var = 0; var < group_ENE_1[g].vars.size(); var++) {
						if (notation == 1) { ENE_match_notation1(t_ps2.psets[var], group_ENE_2[g2].vars[var], NUMINPUTS, mvcube.is_input_phase_assigment); }
						else if (notation == 2) { ENE_match_notation2(t_ps2.psets[var], group_ENE_2[g2].vars[var], mvcube.is_input_phase_assigment); }
					}
					cur_g2 = g2 + 1;
				}
				else
					break;
			}
		}
		for (int var = 0; var < group_ENE_1[g].vars.size(); var++) {
			var2 = group_ENE_1[g].vars[var];
			or_exact_MvVar(ENE_mapping, t_ps2.psets[var], mvcube.mv_var_word_size, mvcube.first_word[var2], mvcube.last_word[var2], mvcube.first_bit[var2], mvcube.last_bit[var2]);
		}
	}


	// X <=> X
	set_clear(t_ps2.psets[0], mvcube.mv_var_size);
	for (int var = 0; var < group_X_2.size(); var++) {
		if (notation == 1) { ENE_match_notation1(t_ps2.psets[0], group_X_2[var], NUMINPUTS, mvcube.is_input_phase_assigment); }
		else if (notation == 2) { ENE_match_notation2(t_ps2.psets[0], group_X_2[var], NUMINPUTS, mvcube.is_input_phase_assigment); }
	}
	//display_pset(t_ps2.psets[0], mvcube.mv_var_size);

	for (int var = 0; var < group_X_1.size(); var++) {
		var2 = group_X_1[var];
		or_exact_MvVar(ENE_mapping, t_ps2.psets[0], mvcube.mv_var_word_size, mvcube.first_word[var2], mvcube.last_word[var2], mvcube.first_bit[var2], mvcube.last_bit[var2]);
	}
	if (print_debug) {
		printf("vvvvvvvvvvvvvvENEvvvvvvvvvvvvvvv\n");
		display_pset_eachvar(ENE_mapping, mvcube.mv_var_size, mvcube.num_mv_var);
		printf("^^^^^^^^^^^^^^ENE^^^^^^^^^^^^^^^\n");
	}
	set_and(totality, totality, ENE_mapping);
}
/*
	KSIG 比對
	KSIG_input 是3維 [output][input].dist[dist]
	在不同output上相交的不同input各自擁有的1~n KSIG距離
	若PLA1 的var與PLA2 的var要相同則兩者在同個output上需要有相同的KSIG距離
*/

#define KSIG_match_notation1(p, pos1, pos2, mv_var_size, num_mv_var, phase_assigment) {set_remove(p, pos1*mv_var_size + pos2); if (phase_assigment) set_remove(p, pos1*mv_var_size + pos2 + num_mv_var);}
#define KSIG_match_notation2(p, pos1, pos2, mv_var_size, phase_assigment) {set_remove(p, pos1*mv_var_size + pos2*2+1); if (phase_assigment) set_remove(p, pos1*mv_var_size + pos2*2);}
void matching_KSIG(pset p, vector<TABLE>& KSIG_input1, vector<TABLE>& KSIG_input2) {
	bool isSame;
	for (int PLA1_input_i = 0; PLA1_input_i < NUMINPUTS; PLA1_input_i++) {
		for (int PLA2_input_i = 0; PLA2_input_i < NUMINPUTS; PLA2_input_i++) {
			isSame = true;
			for (int output_j = 0; output_j < KSIG_input1.size(); output_j++) { //NUMOUPUTS
				if (!(KSIG_input1[output_j][PLA1_input_i] == KSIG_input2[output_j][PLA2_input_i])) {
					
					isSame = false;
					break;
				}
			}
			if (!isSame) {
				if (notation == 1) { KSIG_match_notation1(p, PLA1_input_i, PLA2_input_i, mvcube.mv_var_size, NUMINPUTS, mvcube.is_input_phase_assigment); }
				else if (notation == 2) { KSIG_match_notation2(p, PLA1_input_i, PLA2_input_i, mvcube.mv_var_size, mvcube.is_input_phase_assigment); }
			}
		}
	}
}


bool compare_ksig_output(const pair<int, ksig_value*>& a, pair<int, ksig_value*>& b) {
	for (int i = 0; i < a.second->distx.size(); i++) {
		if (a.second->distx[i] < b.second->distx[i])
			return true;
		else if (a.second->distx[i] > b.second->distx[i])
			return false;
	}
	return false;
}
/*
	KSIG output match
	每個output 的KSIG 內容順序可以不一樣但結果必須一樣
*/
void matching_KSIG_output(pset p_out, TABLE &KSIG_output1, TABLE &KSIG_output2) {
	pset KSIG_p_out = set_new(mvcube.num_output * mvcube.num_output);
	vector<pair<int, ksig_value*>> output1(KSIG_output1.size()), output2(KSIG_output2.size());

	for (int i = 0; i < KSIG_output1.size(); i++) {
		sort(KSIG_output1[i].distx.begin(), KSIG_output1[i].distx.end(), std::greater<cpp_int>());
		output1[i] = make_pair(i, &KSIG_output1[i]);
	}
	for (int i = 0; i < KSIG_output2.size(); i++) {
		sort(KSIG_output2[i].distx.begin(), KSIG_output2[i].distx.end(), std::greater<cpp_int>());
		output2[i] = make_pair(i, &KSIG_output2[i]);
	}

	sort(output1.begin(), output1.end(), compare_ksig_output);
	sort(output2.begin(), output2.end(), compare_ksig_output);

	vector<bool> used1(KSIG_output2.size(), false), used2(KSIG_output2.size(), false);
	vector<int> group1, group2;
	bool same = true;

	for (int i = 0; i < output1.size(); i++) {
		if (used1[i])
			continue;
		used1[i] = true;
		group1.clear();
		group1.push_back(output1[i].first);
		for (int j = i+1; j < output1.size(); j++) {
			if (used1[j])
				continue;
			if (*output1[i].second == *output1[j].second) {
				group1.push_back(output1[j].first);
				used1[j] = true;
			}
		}
		group2.clear();
		for (int j = 0; j < output2.size(); j++) {
			if (used2[j])
				continue;
			if (*output1[i].second == *output2[j].second) {
				group2.push_back(output2[j].first);
				used2[j] = true;
			}
		}
		for (int& g : group1) {
			for (int& g2 : group2) {
				set_insert(KSIG_p_out, g * mvcube.num_output + g2);
			}
		}
	}

	set_and(p_out, p_out, KSIG_p_out);

	
	set_free(KSIG_p_out);
}


/*
	KSIG_input1 3D
	1D: each output (可排序)
	2D: each input (可排序)
	3D: each input distence (不可排序，每一個位置分別代表距離1、2、....、n)

	考慮output有permutation input也有permutation
	須將output 與input由小排到大方便比對

*/
bool compare_ksig_output_with_ksig_input1(const pair<int, vector<pair<int, ksig_value*>>>& a, pair<int, vector<pair<int, ksig_value*>>>& b) {
	for (int i = 0; i < a.second.size(); i++) {
		if (a.second[i].second != b.second[i].second) {
			return *a.second[i].second > *b.second[i].second;
		}
	}
	return false;
}

bool compare_ksig_output_with_ksig_input2(const pair<int, ksig_value*>& a,pair<int, ksig_value*>& b) {
	return *a.second > *b.second;
}

void matching_KSIG_output_with_KSIG_input(pset p_out, vector<TABLE>& KSIG_input1, vector<TABLE>& KSIG_input2) {
	vector<pair<int, vector<pair<int, ksig_value*>>>> output1(KSIG_input1.size()), output2(KSIG_input2.size());

	for (int i = 0; i < KSIG_input1.size(); i++) output1[i].second.resize(KSIG_input1[i].size());
	for (int i = 0; i < KSIG_input2.size(); i++) output2[i].second.resize(KSIG_input2[i].size());

	/*sort input*/
	for (int i = 0; i < KSIG_input1.size(); i++) {
		output1[i].first = i;
		for (int j = 0; j < KSIG_input1[i].size(); j++) {
			output1[i].second[j] = make_pair(j, &KSIG_input1[i][j]);
		}
		sort(output1[i].second.begin(), output1[i].second.end(), compare_ksig_output_with_ksig_input2);
	}
	sort(output1.begin(), output1.end(), compare_ksig_output_with_ksig_input1);

	for (int i = 0; i < KSIG_input2.size(); i++) {
		output2[i].first = i;
		for (int j = 0; j < KSIG_input2[i].size(); j++) {
			output2[i].second[j] = make_pair(j, &KSIG_input2[i][j]);
		}
		sort(output2[i].second.begin(), output2[i].second.end(), compare_ksig_output_with_ksig_input2);
	}
	sort(output2.begin(), output2.end(), compare_ksig_output_with_ksig_input1);

	vector<bool> used1(output1.size(), false), used2(output2.size(), false);
	vector<vector<int>> group1, group2;
	bool same = true;
	int cnt_group = 0;
	// grouping output1
	for (int i = 0; i < output1.size(); i++) {
		if (used1[output1[i].first])
			continue;
		used1[output1[i].first] = true;

		group1.push_back(vector<int>(1, output1[i].first));
		cnt_group++;
		for (int j = i + 1; j < output1.size(); j++) {
			if (used1[output1[j].first])
				continue;
			same = true;
			for (int k = 0; k < output1[i].second.size(); k++) {
				if (*output1[i].second[k].second != *output1[j].second[k].second) {
					same = false;
					break;
				}
			}
			if (same) {
				group1[cnt_group - 1].push_back(output1[j].first);
				used1[output1[j].first] = true;
			}
		}
	}

	cnt_group = 0;
	// grouping output2
	for (int i = 0; i < output2.size(); i++) {
		if (used2[output2[i].first])
			continue;
		used2[output2[i].first] = true;

		group2.push_back(vector<int>(1, output2[i].first));
		cnt_group++;
		for (int j = i + 1; j < output2.size(); j++) {
			if (used2[output2[j].first])
				continue;
			same = true;
			for (int k = 0; k < output2[i].second.size(); k++) {
				if (*output2[i].second[k].second != *output2[j].second[k].second) {
					same = false;
					break;
				}
			}
			if (same) {
				group2[cnt_group - 1].push_back(output2[j].first);
				used2[output2[j].first] = true;
			}
		}
	}


	//分群後 群數必須相同
	if (group1.size() != group2.size()) {
		fprintf(stderr, "KSIG output with input group cnt ERROR %d != %d\n", (int)group1.size(), (int)group2.size());
		exit(EXIT_FAILURE);
	}
	pset KSIG_p_out = set_new(mvcube.num_output * mvcube.num_output);
	for (int i = 0; i < group1.size(); i++) {
		//相同順序的群 群中的值必須相同
		if (group1[i].size() != group2[i].size()) {
			fprintf(stderr, "KSIG output with input group cnt cnt ERROR %d != %d\n", (int)group1.size(), (int)group2.size());
			exit(EXIT_FAILURE);
		}
		//二次驗證有點麻煩跳過


		for (int& g : group1[i]) {
			for (int& g2 : group2[i]) {
				set_insert(KSIG_p_out, g *mvcube.num_output + g2);
			}
		}
	}

	//display_pset_eachvar(KSIG_p_out, mvcube.num_output, mvcube.num_output);
	set_and(p_out, p_out, KSIG_p_out);


	set_free(KSIG_p_out);

	
}


void run_OutputSIG_TEST(pPLA PLA1, pPLA PLA2) {
	clock_t output_begin;
	pset mv_output = set_new(mvcube.num_output * mvcube.num_output);
	// only Unate
	output_begin = clock();
	set_fill(mv_output, mvcube.num_output * mvcube.num_output);
	EXECost(check_Unate_output(mv_output, PLA1, PLA2), "only unate output check");
	printf("cnt_unate %d\n", set_count_ones(mv_output));
	printf("cost_unate %.3lf\n", (double)(clock() - output_begin) / (double)CLOCKS_PER_SEC);

	// only KSO
	output_begin = clock();
	set_fill(mv_output, mvcube.num_output * mvcube.num_output);
	EXECost(check_KSIG(mv_output, PLA1, PLA2, 0, 0), "KSO check");
	printf("cnt_KSO %d\n", set_count_ones(mv_output));
	printf("cost_KSO %.3lf\n", (double)(clock() - output_begin) / (double)CLOCKS_PER_SEC);

	// only KSI
	output_begin = clock();
	set_fill(mv_output, mvcube.num_output * mvcube.num_output);
	EXECost(check_KSIG(mv_output, PLA1, PLA2, NUMINPUTS, 2), "KSI check");
	printf("cnt_KSI %d\n", set_count_ones(mv_output));
	printf("cost_KSI %.3lf\n", (double)(clock() - output_begin) / (double)CLOCKS_PER_SEC);

	// KSO + KSI
	output_begin = clock();
	set_fill(mv_output, mvcube.num_output * mvcube.num_output);
	EXECost(check_KSIG(mv_output, PLA1, PLA2, NUMINPUTS, 3), "KSO_KSI check");
	printf("cnt_KSO_KSI %d\n", set_count_ones(mv_output));
	printf("cost_KSO_KSI %.3lf\n", (double)(clock() - output_begin) / (double)CLOCKS_PER_SEC);


	// Unate + KSO
	output_begin = clock();
	set_fill(mv_output, mvcube.num_output * mvcube.num_output);
	EXECost(check_Unate_output(mv_output, PLA1, PLA2), "unate output check");
	EXECost(check_KSIG(mv_output, PLA1, PLA2, 0, 0), "KSO check");
	printf("cnt_unate_KSO %d\n", set_count_ones(mv_output));
	printf("cost_unate_KSO %.3lf\n", (double)(clock() - output_begin) / (double)CLOCKS_PER_SEC);

	// Unate + KSI
	output_begin = clock();
	set_fill(mv_output, mvcube.num_output * mvcube.num_output);
	EXECost(check_Unate_output(mv_output, PLA1, PLA2), "unate output check");
	EXECost(check_KSIG(mv_output, PLA1, PLA2, NUMINPUTS, 2), "KSI check");
	printf("cnt_unate_KSI %d\n", set_count_ones(mv_output));
	printf("cost_unate_KSI %.3lf\n", (double)(clock() - output_begin) / (double)CLOCKS_PER_SEC);

	// Unate + KSO + KSI
	output_begin = clock();
	set_fill(mv_output, mvcube.num_output * mvcube.num_output);
	EXECost(check_Unate_output(mv_output, PLA1, PLA2), "unate output check");
	EXECost(check_KSIG(mv_output, PLA1, PLA2, NUMINPUTS, 3), "KSO_KSI check");
	printf("cnt_unate_KSO_KSI %d\n", set_count_ones(mv_output));
	printf("cost_unate_KSO_KSI %.3lf\n", (double)(clock() - output_begin) / (double)CLOCKS_PER_SEC);

	set_free(mv_output);
}


/* setp_implies_var -- check if "a-var" implies "b-var" ("b-var" contains "a-var") */
bool setp_implies_var(pset a, pset b, int var)
{
	int i = mycube.last_word[var];
	do {
		if (a[i] & ~b[i] & mycube.var_mask[var][i]) // 檢查b是否包含a, setp_implies(a, b)
			return false;
	} while (--i >= mycube.first_word[var]);
	return true;
}

/* set_and_var -- compute intersection of sets "a" and "b-var" in varibale var */
/* a與b必須要相同size*/
/* true => feasible, false => infeasible*/
bool set_and_var(pset r, pset a, pset b, int var)
{
	/* 先處理要and的部分*/
	int i = mycube.last_word[var];
	bool isEmpty = true;
	do {
		r[i] = a[i] & (b[i] | (~mycube.var_mask[var][i]));
		/* 檢查and的部分是否全為0，只要有一處不為0則整個就不為empty*/
		if (isEmpty && ((r[i] & mycube.var_mask[var][i]) != 0))
			isEmpty = false;
	} while (--i >= mycube.first_word[var]);
	if (isEmpty)
		return false;
	/* 再處理其他單純複製的部分*/
	i = LOOP(a);
	PUTLOOP(r, i); 
	do if (i > mycube.last_word[var] || i < mycube.first_word[var]) r[i] = a[i]; while (--i > 0);
	
	return true;
}





pset set_leftshit_or_disjoint(pset r, pset p, int w_size) {
	for (int w = 1; w < w_size; w++) {
		r[w] = p[w] | (p[w] >> 1) & DISJOINT;
	}
	return r;
}



/*
	mvcube 合法性檢查
	1. 每一個var(row) 至少要對應到一個var(column)
		a. 透過去數每個var(row) 1 的個數
	2. 每一個var(column) 至少要被一個var(row) 對應到

	表示法為aa'bb'cc'時要先將a和a'的部分用shift or在一起
	ex aa'bb'cc' = 101101 => 010101
	
	不用檢查該var是否empty前面and時已檢查過

	return false 表示不合法
*/
#define free_sf_return_false(sf) {sf_free(sf); return false;}
//#define mvcube_feas_check_debug
bool mvcube_feas_check_notation2(pset p) {
	temp_psets tpsets(3);
	/* 儲存(var)temp的結果*/
	tpsets.psets[0] = set_new(mvcube.mv_var_size); pset p_var_t = tpsets.psets[0];
	/* var 的supercube，*/
	tpsets.psets[1] = set_new(mvcube.mv_var_size); pset var_supercube = tpsets.psets[1];
	/* 儲存(var)對應到已知解的結果*/
	tpsets.psets[2] = set_new(mvcube.mv_var_size); pset solution = tpsets.psets[2];
	

	pset_family sf_vars = sf_new(mvcube.num_mv_var, mvcube.mv_var_size);
	

	for (int var = 0; var < mvcube.num_mv_var; var++) {
		/* 將兩個bit合併*/
		get_mvcube_var(p_var_t, p, mvcube.mv_var_word_size, mvcube.first_word[var], mvcube.last_word[var], mvcube.first_bit[var], mvcube.last_bit[var]);
		rshilft_or_disjoint(p_var_t, p_var_t); // (p | (p>>1)) & DISJOINT
		PUTSIZE(p_var_t, set_count_ones(p_var_t));
		sf_addset(sf_vars, p_var_t);
		/* 計算var_supercube*/
		set_or(var_supercube, var_supercube, p_var_t);
	}
	/* supercube 檢查*/
	if (set_count_ones(var_supercube) != mvcube.num_mv_var) {
#ifdef mvcube_feas_check_debug
		printf("var supercube infeasible\n");
		display_pset_eachvar(p, mvcube.mv_var_size, mvcube.num_mv_var);
#endif
		free_sf_return_false(sf_vars);
	}
	/* 小到大排序，一樣的會擺一起*/
	SF_SORT_assend(sf_vars);
	pset last, tp, pre_tp;
	int index_p;
	pre_tp = p_var_t; /* 為了不跳warning 隨便初始化，不會實際用到p_var_t的空間*/

	sf_active(sf_vars);
	int pre_cnt = 0;
	/* 紀錄那些var對應到相同(mvcube)的個數*/
	int cnt_same = 1;
	/* 紀錄那些var已經找到配對*/
	int num_unmatch_var = mvcube.num_mv_var;

	foreachi_set(sf_vars, index_p, tp) {
		/* 若已知解已包含tp全部則tp為無解*/
		if (setp_implies(tp, solution)) {
#ifdef mvcube_feas_check_debug
			printf("已知解已包含tp全部則tp為無解 infeasible\n");
			display_pset_eachvar(p, mvcube.mv_var_size, mvcube.num_mv_var);
#endif
			free_sf_return_false(sf_vars);
		}
		if (SIZE(tp) == 1) {
			/* 1. 一個var只對應到一個解時，該解視為已對應到*/
			set_or(solution, solution, tp);
			RESET(tp, ACTIVE);
			num_unmatch_var--;
		}
		else if (SIZE(tp) < mvcube.num_mv_var) {
			/* 將solution 中的解inverse*/
			set_inverse_disjoint(p_var_t, solution);
			//printf("solution:      "); display_pset(solution, mvcube.mv_var_size);
			//printf("original tp:   "); display_pset(tp, mvcube.mv_var_size);
			//printf("solution inve: "); display_pset(p_var_t, mvcube.mv_var_size);
			
			/* tp去除當前已知解*/
			set_and(tp, p_var_t, tp);
			//printf("removed tp:    "); display_pset(tp, mvcube.mv_var_size);
			if (pre_cnt == SIZE(tp)) {
				if (setp_equal(tp, pre_tp)) {
					cnt_same++;
					/* n個var都有相同的n個解時，解只可能對應到該n個解，故可視為已對應到*/
					if (cnt_same == SIZE(tp)) {
						set_or(solution, solution, tp);
						for (int i = 0; i < cnt_same; i++) {
							RESET(GETSET(sf_vars, index_p - i), ACTIVE);
						}
						num_unmatch_var -= cnt_same;
					}
				}
				else {
					cnt_same = 1;
				}
			}
			else {
				pre_cnt = SIZE(tp);
				pre_tp = tp;
				cnt_same = 1;
			}
		}
		/* 剩餘皆是full cube 跳過不比對*/
		else {
			break;
		}
		/*2. 尚未配對的var比尚未配對的解還多，則表示無解*/
		if (num_unmatch_var > mvcube.num_mv_var - set_count_ones(solution)) {
#ifdef mvcube_feas_check_debug
			printf("尚未配對的var比尚未配對的解還多 infeasible\n");
			display_pset_eachvar(p, mvcube.mv_var_size, mvcube.num_mv_var);
#endif
			free_sf_return_false(sf_vars);
		}
	}

	/* 還有remove的部分可以做*/

	sf_free(sf_vars);
	return true;
}

void sf_removeRedundant(pset_family A, int skip_num) {
	// 如果有被包含則移除
	pset p, p2;
	int i, j;
	sf_active(A);
	foreachi_set(A, i, p) {
		if (!TESTP(p, ACTIVE))
			continue;

		foreachi_set(A, j, p2) {
			if (i < skip_num && j < skip_num) {
				j = skip_num;
				p2 = GETSET(A, j);
			}
			if (j <= i || !TESTP(p2, ACTIVE))
				continue;
			if (setp_implies(p2, p))
				RESET(p2, ACTIVE);
			else if (setp_implies(p, p2))
				RESET(p, ACTIVE);
		}
	}
	sf_inactive(A);
}

/* vvvvvvvvvvvvvvvvvvvvvvvv abcda'b'c'd' 合法性檢查 vvvvvvvvvvvvvvvvvvvvvvvv*/

class v_cube {
public:
	pset s;
	int num_1bit;
	int var;

	v_cube() {
		s = set_new(NUMINPUTS);
		num_1bit = 0;
		var = -1;
	};
	v_cube(const v_cube& vcube) {
		s = set_new(NUMINPUTS);
		set_copy(s, vcube.s);
		num_1bit = vcube.num_1bit;
		var = vcube.var;
	}
	~v_cube() {
		set_free(s);
	}
	void set_pset(pset p, int cnt_1bit, int _var) {
		s = set_new(mvcube.num_mv_var); 
		set_copy(s, p);
		num_1bit = cnt_1bit;
		var = _var;
	}
	void set_num1bit(int cnt_1bit) {
		num_1bit = cnt_1bit;
	}
	v_cube& operator=(const v_cube& vcube) {
		set_free(s);
		s = set_new(NUMINPUTS);
		set_copy(s, vcube.s);
		num_1bit = vcube.num_1bit;
		var = vcube.var;
		return *this;
	}
	bool operator<(const v_cube& vcube) const {
		if (this->num_1bit < vcube.num_1bit)
			return true;
		else if (this->num_1bit == vcube.num_1bit) {
			for (int i = 1; i <= LOOP(s); i++) {
				if (this->s[i] < vcube.s[i])
					return true;
				else if (this->s[i] > vcube.s[i])
					return false;
			}
		}
		return false;
	}
	bool operator==(const v_cube& vcube) const {
		if (this->num_1bit == vcube.num_1bit) {
			for (int i = 1; i <= LOOP(s); i++) {
				if (this->s[i] != vcube.s[i])
					return false;
			}
			return true;
		}
		return false;
	}

};
pset set_Not(pset p, int _size) {
	unsigned int i = 1;
	for (; i < LOOP(p); i++) {
		p[i] = ~p[i];
	}
	p[i] = p[i] ^ (0xFFFFFFFF >> (BPI - (_size % BPI)));

	return p;
}
int cnt1bit(pset p) {
	int cnt = 0;
	for (unsigned int i = 1; i <= LOOP(p); i++) {
		cnt += count_ones(p[i]);
	}
	return cnt;
}

/*
	對mvcube中的某一個特定的var做pset的and
	r 是mvcube, p 是要and的pset
	r 的size > p 的size
		ex: INPUT = 3
			var (0) (1) (2)
			r = 111,111,111
			p = 101 (var 1)

				111 111 111
			and	    101
			----------------
			r = 111 101 111

*/
pset andExactMvVar(pset r, pset p, bool& isEmpty, unsigned int first_word, unsigned int last_word, unsigned int first_bit, unsigned int last_bit) {
	unsigned int val = 0;
	isEmpty = true;
	for (unsigned int M = first_word, P = 1; M <= last_word; M++, P++) {
		val |= (p[P] << first_bit);
		if (M == first_word)
			val |= ~(0xFFFFFFFF << first_bit);
		if (M == last_word)
			val |= (0xFFFFFFFF << (last_bit + 1));
		r[M] &= val;
		if ((BPI - first_bit) != 32)
			val = (p[P] >> (BPI - first_bit)); // p[P]中剩餘的bits
		else
			val = 0;
		if (isEmpty) {
			if (M == first_word && M == last_word) {
				if ((r[M] & (0xFFFFFFFF << first_bit) & (0xFFFFFFFF >> ((BPI - 1) - last_bit))) != 0)
					isEmpty = false;
			}
			else if (M == first_word) {
				if ((r[M] & (0xFFFFFFFF << first_bit)) != 0)
					isEmpty = false;
			}
			else if (M == last_word) {
				if ((r[M] & (0xFFFFFFFF >> ((BPI - 1) - last_bit))) != 0)
					isEmpty = false;
			}
			else {
				if (r[M] != 0)
					isEmpty = false;
			}
		}
	}
	return r;
}

/*  check isFeasible_hard
	1. check if row is empty (each y mapping to each x)
	2. check if column is empty (each x mapping to each y)
*/
bool mvcube_feas_check_notation1(pset A) {
	/* dynamic space for temporary use*/
	temp_psets tpsets(5);
	tpsets.psets[0] = set_new(NUMINPUTS);
	tpsets.psets[1] = set_new(NUMINPUTS);
	tpsets.psets[2] = set_new(NUMINPUTS);
	tpsets.psets[3] = set_new(NUMINPUTS);
	tpsets.psets[4] = set_new(NUMINPUTS);

	pset* left_r = &tpsets.psets[0];	/* get var left part from pset*/
	pset* right_r = &tpsets.psets[1];	/* get var right part from pset*/
	pset* r_supercube = &tpsets.psets[2]; /* supercube of left_r*/
	pset* solution = &tpsets.psets[3];	/*known solution*/
	pset* resolution = &tpsets.psets[4];	/*reverse(known solution)*/

	v_cube vc;
	multiset<v_cube> v_cube_set;
	vector<int> var_status(NUMINPUTS, 0); /* 0 = full, 1 = known solution, 2 = else*/
	int cnt_full_set = 0;

	bool isEmpty = false;

	// solution changed, update resolution
	set_copy(*resolution, *solution);
	set_Not(*resolution, NUMINPUTS);


	for (int var = 0; var < NUMINPUTS; var++) {
		if (!mvcube.is_input_phase_assigment) {
			get_mvcube_var(*left_r, A, mvcube.mv_var_word_size, mvcube.first_word[var], mvcube.last_word[var], mvcube.first_bit[var], mvcube.last_bit[var]);
		}
		else {
			/* break input into two part*/
			/* part 1: positive part*/
			get_mvcube_var(*left_r, A, mvcube.mv_var_word_size, mvcube.first_word[var], mvcube.lhalf_word[var], mvcube.first_bit[var], mvcube.lhalf_bit[var]);
			/* part 2: negative part*/
			get_mvcube_var(*right_r, A, mvcube.mv_var_word_size, mvcube.rhalf_word[var], mvcube.last_word[var], mvcube.rhalf_bit[var], mvcube.last_bit[var]);
			/* OR part 1 and part 2 to check if x match to any y*/
			set_or(*left_r, *left_r, *right_r);
		}
		if (setp_empty(*left_r)) { /*1. check if row is empty (each x mapping to at least one y)*/
			return false;
		}
		if (!setp_full(*left_r, NUMINPUTS)) { // ignore full sets
			set_and(*left_r, *left_r, *resolution);

			set_or(*r_supercube, *r_supercube, *left_r); /* cal supercube*/
			int cnt_1bit = cnt1bit(*left_r);
			if (cnt_1bit == 0)
				return false;
			else if (cnt_1bit == 1) { /* x match only one y*/
				var_status[var] = 1; /* var x matched y*/
				set_or(*solution, *solution, *left_r); /* get solution from those x match only one y*/
				set_copy(*resolution, *solution);
				set_Not(*resolution, NUMINPUTS);
			}
			else { /* sets that maybe be reduced by known solution*/
				var_status[var] = 2;
				vc.set_pset(*left_r, cnt_1bit, var);
				v_cube_set.insert(vc); /* multiset to sort the v_cube (increase by number of 1)*/
			}
		}
		else { // count the full sets (those been ignored) to check rest of the solutions are enough or not
			cnt_full_set++;
		}

	}
	int cnt_dc = cnt_full_set;

	/* reduced known solution from other var*/
	int cnt_same = 1;
	//cout << endl;
	for (multiset<v_cube>::iterator it = v_cube_set.begin(); it != v_cube_set.end(); it++) {

		if (next(it, 1) != v_cube_set.end() && it == next(it, 1)) {
			cnt_same++;
			continue;
		}
		else {
			//printf("before: "); display_pset(it->s, NUMINPUTS);
			set_and(it->s, it->s, *resolution);
			//printf("mask  : "); display_pset(*resolution, NUMINPUTS);
			//printf("after : "); display_pset(it->s, NUMINPUTS);
			int cnt_1bit = cnt1bit(it->s);

			if (cnt_1bit == cnt_same) {
				for (int i = 0; i < cnt_same; i++) {
					int var = prev(it, i)->var;
					var_status[var] = 1;
					if (!mvcube.is_input_phase_assigment) {
						andExactMvVar(A, *resolution, isEmpty, mvcube.first_word[var], mvcube.last_word[var], mvcube.first_bit[var], mvcube.last_bit[var]);
					}
					else {
						/* break input into two part*/
						/* part 1: positive part*/
						andExactMvVar(A, *resolution, isEmpty, mvcube.first_word[var], mvcube.lhalf_word[var], mvcube.first_bit[var], mvcube.lhalf_bit[var]);
						/* part 2: negative part*/
						andExactMvVar(A, *resolution, isEmpty, mvcube.rhalf_word[var], mvcube.last_word[var], mvcube.rhalf_bit[var], mvcube.last_bit[var]);
					}
				}
				// solution changed, update resolution
				set_or(*solution, *solution, it->s);
				set_copy(*resolution, *solution);
				set_Not(*resolution, NUMINPUTS);
			}
			else if (cnt_1bit < cnt_same) {
				return false;
			}
			cnt_same = 1;
		}
	}

	if (cnt1bit(*r_supercube) < NUMINPUTS - cnt_full_set) { /*2. check if column is empty (each y mapping to at least one x)*/
		return false;
	}

	/* remove known soltuion from esle (not full set and not known solution)*/
	//for (int var = 0; var < NUMINPUTS; var++) {
	//	if (var_status[var] != 1) {
	//		if (!mvcube.isInPhase) {
	//			andExactMvVar(mv_cube, *resolution, isEmpty, mvcube.first_word[var], mvcube.last_word[var], mvcube.first_bit[var], mvcube.last_bit[var]);
	//		}
	//		else {
	//			/* break input into two part*/
	//			/* part 1: positive part*/
	//			andExactMvVar(mv_cube, *resolution, isEmpty, mvcube.first_word[var], mvcube.lhalf_word[var], mvcube.first_bit[var], mvcube.lhalf_bit[var]);
	//			/* part 2: negative part*/
	//			andExactMvVar(mv_cube, *resolution, isEmpty, mvcube.rhalf_word[var], mvcube.last_word[var], mvcube.rhalf_bit[var], mvcube.last_bit[var]);
	//		}
	//	}
	//}
	return true;
}

//int cnt_solution(pset_family& totality) {
//	used.resize(NUMINPUTS, 0);
//	tpsets.psets.resize(NUMINPUTS);
//	for (int i = 0; i < NUMINPUTS; i++) {
//		tpsets.psets[i] = set_new(mvcube.num_mv_var);
//	}
//	for (int i = 0; i < totality.mv_list.size(); i++) {
//		for (int j = 0; j < NUMINPUTS; j++)
//			used[j] = 0;
//
//		for (int var = 0; var < NUMINPUTS; var++) {
//			extractMvVar(tpsets.psets[var], totality.mv_list[i].mv_cube, mvcube.first_word[var], mvcube.last_word[var], mvcube.first_bit[var], mvcube.last_bit[var]);
//		}
//
//		find_y(0, 1);
//
//	}
//	printf("solution number: %d\n", cnt_solution_num);
//	return 1;
//}


/* ^^^^^^^^^^^^^^^^^^^^^^^^ abcda'b'c'd' 合法性檢查 ^^^^^^^^^^^^^^^^^^^^^^^^*/