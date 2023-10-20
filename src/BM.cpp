/*
	布林比對的主程式
	所有的功能與指令都須從這使用

	布林比對步驟
	1. 初始設定全域變數與初始化快速計算時所需資料結構
	2. 根據輸入參數決定是CSF比對、ISF比對 等
	3. 對PLA進行初步的轉換，使其更易於且快速計算成MvCube的格式
	4. 計算輸入(出)特徵值，並進行比對，移除不可能的解
	5. 根據輸入參數決定比對的次數與方式(部分平行化、全部平行化等)
	6. 檢查最終計算結果 若有解則比對成功，反之失敗
*/

#include "mainCplusplus.h"
#include "PBM.h"

rand_data DATA; /* 對輸入PLA做打亂/擴展/移除 等處理*/
mv_cube mvcube;	/* 基本運算資訊*/
cube_struct mycube; /* 基本運算資訊 (mask)*/

pm_Detail pm_detail;
string name;
double total_and_time, total_redundant_time;

mutex memory_mtx;


ThreadPool pool;
clock_t AND_time, RR_time;

bool noINPA_AND_Notation2;
bool isNUMINPUTSchange;

void BM(pPLA PLA1, pPLA PLA2, bool INP, bool INPA, bool OUTP, bool OUTPA) {
	AND_time = 0, RR_time = 0;
	name = PLA1->filename;
	//stdout stderr
	noINPA_AND_Notation2 = false;
	isNUMINPUTSchange = false;
	
	if (INPA == false && notation == 2) {
		noINPA_AND_Notation2 = true;
		INPA = true;
	}

	// 是否有擴展輸入個數
	if (extand_input_num != 0) {
		if (extand_input_num > NUMINPUTS) {
			isNUMINPUTSchange = true;
			DATA.extandInput(PLA1, PLA2, extand_input_num);
			//DATA.extandInput(PLA2, extand_input_num);
		}

		// 若輸入個數改變則運算基本資訊重新設定
		if (isNUMINPUTSchange) {
			int numoutputs = NUMOUTPUTS;
			save_cube_struct();
			cube.num_binary_vars = extand_input_num;
			cube.num_vars = extand_input_num + 1;
			cube.part_size = ALLOC(int, cube.num_vars);
			for (int i = 0; i < cube.num_binary_vars; i++) {
				cube.part_size[i] = 2;
			}
			cube.part_size[cube.num_binary_vars] = numoutputs;
			cube_setup();
		}
	}
	
	/* 輸出PLA基本資訊*/
	fprintf(stderr, "%s, #in:%d, #out:%d, #1(on,off):(%d, %d), #2(on,off):(%d, %d), chunk(%d, %d), (thread_num, chunk_size, type): (%d, %d, %s)\n", PLA1->filename, NUMINPUTS, NUMOUTPUTS, PLA1->F->count, PLA1->R->count, PLA2->F->count, PLA2->R->count, row_num, col_num, ppbm_thread_num, ppbm_chunk_size, (ppbm_option == 1) ? "OMP" : (ppbm_option == 2) ? "CPP" : "TBB");

	/* 初始化運算基本資訊*/
	mvcube.init(INP, INPA, OUTP, OUTPA);  //mvcube.display();
	DATA.set(mvcube.num_mv_var, mvcube.num_output);
	/* 亂數測資 設為true*/
	DATA.random(PLA1, PLA2, RIP, RIPA, ROP, ROPA); // 打亂PLA2  亂數PLA2輸入排序, 亂數PLA2輸入相位轉換, 亂數PLA2輸出排序, 亂數PLA2輸出相位轉換
	DATA.display();


	pm_detail.init();

	clock_t start = clock();

	/* 初始化totality, supercube*/
	pset_family Totality = sf_new(1, mvcube.mv_size);
	if (!noINPA_AND_Notation2) {
		sf_addset(Totality, set_full(mvcube.mv_size));
	}
	else {
		sf_addset(Totality, set_full_DISJOINT(mvcube.mv_size));
	}

	/* 檢查所有輸出特徵值的組合*/
	if (EXP_SIG_OUTPUT == 1) {
		cout << "OUTPUT SIG TEST\n";
		run_OutputSIG_TEST(PLA1, PLA2);
		return;
	}

	// 如果移除部分on/off set則為ISF 不可使用輸入特徵值
	if (remove_sf_rate != 0)
		SIG_key = 0;

	/* 檢查所有輸入特徵值的組合*/
	if (EXP_SIG_INPUT == 1) {
		cout << "INPUT SIG TEST\n";
		run_InputSIG_TEST(GETSET(Totality, 0), PLA1, PLA2);
		return;
	}


	if (EXP_ALL_MCQE == 1) {
		/* 計算所有MCQE的組合*/
		MCDB mcdb(PLA1->F, PLA2->R, mvcube.num_mv_var, mvcube.is_input_phase_assigment);
		mcdb.pre_run(mvcube.first_word, mvcube.last_word, mvcube.first_bit, mvcube.last_bit);
		mcdb.run_RULE_TEST();
		return;
	}

	/* 計算並比對輸入特徵值*/
	/* 0b0001: check_unate, 0b0010 = check_ENE, 0b0100 = check_KSIG, 0b1000 = check_cofactor */
	check_character_value(SIG_key, GETSET(Totality, 0), PLA1, PLA2);
	

	/* 移除部分on/off set => ISF*/
	if (remove_sf_rate != 0)
		DATA.removeSF(PLA1, PLA2, remove_sf_rate, false); /* false表示不為亂數移除*/

	/* 初始化supercube 為full cube*/
	pset Supercube;
	Supercube = set_new(mvcube.mv_size);
	set_copy(Supercube, GETSET(Totality, 0));

	/* 計算Totality 1的個數*/
	PUTSIZE(GETSET(Totality, 0), set_count_ones(GETSET(Totality, 0)));

	// 初始化mycube (運算基本資訊(mask))
	mycube.num_binary_vars = 0;
	mycube.num_vars = mvcube.num_mv_var;
	mycube.part_size = ALLOC(int, mycube.num_vars);
	for (int i = mycube.num_binary_vars; i < mycube.num_vars; i++) {
		mycube.part_size[i] = mvcube.mv_var_size;
	}
	cube_setup(mycube);

	total_and_time = 0;
	total_redundant_time = 0;

	
	if (Total_parallel == 1) { /* 全部平行化*/
		Parallel_BM PBM(ppbm_thread_num, 100);
		if (row_num == 0) row_num = PLA1->F->count;
		if (col_num == 0) col_num = PLA2->R->count;
		PBM.setup(PLA1->F, PLA2->R, Totality, Supercube, 0b0111);
		//PBM.setup(PLA1->R, PLA2->F, Totality, Supercube, 0b0111);
		if (ppbm_option == 2)
			pool.setup(ppbm_thread_num);
		EXECost(PBM.run(), "Total parallel Partial mapping");
	}
	else { /* 部分平行化 or 未平行化*/ /* on-set => off-set*/
		if (row_num == 0) row_num = PLA1->F->count;
		if (col_num == 0) col_num = PLA2->R->count;
		/* PLA 前處理 將PLA 轉成 A, A'*/
		MCDB mcdb(PLA1->F, PLA2->R, mvcube.num_mv_var, mvcube.is_input_phase_assigment);
		mcdb.pre_run(mvcube.first_word, mvcube.last_word, mvcube.first_bit, mvcube.last_bit);

		/* 計算MCQE*/
		mcdb.check_Rule(MCQE_key); /* 0bDCBA  A:R1 B:R2 C:R3 D:R4*/
		printf("After MCDB: %d\n", mcdb.cal_TABLE_PQ());
		//mcdb.display();
		//mcdb.display_cost();
	

		/* 非單執行緒才建立執行緒池*/
		if (ppbm_option == 2)
			pool.setup(ppbm_thread_num);

		/* 執行partial mapping on 對 off*/
		EXECost(partial_mapping(Totality, Supercube, mcdb), "Partial mapping");

		/* 若不為CSF則需額外執行off 對on*/
		if (remove_sf_rate != 0) {
			MCDB mcdb2(PLA1->R, PLA2->F, mvcube.num_mv_var, mvcube.is_input_phase_assigment);
			mcdb2.pre_run(mvcube.first_word, mvcube.last_word, mvcube.first_bit, mvcube.last_bit);
			mcdb2.check_Rule(MCQE_key); /* 0bDCBA  A:R1 B:R2 C:R3 D:R4*/
			
			EXECost(partial_mapping(Totality, Supercube, mcdb2), "Partial mapping");
		}
		printf("total and redundant %.3lf %.3lf\n", total_and_time, total_redundant_time);
	}
	/* 全部完成前再次確定所有多餘已被移除*/
	sf_removeRedundant(Totality, 0);
	printf("FINAL sf_totality %d\n", Totality->count);

	/* 如果有解 則輸出其中一個解*/
	if (Totality->count > 0) {
		printf("%d %d\n", mvcube.mv_var_size, mvcube.num_mv_var);
		display_pset_eachvar(GETSET(Totality, 0), mvcube.mv_var_size, mvcube.num_mv_var);
		//display_pset(GETSET(Totality, 0), mvcube.mv_size);
		printf("%s\n", mvcube_feas_check_notation2(GETSET(Totality, 0)) ? "true" : "false");
		//display_sf_eachvar(Totality, mvcube.mv_var_size, mvcube.num_mv_var);
	}

	if (print_detail)
		pm_detail.display();

	printf("AND RR time %.6f %.6f\n", (double)AND_time / (double)CLOCKS_PER_SEC, (double)RR_time / (double)CLOCKS_PER_SEC);
	printf("%s Total: %.6lf\n", PLA1->filename, ((double)clock() - (double)start) / (double)CLOCKS_PER_SEC);
	setdown_cube(mycube);

	if (isNUMINPUTSchange) {
		setdown_cube();
		restore_cube_struct();
	}
}



//#define Supercube_detail
#ifdef Supercube_detail
vector<int> cnt_supercube_pq;
#endif
pset_family partial_mapping(pset_family Totality, pset Supercube, MCDB& mcdb) {
	// 分割資料
	chunk chs(MAX(row_num, 1), MAX(col_num, 1), mcdb.TABLE_PQ->count, mcdb.TABLE_PQ->sf_size);
	chs.order1(); //分割時的資料排序方式
	//chs.display();
	printf("number of chunks: %d\n", chs.chunk_cnt); // 分割的個數

	int sequential_end = chs.chunk_cnt;

	clock_t chunk_start;
	pset_family temp_Totality = sf_new(1, mvcube.mv_size);
	pset_family temp_Totality2 = sf_new(1, mvcube.mv_size);
	
	
	int window = 2; // 這裡window只能是2

	for (int i = 0; i < sequential_end; i++) {
		chunk_start = clock();
		/*以下擇一使用*/

		/*1. 部分平行化，PPBM設為1則為未平行化的版本  (主要使用這個)*/
		if (EXP_WINDOW_SIZE == 0) {
			partial_parallel_mapping_block(chs.chunks[i].f_begin, chs.chunks[i].f_end, chs.chunks[i].g_begin, chs.chunks[i].g_end, Totality, temp_Totality, temp_Totality2, Supercube, mcdb);
		}
		else {
			/*2. cm150a 複數積項配對(未平行化) window size = 2的實驗*/
			partial_parallel_mapping_block_window(chs.chunks[i].f_begin, chs.chunks[i].f_end, chs.chunks[i].g_begin, chs.chunks[i].g_end, Totality, temp_Totality, temp_Totality2, Supercube, mcdb, window);
		}
		/*3. 原始paper的重現實驗 */
		//paper_partial_mapping_block(chs.chunks[i].f_begin, chs.chunks[i].f_end, chs.chunks[i].g_begin, chs.chunks[i].g_end, Totality, temp_Totality, Supercube, mcdb);
		
		/*4. 第一版改良後的PM(完全未平行化版)*/  /*!!!可能有缺少部分加速技巧!!! 未使用此做實驗!!!*/
		//partial_mapping_block(chs.chunks[i].f_begin, chs.chunks[i].f_end, chs.chunks[i].g_begin, chs.chunks[i].g_end, Totality, temp_Totality, Supercube, mcdb);
		
		if (print_detail) printf("Sequential_Totality_%s%d\nSequential_Supercube_%s%d\n", name.c_str(), Totality->count, name.c_str(), set_count_ones(Supercube));
		if (print_detail) printf("%s_chunk cost: %.3f\n", name.c_str(), (float)(clock() - chunk_start) / (float)(CLOCKS_PER_SEC));
		//pm_detail.display();
	}

#ifdef Supercube_detail
	if (pm_detail.outputintersection < 100) {
		for (int i = 0; i < cnt_supercube_pq.size(); i++)
			printf("%sSupercube_cnt %d\n", name.c_str(), cnt_supercube_pq[i]);
	}
	else {
		for (int i = 0; i < cnt_supercube_pq.size(); i++)
			if (i % (cnt_supercube_pq.size() / 100) == 0 || i == cnt_supercube_pq.size() - 1)
				printf("%sSupercube_cnt %d\n", name.c_str(), cnt_supercube_pq[i]);
	}
#endif
	return Totality;
}
//#define S_MAP_debug
bool S_MAP(pset_family Totality, pset_family temp_Totality, pset Supercube, pset literal_pos, pset Q_non_DC, pset p_and) {
	pset p, last;
	int var;
	for (int i = 1; i < SIZE(literal_pos); i++) {
		var = literal_pos[i];
		if (setp_implies_var(Supercube, Q_non_DC, var)) {
			if (print_detail)
				pm_detail.smartMvcube_redundant++;
#ifdef S_MAP_debug
			printf("S_MAP isRedundant\n");
			printf("Supercube\n"); display_pset(Supercube, mvcube.mv_size);
			printf("Q_non_DC\n"); display_pset(Q_non_DC, mvcube.mv_var_size);
			printf("var %d\n", var);
#endif
			return true;
		}
		if (!set_and_var(p_and, Supercube, Q_non_DC, var)) {
			literal_pos[i] = 0xFFFFFFFF;
			if (print_detail)
				pm_detail.smartMvCube_empty_cube_temp++;
#ifdef S_MAP_debug
			printf("S_MAP infeasible\n");
			printf("Supercube\n"); display_pset_eachvar(Supercube, mvcube.mv_var_size, mvcube.num_mv_var);
			printf("Q_non_DC\n"); display_pset(Q_non_DC, mvcube.mv_var_size);
			printf("mask\n"); display_pset_eachvar(mycube.var_mask[var], mvcube.mv_var_size, mvcube.num_mv_var);
			printf("var %d\n", var);
			printf("p_and\n"); display_pset_eachvar(p_and, mvcube.mv_var_size, mvcube.num_mv_var);
#endif
			continue;
		}
		foreach_active_set(Totality, last, p) {
			if (setp_implies_var(p, p_and, var)) {
				sf_addset(temp_Totality, p);
				RESET(p, ACTIVE);
			}
		}
	}
	return false;
}




// 以下為第一版改良後的PM (完全未平行化版)  /*!!!可能有少部分加速技巧!!!*/ /*全部刪除仍可正常運行*/
int partial_mapping_block_cnt_intersect_pair = 0;
clock_t time_begin = clock();
void partial_mapping_block(int f_begin, int f_end, int g_begin, int g_end, pset_family Totality, pset_family temp_Totality, pset Supercube, MCDB& mcdb) {
	bool isRedundant, isEmpty;
	pset p, last, p_lit0, p_lit1, p_and, p_f, Q_non_DC, q_non_DC, temp_Supercube, p_get_var;
	temp_psets t_ps(5);
	t_ps.psets[0] = set_new(mcdb.P->sf_size * 32); p_lit0 = t_ps.psets[0];
	t_ps.psets[1] = set_new(mcdb.P->sf_size * 32); p_lit1 = t_ps.psets[1];
	t_ps.psets[2] = set_new(mvcube.mv_size); p_and = t_ps.psets[2];
	t_ps.psets[3] = set_new(mvcube.mv_size); temp_Supercube = t_ps.psets[3];
	t_ps.psets[4] = set_new(mvcube.mv_var_size); p_get_var = t_ps.psets[4];


	clock_t smartmv_cost, and_cost, redun_cost;
	for (int f_on = f_begin; f_on < f_end; f_on++) {
		p_f = GETSET(mcdb.TABLE_PQ, f_on);

		if (setp_empty(p_f))
			continue;


		for (int g_off = g_begin; g_off < g_end; g_off++) {
			//printf("\n\nPair ( %d %d ) %.3lf\n\n", f_on, g_off, (double)(clock() - time_begin) / (double)(CLOCKS_PER_SEC));


			// 檢查是否符合前處理 (ex output相交...)
			if (is_in_set(p_f, g_off) == 0) {
				continue;
			}
#ifdef Supercube_detail
			cnt_supercube_pq.push_back(set_count_ones(Supercube));
#endif


			pm_detail.outputintersection++;
			temp_Totality->count = 0;

			//printf("ON OFF = (%d, %d)==========================\n", f_on, g_off);

			/* 複製literal postion 變成local*/
			mcdb.lit_pos_copy(p_lit0, GETSET(mcdb.P_lit0_pos, f_on));
			mcdb.lit_pos_copy(p_lit1, GETSET(mcdb.P_lit1_pos, f_on));


			Q_non_DC = GETSET(mcdb.ex_Q_non_DC, g_off);
			q_non_DC = GETSET(mcdb.ex_q_non_DC, g_off);


			/*smartmvcube
			1.將(p_i, q_j)與supercube做檢查，判斷(p_i, q_j)是否可以跳過
			2.同時檢查(p_i, q_j)中的各個pset組合是否feasible
			3.checkredundantbeforeand*/

			sf_active(Totality);

			isRedundant = false;
			if (print_detail) {
				pm_detail.smartMvCube_empty_cube_temp = 0;
			}
			if (S_MAP(Totality, temp_Totality, Supercube, p_lit0, Q_non_DC, p_and)) {
				if (print_detail) pm_detail.smartMvCube_Mvcubes_num_redundant += SIZE(p_lit0) + SIZE(p_lit1);
				continue;
			}

			if (S_MAP(Totality, temp_Totality, Supercube, p_lit1, q_non_DC, p_and)) {
				if (print_detail) pm_detail.smartMvCube_Mvcubes_num_redundant += SIZE(p_lit0) + SIZE(p_lit1);
				continue;
			}

			if (print_detail) {
				pm_detail.smartMvCube_Mvcubes_num += SIZE(p_lit0) + SIZE(p_lit1);
				pm_detail.smartMvCube_empty_cube += pm_detail.smartMvCube_empty_cube_temp;
				if (pm_detail.smartMvCube_empty_cube_temp > 0)
					pm_detail.smartMvCube_empty_cube_pq++;
			}
			sf_inactive(Totality);
			if (Totality->count == 0) {
				if (temp_Totality->count != 0) {
					my_sf_copy(Totality, temp_Totality);
				}
				continue;
			}

			isEmpty = true;
			for (int i = 1; i < SIZE(p_lit0); i++) {
				/* 表示已被supercube移除*/
				if (p_lit0[i] != 0xFFFFFFFF) {
					isEmpty = false;
					break;
				}
			}
			if (!isEmpty) {
				for (int i = 1; i < SIZE(p_lit1); i++) {
					/* 表示已被supercube移除*/
					if (p_lit1[i] != 0xFFFFFFFF) {
						isEmpty = false;
						break;
					}
				}
			}

			if (isEmpty) {
				foreach_set(temp_Totality, last, p) {
					sf_addset(Totality, p);
				}
				continue;
			}

			and_cost = clock();

			//printf("Q_non_DC: "); display_pset(Q_non_DC, mvcube.mv_var_size);
			partital_map_and(Totality, temp_Totality, Supercube, p_lit0, Q_non_DC, p_and);

			//printf("q_non_DC: "); display_pset(q_non_DC, mvcube.mv_var_size);
			partital_map_and(Totality, temp_Totality, Supercube, p_lit1, q_non_DC, p_and);

			total_and_time += (double)(clock() - and_cost) / (double)(CLOCKS_PER_SEC);
			if (temp_Totality->count == 0) {
				printf("temp_Totality EMPTY!!!\n");
				exit(EXIT_FAILURE);
			}
			redun_cost = clock();
			sf_removeRedundant(temp_Totality, 0);
			total_redundant_time += (double)(clock() - redun_cost) / (double)(CLOCKS_PER_SEC);
			sf_copy2(Totality, temp_Totality);


			sf_or2(temp_Supercube, temp_Totality);
			set_copy(Supercube, temp_Supercube);

			//fprintf(stderr, "Sequential_Supercube  %d\n", set_count_ones(Supercube));
			//printf("Sequential_Supercube  %d\n", set_count_ones(Supercube));
		}
		if (temp_Totality->count > 10000000)
			break;
		//printf("(%d, %d) %d %d\n", f_on, mcdb.TABLE_PQ->sf_size, Totality->count, set_count_ones(Supercube));

	}
	if (temp_Totality->count > 10000000)
		printf("\n\n over 10000000\n");
	//printf("%d %d\n", Totality->count, set_count_ones(Supercube));
	//sf_removeRedundant(Totality, 0);
}

void partital_map_and(pset_family Totality, pset_family temp_Totality, pset Supercube, pset literal_pos, pset q, pset p_and) {
	int var;
	pset p, last;
	int index_p;
	/* lit-1 (x->y')*/
	for (int i = 1; i < SIZE(literal_pos); i++) {
		var = literal_pos[i];
		/* 表示已被supercube移除*/
		if (var == 0xFFFFFFFF)
			continue;
		/* 與所有現有的toatlity 做and*/
		foreachi_set(Totality, index_p, p) {
			//printf("map_and %d %d\n", index_p, literal_pos[i]);
			//printf("Before and: "); display_pset(p, mycube.size);
			//printf("and set:    "); display_pset(q, mvcube.mv_var_size);
			if (set_and_var(p_and, p, q, var)) {
				set_and_var(p_and, p_and, Supercube, var);
				//printf("After and:  "); display_pset(p_and, mycube.size);
				if (notation == 1) {
					if (mvcube_feas_check_notation1(p_and)) {
						sf_addset(temp_Totality, p_and);
					}
				}
				else if (notation == 2) {
					if (mvcube_feas_check_notation2(p_and)) {
						sf_addset(temp_Totality, p_and);
					}
				}
				//else printf("INFEASIBLE\n");
			}

			//else printf("EMPTY\n");

		}
	}
}