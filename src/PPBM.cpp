#include "mainCplusplus.h"
#include "PPBM.h"
//#include "BigInteger.h"

int ppbm_cnt_intersect_pair = 0;
clock_t time_begin_ppbm = clock();
//#define PPBM_detail
void partial_parallel_mapping_block(int f_begin, int f_end, int g_begin, int g_end, pset_family Totality, pset_family temp_Totality, pset_family temp_Totality2, pset Supercube, MCDB& mcdb) {
	bool isRedundant, isEmpty;
	pset p, last, p_lit0, p_lit1, p_and, p_f, Q_non_DC, q_non_DC, temp_Supercube;
	temp_psets t_ps(4);
	t_ps.psets[0] = set_new(mcdb.P->sf_size * 32); p_lit0 = t_ps.psets[0];
	t_ps.psets[1] = set_new(mcdb.P->sf_size * 32); p_lit1 = t_ps.psets[1];
	t_ps.psets[2] = set_new(mvcube.mv_size); p_and = t_ps.psets[2];
	t_ps.psets[3] = set_new(mvcube.mv_size); temp_Supercube = t_ps.psets[3];
	int skip_num = 0, pre_totality_size = 0;


	clock_t smartmv_cost, and_cost, redun_cost;
	parallel_and par_AND(ppbm_thread_num, ppbm_chunk_size);
	parallel_remove_redundant par_RR(ppbm_thread_num, ppbm_chunk_size);

	//fprintf(stderr, "partial_parallel_mapping_block begin~\n");

	/* should del*/
	//paper_MCDB oldMCDB(mvcube.mv_size, mvcube.mv_size);
	
	BigInt cnt_AND_Totality = 0;
	BigInt cnt_AND_Skip = 0;

	clock_t cnt_time = clock();

	float total_pp = (f_end - f_begin) * (g_end - g_begin);
	int total_pp_cnt = 0;
	int cnt_total_and = 0;
	for (int f_on = f_begin; f_on < f_end; f_on++) {
		p_f = GETSET(mcdb.TABLE_PQ, f_on);
		skip_num = 0;
		if (setp_empty(p_f))
			continue;

		for (int g_off = g_begin; g_off < g_end; g_off++) {
			//fprintf(stderr, "(on,off) %d %d %d\n", f_on, g_off, Totality->count);
			total_pp_cnt++;
			if (total_pp_cnt % 10000 == 0) {
				//fprintf(stderr, "%d %d %.2f %%\r", f_on, g_off, (float)total_pp_cnt /total_pp);
			}
			//if (total_pp_cnt % 100000 == 0) {
#ifdef PPBM_detail
			fprintf(stderr, "%d %d %.3lf\n", f_on, g_off, (float)(clock() - cnt_time) / (float)CLOCKS_PER_SEC);
#endif
			//}
			//printf("\n\nPair ( %d %d ) %.3lf\n\n", f_on, g_off, (double)(clock() - time_begin_ppbm) / (double)(CLOCKS_PER_SEC));
			// 檢查是否符合前處理 (ex output相交...)
			if (is_in_set(p_f, g_off) == 0) {
				continue;
			}
			
			temp_Totality->count = 0;
			temp_Totality2->count = 0;
			//printf("ON OFF = (%d, %d)==========================\n", f_on, g_off);

			/* 複製literal postion 變成local*/
			mcdb.lit_pos_copy(p_lit0, GETSET(mcdb.P_lit0_pos, f_on));
			mcdb.lit_pos_copy(p_lit1, GETSET(mcdb.P_lit1_pos, f_on));


			Q_non_DC = GETSET(mcdb.ex_Q_non_DC, g_off);
			q_non_DC = GETSET(mcdb.ex_q_non_DC, g_off);
#ifdef PPBM_detail
			printf("before smart %d %d %d\n", f_on, g_off, Totality->count);
#endif
			/*smartmvcube
			1.將(p_i, q_j)與supercube做檢查，判斷(p_i, q_j)是否可以跳過
			2.同時檢查(p_i, q_j)中的各個pset組合是否feasible
			3.checkredundantbeforeand*/
			cnt_AND_Totality += Totality->count;
			sf_active(Totality);


			isRedundant = false;

			if (parallel_S_MAP(Totality, temp_Totality, Supercube, p_lit0, Q_non_DC, p_and)) {
				continue;
			}

			if (parallel_S_MAP(Totality, temp_Totality, Supercube, p_lit1, q_non_DC, p_and)) {
				continue;
			}
			
			sf_inactive(Totality);
			if (Totality->count == 0) {
				if (temp_Totality->count != 0) {
					my_sf_copy(Totality, temp_Totality);
				}
				continue;
			}
			//and_cost = clock();
			skip_num = temp_Totality->count;
			//skip_num = 0;

			
			cnt_AND_Skip += skip_num;
#ifdef PPBM_detail
			printf("before and %d %d %d\n", f_on, g_off, Totality->count);
			int pre_and_call = par_AND.cnt_and_called;
#endif
			
			
			par_AND.init(Totality, temp_Totality, Supercube, p_lit0, Q_non_DC, p_lit1, q_non_DC);

			par_AND.run(ppbm_option); //1:OMP 2:CPP 3:TBB
			cnt_total_and += par_AND.cnt_and;
#ifdef PPBM_detail
			printf("AND call %d %d %d \n", f_on, g_off, par_AND.cnt_and_called- pre_and_call);

			printf("after and %d %d %d\n", f_on, g_off, temp_Totality->count);
#endif

			if (temp_Totality->count == skip_num) {
				my_sf_copy(Totality, temp_Totality);
				continue;
			}

			if (temp_Totality->count == 0) {
				printf("temp_Totality EMPTY!!!\n");
				exit(EXIT_FAILURE);
			}
			
			//if (temp_Totality->count > 10000) {
				redun_cost = clock();
				par_RR.init(temp_Totality, skip_num);
				par_RR.run(ppbm_option); //1:OMP 2:CPP 3:TBB
			//}
#ifdef PPBM_detail
			printf("after redundant %d %d %d\n", f_on, g_off, temp_Totality->count);
#endif
			sf_copy2(Totality, temp_Totality);
			pre_totality_size = Totality->count;

			sf_or2(temp_Supercube, temp_Totality);
			set_copy(Supercube, temp_Supercube);
			
		}
#ifdef PPBM_detail
		printf("f_on %d: %d %d %d\n", f_on, cnt_total_and, Totality->count, set_count_ones(Supercube));
#endif
	}
}




// 僅能測試cm150a 未完全完成
void partial_parallel_mapping_block_window(int f_begin, int f_end, int g_begin, int g_end, pset_family Totality, pset_family temp_Totality, pset_family temp_Totality2, pset Supercube, MCDB& mcdb, int window) {
	bool isRedundant, isEmpty;
	pset p, last, p_lit0, p_lit1, p_and, p_f, Q_non_DC, q_non_DC, temp_Supercube;
	temp_psets t_ps(4);
	t_ps.psets[0] = set_new(mcdb.P->sf_size * 32); p_lit0 = t_ps.psets[0];
	t_ps.psets[1] = set_new(mcdb.P->sf_size * 32); p_lit1 = t_ps.psets[1];
	t_ps.psets[2] = set_new(mvcube.mv_size); p_and = t_ps.psets[2];
	t_ps.psets[3] = set_new(mvcube.mv_size); temp_Supercube = t_ps.psets[3];
	int skip_num = 0, pre_totality_size = 0;

	/*window setting vvvv*/ // 這裡w只能是2
	pset p2, last2;
	vector<pset> window_fs(window); //window 一次看的onset個數
	pset_family window_skip_sf = sf_new(1, mvcube.mv_size);  //window 相同的部分
	temp_psets window_t_ps(window*2);
	for (int w = 0; w < window*2; w++) {
		window_t_ps.psets[w] = set_new(mcdb.P->sf_size * 32); //奇數w個是 p_lit0 //偶數w個是 p_lit1
	}
	temp_psets window_same_t_ps(2);
	for (int w = 0; w < 2; w++) {
		window_same_t_ps.psets[w] = set_new(mcdb.P->sf_size * 32); //奇數w個是 p_lit0 //偶數w個是 p_lit1
	}
	int cnt_temp_same_plit0 = 0, cnt_temp_same_plit1 = 0;
	int cnt_total_same_plit0 = 0, cnt_total_same_plit1 = 0;

	int var;
	pset_family window_omega = sf_new(1, mvcube.mv_size);
	pset_family window_temp_omega = sf_new(1, mvcube.mv_size);
	/*window setting ^^^^ */

	clock_t smartmv_cost, and_cost, redun_cost;
	parallel_and par_AND(ppbm_thread_num, ppbm_chunk_size);
	parallel_remove_redundant par_RR(ppbm_thread_num, ppbm_chunk_size);

	BigInt cnt_AND_Totality = 0;
	BigInt cnt_AND_Skip = 0;

	clock_t cnt_time = clock();

	float total_pp = (f_end - f_begin) * (g_end - g_begin);
	int total_pp_cnt = 0;

	int cnt_total_and = 0;
	for (int f_on = f_begin; f_on < f_end; f_on++) {

		p_f = GETSET(mcdb.TABLE_PQ, f_on);
		
		if (f_on == 0)
			window = 1;
		else
			window = 2;

		if (f_on + window >= f_end) {
			window = f_end - f_on;
		}

		for (int w = 0; w < window; w++) {
			window_fs[w] = GETSET(mcdb.TABLE_PQ, f_on + w);
		}
		
		// 多個output的話要改!!!!
		skip_num = 0;
		if (setp_empty(p_f))
			continue;

		for (int g_off = g_begin; g_off < g_end; g_off++) {
#ifdef PPBM_detail
			fprintf(stderr, "(on,off) %d %d %d\n", f_on, g_off, Totality->count);
#endif

			
			Q_non_DC = GETSET(mcdb.ex_Q_non_DC, g_off);
			q_non_DC = GETSET(mcdb.ex_q_non_DC, g_off);




			//這裡只比對window = 2
			if (window == 2) {
				for (int w = 0; w < window; w++) {
					mcdb.lit_pos_copy(window_t_ps.psets[w * 2], GETSET(mcdb.P_lit0_pos, f_on + w));
					mcdb.lit_pos_copy(window_t_ps.psets[w * 2 + 1], GETSET(mcdb.P_lit1_pos, f_on + w));
				}
				cnt_temp_same_plit0 = 1;
				cnt_temp_same_plit1 = 1;
#ifdef PPBM_detail
				fprintf(stderr, "(on,off) %d %d %d\n", f_on, g_off, Totality->count);
#endif
				for (int i = 1; i < SIZE(window_t_ps.psets[0]); i++) {

					for (int j = 1; j < SIZE(window_t_ps.psets[2]); j++) {
						if (window_t_ps.psets[2][j] == 0xFFFFFFFF) //已經是相同的pos跳過
							continue;
						if (window_t_ps.psets[0][i] < window_t_ps.psets[2][j])
							break;
						if (window_t_ps.psets[0][i] == window_t_ps.psets[2][j]) {
							//printf(" %d plit0 same\n", window_t_ps.psets[0][i]);
							window_same_t_ps.psets[0][cnt_temp_same_plit0++] = window_t_ps.psets[0][i]; //save same
							
							window_t_ps.psets[0][i] = 0xFFFFFFFF; //remove same
							window_t_ps.psets[2][j] = 0xFFFFFFFF; //remove same
						}
					}
				}// printf("\n");
				
				for (int i = 1; i < SIZE(window_t_ps.psets[1]); i++) {

					for (int j = 1; j < SIZE(window_t_ps.psets[3]); j++) {
						if (window_t_ps.psets[3][j] == 0xFFFFFFFF) //已經是相同的pos跳過
							continue;
						if (window_t_ps.psets[1][i] < window_t_ps.psets[3][j])
							break;
						if (window_t_ps.psets[1][i] == window_t_ps.psets[3][j]) {
							//printf(" %d plit1 same\n", window_t_ps.psets[1][i]);
							window_same_t_ps.psets[1][cnt_temp_same_plit1++] = window_t_ps.psets[1][i]; //save same

							window_t_ps.psets[1][i] = 0xFFFFFFFF; //remove same
							window_t_ps.psets[3][j] = 0xFFFFFFFF; //remove same
						}
					}
				} //printf("\n");

				PUTSIZE(window_same_t_ps.psets[0], cnt_temp_same_plit0);
				PUTSIZE(window_same_t_ps.psets[1], cnt_temp_same_plit1);

	
#ifdef PPBM_detail
				printf("cnt lit  %d %d %d %d %d\n", f_on, g_off, SIZE(window_t_ps.psets[0]) + SIZE(window_t_ps.psets[1]) - (SIZE(window_same_t_ps.psets[0]) + SIZE(window_same_t_ps.psets[1])), SIZE(window_t_ps.psets[2]) + SIZE(window_t_ps.psets[3]) - (SIZE(window_same_t_ps.psets[0]) + SIZE(window_same_t_ps.psets[1])), SIZE(window_same_t_ps.psets[0]) + SIZE(window_same_t_ps.psets[1])-2);
#endif
				// 到這裡onset的部分有三個部分 f_0 不相同, f_1 不相同, f_0與f_1相同

				temp_Totality->count = 0;
				sf_active(Totality);
				// f_0與f_1相同 可以用原本的方式進行檢查
				if (parallel_S_MAP(Totality, temp_Totality, Supercube, window_same_t_ps.psets[0], Q_non_DC, p_and, true)) {
					continue;
				}
				if (parallel_S_MAP(Totality, temp_Totality, Supercube, window_same_t_ps.psets[1], q_non_DC, p_and, true)) {
					continue;
				}

				skip_num = temp_Totality->count;
				//skip_num = 0;
#ifdef PPBM_detail
				printf("window totality OAND (%d %d)\n", Totality->count, temp_Totality->count);
				
#endif
				// f_0 不相同 部分檢查不可使用
				if (parallel_S_MAP(Totality, temp_Totality, Supercube, window_t_ps.psets[0], Q_non_DC, p_and, true)) {
					continue;
				}
				if (parallel_S_MAP(Totality, temp_Totality, Supercube, window_t_ps.psets[1], q_non_DC, p_and, true)) {
					continue;
				}
				// f_1 不相同 部分檢查不可使用
				if (parallel_S_MAP(Totality, temp_Totality, Supercube, window_t_ps.psets[2], Q_non_DC, p_and, true)) {
					continue;
				}
				if (parallel_S_MAP(Totality, temp_Totality, Supercube, window_t_ps.psets[3], q_non_DC, p_and, true)) {
					continue;
				}

				sf_inactive(Totality);
				if (Totality->count == 0) {
					if (temp_Totality->count != 0) {
						my_sf_copy(Totality, temp_Totality);
					}
					continue;
				}



				window_omega->count = 0;
				window_temp_omega->count = 0;


				// 將lit 0 f_0 轉成mvcubes
				for (int i = 1; i < SIZE(window_t_ps.psets[0]); i++) {
					var = window_t_ps.psets[0][i];
					if (var == 0xFFFFFFFF)
						continue;
					if (set_and_var(p_and, Supercube, Q_non_DC, var)) {
						sf_addset(window_omega, p_and);
					}
				}
				// 將lit 1 f_0 轉成mvcubes
				for (int i = 1; i < SIZE(window_t_ps.psets[1]); i++) {
					var = window_t_ps.psets[1][i];
					if (var == 0xFFFFFFFF)
						continue;
					if (set_and_var(p_and, Supercube, q_non_DC, var)) {
						sf_addset(window_omega, p_and);
					}
				}

				// TEST vvvvvvvvv
				int cnt_f_1 = 0;
				// 將lit 0 f_0 轉成mvcubes
				for (int i = 1; i < SIZE(window_t_ps.psets[2]); i++) {
					var = window_t_ps.psets[2][i];
					if (var == 0xFFFFFFFF)
						continue;
					cnt_f_1++;
				}
				// 將lit 1 f_0 轉成mvcubes
				for (int i = 1; i < SIZE(window_t_ps.psets[3]); i++) {
					var = window_t_ps.psets[3][i];
					if (var == 0xFFFFFFFF)
						continue;
					cnt_f_1++;
				}
#ifdef PPBM_detail
				printf("cnt_f_1 %d\n", cnt_f_1);
#endif
				// TEST ^^^^^^^^^^^




				// 正常應該還要檢查是否window_t_ps.psets[2]與window_t_ps.psets[3]相加個數是否為0!!!!!!!!!!! 這裡略過
				// f_0 與 f_1 做AND
				
				par_AND.init(window_omega, window_temp_omega, Supercube, window_t_ps.psets[2], Q_non_DC, window_t_ps.psets[3], q_non_DC);
				par_AND.run(ppbm_option); //1:OMP 2:CPP 3:TBB

				cnt_total_and += par_AND.cnt_and;
#ifdef PPBM_detail
				printf("window omega before added same ( %d %d )\n", window_omega->count, window_temp_omega->count);
#endif
				// 將f_same 加入omega中
				for (int i = 1; i < SIZE(window_same_t_ps.psets[0]); i++) {
					var = window_same_t_ps.psets[0][i];
					if (var == 0xFFFFFFFF)
						continue;
					if (set_and_var(p_and, Supercube, Q_non_DC, var)) {
						sf_addset(window_temp_omega, p_and);
					}
				}
				for (int i = 1; i < SIZE(window_same_t_ps.psets[1]); i++) {
					var = window_same_t_ps.psets[1][i];
					if (var == 0xFFFFFFFF)
						continue;
					if (set_and_var(p_and, Supercube, q_non_DC, var)) {
						sf_addset(window_temp_omega, p_and);
					}
				}
#ifdef PPBM_detail
				printf("window omega after added same ( %d %d )\n", window_omega->count, window_temp_omega->count);
#endif
				//printf("window totality AND %d %d %d %d \n", f_on, g_off, Totality->count, window_temp_omega->count);
				// Totality 與omega 做AND
				int cnt_and = 0;
				foreach_set(Totality, last, p) {
					foreach_set(window_temp_omega, last2, p2) {
						cnt_and++;
						set_and(p_and, p, p2);
						PUTSIZE(p_and, set_count_ones(p_and));
						if (notation == 1) {
							printf("UNFINISH !!!!!!!!\n");
							/*if (mvcube_feas_check_notation1(p_and)) {
								sf_addset(temp_Totality, p_and);
							}*/
						}
						else if (notation == 2) {
							if (mvcube_feas_check_notation2(p_and)) {
								sf_addset(temp_Totality, p_and);
							}
						}
					}
				}
#ifdef PPBM_detail
				printf("window totality and omega %d\n", cnt_and);
				printf("window totality ( %d %d )\n", Totality->count, temp_Totality->count);
#endif
				cnt_total_and += cnt_and;
				

				if (temp_Totality->count == 0) {
					printf("window temp totality ERROR!!!\n");
					exit(EXIT_FAILURE);
				}

				par_RR.init(temp_Totality, skip_num);
				par_RR.run(ppbm_option); //1:OMP 2:CPP 3:TBB
				//printf("window totality RR ( %d %d )\n", Totality->count, temp_Totality->count);
				sf_copy2(Totality, temp_Totality);
				pre_totality_size = Totality->count;

				sf_or2(temp_Supercube, temp_Totality);
				set_copy(Supercube, temp_Supercube);
#ifdef PPBM_detail
				printf("f_on %d: %d %d %d\n", f_on, cnt_total_and, Totality->count, set_count_ones(Supercube));
#endif
				continue;
			}
			

			/* 複製literal postion 變成local*/
			mcdb.lit_pos_copy(p_lit0, GETSET(mcdb.P_lit0_pos, f_on));
			mcdb.lit_pos_copy(p_lit1, GETSET(mcdb.P_lit1_pos, f_on));

			total_pp_cnt++;
			if (total_pp_cnt % 10000 == 0) {
				//fprintf(stderr, "%d %d %.2f %%\r", f_on, g_off, (float)total_pp_cnt /total_pp);
			}
			//if (total_pp_cnt % 100000 == 0) {
#ifdef PPBM_detail
			fprintf(stderr, "%d %d %.3lf\n", f_on, g_off, (float)(clock() - cnt_time) / (float)CLOCKS_PER_SEC);
#endif
			//}
			//printf("\n\nPair ( %d %d ) %.3lf\n\n", f_on, g_off, (double)(clock() - time_begin_ppbm) / (double)(CLOCKS_PER_SEC));
			// 檢查是否符合前處理 (ex output相交...)
			if (is_in_set(p_f, g_off) == 0) {
				continue;
			}

			temp_Totality->count = 0;
			temp_Totality2->count = 0;
			//printf("ON OFF = (%d, %d)==========================\n", f_on, g_off);

#ifdef PPBM_detail
			printf("before smart %d %d %d\n", f_on, g_off, Totality->count);
#endif
			/*smartmvcube
			1.將(p_i, q_j)與supercube做檢查，判斷(p_i, q_j)是否可以跳過
			2.同時檢查(p_i, q_j)中的各個pset組合是否feasible
			3.checkredundantbeforeand*/
			cnt_AND_Totality += Totality->count;
			sf_active(Totality);


			isRedundant = false;
			if (parallel_S_MAP(Totality, temp_Totality, Supercube, p_lit0, Q_non_DC, p_and)) {
				continue;
			}

			if (parallel_S_MAP(Totality, temp_Totality, Supercube, p_lit1, q_non_DC, p_and)) {
				continue;
			}

			sf_inactive(Totality);
			if (Totality->count == 0) {
				if (temp_Totality->count != 0) {
					my_sf_copy(Totality, temp_Totality);
				}
				
				continue;
			}
			//and_cost = clock();
			skip_num = temp_Totality->count;
			//skip_num = 0;


			cnt_AND_Skip += skip_num;
#ifdef PPBM_detail
			printf("before and %d %d %d\n", f_on, g_off, Totality->count);
#endif
			
			par_AND.init(Totality, temp_Totality, Supercube, p_lit0, Q_non_DC, p_lit1, q_non_DC);
			par_AND.run(ppbm_option); //1:OMP 2:CPP 3:TBB
			cnt_total_and += par_AND.cnt_and;
#ifdef PPBM_detail
			printf("after and %d %d %d\n", f_on, g_off, temp_Totality->count);
#endif
			

			if (temp_Totality->count == skip_num) {
				my_sf_copy(Totality, temp_Totality);
				continue;
			}

			/*foreach_set(temp_Totality, last, p) {
				printf("%d ", SIZE(p));
			}printf("\n\n\n");*/

			
			if (temp_Totality->count == 0) {
				printf("temp_Totality EMPTY!!!\n");
				exit(EXIT_FAILURE);
			}

			//if (temp_Totality->count > 10000) {
			redun_cost = clock();
			par_RR.init(temp_Totality, 0);
			par_RR.run(ppbm_option); //1:OMP 2:CPP 3:TBB
		//}
#ifdef PPBM_detail
			printf("after redundant %d %d %d\n", f_on, g_off, temp_Totality->count);
#endif
			//parallel_sf_removeRedundant(temp_Totality);
			//sf_removeRedundant(temp_Totality);
			sf_copy2(Totality, temp_Totality);
			pre_totality_size = Totality->count;

			sf_or2(temp_Supercube, temp_Totality);
			set_copy(Supercube, temp_Supercube);
		}
		// window
		f_on+= window-1;
#ifdef PPBM_detail
		printf("f_on %d: %d %d %d\n", f_on, cnt_total_and, Totality->count, set_count_ones(Supercube));
#endif
	}
}


/* sf_sort -- sort the sets of A */
pset* parallel_sf_sort(pset_family A, qsort_compare_func compare)
{
	pset p, last, * pdest, * A1;

	/* Create a single array pointing to each cube of A */
	pdest = A1 = ALLOC(pset, A->count + 1);
	foreach_set(A, last, p) {
		PUTSIZE(p, set_ord(p));         /* compute the set size */
		*pdest++ = p;                   /* save the pointer */
	}
	*pdest = NULL;                      /* Sentinel -- never seen by sort */

	/* Sort cubes by size */
	qsort((char*)A1, A->count, sizeof(pset), compare);
	return A1;
}

/* sf_unlist -- make a set family out of a list of pointers to sets */
pset_family parallel_sf_unlist(pset_family R, pset* A1)
{
	pset pr, p, * pa;
	R->count = R->capacity;
	for (pr = R->data, pa = A1; (p = *pa++) != NULL; pr += R->wsize)
		INLINEset_copy(pr, p);
	FREE(A1);
	return R;
}


void parallel_removeRedundant_work(pset_family A, vector<atomic<bool>>& isRedundants, atomic<int>& chunk_id, int chunk_size, int chunk_cnt) {
	pset p, p2;
	int i, j;
	int id, end_i;
	while ((id = chunk_id.fetch_add(1)) < chunk_cnt) {
		end_i = (id == chunk_cnt-1) ? A->count : id * chunk_size + chunk_size;
		for (p = GETSET(A, id * chunk_size), i = id * chunk_size; i < end_i; p += A->wsize, i++) {
			if (isRedundants[i].load())
				continue;
			for (p2 = GETSET(A, i + 1), j = i + 1; j < A->count; p2 += A->wsize, j++) {
				if (isRedundants[j].load())
					continue;
				if (setp_implies(p2, p))
					isRedundants[j].store(true);
				else if (setp_implies(p, p2))
					isRedundants[i].store(true);
			}
		}
	}
}



void parallel_sf_removeRedundant(pset_family A) {
	int thread_num = 24;
	int chunk_size = 50;
	int chunk_cnt = ceil((double)A->count / (double)chunk_size);
	atomic<int> chunk_id(0);
	vector<atomic<bool>> isRedundants(A->count);
	for (atomic<bool>& isRedundant : isRedundants) {
		isRedundant.store(false);
	}
	vector<thread> threads;
	for (int i = 0; i < thread_num; i++) {
		if (chunk_id.load() >= chunk_cnt)
			break;
		threads.emplace_back(parallel_removeRedundant_work, ref(A), ref(isRedundants), ref(chunk_id), chunk_size, chunk_cnt);
	}

	for (int i = 0; i < threads.size(); i++) {
		threads[i].join();
	}

	pset p, p2, last, pdest;
	int i = 0;
	pdest = A->data;
	foreach_set(A, last, p) {
		if (!isRedundants[i]) {
			if (pdest != p) {
				INLINEset_copy(pdest, p);
			}
			pdest += A->wsize;
		}
		else {
			A->count--;
		}
		i++;
	}


	//pset_family_mutex sf_mutex(A);

	// 如果有被包含則移除
}


/* sf_active -- make all members of the set family active */
pset_family test_sf_active(pset_family A)
{
	pset p, last;
	foreach_set(A, last, p) {
		SET(p, ACTIVE);
	}
	A->active_count = A->count;
	return A;
}




/* sf_inactive -- remove all inactive cubes in a set family */
pset_family test_sf_inactive(pset_family A)
{
	pset p, last, pdest;

	pdest = A->data;
	foreach_set(A, last, p) {
		if (TESTP(p, ACTIVE)) {
			if (pdest != p) {
				INLINEset_copy(pdest, p);
			}
			pdest += A->wsize;
		}
		else {
			A->count--;
		}
	}
	return A;
}

//#define S_MAP_debug
bool parallel_S_MAP(pset_family Totality, pset_family temp_Totality, pset Supercube, pset literal_pos, pset Q_non_DC, pset p_and, bool window) {
	pset p, last;
	int var;
	for (int i = 1; i < SIZE(literal_pos); i++) {
		//printf("\n\n\n");
		//printf(" %d", literal_pos[i]);
		var = literal_pos[i];
		if (var == 0xFFFFFFFF)
			continue;
		//display_pset(Supercube, mvcube.mv_size);
		//display_pset(p_g, mvcube.mv_size);
		if (!window && setp_implies_var(Supercube, Q_non_DC, var)) {
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
		if (!window) {
			foreach_active_set(Totality, last, p) {
				if (setp_implies_var(p, p_and, var)) {
					sf_addset(temp_Totality, p);
					RESET(p, ACTIVE);
				}
			}
		}
	}//printf("\n");
	return false;
}