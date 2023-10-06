/*
	Parallel Boolean Matching
	全部平行化的布林比對 (CPP Only)


	1. 將onset 切塊
	2. 同時執行n 塊的onset給n個thread執行 Boolean Matching的計算，得到n個計算結果 (計算過程中supercube為全域變數 所有執行緒共用 故只會越來越好)
	3. 等待所有執行緒執行完畢 合併n個計算結果
	4. 合併過程分為三步 (每個步驟皆有CPP平行化)
		a. 分別用supercube檢查n個計算結果 ，因supercube是在逐漸變好 可能有些較早算的未使用到supercube的好處
		b. 將n個計算結果互相檢查是否有重複多餘的並移除
		c. 將n個計算結果中的所有內容互相做AND來合併
*/

#ifndef PBM_H
#define PBM_H

#include <iostream>
#include <mutex>
#include <thread>
#include <atomic>
#include <tbb\parallel_reduce.h>
#include "parallel_tool.h"
#include "PPBM.h"

using namespace std;



class Parallel_BM {
private:
	/* 基本的參數*/
	int thread_num, chunk_size, chunk_cnt;
	pset_family Totality, temp_Totality;
	pset Supercube;
	chunk chs;
	MCDB mcdb;

	/* 平行化需要的參數*/
	/* main*/
	atomic<int> chunk_id;
	shared_mutex sc_mtx; // Supercube mutex
	mutex totality_mtx;
	vector<pset_family> Totalitys, temp_Totalitys;
	vector<temp_psets> t_pss;
	Barrier BAR;
	



	/* and*/
	vector<pset_family> sf_vars, sf_sort;

	/* redundant*/
	
	/* combine*/
	int combine_chunk_size, combine_chunk_cnt_step1A, combine_chunk_cnt_step1B, combine_chunk_cnt_step2, combine_chunk_size_step2, combine_chunk_cnt_step3, combine_chunk_size_step3;

public:
	Parallel_BM(int thread_num_, int chunk_size_) :thread_num(thread_num_), chunk_size(chunk_size_), chunk_cnt(0) {
		/* main*/
		Totalitys.resize(thread_num);
		temp_Totalitys.resize(thread_num);
		t_pss.resize(thread_num);
		/* and*/
		sf_vars.resize(thread_num);
		sf_sort.resize(thread_num);
		
		temp_Totality = sf_new(1, mvcube.mv_size);

		for (int i = 0; i < thread_num; i++) {
			/* main*/
			Totalitys[i] = sf_new(1, mvcube.mv_size);
			temp_Totalitys[i] = sf_new(1, mvcube.mv_size);
			t_pss[i].psets.resize(8);
			t_pss[i].psets[0] = set_new(mvcube.num_mv_var*2 * 32); //p_lit0
			t_pss[i].psets[1] = set_new(mvcube.num_mv_var*2 * 32); //p_lit1
			t_pss[i].psets[2] = set_new(mvcube.mv_size); //p_and
			t_pss[i].psets[3] = set_new(mvcube.mv_size); //temp_Supercube
			/* and*/
			sf_vars[i] = sf_new(mvcube.num_mv_var, mvcube.mv_var_size);
			sf_sort[i] = sf_new(mvcube.num_mv_var, mvcube.mv_var_size);
			t_pss[i].psets[4] = set_new(mvcube.mv_var_size); //p_var_t
			t_pss[i].psets[5] = set_new(mvcube.mv_var_size); //var_supercube
			t_pss[i].psets[6] = set_new(mvcube.mv_var_size); //solution
			/* merge*/
			t_pss[i].psets[7] = set_new(mvcube.mv_size); //pset_combine

		}
		

		chunk_id.store(0);


		BAR.setup(thread_num);
	}

	~Parallel_BM() {
		sf_free(temp_Totality);
		for (int i = 0; i < thread_num; i++) {
			sf_free(Totalitys[i]);
			sf_free(temp_Totalitys[i]);
			sf_free(sf_vars[i]);
			sf_free(sf_sort[i]);
		}
		//printf("Parallel_BM deconstruct finish\n");
	}

	void setup(pset_family F, pset_family G, pset_family Totality_, pset Supercube_, int key) {
		Totality = Totality_;
		Supercube = Supercube_;

		mcdb.setup(F, G, mvcube.num_mv_var, mvcube.is_input_phase_assigment);
		mcdb.pre_run(mvcube.first_word, mvcube.last_word, mvcube.first_bit, mvcube.last_bit);
		mcdb.check_Rule(key); /* 0bDCBA  A:R1 B:R2 C:R3 D:R4*/

		//if (row_size == 0) row_size = 1; else row_size = mcdb.TABLE_PQ->count / row_size;
		//if (col_size == 0) col_size = 1; else col_size = mcdb.TABLE_PQ->sf_size / col_size;

		chs.setup(MAX(row_num, 1), MAX(col_num, 1), mcdb.TABLE_PQ->count, mcdb.TABLE_PQ->sf_size);
		chs.order1();

		chunk_cnt = chs.chunk_cnt;
	}

	void run() {

		printf("chunk cnt %d\n", chunk_cnt);
		run_cpp();
	}

	void run_cpp() {
		
		pset_family sf_temp = sf_new(1, mvcube.mv_size);
		int id = 0;

		parallel_remove_redundant par_RR(ppbm_thread_num, ppbm_chunk_size);
		while (id < chunk_cnt) {
			for (int i = 0; i < thread_num; i++) {
				Totalitys[i]->count = 0;
				temp_Totalitys[i]->count = 0;
			}


			for (int cnt = 0;cnt < thread_num && id < chunk_cnt; id++, cnt++) {
				pool.task_begin();
				pool.enqueue(&Parallel_BM::work2, this, id);
			}

			pool.wait();

			clock_t c = clock();
			vector<pair<int, int>> Totality_sizes(thread_num);
			for (int i = 0; i < thread_num; i++) {
				Totality_sizes[i] = make_pair(i, Totalitys[i]->count);
			}
			sort(Totality_sizes.begin(), Totality_sizes.end(), [](const pair<int, int>& a, const pair<int, int>& b) {
				return a.second < b.second;
				});

			int Totalitys_id;
			
			Totality->count = 0;

			for (int i = 0; i < thread_num; i++) {
				Totalitys_id = Totality_sizes[i].first;
				if (Totalitys[Totalitys_id]->count != 0) {
					if (Totality->count == 0) {
						sf_copy3(Totality, Totalitys[Totalitys_id]);
						fprintf(stderr, "skip parallel_combine\n");
					}
					else {
						sf_copy3(sf_temp, Totality);
						Totality->count = 0;
						//combine(Totality, sf_temp, Totalitys[i]); // 單執行緒的merge
						fprintf(stderr, "parallel_combine start\n");
						parallel_combine(Totality, sf_temp, Totalitys[Totalitys_id]);
						fprintf(stderr, "parallel_combine end\n");
						printf("(%d %d %d %d)\n", Totalitys_id, sf_temp->count, Totalitys[Totalitys_id]->count, Totality->count);
					}
					fprintf(stderr, "reductipon start\n");

					//sf_removeRedundant(Totality, 0);
					par_RR.init(Totality, 0);
					par_RR.run(ppbm_option); //1:OMP 2:CPP 3:TBB

					fprintf(stderr, "reductipon end\n");
					sf_or2(t_pss[0].psets[3], Totality);
					{
						lock_guard<shared_mutex> lock(sc_mtx);
						set_and(Supercube, Supercube, t_pss[0].psets[3]);
					}
				}
			}
			printf("merge cost: %.6f\n", (float)(clock() - c) / (float)CLOCKS_PER_SEC);
			//printf("combine: %d\n", Totality->count);
		}
		
		sf_free(sf_temp);


		
	}
	
	/* 用barrier*/
	void work(int thread_id) {
		int id;
		bool isRedundant, isEmpty;
		pset p, last, p_lit0, p_lit1, p_and, p_f, Q_non_DC, q_non_DC, temp_Supercube, p_get_var;

		p_lit0 = t_pss[thread_id].psets[0];
		p_lit1 = t_pss[thread_id].psets[1];
		p_and = t_pss[thread_id].psets[2];
		temp_Supercube = t_pss[thread_id].psets[3];
		
		
		sf_copy3(Totalitys[thread_id], Totality);

		
		while ((id = chunk_id.fetch_add(1)) < chunk_cnt) {
			printf("thread_id %d: %d %d\n", thread_id, id, Totalitys[thread_id]->count);
			for (int f_on = chs.chunks[id].f_begin; f_on < chs.chunks[id].f_end; f_on++) {
				p_f = GETSET(mcdb.TABLE_PQ, f_on);

				if (setp_empty(p_f))
					continue;
				for (int g_off = chs.chunks[id].g_begin; g_off < chs.chunks[id].g_end; g_off++) {
					// 檢查是否符合前處理 (ex output相交...)
					if (is_in_set(p_f, g_off) == 0) {
						continue;
					}
					
					temp_Totalitys[thread_id]->count = 0;

					/* 複製literal postion 變成local*/
					mcdb.lit_pos_copy(p_lit0, GETSET(mcdb.P_lit0_pos, f_on));
					mcdb.lit_pos_copy(p_lit1, GETSET(mcdb.P_lit1_pos, f_on));

					Q_non_DC = GETSET(mcdb.ex_Q_non_DC, g_off);
					q_non_DC = GETSET(mcdb.ex_q_non_DC, g_off);


					/*smartmvcube
					1.將(p_i, q_j)與supercube做檢查，判斷(p_i, q_j)是否可以跳過
					2.同時檢查(p_i, q_j)中的各個pset組合是否feasible
					3.checkredundantbeforeand*/
					sf_active(Totalitys[thread_id]);

					{
						lock_guard<shared_mutex> lock(sc_mtx);
						set_copy(temp_Supercube, Supercube);
					}
					isRedundant = false;
					if (S_MAP(Totalitys[thread_id], temp_Totalitys[thread_id], temp_Supercube, p_lit0, Q_non_DC, p_and)) {
						continue;
					}
					if (S_MAP(Totalitys[thread_id], temp_Totalitys[thread_id], temp_Supercube, p_lit1, q_non_DC, p_and)) {
						continue;
					}
					sf_inactive(Totalitys[thread_id]);

					if (Totalitys[thread_id]->count == 0) {
						if (temp_Totalitys[thread_id]->count != 0) {
							my_sf_copy(Totalitys[thread_id], temp_Totalitys[thread_id]);
						}
						continue;
					}
					sf_and_mvcubes(thread_id, Totalitys[thread_id], temp_Totalitys[thread_id], p_lit0, Q_non_DC, p_lit1, q_non_DC, p_and);

					if (temp_Totalitys[thread_id]->count == 0) {
						printf("temp_Totality EMPTY!!!\n");
						exit(EXIT_FAILURE);
					}

					printf("%d\n", temp_Totalitys[thread_id]->count);

					sf_removeRedundant(temp_Totalitys[thread_id], 0);
					sf_copy2(Totalitys[thread_id], temp_Totalitys[thread_id]);


					sf_or2(temp_Supercube, temp_Totalitys[thread_id]);
					{
						lock_guard<shared_mutex> lock(sc_mtx);
						set_copy(Supercube, temp_Supercube);
					}
				}
			}
			BAR.Wait();
		}
		BAR.Wait();
	}

	/* thread pool*/
	void work2(int id) {
		clock_t c = clock();
		int thread_id = pool.get_thread_ids(this_thread::get_id());
		bool isRedundant, isEmpty;
		pset p, last, p_lit0, p_lit1, p_and, p_f, Q_non_DC, q_non_DC, temp_Supercube, p_get_var;

		p_lit0 = t_pss[thread_id].psets[0];
		p_lit1 = t_pss[thread_id].psets[1];
		p_and = t_pss[thread_id].psets[2];
		temp_Supercube = t_pss[thread_id].psets[3];
		int skip_num = 0, pre_totality_size = 0;
		//if (Totalitys[thread_id]->count == 0)
			sf_copy3(Totalitys[thread_id], Totality);


		printf("thread_id: %d %d\n",thread_id, set_count_ones(Supercube));
		//printf("thread_id %d: %d %d\n", thread_id, id, Totalitys[thread_id]->count);
		for (int f_on = chs.chunks[id].f_begin; f_on < chs.chunks[id].f_end; f_on++) {
			p_f = GETSET(mcdb.TABLE_PQ, f_on);

			if (setp_empty(p_f))
				continue;
			for (int g_off = chs.chunks[id].g_begin; g_off < chs.chunks[id].g_end; g_off++) {
				// 檢查是否符合前處理 (ex output相交...)
				if (is_in_set(p_f, g_off) == 0) {
					continue;
				}

				temp_Totalitys[thread_id]->count = 0;

				/* 複製literal postion 變成local*/
				mcdb.lit_pos_copy(p_lit0, GETSET(mcdb.P_lit0_pos, f_on));
				mcdb.lit_pos_copy(p_lit1, GETSET(mcdb.P_lit1_pos, f_on));

				Q_non_DC = GETSET(mcdb.ex_Q_non_DC, g_off);
				q_non_DC = GETSET(mcdb.ex_q_non_DC, g_off);


				/*smartmvcube
				1.將(p_i, q_j)與supercube做檢查，判斷(p_i, q_j)是否可以跳過
				2.同時檢查(p_i, q_j)中的各個pset組合是否feasible
				3.checkredundantbeforeand*/
				sf_active(Totalitys[thread_id]);

				{
					lock_guard<shared_mutex> lock(sc_mtx);
					set_copy(temp_Supercube, Supercube);
				}


				isRedundant = false;
				if (S_MAP(Totalitys[thread_id], temp_Totalitys[thread_id], temp_Supercube, p_lit0, Q_non_DC, p_and)) {
					continue;
				}
				if (S_MAP(Totalitys[thread_id], temp_Totalitys[thread_id], temp_Supercube, p_lit1, q_non_DC, p_and)) {
					continue;
				}

				sf_inactive(Totalitys[thread_id]);

				if (Totalitys[thread_id]->count == 0) {
					if (temp_Totalitys[thread_id]->count != 0) {
						my_sf_copy(Totalitys[thread_id], temp_Totalitys[thread_id]);
					}
					continue;
				}

				skip_num = temp_Totalitys[thread_id]->count;

				sf_and_mvcubes(thread_id, Totalitys[thread_id], temp_Totalitys[thread_id], p_lit0, Q_non_DC, p_lit1, q_non_DC, p_and);

				if (temp_Totalitys[thread_id]->count == skip_num) {
					my_sf_copy(Totalitys[thread_id], temp_Totalitys[thread_id]);
					continue;
				}


				if (temp_Totalitys[thread_id]->count == 0) {
					printf("temp_Totality EMPTY!!!\n");
					exit(EXIT_FAILURE);
				}



				printf("%d %d\n", thread_id, temp_Totalitys[thread_id]->count);
				
				if (pre_totality_size < temp_Totalitys[thread_id]->count) {
					sf_removeRedundant(temp_Totalitys[thread_id], skip_num);
				}

				
				sf_copy2(Totalitys[thread_id], temp_Totalitys[thread_id]);
				pre_totality_size = Totalitys[thread_id]->count;

				sf_or2(temp_Supercube, temp_Totalitys[thread_id]);
				{
					lock_guard<shared_mutex> lock(sc_mtx);
					set_and(Supercube, Supercube, temp_Supercube);
					//printf("thread_id: %d %d\n",thread_id, set_count_ones(Supercube));
				}
			}
		}
		//printf("merge cose: %.6f\n", (float)(c - clock()) / (float)CLOCKS_PER_SEC);
		printf("IN WORK2 %d %d %d %d %.6f\n", thread_id, id, Totalitys[thread_id]->count, set_count_ones(temp_Supercube), (float)(clock() - c) / (float)CLOCKS_PER_SEC);
		pool.task_end();

	}

	/* thread pool 大量lock*/
	void work3(int id) {
		int thread_id = pool.get_thread_ids(this_thread::get_id());
		bool isRedundant, isEmpty;
		pset p, last, p_lit0, p_lit1, p_and, p_f, Q_non_DC, q_non_DC, temp_Supercube, p_get_var;

		p_lit0 = t_pss[thread_id].psets[0];
		p_lit1 = t_pss[thread_id].psets[1];
		p_and = t_pss[thread_id].psets[2];
		//temp_Supercube = t_pss[thread_id].psets[3];

		

		//printf("thread_id %d: %d %d\n", thread_id, id, Totalitys[thread_id]->count);
		for (int f_on = chs.chunks[id].f_begin; f_on < chs.chunks[id].f_end; f_on++) {
			p_f = GETSET(mcdb.TABLE_PQ, f_on);

			if (setp_empty(p_f))
				continue;
			for (int g_off = chs.chunks[id].g_begin; g_off < chs.chunks[id].g_end; g_off++) {
				// 檢查是否符合前處理 (ex output相交...)
				if (is_in_set(p_f, g_off) == 0) {
					continue;
				}

				

				/* 複製literal postion 變成local*/
				mcdb.lit_pos_copy(p_lit0, GETSET(mcdb.P_lit0_pos, f_on));
				mcdb.lit_pos_copy(p_lit1, GETSET(mcdb.P_lit1_pos, f_on));

				Q_non_DC = GETSET(mcdb.ex_Q_non_DC, g_off);
				q_non_DC = GETSET(mcdb.ex_q_non_DC, g_off);


				/*smartmvcube
				1.將(p_i, q_j)與supercube做檢查，判斷(p_i, q_j)是否可以跳過
				2.同時檢查(p_i, q_j)中的各個pset組合是否feasible
				3.checkredundantbeforeand*/

				{
					lock_guard<mutex> lk(totality_mtx);
					temp_Totality->count = 0;
					sf_active(Totality);
					isRedundant = false;
					if (S_MAP(Totality, temp_Totality, Supercube, p_lit0, Q_non_DC, p_and)) {
						continue;
					}
					if (S_MAP(Totality, temp_Totality, Supercube, p_lit1, q_non_DC, p_and)) {
						continue;
					}
					
					sf_inactive(Totality);

					if (Totality->count == 0) {
						if (temp_Totality->count != 0) {
							my_sf_copy(Totality, temp_Totality);
						}
						continue;
					}
					sf_and_mvcubes(thread_id, Totality, temp_Totality, p_lit0, Q_non_DC, p_lit1, q_non_DC, p_and);
					
					//printf("%d %d\n", thread_id, temp_Totality->count);
					
					if (temp_Totality->count == 0) {
						printf("temp_Totality EMPTY!!!\n");
						exit(EXIT_FAILURE);
					}

					



					sf_removeRedundant(temp_Totality, 0);
					sf_copy2(Totality, temp_Totality);
					sf_or2(Supercube, temp_Totality);

					/*sf_or2(temp_Supercube, temp_Totality);
					{
						lock_guard<shared_mutex> lock(sc_mtx);
						set_and(Supercube, Supercube, temp_Supercube);
						printf("%d\n", set_count_ones(Supercube));
					}*/
				}
			}
		}
		pool.task_end();

	}

	//void barrier(std::atomic<int>& counter, int num_threads) {
	//	counter.fetch_add(1);
	//	while (counter < num_threads) { // Wait for all threads to reach the barrier
	//		std::this_thread::yield(); // Yield the CPU to avoid busy waiting
	//	}
	//}

	void sf_and_mvcubes(int thread_id, pset_family and_Totality, pset_family and_temp_totality, pset literal_0_pos, pset Q, pset literal_1_pos, pset q, pset p_and) {
		pset p, last, using_q, p_var_t, var_supercube, solution;
		int var, id, chunk_end;
		
		/* 儲存(var)temp的結果*/
		p_var_t = t_pss[thread_id].psets[4];
		/* var 的supercube，*/
		var_supercube = t_pss[thread_id].psets[5];
		/* 儲存(var)對應到已知解的結果*/
		solution = t_pss[thread_id].psets[6];

		foreach_set(and_Totality, last, p) {
			for (int i = 1; i < SIZE(literal_0_pos) + SIZE(literal_1_pos) - 1; i++) {
				if (i < SIZE(literal_0_pos)) {
					var = literal_0_pos[i];
					using_q = Q;
				}
				else {
					var = literal_1_pos[1 + i - SIZE(literal_0_pos)];
					using_q = q;
				}
				/* 表示已被supercube移除*/
				if (var == 0xFFFFFFFF)
					continue;

				/* 與所有現有的toatlity 做and*/
				if (set_and_var(p_and, p, using_q, var)) {
					if (notation == 1) {
						printf("UNFINISH !!!!!!!!\n");
						/*if (mvcube_feas_check_notation1(p_and)) {
							sf_addset(temp_Totality, p_and);
						}*/
					}
					else if (notation == 2) {
						if (parallel_mvcube_feas_check_notation2(p_and, sf_vars[thread_id], sf_sort[thread_id], p_var_t, var_supercube, solution)) {
							sf_addset(and_temp_totality, p_and);
						}
					}
				}
			}
		}
	}

	/*
		merge 2 pset_family
		step1: 先用supercube分別對A和B進行過濾
		step2: 若有包含關係則跳過merge
		step3: 不相同則兩兩相乘做AND

		//sf A B C 皆沒有free
		*********** AND 和 merge無法同時執行，有部分暫存變數重疊*************
	*/
	void combine(pset_family sf_C, pset_family sf_A, pset_family sf_B) {
		
		sf_active(sf_A);
		sf_active(sf_B);
		pset p1, p2, last1, last2;
		int thread_id = 0;
		pset combine = t_pss[thread_id].psets[7];
		/* 儲存(var)temp的結果*/
		pset p_var_t = t_pss[thread_id].psets[4];
		/* var 的supercube，*/
		pset var_supercube = t_pss[thread_id].psets[5];
		/* 儲存(var)對應到已知解的結果*/


		pset solution = t_pss[thread_id].psets[6];
		{
			lock_guard<shared_mutex> lock(sc_mtx);
			set_copy(combine, Supercube);
		}

		// step 1 先用supercube分別對A和B進行過濾
		foreach_set(sf_A, last1, p1) {
			{
				lock_guard<shared_mutex> lock(sc_mtx);
				set_and(combine, p1, Supercube);
			}
			if (!parallel_mvcube_feas_check_notation2(combine, sf_vars[thread_id], sf_sort[thread_id], p_var_t, var_supercube, solution)) {
				RESET(p1, ACTIVE);
				sf_A->active_count--;
			}
		}
		if (sf_A->active_count == 0) {
			printf("combine sf_A = 0\n");
			exit(EXIT_FAILURE);
		}

		foreach_set(sf_B, last2, p2) {
			{
				lock_guard<shared_mutex> lock(sc_mtx);
				set_and(combine, p2, Supercube);
			}
			if (!parallel_mvcube_feas_check_notation2(combine, sf_vars[thread_id], sf_sort[thread_id], p_var_t, var_supercube, solution)) {
				RESET(p2, ACTIVE);
				sf_B->active_count--;
			}
		}
		if (sf_B->active_count == 0) {
			printf("combine sf_B = 0\n");
			exit(EXIT_FAILURE);
		}

		//printf("step1 activate cnt A B (%d %d)\n", sf_A->active_count, sf_B->active_count);

		// step 2 找相同的
		foreach_active_set(sf_A, last1, p1) {
			foreach_active_set(sf_B, last2, p2) {
				/* p2 包含 p1*/
				if (setp_implies(p1, p2)) {
					RESET(p1, ACTIVE);
					sf_addset(sf_C, p1);
					sf_A->active_count--;
				}
				/* p1 包含 p2*/
				if (setp_implies(p2, p1)) {
					RESET(p2, ACTIVE);
					sf_addset(sf_C, p2);
					sf_B->active_count--;
				}
			}
		}
		//printf("step2 activate cnt A B (%d %d)\n", sf_A->active_count, sf_B->active_count);

		// step 3
		foreach_active_set(sf_A, last1, p1) {
			/*printf("\n\n\n");
			display_pset(p1, mvcube.mv_size);*/
			foreach_active_set(sf_B, last2, p2) {
				/*printf("p2\n");
				display_pset(p2, mvcube.mv_size);*/
				set_and(combine, p1, p2);
				/*printf("AND\n");
				display_pset(combine, mvcube.mv_size);*/
				if (notation == 1) {
					printf("UNFINISH !!!!!!!!\n");
					/*if (mvcube_feas_check_notation1(p_and)) {
						sf_addset(temp_Totality, p_and);
					}*/
				}
				else if (notation == 2) {
					if (parallel_mvcube_feas_check_notation2(combine, sf_vars[thread_id], sf_sort[thread_id], p_var_t, var_supercube, solution)) {
						
						sf_addset(sf_C, combine);
					}
				}
			}
		}
	}


	void parallel_combine(pset_family sf_C, pset_family sf_A, pset_family sf_B) {
		sf_active(sf_A);
		sf_active(sf_B);
		mutex mutex_C;
		vector<mutex> mutexs_A(sf_A->count), mutexs_B(sf_B->count);
		//printf("mutex: (%d, %d)\n", mutexs_A.size(), mutexs_B.size());
		combine_chunk_size = 10;
		/* step1*/ /* */
		combine_chunk_cnt_step1A = ceil((double)sf_A->count / (double)combine_chunk_size);
		combine_chunk_cnt_step1B = ceil((double)sf_B->count / (double)combine_chunk_size);
		run_CPP_step1(sf_A, sf_B, ref(mutexs_A), ref(mutexs_B));
		/* step2*/
		combine_chunk_cnt_step2 = ceil((double)sf_A->count / (double)combine_chunk_size);
		run_CPP_step2(sf_A, sf_B, sf_C, ref(mutexs_A), ref(mutexs_B), ref(mutex_C));
		
		/* step3*/
		combine_chunk_cnt_step3 = ceil((double)sf_A->count / (double)combine_chunk_size);
		run_CPP_step3(sf_A, sf_B, sf_C, ref(mutex_C));
	}

	void run_CPP_step1(pset_family sf_A, pset_family sf_B, vector<mutex>& mutexs_A, vector<mutex>& mutexs_B) {
		//printf("STEP1: (%d %d) (%d %d)\n", sf_A->count, sf_B->count, combine_chunk_cnt_step1A, combine_chunk_cnt_step1B);
		//printf("shared_mutex2: (%d)\n", mutexs_A.size());
		for (int id = 0; id < combine_chunk_cnt_step1A; id++) {
			pool.task_begin();
			pool.enqueue(&Parallel_BM::parallel_combine_step1, this, id, combine_chunk_cnt_step1A, sf_A, ref(mutexs_A));
		}

		for (int id = 0; id < combine_chunk_cnt_step1B; id++) {
			pool.task_begin();
			pool.enqueue(&Parallel_BM::parallel_combine_step1, this, id, combine_chunk_cnt_step1B, sf_B, ref(mutexs_B));
		}

		pool.wait();
		sf_inactive(sf_A);
		sf_inactive(sf_B);
	}

	void run_CPP_step2(pset_family sf_A, pset_family sf_B, pset_family sf_C, vector<mutex>& mutexs_A, vector<mutex>& mutexs_B, mutex& mutex_C) {
		for (int id = 0; id < combine_chunk_cnt_step2; id++) {
			pool.task_begin();
			pool.enqueue(&Parallel_BM::parallel_combine_step2, this, id, sf_A, sf_B, sf_C, ref(mutexs_A), ref(mutexs_B), ref(mutex_C));
		}
		pool.wait();
		sf_inactive(sf_A);
		sf_inactive(sf_B);
	}

	void run_CPP_step3(pset_family sf_A, pset_family sf_B, pset_family sf_C, mutex& mutex_C) {
		for (int id = 0; id < combine_chunk_cnt_step3; id++) {
			pool.task_begin();
			pool.enqueue(&Parallel_BM::parallel_combine_step3, this, id, sf_A, sf_B, sf_C, ref(mutex_C));
		}
		pool.wait();
	}

	void parallel_combine_step1(int id, int cnt, pset_family sf, vector<mutex>& sf_mutexs) {
		int thread_id = pool.get_thread_ids(this_thread::get_id());
		pset p;
		/* 儲存(var)temp的結果*/
		pset p_var_t = t_pss[thread_id].psets[4];
		/* var 的supercube，*/
		pset var_supercube = t_pss[thread_id].psets[5];
		/* 儲存(var)對應到已知解的結果*/
		pset solution = t_pss[thread_id].psets[6];

		pset combine = t_pss[thread_id].psets[7];

		int i, end_i;
		end_i = (id == cnt - 1) ? sf->count : (id+1) * combine_chunk_size;
		//printf("combine_step1 %d %d %d %d\n", thread_id, cnt, id, end_i);
		//printf("shared_mutex3: (%d)\n", sf_mutexs.size());
		for (p = GETSET(sf, id * combine_chunk_size), i = id * combine_chunk_size; i < end_i; p += sf->wsize, i++) {
			{
				lock_guard<mutex> lock(sf_mutexs[i]);
				if (!TESTP(p, ACTIVE))
					continue;
			}
			{
				lock_guard<shared_mutex> lock(sc_mtx);
				set_and(combine, p, Supercube);
			}
			if (!parallel_mvcube_feas_check_notation2(combine, sf_vars[thread_id], sf_sort[thread_id], p_var_t, var_supercube, solution)) {
				lock_guard<mutex> lk(sf_mutexs[i]);
				RESET(p, ACTIVE);
				//sf->active_count--;
			}
		}
		pool.task_end();
	}

	void parallel_combine_step2(int id, pset_family sf_A, pset_family sf_B, pset_family sf_C, vector<mutex>& mutexs_A, vector<mutex>& mutexs_B, mutex& mutex_C) {
		int i, j, end_i;
		pset p1, p2;
		end_i = (id == combine_chunk_cnt_step2 - 1) ? sf_A->count : (id + 1) * combine_chunk_size;
		for (p1 = GETSET(sf_A, id * combine_chunk_size), i = id * combine_chunk_size; i < end_i; p1 += sf_A->wsize, i++) {
			if (!TESTP(p1, ACTIVE))
				continue;

			foreachi_active_set(sf_B, j, p2) {
				/* p2 包含 p1*/
				if (setp_implies(p1, p2)) {
					{
						lock_guard<mutex> lk(mutexs_A[i]);
						RESET(p1, ACTIVE);
					}
					{
						lock_guard<mutex> lk(mutex_C);
						sf_addset(sf_C, p1);
					}
					//sf_A->active_count--;
				}
				/* p1 包含 p2*/
				if (setp_implies(p2, p1)) {
					{
						lock_guard<mutex> lk(mutexs_B[j]);
						RESET(p2, ACTIVE);
					}
					{
						lock_guard<mutex> lk(mutex_C);
						sf_addset(sf_C, p2);
					}
					//sf_B->active_count--;
				}
			}
		}
		pool.task_end();
	}
	
	

	void parallel_combine_step3(int id, pset_family sf_A, pset_family sf_B, pset_family sf_C, mutex& mutex_C) {
		int thread_id = pool.get_thread_ids(this_thread::get_id());
		pset p;
		/* 儲存(var)temp的結果*/
		pset p_var_t = t_pss[thread_id].psets[4];
		/* var 的supercube，*/
		pset var_supercube = t_pss[thread_id].psets[5];
		/* 儲存(var)對應到已知解的結果*/
		pset solution = t_pss[thread_id].psets[6];

		pset combine = t_pss[thread_id].psets[7];
		int i, j, end_i;
		pset p1, p2;
		end_i = (id == combine_chunk_cnt_step3 - 1) ? sf_A->count : (id + 1) * combine_chunk_size;

		for (p1 = GETSET(sf_A, id * combine_chunk_size), i = id * combine_chunk_size; i < end_i; p1 += sf_A->wsize, i++) {
			if (!TESTP(p1, ACTIVE))
				continue;

			foreachi_active_set(sf_B, j, p2) {
				//printf("%d (%d, %d)\n", thread_id, i, j);
				set_and(combine, p1, p2);
				/*printf("AND\n");
				display_pset(combine, mvcube.mv_size);*/
				if (notation == 1) {
					printf("UNFINISH !!!!!!!!\n");
					/*if (mvcube_feas_check_notation1(p_and)) {
						sf_addset(temp_Totality, p_and);
					}*/
				}
				else if (notation == 2) {
					if (parallel_mvcube_feas_check_notation2(combine, sf_vars[thread_id], sf_sort[thread_id], p_var_t, var_supercube, solution)) {
						lock_guard<mutex> lk(mutex_C);
						sf_addset(sf_C, combine);
					}
				}
			}
		}
		pool.task_end();
	}
};

#endif 