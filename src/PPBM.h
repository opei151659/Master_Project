/*
	Parallel Partial Boolean Mapping
	部分平行化的布林比對
	分為AND與Redundant兩個部分
	兩個部分各自都有三種平行化函式庫的實作方式
	1. OMP (openMP)
	2. CPP (c++ 11 concurrency)
	3. TBB (intel TBB)

	其中OMP 與TBB 可再細分為原paper作法與使用了一些加速方法的版本(較早期 故不完整 不包含在實驗中)
*/
#ifndef PPBM_H
#define PPBM_H
#include <iostream>
#include <thread>
#include <atomic>
#include <vector>
#include <deque>
#include <omp.h>
#include <tbb/parallel_for.h>
#include "tool.h"
#include "parallel_tool.h"



using namespace std;

bool parallel_mvcube_feas_check_notation2(pset p, pset_family sf_vars, pset_family sf_sort, pset p_var_t, pset var_supercube, pset solution);
void parallel_sf_removeRedundant(pset_family A);
bool parallel_S_MAP(pset_family Totality, pset_family temp_Totality, pset Supercube, pset literal_pos, pset Q_non_DC, pset p_and, bool window = false);



class parallel_and {
private:
	int thread_num, current_thread_num, chunk_size, chunk_cnt;
	pset_family Totality, temp_Totality;
	pset Supercube, literal_0_pos, literal_1_pos, Q, q;
	vector<pset_family> temp_Totalitys, sf_vars, sf_sort;
	vector<temp_psets> t_pss;
	atomic<int> chunk_id;

public:
	/* collect data*/
	double average_thread_used;
	int cnt_and_called;
	int cnt_and;

	parallel_and(int thread_num_, int chunk_size_) : thread_num(thread_num_), chunk_size(chunk_size_), current_thread_num(0), chunk_cnt(0), average_thread_used(0.0), cnt_and_called(0) {
		setup();
	};
	~parallel_and() {
		for (int i = 0; i < thread_num; i++) {
			sf_free(temp_Totalitys[i]);
			sf_free(sf_vars[i]);
			sf_free(sf_sort[i]);
		}
	}
	/* 只會執行一次*/
	void setup() {
		temp_Totalitys.resize(thread_num);
		sf_vars.resize(thread_num);
		sf_sort.resize(thread_num);
		t_pss.resize(thread_num);

		for (int i = 0; i < thread_num; i++) {
			temp_Totalitys[i] = sf_new(1, mvcube.mv_size);
			sf_vars[i] = sf_new(mvcube.num_mv_var, mvcube.mv_var_size);
			sf_sort[i] = sf_new(mvcube.num_mv_var, mvcube.mv_var_size);

			/* feasible notation 2*/
			t_pss[i].psets.resize(4);
			t_pss[i].psets[0] = set_new(mvcube.mv_var_size);
			t_pss[i].psets[1] = set_new(mvcube.mv_var_size);
			t_pss[i].psets[2] = set_new(mvcube.mv_var_size);
			t_pss[i].psets[3] = set_new(mvcube.mv_size);
		}
	}
	/* 每次使用前執行一次*/
	void init(pset_family Totality_, pset_family temp_Totality_, pset Supercube_, pset literal_0_pos_, pset Q_, pset literal_1_pos_, pset q_) {
		Totality = Totality_;
		temp_Totality = temp_Totality_;
		Supercube = Supercube_;
		literal_0_pos = literal_0_pos_;
		literal_1_pos = literal_1_pos_;
		Q = Q_;
		q = q_;

		if (thread_num < 0)
			fprintf(stderr, "parallel_and wrong initialization\n");

		chunk_cnt = ceil((double)Totality->count / (double)chunk_size);
		for (int i = 0; i < thread_num; i++) {
			temp_Totalitys[i]->count = 0;
		}
		chunk_id.store(0);
	}
	void run(int option) {
		cnt_and = 0;
		auto start = clock();
		if (option == 1) {
			//run_OMP();
			run_OMP_origin();
		}
		else if (option == 2) {
			run_CPP();
		}
		else if (option == 3) {
			//run_TBB();
			run_TBB_origin();
		}
		else {
			fprintf(stderr, "wrong parallel and option\n");
		}
		AND_time += clock() - start;

		merge();
	}
	/* 加速後的OMP平行化*/
	void run_OMP() {
		current_thread_num = 0;
#pragma omp for
		for (int thread_id = 0; thread_id < thread_num; thread_id++) {
			if (chunk_id.load() >= chunk_cnt)
				break;
			parallel_and::work(thread_id);
#pragma omp critical
			current_thread_num += 1;
		}
	}
	/* 未加速後的OMP平行化*/
	void run_OMP_origin() {
		omp_set_num_threads(thread_num);
		current_thread_num = thread_num;
		omp_lock_t lock;
		omp_init_lock(&lock); // Initialize the lock

#pragma omp parallel for
		for (int p_i = 0; p_i < Totality->count; p_i++) {
			//int thread_id = omp_get_thread_num(); printf("%d %d\n", Totality->count, thread_id);
			pset p, last, using_q, p_var_t, var_supercube, solution, p_and;
			int var;
			p_var_t = set_new(mvcube.mv_var_size);
			var_supercube = set_new(mvcube.mv_var_size);
			solution = set_new(mvcube.mv_var_size);
			p_and = set_new(mvcube.mv_size);
			p = GETSET(Totality, p_i);
			pset_family sf_var_t = sf_new2(mvcube.num_mv_var, mvcube.mv_var_size);
			pset_family sf_sort_t = sf_new2(mvcube.num_mv_var, mvcube.mv_var_size);
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
						if (parallel_mvcube_feas_check_notation2(p_and, ref(sf_var_t), ref(sf_sort_t), p_var_t, var_supercube, solution)) {
							PUTSIZE(p_and, set_count_ones(p_and));
							omp_set_lock(&lock);
							{
								//sf_addset(temp_Totalitys[thread_id], p_and);
								sf_addset(temp_Totality, p_and);
							}
							omp_unset_lock(&lock);
						}
					}
				}

			}
			set_free(p_var_t);
			set_free(var_supercube);
			set_free(solution);
			set_free(p_and);
			sf_free2(sf_var_t);
			sf_free2(sf_sort_t);
		}
	}

	/* 加速後的CPP平行化*/
	void run_CPP() {
		/* thread pool*/
		if (use_with_thread_pool) {
			pool.reset_Thread_task_used_cnt();
			current_thread_num = thread_num;
			for (int id = 0; id < chunk_cnt; id++) {
				pool.task_begin();
				pool.enqueue(&parallel_and::work2, this, id);
			}

			pool.wait();
			cnt_and_called++;
			//printf("AND chunk_cnt %d %d %d\n", chunk_cnt, pool.display_Thread_task_used_cnt(), (double)pool.display_Thread_task_used_cnt());
			average_thread_used += (double)pool.display_Thread_task_used_cnt();//(double)chunk_cnt / (double)pool.display_Thread_task_used_cnt();
		}
		
		/* bad way without thread pool*/
		else {
			current_thread_num = MIN(thread_num, chunk_cnt);
			for (int id_chunk = 0; id_chunk < chunk_cnt; ) {
				vector<thread> threads;
				for (int thread_id = 0; thread_id < thread_num && id_chunk < chunk_cnt; thread_id++, id_chunk++) {
					threads.emplace_back(&parallel_and::work_no_thread_pool, this, thread_id, id_chunk, id_chunk + 1);
				}
				for (int thread_id = 0; thread_id < threads.size(); thread_id++) {
					threads[thread_id].join();
				}

			}
		}
	}

	/* 加速後的TBB平行化*/
	void run_TBB() {
		current_thread_num = 0;
		tbb::parallel_for(0, thread_num, 1, [&](unsigned long long thread_id) {
				parallel_and::work(thread_id);
				current_thread_num += 1;
			});
	}

	/* 未加速後的TBB平行化*/
	void run_TBB_origin() {
		tbb::spin_mutex mutex;
		current_thread_num = thread_num;
		tbb::parallel_for(tbb::blocked_range<int>(0, Totality->count),
			[&](tbb::blocked_range<int> r) {
			for (int p_i = r.begin(); p_i < r.end(); ++p_i) {
				//int thread_id = std::this_thread::get_id(); printf("%d %d\n", Totality->count, thread_id);
				pset p, last, using_q, p_var_t, var_supercube, solution, p_and;
				int var;
				p_var_t = set_new(mvcube.mv_var_size);
				var_supercube = set_new(mvcube.mv_var_size);
				solution = set_new(mvcube.mv_var_size);
				p_and = set_new(mvcube.mv_size);
				pset_family sf_var_t = sf_new2(mvcube.num_mv_var, mvcube.mv_var_size);
				pset_family sf_sort_t = sf_new2(mvcube.num_mv_var, mvcube.mv_var_size);
				p = GETSET(Totality, p_i);
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
						if (parallel_mvcube_feas_check_notation2(p_and, ref(sf_var_t), ref(sf_sort_t), p_var_t, var_supercube, solution)) {
							PUTSIZE(p_and, set_count_ones(p_and));
							{
								tbb::spin_mutex::scoped_lock lock(mutex);
								//sf_addset(temp_Totalitys[thread_id], p_and);
								sf_addset(temp_Totality, p_and);
							}
						}
					}
				}

			}
			set_free(p_var_t);
			set_free(var_supercube);
			set_free(solution);
			set_free(p_and);
			sf_free2(sf_var_t);
			sf_free2(sf_sort_t);
			}
			});
	}
	
	/* 未加速後的Parital mapping*/
	void work(int thread_id) {
		pset p, last, using_q, p_var_t, var_supercube, solution, p_and;
		int var, id, chunk_end;

		/* 儲存(var)temp的結果*/
		p_var_t = t_pss[thread_id].psets[0];
		/* var 的supercube，*/
		var_supercube = t_pss[thread_id].psets[1];
		/* 儲存(var)對應到已知解的結果*/
		solution = t_pss[thread_id].psets[2];
		/* 儲存and的結果*/
		p_and = t_pss[thread_id].psets[3];
		
		while ((id = chunk_id.fetch_add(1)) < chunk_cnt) {
			chunk_end = (id == chunk_cnt - 1) ? Totality->count : (id + 1) * chunk_size;
			for (p = GETSET(Totality, id * chunk_size), last = GETSET(Totality, chunk_end); p < last; p += Totality->wsize) {
				
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
							if (parallel_mvcube_feas_check_notation2(p_and, ref(sf_vars[thread_id]), ref(sf_sort[thread_id]), p_var_t, var_supercube, solution)) {
								if (mvcube_feas_check_notation2(p_and)) {
									PUTSIZE(p_and, set_count_ones(p_and));
									sf_addset(temp_Totalitys[thread_id], p_and);
								}
							}
						}
					}
				}
			}
		}
	}

	/* 加速後的Parital mapping使用類似執行緒池的簡化實作方式 但無實際的執行緒池*/
	void work_no_thread_pool(int thread_id, int begin_chunk_id, int end_chunk_id) {
		pset p, last, using_q, p_var_t, var_supercube, solution, p_and;
		int var, chunk_end;

		/* 儲存(var)temp的結果*/
		p_var_t = t_pss[thread_id].psets[0];
		/* var 的supercube，*/
		var_supercube = t_pss[thread_id].psets[1];
		/* 儲存(var)對應到已知解的結果*/
		solution = t_pss[thread_id].psets[2];
		/* 儲存and的結果*/
		p_and = t_pss[thread_id].psets[3];

		for (int id = begin_chunk_id; id < end_chunk_id; id++) {
			chunk_end = (id == chunk_cnt - 1) ? Totality->count : (id + 1) * chunk_size;
			for (p = GETSET(Totality, id * chunk_size), last = GETSET(Totality, chunk_end); p < last; p += Totality->wsize) {

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
							if (parallel_mvcube_feas_check_notation2(p_and, ref(sf_vars[thread_id]), ref(sf_sort[thread_id]), p_var_t, var_supercube, solution)) {
								PUTSIZE(p_and, set_count_ones(p_and));
								sf_addset(temp_Totalitys[thread_id], p_and);
							}
						}
					}

				}

			}
		}
	}
	/* 加速後的Parital mapping使用執行緒池*/
	void work2(int id) {
		pset p, last, using_q, p_var_t, var_supercube, solution, p_and;
		int var, chunk_end, thread_id = pool.get_thread_ids(this_thread::get_id());

		/* 儲存(var)temp的結果*/
		p_var_t = t_pss[thread_id].psets[0];
		/* var 的supercube，*/
		var_supercube = t_pss[thread_id].psets[1];
		/* 儲存(var)對應到已知解的結果*/
		solution = t_pss[thread_id].psets[2];
		/* 儲存and的結果*/
		p_and = t_pss[thread_id].psets[3];

		
		chunk_end = (id == chunk_cnt - 1) ? Totality->count : (id + 1) * chunk_size;
		for (p = GETSET(Totality, id * chunk_size), last = GETSET(Totality, chunk_end); p < last; p += Totality->wsize) {
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
				cnt_and++;
				
				/* 與所有現有的toatlity 做and*/
				if (set_and_var(p_and, p, using_q, var)) {
					if (notation == 1) {
						printf("UNFINISH !!!!!!!!\n");
						/*if (mvcube_feas_check_notation1(p_and)) {
							sf_addset(temp_Totality, p_and);
						}*/
					}
					else if (notation == 2) {
						if (parallel_mvcube_feas_check_notation2(p_and, ref(sf_vars[thread_id]), ref(sf_sort[thread_id]), p_var_t, var_supercube, solution)) {
							PUTSIZE(p_and, set_count_ones(p_and));
							sf_addset(temp_Totalitys[thread_id], p_and);
						}
					}
				}

			}
		}
		pool.task_end();
	}




	/* 整理合併每個執行緒的計算結果*/
	void merge() {
		int total_size = 0;
		for (int i = 0; i < current_thread_num; i++) {
			total_size += temp_Totalitys[i]->count;
		}
		if (temp_Totality->capacity < temp_Totality->count + total_size) {
			temp_Totality->capacity = temp_Totality->count + total_size;
			temp_Totality->data = REALLOC(unsigned int, temp_Totality->data, (long)temp_Totality->capacity * temp_Totality->wsize);
		}
		unsigned int* begin = temp_Totality->data + temp_Totality->wsize * temp_Totality->count;
		for (int i = 0; i < current_thread_num; i++) {
			local_intcpy(begin, temp_Totalitys[i]->data, temp_Totalitys[i]->wsize * temp_Totalitys[i]->count);
			begin += temp_Totalitys[i]->wsize * temp_Totalitys[i]->count;
		}
		temp_Totality->count += total_size;
	}
};

class parallel_remove_redundant {
private:
	int thread_num, chunk_size, chunk_cnt;
	atomic<int> chunk_id;
	int skip_num;
	pset_family sf;
public:
	/* collect data*/
	double average_thread_used;
	int cnt_and_called;

	parallel_remove_redundant(int thread_num_, int chunk_size_) : thread_num(thread_num_), chunk_size(chunk_size_), chunk_cnt(0), average_thread_used(0.0), cnt_and_called(0) {};
	~parallel_remove_redundant() {};
	void init(pset_family sf_, int _skip_num) {
		sf = sf_;
		skip_num = _skip_num;
		if (thread_num < 0)
			fprintf(stderr, "parallel_and wrong initialization\n");

		chunk_cnt = ceil((double)sf->count / (double)chunk_size);
		chunk_id.store(0);
	}
	void run(int option) {
		vector<mutex> mutexs(sf->count);
		sf_active(sf);
		auto start = clock();
		if (option == 1) {
			//run_OMP(mutexs);
			run_OMP_origin(mutexs);
		}
		else if (option == 2) {
			run_CPP(mutexs);
		}
		else if (option == 3) {
			//run_TBB(mutexs);
			run_TBB_origin(mutexs);
		}
		else {
			fprintf(stderr, "wrong parallel and option\n");
		}
		RR_time += clock() - start;
		sf_inactive(sf);
	}
	void run_OMP(vector<mutex>& mutexs) {
		omp_set_num_threads(thread_num);
#pragma omp parallel
		{
#pragma omp for
			for (int thread_id = 0; thread_id < thread_num; thread_id++) {
				if (chunk_id.load() >= chunk_cnt)
					break;
				work2(ref(mutexs));
			}
		}
	}

	void run_OMP_origin(vector<mutex>& mutexs) {
		omp_lock_t lock; 
		omp_init_lock(&lock); // Initialize the lock
		omp_set_num_threads(thread_num);
#pragma omp parallel for
		for (int i = 0; i < sf->count; i++) {
			//int thread_id = omp_get_thread_num(); printf("%d\n", thread_id);
			pset p = GETSET(sf, i);
			if (!TESTP(p, ACTIVE))
				continue;
			for (int j = i + 1; j < sf->count; j++) {
				pset p2 = GETSET(sf, j);
				if (i < skip_num && j < skip_num) {
					j = skip_num;
					p2 = GETSET(sf, j);
				}
				if (!TESTP(p2, ACTIVE))
					continue;
				if (SIZE(p) >= SIZE(p2) && setp_implies(p2, p)) {
					lock_guard<mutex> lock(mutexs[j]);
					RESET(p2, ACTIVE);
				}
				else if (SIZE(p) <= SIZE(p2) && setp_implies(p, p2)) {
					lock_guard<mutex> lock(mutexs[j]);
					RESET(p, ACTIVE);
				}
			}
		}
	}

	void run_CPP(vector<mutex>& mutexs) { 
		/* thread_pool*/
		if (use_with_thread_pool) {
			pool.reset_Thread_task_used_cnt();
			for (int id = 0; id < chunk_cnt; id++) {
				pool.task_begin();
				pool.enqueue(&parallel_remove_redundant::work2_2, this, id, ref(mutexs));
			}
			pool.wait();
			cnt_and_called++;
			//printf("RR chunk_cnt %d %d %d\n", chunk_cnt, pool.display_Thread_task_used_cnt(), (double)pool.display_Thread_task_used_cnt());
			average_thread_used += (double)pool.display_Thread_task_used_cnt();//(double)chunk_cnt / (double)pool.display_Thread_task_used_cnt();
			/*if (chunk_cnt > 1) {


			}*/
		}
		/* without thread pool*/
		else {
			vector<thread> threads;
			int chunk_cnt_size = chunk_cnt / thread_num;
			for (int thread_id = 0; thread_id < thread_num; thread_id++) {
				if ((thread_id + 1) * chunk_cnt_size >= chunk_cnt) {
					threads.emplace_back(&parallel_remove_redundant::work2_without_thread_pool, this, ref(mutexs), thread_id * chunk_cnt_size, chunk_cnt);
					break;
				}
				else {
					threads.emplace_back(&parallel_remove_redundant::work2_without_thread_pool, this, ref(mutexs), thread_id * chunk_cnt_size, (thread_id + 1) * chunk_cnt_size);
				}
			}
			for (int thread_id = 0; thread_id < threads.size(); thread_id++) {
				threads[thread_id].join();
			}
		}
	}
	void run_TBB(vector<mutex>& mutexs) {
		tbb::parallel_for(0, thread_num, 1, [&](unsigned long long thread_id) {
			parallel_remove_redundant::work2(ref(mutexs));
			});
	}

	void run_TBB_origin(vector<mutex>& mutexs) {
		tbb::spin_mutex mutex;		
		tbb::parallel_for(tbb::blocked_range<int>(0, sf->count),
			[&](tbb::blocked_range<int> r) {
				for (int i = r.begin(); i < r.end(); ++i) {
					pset p = GETSET(sf, i);
					if (!TESTP(p, ACTIVE))
						continue;
					for (int j = i + 1; j < sf->count; j++) {
						pset p2 = GETSET(sf, j);
						if (i < skip_num && j < skip_num) {
							j = skip_num;
							p2 = GETSET(sf, j);
						}
						if (!TESTP(p2, ACTIVE))
							continue;
						if (SIZE(p) >= SIZE(p2) && setp_implies(p2, p)) {
							tbb::spin_mutex::scoped_lock lock(mutex);
							RESET(p2, ACTIVE);
						}
						else if (SIZE(p) <= SIZE(p2) && setp_implies(p, p2)) {
							tbb::spin_mutex::scoped_lock lock(mutex);
							RESET(p, ACTIVE);
						}
					}
				}
			});
	}


	void work2(vector<mutex>& mutexs) {
		
		pset p, p2;
		int i, j;
		int id, end_i;
		while ((id = chunk_id.fetch_add(1)) < chunk_cnt) {
			
			end_i = (id == chunk_cnt - 1) ? sf->count : id * chunk_size + chunk_size;
			for (p = GETSET(sf, id * chunk_size), i = id * chunk_size; i < end_i; p += sf->wsize, i++) {
				if (!TESTP(p, ACTIVE))
					continue;
				
				for (p2 = GETSET(sf, i + 1), j = i + 1; j < sf->count; p2 += sf->wsize, j++) {
					if (i < skip_num && j < skip_num) {
						j = skip_num;
						p2 = GETSET(sf, j);
					}
					if (!TESTP(p2, ACTIVE))
						continue;
					if (SIZE(p) >= SIZE(p2) && setp_implies(p2, p)) {
						lock_guard<mutex> lock(mutexs[j]);
						RESET(p2, ACTIVE);
					}
					else if(SIZE(p) <= SIZE(p2) && setp_implies(p, p2)) {
						lock_guard<mutex> lock(mutexs[i]);
						RESET(p, ACTIVE);
					}
				}
			}
		}
	}
	void work2_without_thread_pool(vector<mutex>& mutexs, int begin_chunk_id, int end_chunk_id) {
		pset p, p2;
		int i, j;
		int end_i;
		for (int id = begin_chunk_id; id < end_chunk_id; id++) {
			end_i = (id == chunk_cnt - 1) ? sf->count : id * chunk_size + chunk_size;
			for (p = GETSET(sf, id * chunk_size), i = id * chunk_size; i < end_i; p += sf->wsize, i++) {
				if (!TESTP(p, ACTIVE))
					continue;

				for (p2 = GETSET(sf, i + 1), j = i + 1; j < sf->count; p2 += sf->wsize, j++) {
					if (i < skip_num && j < skip_num) {
						j = skip_num;
						p2 = GETSET(sf, j);
					}
					if (!TESTP(p2, ACTIVE))
						continue;
					if (SIZE(p) >= SIZE(p2) && setp_implies(p2, p)) {
						lock_guard<mutex> lock(mutexs[j]);
						RESET(p2, ACTIVE);
					}
					else if (SIZE(p) <= SIZE(p2) && setp_implies(p, p2)) {
						lock_guard<mutex> lock(mutexs[i]);
						RESET(p, ACTIVE);
					}
				}
			}
		}
	}
	/* with thread_pool*/
	void work2_2(int id , vector<mutex>& mutexs) {
		pset p, p2;
		int i, j;
		int end_i;
		

		end_i = (id == chunk_cnt - 1) ? sf->count : id * chunk_size + chunk_size;
		for (p = GETSET(sf, id * chunk_size), i = id * chunk_size; i < end_i; p += sf->wsize, i++) {
			if (!TESTP(p, ACTIVE))
				continue;

			for (p2 = GETSET(sf, i + 1), j = i + 1; j < sf->count; p2 += sf->wsize, j++) {
				if (i < skip_num && j < skip_num) {
					j = skip_num;
					p2 = GETSET(sf, j);
				}
				if (!TESTP(p2, ACTIVE))
					continue;
				if (SIZE(p) >= SIZE(p2) && setp_implies(p2, p)) {
					lock_guard<mutex> lock(mutexs[j]);
					RESET(p2, ACTIVE);
				}
				else if (SIZE(p) <= SIZE(p2) && setp_implies(p, p2)) {
					lock_guard<mutex> lock(mutexs[i]);
					RESET(p, ACTIVE);
				}
			}
		}
		pool.task_end();
	}
};

#endif
