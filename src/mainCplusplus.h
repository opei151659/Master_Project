/*
	所有header檔的最上層，如需使用其他header 也務必將其在此header進行include確保為全域的include
	包含所有功能的全域呼叫函示與所有的全域變數
*/

#ifndef MAINCPLUSPLUS_H
#define MAINCPLUSPLUS_H

/* 使用boost的bigInteger*/
#include <boost/assign.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <algorithm>
#include <numeric>
#include <set>
#include <map>
#include <time.h>
#include <iostream>
#include <fstream>
#include <direct.h>
#define NDEBUG 
#include <assert.h>
#include <cstdlib> /* 亂數相關函數 */
#include <cstddef>
#include <ctime>   /* 時間相關函數 */
#include <chrono>

/* 平行化函式庫*/
#include <future>
#include <thread>
#include <shared_mutex>
#include <tbb/tbb.h>
#include <omp.h>

extern "C" {
	#include "espresso.h"   /* ESPRESSO.lib*/
}

#include "tool.h" //常用PLA運算
#include "rand_data.h" // 打亂資料


//特徵值
#include "unate.h"
#include "ene.h"
#include "ksig_new.h"
#include "cofactor.h"
#include "support.h"


#include "mcdb.h" //PLA 前置處理與MCQE
#include "chunk.h" //資料分割



using namespace std;

#define print_debug 0
#define print_cost 1
#define print_detail 0
#define notation 2			/* 1:表示abcda'b'c'd' , 2:表示aa'bb'cc'dd'*/ /*notation 1 的產生部分未完成*/
#define EXECost(fct, str) {clock_t t = clock();fct;if (print_cost) printf("%s: %.3lf\n", str, ((double)clock() - (double)t) / (double)CLOCKS_PER_SEC);}
#define EXECost_W(fct, cost) {cost = clock();fct; cost = clock()-cost;}
#define SF_SORT_assend(P) {P = sf_unlist(sf_sort(P, (qsort_compare_func)ascend), P->count, P->sf_size);}



class mv_cube {
public:
	int num_mv_var;
	int mv_var_size;
	int mv_var_word_size;
	int mv_size;
	int num_output;
	
	bool is_input_permutation;
	bool is_input_phase_assigment;
	bool is_output_permutation;
	bool is_output_phase_assigment;

	unsigned int* first_word;	/*each var*/
	unsigned int* lhalf_word;	/*each var*/
	unsigned int* rhalf_word;	/*each var*/
	unsigned int* last_word;	/*each var*/
	unsigned int* first_bit;	/*each var*/
	unsigned int* lhalf_bit;	/*each var*/
	unsigned int* rhalf_bit;	/*each var*/
	unsigned int* last_bit;		/*each var*/

	mv_cube() {
		num_mv_var = 0;
		mv_var_size = 0;
		mv_var_word_size = 0;
		mv_size = 0;
		num_output = 0;
		is_input_permutation = false;
		is_input_phase_assigment = false;
		is_output_permutation = false;
		is_output_phase_assigment = false;
		// initial space to prevent warning
		first_word = new unsigned int[1]; // each vars first word index
		last_word = new unsigned int[1];  // each vars last word index
		first_bit = new unsigned int[1];  // each vars first bit index
		last_bit = new unsigned int[1];	  // each vars last bit index
	};
	mv_cube(bool inpermu, bool inphase, bool outpermu, bool outphase) {
		init(inpermu, inphase, outpermu, outphase);
	};
	~mv_cube() {
		delete first_word;
		delete last_word;
		delete first_bit;
		delete last_bit;
	};
	void init(bool inpermu, bool inphase, bool outpermu, bool outphase) {
		is_input_permutation = inpermu;
		is_input_phase_assigment = inphase;
		is_output_permutation = outpermu;
		is_output_phase_assigment = outphase;

		num_mv_var = NUMINPUTS; // number of var x matching to vars y
		mv_var_size = NUMINPUTS;
		num_output = NUMOUTPUTS;
		if (is_input_phase_assigment)
			mv_var_size += NUMINPUTS;
		mv_var_word_size = floor(mv_var_size / 32.0) + 1;
		mv_size = num_mv_var * mv_var_size;

		// resize the initial dynamic space
		first_word = new unsigned int[NUMINPUTS];
		lhalf_word = new unsigned int[NUMINPUTS];
		rhalf_word = new unsigned int[NUMINPUTS];
		last_word = new unsigned int[NUMINPUTS];
		first_bit = new unsigned int[NUMINPUTS];
		lhalf_bit = new unsigned int[NUMINPUTS];
		rhalf_bit = new unsigned int[NUMINPUTS];
		last_bit = new unsigned int[NUMINPUTS];

		/* calculate each part index*/
		/*
			inputs with no phase:  (first_word, first_bit) --- mv_var_size --- (last_word, last_bit)
			inputs with phse: (first_word, first_bit) --- mv_var_size --- (lhalf_word, lhalf_bit) +
								(rhalf_word, rhalf_bit) --- mv_var_size --- (last_word, last_bit)
		*/
		for (int i = 0; i < NUMINPUTS; i++) {
			first_word[i] = (i * mv_var_size) / BPI + 1;
			last_word[i] = ((i + 1) * mv_var_size - 1) / BPI + 1;
			first_bit[i] = (i * mv_var_size) % BPI;
			last_bit[i] = ((i + 1) * mv_var_size - 1) % BPI;
		}
	};
	void display() {
		printf("num_mv_var: %d\n", num_mv_var);
		// word index
		printf("first_word: "); for (int i = 0; i < NUMINPUTS; i++) printf("%4d", first_word[i]); printf("\n");
		printf("last_word:  "); for (int i = 0; i < NUMINPUTS; i++) printf("%4d", last_word[i]); printf("\n");
		// bit index
		printf("first_bit:  "); for (int i = 0; i < NUMINPUTS; i++) printf("%4d", first_bit[i]); printf("\n");
		printf("last_bit:   "); for (int i = 0; i < NUMINPUTS; i++) printf("%4d", last_bit[i]); printf("\n");
	};
};


struct pm_Detail {
	mutex mtx;
	unsigned long long int smartMvCube_Mvcubes_num;				/* 每個(p,q) 所產生的所有mvcubes個數*/
	unsigned long long int smartMvCube_Mvcubes_num_redundant;	/* 當前(p,q) 所產生的所有mvcubes個數且該p,q 為redundant*/
	unsigned long long int smartMvCube_empty_cube;				/* 在S_MAP中被Supercube AND後為empty cube的mvcube個數*/
	unsigned long long int smartMvCube_empty_cube_temp;			/* 在當前S_MAP中被Supercube AND後為empty cube的mvcube個數*/
	unsigned long long int smartMvCube_empty_cube_pq;			/* 以(p, q)為單位，該S_MAP中有出現empty的個數*/
	unsigned long long int smartMvcube_redundant;				/* 以(p, q)為單位，該S_MAP的mvcubes全部為redundant的個數*/
	unsigned long long int outputintersection;					/* 以(p, q)為單位，output有相交的個數*/
	void init() {
		smartMvCube_Mvcubes_num = 0;
		smartMvCube_Mvcubes_num_redundant = 0;
		smartMvCube_empty_cube = 0;
		smartMvCube_empty_cube_temp = 0;
		smartMvCube_empty_cube_pq = 0;
		smartMvcube_redundant = 0;
		outputintersection = 0;
	}
	void display() {
		printf("Cube NUM:\n");
		printf("total_cube_num: %lld\n", smartMvCube_Mvcubes_num);
		printf("empty_cube_num: %lld\n", smartMvCube_empty_cube);
		printf("redundant_cube_num: %lld\n", smartMvCube_Mvcubes_num_redundant);
		printf("P,Q NUM:\n");
		printf("outputintersection_P_Q_num: %lld\n", outputintersection);
		printf("empty_cube_P_Q_num: %lld\n", smartMvCube_empty_cube_pq);
		printf("redundant_P_Q_num: %lld\n", smartMvcube_redundant);
	}
};

/**********************/
/***   global var	***/
/**********************/
extern int CAL_SIG_OUTPUT;	/* 是否計算輸出特徵值  (1/0)*/
extern int CAL_SIG_INPUT;
extern int CAL_ALL_MCQE;
extern mv_cube mvcube; 
extern cube_struct mycube;
extern pm_Detail pm_detail;
extern int row_num;
extern int col_num;
extern mutex memory_mtx;
extern int ppbm_thread_num;
extern int ppbm_chunk_size;
extern int ppbm_option;
extern double total_and_time;
extern double total_redundant_time;
extern float remove_sf_rate;
extern int extand_input_num;
extern int SIG_key;
extern int MCQE_key;
extern int Total_parallel;

extern clock_t AND_time, RR_time;
/**********************/
/***	function	***/
/**********************/
/* BM.cpp*/			extern void BM(pPLA PLA, pPLA PLA1, bool INP, bool INPA, bool OUTP, bool OUTPA);
/* BM.cpp */		extern bool S_MAP(pset_family Totality, pset_family temp_Totality, pset Supercube, pset literal_pos, pset Q_non_DC, pset p_and);

/* 4種 partial_mapping實作方式*/
/* BM.cpp */		extern pset_family partial_mapping(pset_family Totality, pset Supercube, MCDB& mcdb);

	/*1. 部分平行化，PPBM設為1則為未平行化的版本  (主要使用這個)*/
/* PPBM.cpp */ extern void partial_parallel_mapping_block(int f_begin, int f_end, int g_begin, int g_end, pset_family Totality, pset_family temp_Totality, pset_family temp_Totality2, pset Supercube, MCDB& mcdb);

	/*2. cm150a 複數積項配對(未平行化) window size = 2的實驗*/
/* PPBM.cpp */ extern void partial_parallel_mapping_block_window(int f_begin, int f_end, int g_begin, int g_end, pset_family Totality, pset_family temp_Totality, pset_family temp_Totality2, pset Supercube, MCDB& mcdb, int window);
	
	/*3. 原paper的實作重現*/
/* paper.cpp */ extern void paper_partial_mapping_block(int f_begin, int f_end, int g_begin, int g_end, pset_family Totality, pset_family temp_Totality, pset Supercube, MCDB& mcdb);


	/*4. 第一版改良後的PM(完全未平行化版)*/  /*!!!可能有缺少部分加速技巧!!!*/
/* BM.cpp */		extern void partital_map_and(pset_family Totality, pset_family temp_Totality, pset Supercube, pset literal_pos, pset Q_non_DC, pset p_and);
/* BM.cpp */		extern void partial_mapping_block(int f_begin, int f_end, int g_begin, int g_end, pset_family Totality, pset_family temp_Totality, pset Supercube, MCDB& mcdb);


/* 輸入特徵值*/
/* cplusplus.cpp */ extern void check_unate(pset p, pPLA PLA1, pPLA PLA2);
/* cplusplus.cpp */ extern void check_ENE(pset p, pPLA PLA1, pPLA PLA2, int option);
/* cplusplus.cpp */ extern void grouping_ENE(pset_family sf, vector<group_vars>& group_E, vector<group_vars>& group_NE, vector<group_vars>& group_ENE, vector<int>& group_X);
/* cplusplus.cpp */ extern void matching_ENE(vector<group_vars>& group_E_1, vector<group_vars>& group_NE_1, vector<group_vars>& group_ENE_1, vector<int>& group_X_1
											, vector<group_vars>& group_E_2, vector<group_vars>& group_NE_2, vector<group_vars>& group_ENE_2, vector<int>& group_X_2, pset totality);
/* cplusplus.cpp */ extern void check_KSIG(pset p, pPLA PLA1, pPLA PLA2, int distN, int mode);
/* cplusplus.cpp */ extern void matching_KSIG(pset p, vector<TABLE>& KSIG_input1, vector<TABLE>&  KSIG_input2);
/* cplusplus.cpp */ extern void check_cofactor(pset p, pPLA PLA1, pPLA PLA2);
/* cplusplus.cpp */ extern void check_character_value(int key, pset totality, pPLA PLA1, pPLA PLA2);
/* cplusplus.cpp */ extern void run_InputSIG_TEST(pset totality, pPLA PLA1, pPLA PLA2);

/* 輸出特徵值*/
/* cplusplus.cpp */ extern void matching_KSIG_output(pset p_out, TABLE& KSIG_output1, TABLE& KSIG_output2);
/* cplusplus.cpp */ extern void matching_KSIG_output_with_KSIG_input(pset p_out, vector<TABLE>& KSIG_output1, vector<TABLE>& KSIG_output2);
/* cplusplus.cpp */ extern void check_Unate_output(pset p_out, pPLA PLA1, pPLA PLA2);
/* cplusplus.cpp */ extern void run_OutputSIG_TEST(pPLA PLA1, pPLA PLA2);

/* pset 運算*/
/* cplusplus.cpp */ extern bool setp_implies_var(pset a, pset b, int var);
/* cplusplus.cpp */ extern bool set_and_var(pset r, pset a, pset b, int var);

/* mvcube檢查合法性*/
/* cplusplus.cpp */ extern bool mvcube_feas_check_notation1(pset p); /* 可能有問題*/
/* cplusplus.cpp */ extern bool mvcube_feas_check_notation2(pset p);

/* 移除多餘的mvcube*/
/* cplusplus.cpp */ extern void sf_removeRedundant(pset_family A, int skip_num);

#endif //MAINCPLUSPLUS_H