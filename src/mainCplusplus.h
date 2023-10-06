/*
	�Ҧ�header�ɪ��̤W�h�A�p�ݨϥΨ�Lheader �]�ȥ��N��b��header�i��include�T�O�����쪺include
	�]�t�Ҧ��\�઺����I�s��ܻP�Ҧ��������ܼ�
*/

#ifndef MAINCPLUSPLUS_H
#define MAINCPLUSPLUS_H

/* �ϥ�boost��bigInteger*/
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
#include <cstdlib> /* �üƬ������ */
#include <cstddef>
#include <ctime>   /* �ɶ�������� */
#include <chrono>

/* ����ƨ禡�w*/
#include <future>
#include <thread>
#include <shared_mutex>
#include <tbb/tbb.h>
#include <omp.h>

extern "C" {
	#include "espresso.h"   /* ESPRESSO.lib*/
}

#include "tool.h" //�`��PLA�B��
#include "rand_data.h" // ���ø��


//�S�x��
#include "unate.h"
#include "ene.h"
#include "ksig_new.h"
#include "cofactor.h"
#include "support.h"


#include "mcdb.h" //PLA �e�m�B�z�PMCQE
#include "chunk.h" //��Ƥ���



using namespace std;

#define print_debug 0
#define print_cost 1
#define print_detail 0
#define notation 2			/* 1:���abcda'b'c'd' , 2:���aa'bb'cc'dd'*/ /*notation 1 �����ͳ���������*/
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
	unsigned long long int smartMvCube_Mvcubes_num;				/* �C��(p,q) �Ҳ��ͪ��Ҧ�mvcubes�Ӽ�*/
	unsigned long long int smartMvCube_Mvcubes_num_redundant;	/* ��e(p,q) �Ҳ��ͪ��Ҧ�mvcubes�ӼƥB��p,q ��redundant*/
	unsigned long long int smartMvCube_empty_cube;				/* �bS_MAP���QSupercube AND�ᬰempty cube��mvcube�Ӽ�*/
	unsigned long long int smartMvCube_empty_cube_temp;			/* �b��eS_MAP���QSupercube AND�ᬰempty cube��mvcube�Ӽ�*/
	unsigned long long int smartMvCube_empty_cube_pq;			/* �H(p, q)�����A��S_MAP�����X�{empty���Ӽ�*/
	unsigned long long int smartMvcube_redundant;				/* �H(p, q)�����A��S_MAP��mvcubes������redundant���Ӽ�*/
	unsigned long long int outputintersection;					/* �H(p, q)�����Aoutput���ۥ檺�Ӽ�*/
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
extern int CAL_SIG_OUTPUT;	/* �O�_�p���X�S�x��  (1/0)*/
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

/* 4�� partial_mapping��@�覡*/
/* BM.cpp */		extern pset_family partial_mapping(pset_family Totality, pset Supercube, MCDB& mcdb);

	/*1. ��������ơAPPBM�]��1�h��������ƪ�����  (�D�n�ϥγo��)*/
/* PPBM.cpp */ extern void partial_parallel_mapping_block(int f_begin, int f_end, int g_begin, int g_end, pset_family Totality, pset_family temp_Totality, pset_family temp_Totality2, pset Supercube, MCDB& mcdb);

	/*2. cm150a �Ƽƿn���t��(�������) window size = 2������*/
/* PPBM.cpp */ extern void partial_parallel_mapping_block_window(int f_begin, int f_end, int g_begin, int g_end, pset_family Totality, pset_family temp_Totality, pset_family temp_Totality2, pset Supercube, MCDB& mcdb, int window);
	
	/*3. ��paper����@���{*/
/* paper.cpp */ extern void paper_partial_mapping_block(int f_begin, int f_end, int g_begin, int g_end, pset_family Totality, pset_family temp_Totality, pset Supercube, MCDB& mcdb);


	/*4. �Ĥ@����}�᪺PM(����������ƪ�)*/  /*!!!�i�঳�ʤֳ����[�t�ޥ�!!!*/
/* BM.cpp */		extern void partital_map_and(pset_family Totality, pset_family temp_Totality, pset Supercube, pset literal_pos, pset Q_non_DC, pset p_and);
/* BM.cpp */		extern void partial_mapping_block(int f_begin, int f_end, int g_begin, int g_end, pset_family Totality, pset_family temp_Totality, pset Supercube, MCDB& mcdb);


/* ��J�S�x��*/
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

/* ��X�S�x��*/
/* cplusplus.cpp */ extern void matching_KSIG_output(pset p_out, TABLE& KSIG_output1, TABLE& KSIG_output2);
/* cplusplus.cpp */ extern void matching_KSIG_output_with_KSIG_input(pset p_out, vector<TABLE>& KSIG_output1, vector<TABLE>& KSIG_output2);
/* cplusplus.cpp */ extern void check_Unate_output(pset p_out, pPLA PLA1, pPLA PLA2);
/* cplusplus.cpp */ extern void run_OutputSIG_TEST(pPLA PLA1, pPLA PLA2);

/* pset �B��*/
/* cplusplus.cpp */ extern bool setp_implies_var(pset a, pset b, int var);
/* cplusplus.cpp */ extern bool set_and_var(pset r, pset a, pset b, int var);

/* mvcube�ˬd�X�k��*/
/* cplusplus.cpp */ extern bool mvcube_feas_check_notation1(pset p); /* �i�঳���D*/
/* cplusplus.cpp */ extern bool mvcube_feas_check_notation2(pset p);

/* �����h�l��mvcube*/
/* cplusplus.cpp */ extern void sf_removeRedundant(pset_family A, int skip_num);

#endif //MAINCPLUSPLUS_H