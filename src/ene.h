/*
	��ٯS�����ˬd���p��P���L���ɪ�����
*/
#ifndef ENE_H
#define ENE_H
#include <iostream>
#include <vector>
#include <tbb/spin_mutex.h>
#include <tbb/parallel_for.h>
#include <omp.h>
#include <shared_mutex>
#include <time.h>
#include <functional> //mem_fn

#include "tool.h"

extern "C" {
#include "espresso.h"   /* ESPRESSO.lib*/
}

// �^��pset_family
// sym table �ϥ� vector<unique_ptr<data_mutex>>>  data_mutex�i�H���cpp�Btbb�άOopenmp������
// cpp_mutex tbb_mutex mp_mutex �@���T��mutex ��l��cpp��shared_mutex
// sym table �|�s�� pset_family(����P�@point)
// ���פ����
// �S���e�B�z

#define THREADNUM 24
#define Bat 500
#define print_sym 0

class ENECLASS { //using vdm(vector data_mutex) for symmetry in CPP
private:
	class cpp_mutex {
	private:
		mutable std::shared_mutex dm;
		pset p;
	public:
		cpp_mutex(pset p_) {
			p = p_;
		}
		void to_and(pset k) {
			std::lock_guard<std::shared_mutex> lock(dm);
			set_and(p, p, k);
		}
		void to_remove(int index) {
			std::lock_guard<std::shared_mutex> lock(dm);
			set_remove(p, index);
		}
		pset get_pset() {
			return p;
		}
		void print_set() {
			std::cout << pbv1(p, NUMINPUTS * 2) << std::endl;
		}
	};
	class tbb_mutex {
	private:
		tbb::spin_mutex dm;
		pset p;
	public:
		tbb_mutex(pset p_) {
			p = p_;
		}
		void to_and(pset k) {
			tbb::spin_mutex::scoped_lock lock(dm);
			set_and(p, p, k);
		}
		void to_remove(int index) {
			tbb::spin_mutex::scoped_lock lock(dm);
			set_remove(p, index);
		}
		pset get_pset() {
			return p;
		}
		void print_set() {
			std::cout << pbv1(p, NUMINPUTS * 2) << std::endl;
		}
	};
	class mp_mutex {
	private:
		omp_lock_t dm;
		pset p;
	public:
		mp_mutex(pset p_) {
			omp_init_lock(&dm);
			p = p_;
		}
		void to_and(pset k) {
			omp_set_lock(&dm);
			set_and(p, p, k);
			omp_unset_lock(&dm);
		}
		void to_remove(int index) {
			omp_set_lock(&dm);
			set_remove(p, index);
			omp_unset_lock(&dm);
		}
		pset get_pset() {
			return p;
		}
		void print_set() {
			std::cout << pbv1(p, NUMINPUTS * 2) << std::endl;
		}
	};
	pset_family F;
	pset_family R;
	pset_family sym;
	std::vector<std::unique_ptr<cpp_mutex>> ddd; //��ƩMmutex��@�_
	int is_not_sym;
	unsigned long long maxpset;
	std::atomic<unsigned long long> aui;
	int chunk_size;
	unsigned long long col, row;
	std::vector<unsigned long long> chunk_table;

public:
	ENECLASS(pset_family F_, pset_family R_) : ddd(NUMINPUTS) {
		F = F_;
		R = R_;
		sym = sf_new(NUMINPUTS, NUMINPUTS * 2);
		for (int i = 0; i < NUMINPUTS; i++) {
			sf_addset(sym, set_full(NUMINPUTS * 2));
			ddd[i].reset(new cpp_mutex(sym->data + sym->wsize * i));
		}
		terminal_cond = set_new(NUMINPUTS);
		is_not_sym = false;
		chunk_size = 500;
		col = ceil((double)R->count / chunk_size);
		row = ceil((double)F->count / chunk_size);
		maxpset = col * row;
		chunk_table.resize(maxpset);
		aui.store(0);
	}
	pset terminal_cond;// �Ω�O�svar�O�_�٦�sym(�פ����)
	int chunk;
	pset_family cpp_job();
	pset_family tbb_job();
	pset_family openmp_job();
	void for_symmetry1();
	void for_symmetry2(unsigned long long pairid);
	void outputSym();
	void flip();
	void allflip();
	void printSym();
	void chunk_num();
	bool check_symmetry();
	~ENECLASS() {
		FREE(terminal_cond);
		sf_free(sym);
	}
};

/*
	�������դ���var ������
*/
struct group_vars {
	int group_id;
	int group_size;
	std::vector<int> vars;
	group_vars() : group_id(-1), group_size(0) {};
	bool operator>(const group_vars& gvs) const {
		if (this->vars.size() > gvs.vars.size())
			return true;
		return false;
	}
	void display() {
		printf("(");
		for (int i = 0; i < vars.size(); i++) {
			printf(" %d", vars[i]);
		}
		printf(")\n");
	}
};

// group = -1 ��ܥ��Q���s
struct E_NE_group_id {
	int E_group;
	int NE_group;
	int ENE_group;
	bool X_group;
	E_NE_group_id() : E_group(-1), NE_group(-1), ENE_group(-1), X_group(true) {};

};


#define E_symm 2
#define NE_symm 1
#define ENE_symm 3
#define X_symm 0
class E_NE_group {
	std::vector<E_NE_group_id> var_group; //�����C��var��e�O�Q���b E/NE/ENE���ӧO�����@�s
	int E_group_cnt;
	int NE_group_cnt;
	int ENE_group_cnt;
public:
	E_NE_group() : E_group_cnt(0), NE_group_cnt(0), ENE_group_cnt(0) {
		var_group.resize(NUMINPUTS);
	}
	void set_group_id(int& A_id, int& B_id, int& cnt, int type) {
		if (A_id == -1) { //A_id(�L)
			if (B_id == -1) { //A_id(�L)B_id(�L)
				A_id = cnt++;
				B_id = A_id;
			}
			else { //A_id(�L)B_id(��)
				A_id = B_id;
			}
		}
		else { //A_id(��)
			if (B_id == -1) { //A_id(��)B_id(�L)
				B_id = A_id;
			}
			else { //A_id(��)B_id(��) ��Ҧ��PB_id�P�s���X�iA_id
				if (type == NE_symm)
					for (int k = 0; k < NUMINPUTS; k++) {
						if (var_group[k].NE_group == B_id) {
							var_group[k].NE_group = A_id;
						}
					}
				else if (type == ENE_symm)
					for (int k = 0; k < NUMINPUTS; k++) {
						if (var_group[k].ENE_group == B_id) {
							var_group[k].ENE_group = A_id;
						}
					}
				else {
					printf("ERROR!!!!!\n");
				}
			}
		}
	}
	void grouping(int i, int j, int type, std::vector<group_vars>& group_E) {
		if (type == X_symm) // �YE symmetry �MNE symmertry�ҨS���h�����k���P�@�s
			return;

		var_group[i].X_group = false; var_group[j].X_group = false;
		if (type == E_symm) { // E symmetry �@�w�O���@�s
			//printf("TEST %d %d\n", group_E.size(), E_group_cnt);
			group_E[E_group_cnt].vars[0] = i;
			group_E[E_group_cnt++].vars[1] = j;
		}
		else if (type == NE_symm) {
			set_group_id(var_group[i].NE_group, var_group[j].NE_group, NE_group_cnt, NE_symm);
		}
		else if (type == ENE_symm) {
			set_group_id(var_group[i].ENE_group, var_group[j].ENE_group, ENE_group_cnt, ENE_symm);
		}

		if (E_group_cnt > NUMINPUTS * (NUMINPUTS - 1) / 2 || NE_group_cnt > NUMINPUTS) {
			printf("(E, NE, ENE) = (%d, %d, %d)\n", E_group_cnt, NE_group_cnt, ENE_group_cnt);
			printf("grouping ERROR!!!!!!\n");
			exit(EXIT_FAILURE);
		}
	}
	void classify(std::vector<group_vars>& group_E, std::vector<group_vars>& group_NE, std::vector<group_vars>& group_ENE, std::vector<int>& group_X) {
		group_E.resize(E_group_cnt);
		group_E.shrink_to_fit();
		group_NE.resize(NE_group_cnt);
		group_NE.shrink_to_fit();
		group_ENE.resize(ENE_group_cnt);
		group_ENE.shrink_to_fit();
		// ��z�C��var�����쪺X /NE /ENE�s
		for (int i = 0; i < NUMINPUTS; i++) {
			if (var_group[i].X_group)
				group_X.push_back(i);
			/*else if (var_group[i].E_group >= 0)
				group_E[var_group[i].E_group].vars.push_back(i);*/
			else if (var_group[i].NE_group >= 0)
				group_NE[var_group[i].NE_group].vars.push_back(i);
			else if (var_group[i].ENE_group >= 0)
				group_ENE[var_group[i].ENE_group].vars.push_back(i);
			/*else
				printf("Classify ERROR!!!!!!!");*/

		}
		// �M�����ϥΪ��sid (�M�����ϥΪ�id�A�T�On�Ӹs�� 1~n�s���O���Ϊ�)
		for (int i = group_NE.size() - 1; i >= 0; i--) {
			if (group_NE[i].vars.size() == 0) {
				auto it = group_NE.begin() + i;
				*it = std::move(group_NE.back());
				group_NE.pop_back();
			}
		}
		// �M�����ϥΪ��sid (�M�����ϥΪ�id�A�T�On�Ӹs�� 1~n�s���O���Ϊ�)
		for (int i = group_ENE.size() - 1; i >= 0; i--) {
			if (group_ENE[i].vars.size() == 0) {
				auto it = group_ENE.begin() + i;
				*it = std::move(group_ENE.back());
				group_ENE.pop_back();
			}
		}
	}
	/* �YNE����E�����Y�hNE�X��*/
	void NE_E_NE_grouping(pset_family sf, std::vector<group_vars>& group_E, std::vector<group_vars>& group_NE) {
		temp_psets t_ps(1);
		t_ps.psets[0] = set_new(NUMINPUTS * 2);
		pset p = t_ps.psets[0];
		std::vector<bool> var_grouped(NUMINPUTS, false);
		std::vector<bool> remove_group(group_NE.size(), false);
		// �ˬd�P��NEgroup���O�_�P���Nvar����E�����Y
		for (int i = 0; i < group_NE.size(); i++) {
			set_fill(p, NUMINPUTS * 2);
			for (int j = 0; j < group_NE[i].vars.size(); j++) {
				if (var_grouped[group_NE[i].vars[j]]) {
					// group�������@��var�w�Q�X�֡A��ܾ��group���w�Q�X��
					remove_group[i] = true;
					break;
				}
				set_and(p, p, GETSET(sf, group_NE[i].vars[j]));
			}
			if (remove_group[i])
				continue;
			for (int j = 0; j < NUMINPUTS; j++) {
				// �Y�Ҧ����Pvar��E�����Y�A�h�Nvar�X�֨�P��NE_group;
				if (GETINPUT(p, j) == E_symm || GETINPUT(p, j) == ENE_symm) {
					group_NE[i].vars.push_back(j);
					var_grouped[j] = true;
				}
			}
		}
		// �����Q�X�֪�NEgroup
		for (int i = remove_group.size() - 1; i >= 0; i--) {
			if (remove_group[i]) {
				auto it = group_NE.begin() + i;
				*it = std::move(group_NE.back());
				group_NE.pop_back();
			}
		}
		// �����Q�X�֪�Egroup
		bool g_remove;
		for (int i = group_E.size() - 1; i >= 0; i--) {
			g_remove = false;
			for (int j = 0; j < group_E[i].vars.size(); j++) {
				if (var_grouped[group_E[i].vars[j]]) {
					g_remove = true;
					break;
				}
			}
			if (g_remove) {
				auto it = group_E.begin() + i;
				*it = std::move(group_E.back());
				group_E.pop_back();
			}
		}
	}
	void display() {
		for (int i = 0; i < NUMINPUTS; i++) {
			printf("var %d:", i);
			if (var_group[i].X_group)
				printf("(X)\n");
			else
				printf("(%3d, %3d, %3d)\n", var_group[i].E_group, var_group[i].NE_group, var_group[i].ENE_group);
		}
	}
};


void ene(pset_family ENE, pset_family F, pset_family R, int option);
pset set_dc(pset r, pset a);// �M��pset a��dcpart�A���x�s��r
int cdist123(pset a, pset b, int dist_on[2]);// �p���pset�Z���A�H�ζZ�����ͦb����var
pset remove_symmetry(pset r, pset one, pset two); // dist1������k
void find_Symmetry(pset_family A); //�Ω��var symmetry���Y
bool check_symmetry(pset tc, pset_family A);// �d��symmetry table�O�_��var�w�g�S��symmetry
// debug��
void binary_out(unsigned int a); //�Ω�d��pset����word����
void printSym(pset_family A);// ��X�a���s����symmetry table
void printftriangle(pset_family A);// �u��X�W�T����
// �u��single thread���
void flip_Sym(pset_family A);// ½��symmetry table��W�T����
#endif ENE_H