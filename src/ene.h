/*
	對稱特爭值檢查的計算與布林比對時的應用
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

// 回傳pset_family
// sym table 使用 vector<unique_ptr<data_mutex>>>  data_mutex可以選擇cpp、tbb或是openmp的版本
// cpp_mutex tbb_mutex mp_mutex 共有三種mutex 初始為cpp的shared_mutex
// sym table 會連接 pset_family(指到同一point)
// 有終止條件
// 沒有前處理

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
	std::vector<std::unique_ptr<cpp_mutex>> ddd; //資料和mutex放一起
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
	pset terminal_cond;// 用於保存var是否還有sym(終止條件)
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
	紀錄分組內的var 有哪些
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

// group = -1 表示未被分群
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
	std::vector<E_NE_group_id> var_group; //紀錄每個var當前是被分在 E/NE/ENE中個別的哪一群
	int E_group_cnt;
	int NE_group_cnt;
	int ENE_group_cnt;
public:
	E_NE_group() : E_group_cnt(0), NE_group_cnt(0), ENE_group_cnt(0) {
		var_group.resize(NUMINPUTS);
	}
	void set_group_id(int& A_id, int& B_id, int& cnt, int type) {
		if (A_id == -1) { //A_id(無)
			if (B_id == -1) { //A_id(無)B_id(無)
				A_id = cnt++;
				B_id = A_id;
			}
			else { //A_id(無)B_id(有)
				A_id = B_id;
			}
		}
		else { //A_id(有)
			if (B_id == -1) { //A_id(有)B_id(無)
				B_id = A_id;
			}
			else { //A_id(有)B_id(有) 把所有與B_id同群的合進A_id
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
		if (type == X_symm) // 若E symmetry 和NE symmertry皆沒有則全部歸為同一群
			return;

		var_group[i].X_group = false; var_group[j].X_group = false;
		if (type == E_symm) { // E symmetry 一定是兩兩一群
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
		// 整理每個var對應到的X /NE /ENE群
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
		// 清除未使用的群id (清除未使用的id，確保n個群中 1~n群都是有用的)
		for (int i = group_NE.size() - 1; i >= 0; i--) {
			if (group_NE[i].vars.size() == 0) {
				auto it = group_NE.begin() + i;
				*it = std::move(group_NE.back());
				group_NE.pop_back();
			}
		}
		// 清除未使用的群id (清除未使用的id，確保n個群中 1~n群都是有用的)
		for (int i = group_ENE.size() - 1; i >= 0; i--) {
			if (group_ENE[i].vars.size() == 0) {
				auto it = group_ENE.begin() + i;
				*it = std::move(group_ENE.back());
				group_ENE.pop_back();
			}
		}
	}
	/* 若NE間有E的關係則NE合併*/
	void NE_E_NE_grouping(pset_family sf, std::vector<group_vars>& group_E, std::vector<group_vars>& group_NE) {
		temp_psets t_ps(1);
		t_ps.psets[0] = set_new(NUMINPUTS * 2);
		pset p = t_ps.psets[0];
		std::vector<bool> var_grouped(NUMINPUTS, false);
		std::vector<bool> remove_group(group_NE.size(), false);
		// 檢查同個NEgroup中是否與任意var都有E的關係
		for (int i = 0; i < group_NE.size(); i++) {
			set_fill(p, NUMINPUTS * 2);
			for (int j = 0; j < group_NE[i].vars.size(); j++) {
				if (var_grouped[group_NE[i].vars[j]]) {
					// group中有任一個var已被合併，表示整個group都已被合併
					remove_group[i] = true;
					break;
				}
				set_and(p, p, GETSET(sf, group_NE[i].vars[j]));
			}
			if (remove_group[i])
				continue;
			for (int j = 0; j < NUMINPUTS; j++) {
				// 若所有都與var有E的關係，則將var合併到同的NE_group;
				if (GETINPUT(p, j) == E_symm || GETINPUT(p, j) == ENE_symm) {
					group_NE[i].vars.push_back(j);
					var_grouped[j] = true;
				}
			}
		}
		// 移除被合併的NEgroup
		for (int i = remove_group.size() - 1; i >= 0; i--) {
			if (remove_group[i]) {
				auto it = group_NE.begin() + i;
				*it = std::move(group_NE.back());
				group_NE.pop_back();
			}
		}
		// 移除被合併的Egroup
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
pset set_dc(pset r, pset a);// 尋找pset a的dcpart，並儲存到r
int cdist123(pset a, pset b, int dist_on[2]);// 計算兩pset距離，以及距離產生在哪個var
pset remove_symmetry(pset r, pset one, pset two); // dist1移除辦法
void find_Symmetry(pset_family A); //用於看var symmetry關係
bool check_symmetry(pset tc, pset_family A);// 查詢symmetry table是否有var已經沒有symmetry
// debug用
void binary_out(unsigned int a); //用於查看pset中的word為何
void printSym(pset_family A);// 輸出帶有編號的symmetry table
void printftriangle(pset_family A);// 只輸出上三角形
// 只有single thread能用
void flip_Sym(pset_family A);// 翻轉symmetry table到上三角形
#endif ENE_H