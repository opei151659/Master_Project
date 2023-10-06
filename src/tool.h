/*
	常用的pset計算的函式，類似於epsresso 中的set.c
	包含顯示pset與pset family (sf)的函式
*/
#ifndef TOOL_H
#define TOOL_H
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
//#include <utility> //pair

//#include "support.h"
#define varnum BPI/2

extern "C" {
#include "espresso.h"   /* ESPRESSO.lib*/
}

#define set_full_DISJOINT(size)	set_fill_disjoint(ALLOC(unsigned int, SET_SIZE(size)), size)

/*
	pset 的暫存，結束使用會自動釋放空間
*/
class temp_psets {
public:
	std::vector<pset> psets;
	// Default constructor
	temp_psets() {};
	temp_psets(int _size) {
		psets.resize(_size);
	};
	// Destructor
	~temp_psets() {
		//printf("temp_psets deconstruct start \n");
		for (int i = 0; i < psets.size(); i++) {
			//printf("%d %d\n", i, psets.size());
			if (psets[i] != nullptr) {
				//printf("release memory %d\n", i);
				set_free(psets[i]);
				
			}
		}
		psets.clear();
		psets.shrink_to_fit();
		//printf("temp_psets deconstruct finish \n");
	};
	// copy constructor
	temp_psets(const temp_psets& tps) {
		psets = tps.psets;
	}
	// move constructor
	temp_psets(temp_psets&& other) noexcept : psets(std::move(other.psets)) {};
	temp_psets& operator=(const temp_psets& tps) {
		psets.clear();
		psets.shrink_to_fit();
		psets = tps.psets;
		return *this;
	}
};
	
	int cdist1(pset a, pset b);
	char getvar(pset a, int position);// 00 -> '?' 10 -> '0' 01 -> '1' 11 -> '2'
	bool OutputIntersection(pset A, pset B);
	pset_family sf_inverse(pset_family sf);
	pset_family my_sf_copy(pset_family R, pset_family A);
	pset_family sf_copy2(pset_family R, pset_family A);
	pset_family sf_copy3(pset_family R, pset_family A);
	pset_family my_sf_append(pset_family A, pset_family B);
	pset sf_or2(pset or , pset_family A);
	pset or_exact_MvVar(pset r, pset p, int p_w_size, unsigned int first_word, unsigned int last_word, unsigned int first_bit, unsigned int last_bit);
	pset get_mvcube_var(pset r, pset p, int p_w_size, unsigned int first_word, unsigned int last_word, unsigned int first_bit, unsigned int last_bit);
	pset rshilft_or_disjoint(pset r, pset p);
	pset set_not(pset r, pset a);
	pset set_inverse_disjoint(pset r, pset p);
	int set_count_ones(pset p);
	void pset_grouping(pset p, int var_size, int var_num, bool display);
	pset string2pset(std::string str);
	pset_family sf_new2(int num, int size);
	void sf_free2(pset_family A);
	pset set_fill_disjoint(pset r, int size);


	void cube_setup(cube_struct& cube);
	void setdown_cube(cube_struct& cube);

	void display_pset(pset p, int _size);
	void display_pset_eachvar(pset p, int var_size, int var_num);
	void display_binary(int p);
	void display_sf(pset_family sf);
	void display_sf_eachvar(pset_family sf, int _size, int _num);
	void display_pset_01X(pset p, int input_num, int output_num);


	
#endif