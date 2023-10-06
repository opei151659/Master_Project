/*
	餘因子特徵值檢查的實作與比對， 實作方式為計算所有的onset， 較慢， 若要加速則可以使用遞迴來實作
*/
#ifndef COFACTOR_H
#define COFACTOR_H
#include <iostream>
#include <vector>
#include <time.h>
#include <algorithm>
#include "tool.h"
#include "BigInteger.h"

extern "C" {
#include "espresso.h"   /* ESPRESSO.lib*/
}

struct var_cofactor {
	int var;
	
	vector <BigInt> positive;
	vector <BigInt> negative;
	var_cofactor() : var(-1) {
	}
	void init(int numoutput){
		positive.resize(numoutput);
		negative.resize(numoutput);
	}
	bool operator==(const var_cofactor& other) const {
		for (int i = 0; i < positive.size(); i++) {
			if (positive[i] != other.positive[i])
				return false;
			if (negative[i] != other.negative[i])
				return false;
		}
		return true;
	}
	bool compareSamePos(const var_cofactor& other) const{
		for (int i = 0; i < positive.size(); i++) {
			if (positive[i] != other.positive[i])
				return false;
			if (negative[i] != other.negative[i])
				return false;
		}
		return true;
	}
	bool compareSameNeg(const var_cofactor& other) const {
		for (int i = 0; i < positive.size(); i++) {
			if (positive[i] != other.negative[i])
				return false;
			if (negative[i] != other.positive[i])
				return false;
		}
		return true;
	}

	var_cofactor& operator=(const var_cofactor& rhs) {
		var = rhs.var;
		positive = rhs.positive;
		negative = rhs.negative;
		return *this;
	}
};

//bool compare_positive(const var_cofactor& p1, const var_cofactor& p2) {
//	if (p1.positive < p2.positive) {
//		return true;
//	}
//	else if (p1.positive == p2.positive) {
//		return p1.negative < p2.negative;
//	}
//	return false;
//}

class Cofactor {
private:
	pset p, last;
	pset_family F;
	/* count DC num*/
	int set_count_inputs_DC(pset p_) {
		int cnt = 0, i = cube.last_word[cube.num_binary_vars - 1];
		int w = 1;
		unsigned int a;
		for (; w < cube.last_word[cube.num_binary_vars]; w++) {
			a = ((p_[w] >> 1 & p_[w]) & DISJOINT);
			cnt += count_ones(a);
		}
		a = ((p_[w] >> 1 & p_[w]) & DISJOINT) & cube.binary_mask[w];
		cnt += count_ones(a);

		return cnt;
	};
public:
	std::vector<var_cofactor> var_cof; /* 記錄各個variable 在cofactor下在onset 中的minterm 個數*/

	Cofactor(pset_family F_) {
		F = sf_new(1, F_->sf_size); 
		sf_copy3(F, F_); 
		F = make_disjoint(F); 
		F = sf_contain(F);
		foreach_set(F, last, p) {
			PUTSIZE(p, set_count_inputs_DC(p));
			//printf("%d ", SIZE(p));  display_pset_01X(p, NUMINPUTS , NUMOUTPUTS);
		}
		var_cof.resize(NUMINPUTS);
		for (int i = 0; i < NUMINPUTS; i++) {
			var_cof[i].var = i;
			var_cof[i].init(NUMOUTPUTS);
		}
	};
	void cal_cofactor() {
		int var;
		foreach_set(F, last, p) {
			for (int i = 0; i < NUMINPUTS; i++) {
				var = GETINPUT(p, i);
				for (int j = 0; j < NUMOUTPUTS; j++) {
					if (GETOUTPUT(p, j)) {
						if (var == ONE) {
							var_cof[i].positive[j] += pow(2, SIZE(p));
						}
						else if (var == ZERO) {
							var_cof[i].negative[j] += pow(2, SIZE(p));
						}
						else {
							var_cof[i].positive[j] += pow(2, (SIZE(p) - 1));
							var_cof[i].negative[j] += pow(2, (SIZE(p) - 1));
						}
					}
				}
			}
		}
	}

	void matching_cofactor(pset p, std::vector<var_cofactor>&var_cof2, int INPA, int notation) {
		for (int i = 0; i < var_cof.size(); i++) {
			for (int j = 0; j < var_cof2.size(); j++) {
				/* (pos, neg) == (pos2, neg2)*/
				if (var_cof[i].compareSamePos(var_cof2[j])) {
					if (notation == 1) {
						if (!INPA) set_insert(p, i * NUMINPUTS + j);
						else set_insert(p, i * (NUMINPUTS * 2) + j);
					}
					else if (notation == 2) {
						set_insert(p, i * (NUMINPUTS * 2) + (j * 2 + 1));
					}
				}
				if (INPA && var_cof[i].compareSameNeg(var_cof2[j])) {
					if (notation == 1) {
						set_insert(p, i * (NUMINPUTS * 2) + j + NUMINPUTS);
					}
					else if (notation == 2) {
						set_insert(p, i * (NUMINPUTS * 2) + (j * 2));
					}
				}
			}
		}
	}

	void display() {
		for (int i = 0; i < NUMINPUTS; i++) {
			for (int j = 0; j < NUMOUTPUTS; j++) {
				var_cof[i].positive[j].print_value(); printf(" ");
			}printf("\n");
			for (int j = 0; j < NUMOUTPUTS; j++) {
				var_cof[i].negative[j].print_value(); printf(" ");
			}printf("\n");
			
		}printf("\n");

	}

	~Cofactor() {};
};




#endif
