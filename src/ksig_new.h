/*
	計算PLA的K特徵值
	分為輸入與輸出兩種
	輸入的特徵值為輸出的詳細版K特徵值
	使用Boots 中的BigInteger 來儲存K特徵值的計算結果
*/

#ifndef KSIG_NEW_H
#define KSIG_NEW_H


#include <iostream>
#include <vector>
#include <time.h>
#include <string>
#include <thread>
#include <tbb/tbb.h>
#include <boost/multiprecision/cpp_int.hpp>

extern "C" {
#include "espresso.h"
}

#include "support.h"

using namespace std;
using namespace boost::multiprecision;
namespace mp = boost::multiprecision;
/*基本定義與結構*/
#define ullg unsigned long long 

#define varnum 16 //在binary-input的情況下，一個word會有16個input
#define low cube.last_word[cube.output] // Last output word 最後一個output的word
#define fow cube.first_word[cube.output] // First output word 第一個output的word
#define var_insert(set, e, value) (set[WHICH_WORD(e * 2)] |= (value+1) << WHICH_BIT(e * 2)) //在pset中第e個變數插入value, 2 -> 11, 1 -> 01, 0 -> 10
#define var_remove(set, e) (set[WHICH_WORD(e * 2)] &= ~(3 << WHICH_BIT(e * 2))) //移除pset中第e個變數


#define LastInputWord WHICH_WORD(NUMINPUTS * 2)

#define Bat 500 //用於對稱檢查的檢查點，用於檢查Sym-Table的對稱情況 sym使用

/*ksig使用*/
#define TABLE vector<ksig_value>

struct ksig_value {
	vector<cpp_int> distx;
	ksig_value(int n) : distx(n) { }
	ksig_value() : distx(1) { }

	ksig_value& operator= (const ksig_value a) {
		if (distx.size() == a.distx.size()) {
			for (int i = 0; i < distx.size(); i++) {
				distx[i] = a.distx[i];
			}
		}
		return *this;
	}

	ksig_value operator+= (const ksig_value a) {
		if (distx.size() == a.distx.size()) {
			for (int i = 0; i < distx.size(); i++) {
				distx[i] += a.distx[i];
			}
		}
		return *this;
	}


	bool operator== (const ksig_value& a) const {
		if (distx.size() != a.distx.size())
			return false;

		for (int i = 0; i < distx.size(); i++) {
			if (distx[i] != a.distx[i]) {
				return false;
			}
		}
		return true;
	}
	bool operator!= (const ksig_value& a) const {
		if (distx.size() != a.distx.size())
			return true;

		for (int i = 0; i < distx.size(); i++) {
			if (distx[i] != a.distx[i]) {
				return true;
			}
		}
		return false;
	}


	bool operator<(const ksig_value& a) const {
		if (distx.size() != a.distx.size()) {
			fprintf(stderr, "ksig_value compare error %d != %d\n", distx.size(), a.distx.size());
			exit(EXIT_FAILURE);
		}
		for (int i = 0; i < distx.size(); i++) {
			if (distx[i] != a.distx[i]) {
				return distx[i] < a.distx[i];
			}
		}
		return false;
	}

	bool operator>(const ksig_value& a) const {
		if (distx.size() != a.distx.size()) {
			fprintf(stderr, "ksig_value compare error %d != %d\n", distx.size(), a.distx.size());
			exit(EXIT_FAILURE);
		}
		for (int i = 0; i < distx.size(); i++) {
			if (distx[i] != a.distx[i]) {
				return distx[i] > a.distx[i];
			}
		}
		return false;
	}
};


void cal_ksignature(pPLA PLA, TABLE& KSO, vector<TABLE>& KSI, int N, int Mode);

#endif //KSIG_NEW_H