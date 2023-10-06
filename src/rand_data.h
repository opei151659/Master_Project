#ifndef RAND_DATA_H
#define RAND_DATA_H
/*
	打亂輸入測資，擴展input或是移除on-set off-set
*/

#include <iostream>
#include <vector>
#include <time.h>
#include <string>
extern "C" {
#include "espresso.h"   /* ESPRESSO.lib*/
}

class rand_data {
private:
	std::vector<int> input_order;
	std::vector<int> output_order;
public:
	bool isINP;		//input permutation
	bool isINPA;	//input phase assigment
	bool isOUTP;	//output permutation
	bool isOUTPA;	//output phase assigment


	rand_data() {
		isINP = false;
		isINPA = false;
		isOUTP = false;
		isOUTPA = false;
	}
	void set(int num_inputs, int num_outputs) {
		input_order.resize(num_inputs);
		for (int i = 0; i < num_inputs; i++)
			input_order[i] = i;
		output_order.resize(num_outputs);
		for (int i = 0; i < num_outputs; i++)
			output_order[i] = i;
	}
	void input_order_permutation() {
		srand((unsigned int)time(NULL));
		std::vector<int> temp_order = input_order;
		int index;
		for (int i = 0; i < input_order.size(); i++) {
			index = rand() % (input_order.size() - i);
			input_order[i] = temp_order[index];
			auto it = temp_order.begin() + index;
			*it = std::move(temp_order.back());
			temp_order.pop_back();
		}
	}
	void output_order_permutation() {
		srand((unsigned int)time(NULL));
		std::vector<int> temp_order = output_order;
		int index;
		for (int i = 0; i < output_order.size(); i++) {
			index = rand() % (output_order.size() - i);
			output_order[i] = temp_order[index];
			auto it = temp_order.begin() + index;
			*it = std::move(temp_order.back());
			temp_order.pop_back();
		}
	}
	void input_order_phase_assigment() {
		srand((unsigned int)time(NULL));
		for (int i = 0; i < input_order.size(); i++) {
			if (rand() % 2) {
				input_order[i] = -input_order[i];
			}
		}
	}
	void output_order_phase_assigment() {
		srand((unsigned int)time(NULL));
		for (int i = 0; i < output_order.size(); i++) {
			if (rand() % 2) {
				output_order[i] = -output_order[i];
			}
		}
	}


	void pla_input_permutation(pPLA PLA1, pPLA PLA2) {
		isINP = true;
		pset p, p2;
		int index_p;
		foreachi_set(PLA2->F, index_p, p) {
			p2 = GETSET(PLA1->F, index_p);
			for (int i = 0; i < NUMINPUTS; i++) {
				PUTINPUT(p, i, GETINPUT(p2, input_order[i]));
			}
		}
		foreachi_set(PLA2->R, index_p, p) {
			p2 = GETSET(PLA1->R, index_p);
			for (int i = 0; i < NUMINPUTS; i++) {
				PUTINPUT(p, i, GETINPUT(p2, input_order[i]));
			}
		}
	}
	void pla_output_permutation(pPLA PLA1, pPLA PLA2) {
		isOUTP = true;
		pset p, p2;
		int index_p;
		foreachi_set(PLA2->F, index_p, p) {
			p2 = GETSET(PLA1->F, index_p);
			for (int i = 0; i < NUMOUTPUTS; i++) {
				if (is_in_set(p2, NUMINPUTS * 2 + output_order[i]))
					set_insert(p, NUMINPUTS * 2 + i);
				else
					set_remove(p, NUMINPUTS * 2 + i);
				//PUTOUTPUT(p, i, 0); //GETOUTPUT(p2, output_order[i]));
			}
		}
		foreachi_set(PLA2->R, index_p, p) {
			p2 = GETSET(PLA1->R, index_p);
			for (int i = 0; i < NUMOUTPUTS; i++) {
				if (is_in_set(p2, NUMINPUTS * 2 + output_order[i]))
					set_insert(p, NUMINPUTS * 2 + i);
				else
					set_remove(p, NUMINPUTS * 2 + i);
				//PUTOUTPUT(p, i, 0); //GETOUTPUT(p2, output_order[i]));
			}
		}
	}
	void pla_input_phase_assigment(pPLA PLA) {
		isINPA = true;
		pset p, last;
		int temp;
		foreach_set(PLA->F, last, p) {
			for (int i = 0; i < NUMINPUTS; i++) {
				if (input_order[i] < 0) {
					temp = GETINPUT(p, i);
					if (temp == ONE) {
						PUTINPUT(p, i, ZERO);
					}
					else if (temp == ZERO) {
						PUTINPUT(p, i, ONE);
					}
				}
			}
		}
		foreach_set(PLA->R, last, p) {
			for (int i = 0; i < NUMINPUTS; i++) {
				if (input_order[i] < 0) {
					temp = GETINPUT(p, i);
					//printf("%d\n", temp);
					if (temp == ONE) {
						PUTINPUT(p, i, ZERO);
					}
					else if (temp == ZERO) {
						PUTINPUT(p, i, ONE);
					}
				}
			}
		}
	}
	/* 未經驗證*/
	void pla_output_phase_assigment(pPLA PLA) {
		isOUTPA = true;
		pset p, last;
		int temp;
		foreach_set(PLA->F, last, p) {
			for (int i = 0; i < NUMOUTPUTS; i++) {
				if (input_order[i] < 0) {
					temp = GETOUTPUT(p, i);
					if (temp == 1) {
						PUTOUTPUT(p, i, 0);
					}
					else if (temp == 0) {
						PUTOUTPUT(p, i, 1);
					}
				}
			}
		}
		foreach_set(PLA->R, last, p) {
			for (int i = 0; i < NUMOUTPUTS; i++) {
				if (input_order[i] < 0) {
					temp = GETOUTPUT(p, i);
					//printf("%d\n", temp);
					if (temp == 1) {
						PUTOUTPUT(p, i, 0);
					}
					else if (temp == 0) {
						PUTOUTPUT(p, i, 1);
					}
				}
			}
		}
	}

	void random(pPLA PLA1, pPLA PLA2, bool random_input_permutation, bool random_input_phase_assigment, bool random_output_permutation, bool random_output_phase_assigment) {
		if (random_input_permutation) {
			input_order_permutation();
			pla_input_permutation(PLA1, PLA2);
		}
		if (random_input_phase_assigment) {
			input_order_phase_assigment();
			pla_input_phase_assigment(PLA2);
		}
		if (random_output_permutation) {
			output_order_permutation();
			pla_output_permutation(PLA1, PLA2);
		}
		if (random_output_phase_assigment) {
			output_order_phase_assigment();
			pla_output_phase_assigment(PLA2);
		}
		

	}

	void extandInput(pPLA PLA1, pPLA PLA2, int extand_input_num) {
		int extand_input_size = extand_input_num - (NUMINPUTS * 2 + NUMOUTPUTS);
		//pPLA pla = new_PLA();
		//pla->F = sf_new(PLA->F->count, extand_input_num * 2 + NUMOUTPUTS);
		//pla->R = sf_new(PLA->R->count, extand_input_num * 2 + NUMOUTPUTS);
		//free_PLA(pla);
		pset_family F = sf_new(PLA1->F->count, extand_input_num * 2 + NUMOUTPUTS);
		pset_family R = sf_new(PLA1->R->count, extand_input_num * 2 + NUMOUTPUTS);
		pset_family D = sf_new(PLA1->F->count, extand_input_num * 2 + NUMOUTPUTS);

		/* 完全亂數*/
		//extandSF(F, PLA1->F, extand_input_num);
		//extandSF(R, PLA1->R, extand_input_num);
		
		/* 非亂數*/
		extandNotRandomSF(F, PLA1->F, extand_input_num);
	
		/* 非亂數*/
		extandNotRandomSF(R, PLA1->R, extand_input_num);

		
		sf_free(PLA1->F);
		sf_free(PLA1->R);
		sf_free(PLA2->F);
		sf_free(PLA2->R);


		PLA1->F = F;
		PLA1->R = R;
		PLA2->F = F;
		PLA2->R = R;
	}

	/* 完全亂數*/ /*!!!CSF對F作亂數則R要重新產生!!!*/
	void extandSF(pset_family new_sf, pset_family old_sf, int extand_input_num) {
		srand((unsigned int)time(NULL));
		pset p, new_p;
		new_p = set_new(new_sf->sf_size);
		int index_p, i, j;

		foreachi_set(old_sf, index_p, p) {
			set_clear(new_p, new_sf->sf_size);
			for (i = 0; i < NUMINPUTS; i++) {
				PUTINPUT(new_p, i, GETINPUT(p, i));
			}
			for (; i < extand_input_num; i++) {
				PUTINPUT(new_p, i, (rand()%3 + 1));
			}
			for (j = 0; j < NUMOUTPUTS; j++) {
				if (GETOUTPUT(p, j) > 0) {
					set_insert(new_p, extand_input_num*2 + j);
				}
			}
			sf_addset(new_sf, new_p);
		}
		set_free(new_p);
	}

	/* 非亂數*/ /*!!!CSF對F作亂數則R要重新產生!!!*/
	void extandNotRandomSF(pset_family new_sf, pset_family old_sf, int extand_input_num) {
		srand((unsigned int)time(NULL));
		pset p, new_p;
		new_p = set_new(new_sf->sf_size);
		int index_p, i, j;
		/* 手動新增要擴展的變數*/
		std::string extand_bits = "012221012221012221012221012221012221012221012221012221012221012221012221012221012221012221012221012221012221012221012221012221";


		foreachi_set(old_sf, index_p, p) {
			set_clear(new_p, new_sf->sf_size);
			for (i = 0; i < NUMINPUTS; i++) {
				PUTINPUT(new_p, i, GETINPUT(p, i));
			}
			for (j = 0; i < extand_input_num; i++, j++) {
				PUTINPUT(new_p, i, extand_bits[j]-'0');
			}
			for (j = 0; j < NUMOUTPUTS; j++) {
				if (GETOUTPUT(p, j) > 0) {
					set_insert(new_p, extand_input_num * 2 + j);
				}
			}
			sf_addset(new_sf, new_p);
		}
		set_free(new_p);
	}

	/* 輸出成PLA格式*/
	void writePLA(pPLA PLA, int numinputs, int numoutputs, std::string filename) {
		FILE* output_file = fopen(filename.c_str(), "w");
		if (output_file == NULL) {
			printf("can't open file %s\n", filename);
		}
		
		fprintf(output_file, ".i %d\n", numinputs);
		fprintf(output_file, ".o %d\n", numoutputs);
		fprintf(output_file, ".p %d\n", PLA->F->count + PLA->R->count);
		pset p, last;
		int i, j;
		foreach_set(PLA->F, last, p) {
			for (i = 0; i < numinputs; i++) {
				fprintf(output_file, "%c", "?01-"[GETINPUT(p, i)]);
			}
			fprintf(output_file, " ");
			for (j = 0; j < numoutputs; j++) {
				if (is_in_set(p, i*2+j))
					fprintf(output_file, "1");
				else
					fprintf(output_file, "0");
			}
			fprintf(output_file, "\n");
		}
		foreach_set(PLA->R, last, p) {
			for (i = 0; i < numinputs; i++) {
				fprintf(output_file, "%c", "?01-"[GETINPUT(p, i)]);
			}
			fprintf(output_file, " ");
			for (j = 0; j < numoutputs; j++) {
				if (is_in_set(p, i * 2 + j))
					fprintf(output_file, "1");
				else
					fprintf(output_file, "0");
			}
			fprintf(output_file, "\n");
		}
		fprintf(output_file, ".e\n");
		fclose(output_file);
	}

	/* 移除一定比例的on-set 與off-set*/ /* 這裡非亂數*/
	void removeSF(pPLA PLA1, pPLA PLA2, float remove_SF_rate, bool random) {
		if (remove_SF_rate != 0) {
			printf("before remove sf %d %d\n", PLA2->F->count, PLA2->R->count);
			if (random == true) {
				randomRemoveRate(PLA2->F, remove_SF_rate);
				randomRemoveRate(PLA2->R, remove_SF_rate);
			}
			else {
				removeRate(PLA2->F, remove_SF_rate);
				removeRate(PLA2->R, remove_SF_rate);
			}
			printf("after remove sf %d %d\n", PLA2->F->count, PLA2->R->count);
		}
	}

	void removeRate(pset_family sf, float rate) {
		int i, begin_remove = sf->count * (1.0 - rate);

		pset p;
		sf_active(sf);
		foreachi_set(sf, i, p) {
			if (i > begin_remove) {
				RESET(p, ACTIVE);
			}
		}
		sf_inactive(sf);
	}

	void randomRemoveRate(pset_family sf, float rate) {
		srand((unsigned int)time(NULL));
		pset p, last;
		sf_active(sf);
		int removecnt = sf->count * rate;
		int removebound = 100 * rate;
		if (removecnt >= sf->count) {
			printf("ERROR randomRemoveRate %d >= %d\n", removecnt, sf->count);
			exit(EXIT_FAILURE);
		}
		while (removecnt > 0) {
			foreach_active_set(sf, last, p) {
				if (rand() % 100 < removebound) {
					RESET(p, ACTIVE);
					removecnt--;
					if (removecnt <= 0)
						break;
				}
			}
		}
		sf_inactive(sf);
	}


	void display() {
		printf("input random permutation: %s\n", isINP ? "Y" : "N");
		printf("input random phase assigment: %s\n", isINPA ? "Y" : "N");
		printf("output random permutation: %s\n", isOUTP ? "Y" : "N");
		printf("output random phase assigment: %s\n", isOUTPA ? "Y" : "N");
		if (isINP || isINPA) {
			std::cout << "input order:";
			for (const auto& order : input_order) {
				std::cout << " " << order;
			}std::cout << std::endl;
		}
		else {
			printf("Original input order\n");
		}
		if (isOUTP || isOUTPA) {
			std::cout << "output order:";
			for (const auto& order : output_order) {
				std::cout << " " << order;
			}std::cout << std::endl;
		}
		else {
			printf("Original output order\n");
		}
	}
};

#endif //RAND_DATA_H