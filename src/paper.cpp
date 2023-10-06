/*
	根據老師的paper實作的原版布林比對演算法
*/

#include "mainCplusplus.h"

pset PaperandExactMvVar(pset r, pset p, bool& isEmpty, unsigned int first_word, unsigned int last_word, unsigned int first_bit, unsigned int last_bit);
pset PaperorExactMvVar(pset r, pset p, unsigned int first_word, unsigned int last_word, unsigned int first_bit, unsigned int last_bit);
pset PaperextractMvVar(pset r, pset p, unsigned int first_word, unsigned int last_word, unsigned int first_bit, unsigned int last_bit);
void Paperexpand_sf_Q(pset r, pset a);


//#define test_onlyMCDB

void paper_partial_mapping_block(int f_begin, int f_end, int g_begin, int g_end, pset_family Totality, pset_family temp_Totality, pset Supercube, MCDB& mcdb) {
	clock_t memory = clock();
	bool isRedundant, isEmpty;
	pset_family mvcubes = sf_new(1, mvcube.mv_size);
	pset p, p2, last, last2, p_lit0, p_lit1, p_and, p_f, Q_non_DC, q_non_DC, temp_Supercube, p_get_var, lit01, lit10, litQ, litq, onset, offset;
	temp_psets t_ps(13);
	t_ps.psets[0] = set_new(mcdb.P->sf_size * 32); p_lit0 = t_ps.psets[0];
	t_ps.psets[1] = set_new(mcdb.P->sf_size * 32); p_lit1 = t_ps.psets[1];
	t_ps.psets[2] = set_new(mvcube.mv_size); p_and = t_ps.psets[2];
	t_ps.psets[3] = set_new(mvcube.mv_size); temp_Supercube = t_ps.psets[3];
	t_ps.psets[4] = set_new(mvcube.mv_var_size); p_get_var = t_ps.psets[4];
	t_ps.psets[5] = set_new(mvcube.num_mv_var * 2); lit01 = t_ps.psets[5];
	t_ps.psets[6] = set_new(mvcube.num_mv_var * 2); lit10 = t_ps.psets[6];
	t_ps.psets[7] = set_new(mvcube.mv_size); Q_non_DC = t_ps.psets[7];
	t_ps.psets[8] = set_new(mvcube.mv_size); q_non_DC = t_ps.psets[8];

	t_ps.psets[9] = set_new(mvcube.mv_size); pset TEST_Q = t_ps.psets[9];
	t_ps.psets[10] = set_new(mvcube.mv_size); pset TEST_q = t_ps.psets[10];
	t_ps.psets[11] = set_new(mvcube.num_mv_var * 2); litQ = t_ps.psets[11];
	t_ps.psets[12] = set_new(mvcube.num_mv_var * 2); litq = t_ps.psets[12];
	unsigned long long cnt_set_and = 0;
	
	clock_t onlyOLDMCDB_begin = clock(), mcdb_begin, and_begin, redundant_begin, smart_begin, OutInSec_begin;
	double total_and = 0, total_redundant = 0, total_mcdb = 0, total_smart = 0, total_OutInSec = 0;
	printf("memory %.3lf\n", (double)(clock() - memory) / (double)CLOCKS_PER_SEC);
	memory = clock();
	paper_MCDB oldMCDB(mvcube.mv_size, mvcube.mv_size);

	printf("memory mcdb %.3lf\n", (double)(clock() - memory) / (double)CLOCKS_PER_SEC);

	for (int f_on = f_begin; f_on < f_end; f_on++) {
		//mcdb.non_DC(lit01, lit10, mcdb.getQset(f_on));
		onset = mcdb.getPset(f_on);
		fprintf(stderr, "f_on: %d\r", f_on); 

			
		for (int g_off = g_begin; g_off < g_end; g_off++) {
			OutInSec_begin = clock();

			if (!mcdb.check_output_intersection_pair(f_on, g_off)) {
				total_OutInSec += (double)(clock() - OutInSec_begin) / (double)CLOCKS_PER_SEC;
				continue;
			}
			total_OutInSec += (double)(clock() - OutInSec_begin) / (double)CLOCKS_PER_SEC;
			offset = mcdb.getQset(g_off);
		
			mvcubes->count = 0;

			smart_begin = clock();	//smart mvcube;
			bool isEmpty;
			for (int var = 0; var < mvcube.num_mv_var; var++) {
				mcdb.non_DC(litQ, litq, offset);
				set_fill(p_and, mvcube.mv_size);
				pset Q_non_DC = GETSET(mcdb.ex_Q_non_DC, g_off);
				pset q_non_DC = GETSET(mcdb.ex_q_non_DC, g_off);
				if (GETINPUT(onset, var) == ONE) {
					//display_pset(litq, mvcube.mv_var_size);
					//PaperandExactMvVar(p_and, litq, isEmpty, mvcube.first_word[var], mvcube.last_word[var], mvcube.first_bit[var], mvcube.last_bit[var]);
					set_and_var(p_and, Supercube, Q_non_DC, var);
				}
				else if (GETINPUT(onset, var) == ZERO) {
					//display_pset(litQ, mvcube.mv_var_size);
					//PaperandExactMvVar(p_and, litQ, isEmpty, mvcube.first_word[var], mvcube.last_word[var], mvcube.first_bit[var], mvcube.last_bit[var]);
					set_and_var(p_and, Supercube, q_non_DC, var);
				}
				else {
					continue;
				}
				//set_and(p_and, p_and, Supercube);
				//if (mvcube_feas_check_notation2(p_and)) {
				set_and(p_and, p_and, Supercube);
				//display_pset(p_and, mvcube.mv_size);
				if (mvcube_feas_check_notation2(p_and)) {
					sf_addset(mvcubes, p_and);
				}
				//}
			}
			//printf("mvcubes %d\n", mvcubes->count);
			total_smart += (double)(clock() - smart_begin) / (double)CLOCKS_PER_SEC;
			
			if (mvcubes->count == 0) {
				continue;
			}

			continue;
			//display_sf_eachvar(mvcubes, mvcube.mv_var_size, mvcube.num_mv_var);
			//printf("\n\n");
			mcdb_begin = clock();
			if (oldMCDB.isSFAbsorbed(mvcubes)) {
				//printf("oldMCDB continue\n");
				//printf("( %d, %d)\n", f_on, g_off);
				oldMCDB.cnt_skip++;
				total_mcdb += (double)(clock() - mcdb_begin) / (double)CLOCKS_PER_SEC;
				continue;
			}
			total_mcdb += (double)(clock() - mcdb_begin) / (double)CLOCKS_PER_SEC;

			temp_Totality->count = 0;
			and_begin = clock();
			foreach_set(Totality, last, p) {
				foreach_set(mvcubes, last2, p2) {
					set_and(p_and, p, p2);
					if (mvcube_feas_check_notation2(p_and)) {
						sf_addset(temp_Totality, p_and);
					}
				}
			}
			total_and += (double)(clock() - and_begin) / (double)CLOCKS_PER_SEC;
			if (temp_Totality->count == 0) {
				printf("ERROR\n");
				exit(EXIT_FAILURE);
			}



			
#ifndef test_onlyMCDB
			/*printf("BEFORE\n");
			display_sf_eachvar(temp_Totality, mvcube.mv_var_size, mvcube.num_mv_var);*/
			redundant_begin = clock();

			if (temp_Totality->count > 10000) {
				fprintf(stderr, "%d\n", temp_Totality->count);
				sf_removeRedundant(temp_Totality, 0);
			}
			total_redundant += (double)(clock() - redundant_begin) / (double)CLOCKS_PER_SEC;
			/*printf("AFTER\n");
			display_sf_eachvar(temp_Totality, mvcube.mv_var_size, mvcube.num_mv_var);*/
			sf_copy2(Totality, temp_Totality);


			sf_or2(temp_Supercube, temp_Totality);
			set_copy(Supercube, temp_Supercube);
#endif

			

			

		}
		if (temp_Totality->count > 10000) {
			fprintf(stderr, "%d\n", temp_Totality->count);
			sf_removeRedundant(temp_Totality, 0);
		}
		/*if (temp_Totality->count > 30000)
			break;*/
			//printf("(%d, %d) %d %d\n", f_on, mcdb.TABLE_PQ->sf_size, Totality->count, set_count_ones(Supercube));

	}

	printf("OLD MCDB %d\n", oldMCDB.cnt_skip);
	printf("OLDMCDB cost %.3lf\n", (double)(clock() - onlyOLDMCDB_begin) / (double)CLOCKS_PER_SEC);

	printf("TOTAL_SMART_MCDB_AND_RR  %.3lf %.3lf  %.3lf  %.3lf %.3lf\n", total_OutInSec, total_smart, total_mcdb, total_and, total_redundant);

	//printf("cnt_set_and: %lld\n", cnt_set_and);
	//printf("%d %d\n", Totality->count, set_count_ones(Supercube));
	//sf_removeRedundant(Totality, 0);
}


pset PaperandExactMvVar(pset r, pset p, bool& isEmpty, unsigned int first_word, unsigned int last_word, unsigned int first_bit, unsigned int last_bit) {
	unsigned int val = 0;
	isEmpty = true;
	for (unsigned int M = first_word, P = 1; M <= last_word; M++, P++) {
		val |= (p[P] << first_bit);
		if (M == first_word)
			val |= ~(0xFFFFFFFF << first_bit);
		if (M == last_word)
			val |= (0xFFFFFFFF << (last_bit + 1));
		r[M] &= val;
		if ((BPI - first_bit) != 32)
			val = (p[P] >> (BPI - first_bit)); // p[P]中剩餘的bits
		else
			val = 0;
		if (isEmpty) {
			if (M == first_word && M == last_word) {
				if ((r[M] & (0xFFFFFFFF << first_bit) & (0xFFFFFFFF >> ((BPI - 1) - last_bit))) != 0)
					isEmpty = false;
			}
			else if (M == first_word) {
				if ((r[M] & (0xFFFFFFFF << first_bit)) != 0)
					isEmpty = false;
			}
			else if (M == last_word) {
				if ((r[M] & (0xFFFFFFFF >> ((BPI - 1) - last_bit))) != 0)
					isEmpty = false;
			}
			else {
				if (r[M] != 0)
					isEmpty = false;
			}
		}
	}

	return r;
}

pset PaperorExactMvVar(pset r, pset p, unsigned int first_word, unsigned int last_word, unsigned int first_bit, unsigned int last_bit) {
	unsigned int val = 0;
	for (unsigned int M = first_word, P = 1; M <= last_word; M++, P++) {
		if (P <= mvcube.num_mv_var / 32 + 1)
			val |= (p[P] << first_bit);
		r[M] |= val;
		if (BPI - first_bit != 32)
			val = ((unsigned int)(p[P]) >> (BPI - first_bit)); // p[P]中剩餘的bits
		else
			val = 0;
	}
	return r;
}

pset PaperextractMvVar(pset r, pset p, unsigned int first_word, unsigned int last_word, unsigned int first_bit, unsigned int last_bit) {
	unsigned int P = first_word, R = 1;
	unsigned int val;
	/*split words without last*/
	//printf("\n\n%d %d %d %d\n", first_word, last_word, first_bit, last_bit);
	for (; P <= last_word && R <= LOOP(r); P++, R++) {
		val = p[P] >> first_bit;
		if (P != last_word && first_bit != 0)
			val |= (p[P + 1] << (BPI - first_bit));
		if (R == LOOP(r))
			val &= ~(0xFFFFFFFF << (((BPI - first_bit) + (last_bit + 1)) % 32));
		r[R] = val;
		/*printf("P:"); display_binary(p[P]);
		printf("R:"); display_binary(r[R]);*/
	}
	return r;
}


/* 將sf 中的 pset a 擴展成size倍數的pset r  ex: a=123, 3 => r=123123123*/
void Paperexpand_sf_Q(pset r, pset a) {
	
	set_clear(r, mvcube.mv_size);
	for (int var = 0; var < mvcube.num_mv_var; var++) {
		//printf("first_word: %d\nlast_word: %d\nfirst_bit: %d\nlast_bit: %d\n", first_word[var], last_word[var], first_bit[var], last_bit[var]);
		PaperorExactMvVar(r, a, mvcube.first_word[var], mvcube.last_word[var], mvcube.first_bit[var], mvcube.last_bit[var]);
	}
}