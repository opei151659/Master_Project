#include "parallel_tool.h"


void local_intcpy(unsigned int* d, unsigned int* s, long int n)
{
	int i;
	for (i = 0; i < n; i++) {
		*d++ = *s++;
	}
}

/*
	mvcube �X�k���ˬd
	1. �C�@��var(row) �ܤ֭n������@��var(column)
		a. �z�L�h�ƨC��var(row) 1 ���Ӽ�
	2. �C�@��var(column) �ܤ֭n�Q�@��var(row) ������

	��ܪk��aa'bb'cc'�ɭn���Na�Ma'��������shift or�b�@�_
	ex aa'bb'cc' = 101101 => 010101

	�����ˬd��var�O�_empty�e��and�ɤw�ˬd�L

	return false ��ܤ��X�k
*/
//#define free_sf_return_false(sf) {sf_free(sf); return false;}
//#define free_sf_return_false(sf) {return false;}
//#define mvcube_feas_check_debug
pset* parallel_sf_sort(pset_family A, qsort_compare_func compare);
pset_family parallel_sf_unlist(pset_family R, pset* A1);
bool parallel_mvcube_feas_check_notation2(pset p, pset_family sf_vars, pset_family sf_sort, pset p_var_t, pset var_supercube, pset solution) {
	//parallel_temp_psets tpsets(3, memory_mtx);

	//pset_family sf_vars = sf_new(mvcube.num_mv_var, mvcube.mv_var_size);
	sf_vars->count = 0;
	set_clear(var_supercube, mvcube.mv_var_size);
	set_clear(solution, mvcube.mv_var_size);
	for (int var = 0; var < mvcube.num_mv_var; var++) {
		/* �N���bit�X��*/
		get_mvcube_var(p_var_t, p, mvcube.mv_var_word_size, mvcube.first_word[var], mvcube.last_word[var], mvcube.first_bit[var], mvcube.last_bit[var]);
		rshilft_or_disjoint(p_var_t, p_var_t); // (p | (p>>1)) & DISJOINT
		PUTSIZE(p_var_t, set_count_ones(p_var_t));
		sf_addset(sf_vars, p_var_t);
		/* �p��var_supercube*/
		set_or(var_supercube, var_supercube, p_var_t);
	}
	/* supercube �ˬd*/
	if (set_count_ones(var_supercube) != mvcube.num_mv_var) {
#ifdef mvcube_feas_check_debug
		printf("var supercube infeasible\n");
		display_pset_eachvar(p, mvcube.mv_var_size, mvcube.num_mv_var);
#endif
		//free_sf_return_false(sf_vars);
		return false;
	}
	//#define SF_SORT_assend(P) {P = sf_unlist(sf_sort(P, (qsort_compare_func)ascend), P->count, P->sf_size);}
		/* �p��j�ƧǡA�@�˪��|�\�@�_*/
		//sf_sort(sf_vars, (qsort_compare_func)ascend);
		//SF_SORT_assend(sf_vars);
	sf_vars = parallel_sf_unlist(sf_sort, parallel_sf_sort(sf_vars, (qsort_compare_func)ascend));
	pset last, tp, pre_tp;
	int index_p;
	pre_tp = p_var_t; /* ���F����warning �H�K��l�ơA���|��ڥΨ�p_var_t���Ŷ�*/

	sf_active(sf_vars);
	int pre_cnt = 0;
	/* ��������var������ۦP(mvcube)���Ӽ�*/
	int cnt_same = 1;
	/* ��������var�w�g���t��*/
	int num_unmatch_var = mvcube.num_mv_var;

	foreachi_set(sf_vars, index_p, tp) {
		/* �Y�w���Ѥw�]�ttp�����htp���L��*/
		if (setp_implies(tp, solution)) {
#ifdef mvcube_feas_check_debug
			printf("�w���Ѥw�]�ttp�����htp���L�� infeasible\n");
			display_pset_eachvar(p, mvcube.mv_var_size, mvcube.num_mv_var);
#endif
			//free_sf_return_false(sf_vars);
			return false;
		}
		if (SIZE(tp) == 1) {
			/* 1. �@��var�u������@�ӸѮɡA�Ӹѵ����w������*/
			set_or(solution, solution, tp);
			RESET(tp, ACTIVE);
			num_unmatch_var--;
		}
		else if (SIZE(tp) < mvcube.num_mv_var) {
			/* �Nsolution ������inverse*/
			set_inverse_disjoint(p_var_t, solution);
			//printf("solution:      "); display_pset(solution, mvcube.mv_var_size);
			//printf("original tp:   "); display_pset(tp, mvcube.mv_var_size);
			//printf("solution inve: "); display_pset(p_var_t, mvcube.mv_var_size);

			/* tp�h����e�w����*/
			set_and(tp, p_var_t, tp);
			//printf("removed tp:    "); display_pset(tp, mvcube.mv_var_size);
			if (pre_cnt == SIZE(tp)) {
				if (setp_equal(tp, pre_tp)) {
					cnt_same++;
					/* n��var�����ۦP��n�ӸѮɡA�ѥu�i��������n�ӸѡA�G�i�����w������*/
					if (cnt_same == SIZE(tp)) {
						set_or(solution, solution, tp);
						for (int i = 0; i < cnt_same; i++) {
							RESET(GETSET(sf_vars, index_p - i), ACTIVE);
						}
						num_unmatch_var -= cnt_same;
					}
				}
				else {
					cnt_same = 1;
				}
			}
			else {
				pre_cnt = SIZE(tp);
				pre_tp = tp;
				cnt_same = 1;
			}
		}
		/* �Ѿl�ҬOfull cube ���L�����*/
		else {
			break;
		}
		/*2. �|���t�諸var��|���t�諸���٦h�A�h��ܵL��*/
		if (num_unmatch_var > mvcube.num_mv_var - set_count_ones(solution)) {
#ifdef mvcube_feas_check_debug
			printf("�|���t�諸var��|���t�諸���٦h infeasible\n");
			display_pset_eachvar(p, mvcube.mv_var_size, mvcube.num_mv_var);
#endif
			//free_sf_return_false(sf_vars);
			return false;
		}
	}

	/* �٦�remove�������i�H��*/

	//sf_free(sf_vars);
	return true;
}