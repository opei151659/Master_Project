#include "parallel_tool.h"


void local_intcpy(unsigned int* d, unsigned int* s, long int n)
{
	int i;
	for (i = 0; i < n; i++) {
		*d++ = *s++;
	}
}

/*
	mvcube 合法性檢查
	1. 每一個var(row) 至少要對應到一個var(column)
		a. 透過去數每個var(row) 1 的個數
	2. 每一個var(column) 至少要被一個var(row) 對應到

	表示法為aa'bb'cc'時要先將a和a'的部分用shift or在一起
	ex aa'bb'cc' = 101101 => 010101

	不用檢查該var是否empty前面and時已檢查過

	return false 表示不合法
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
		/* 將兩個bit合併*/
		get_mvcube_var(p_var_t, p, mvcube.mv_var_word_size, mvcube.first_word[var], mvcube.last_word[var], mvcube.first_bit[var], mvcube.last_bit[var]);
		rshilft_or_disjoint(p_var_t, p_var_t); // (p | (p>>1)) & DISJOINT
		PUTSIZE(p_var_t, set_count_ones(p_var_t));
		sf_addset(sf_vars, p_var_t);
		/* 計算var_supercube*/
		set_or(var_supercube, var_supercube, p_var_t);
	}
	/* supercube 檢查*/
	if (set_count_ones(var_supercube) != mvcube.num_mv_var) {
#ifdef mvcube_feas_check_debug
		printf("var supercube infeasible\n");
		display_pset_eachvar(p, mvcube.mv_var_size, mvcube.num_mv_var);
#endif
		//free_sf_return_false(sf_vars);
		return false;
	}
	//#define SF_SORT_assend(P) {P = sf_unlist(sf_sort(P, (qsort_compare_func)ascend), P->count, P->sf_size);}
		/* 小到大排序，一樣的會擺一起*/
		//sf_sort(sf_vars, (qsort_compare_func)ascend);
		//SF_SORT_assend(sf_vars);
	sf_vars = parallel_sf_unlist(sf_sort, parallel_sf_sort(sf_vars, (qsort_compare_func)ascend));
	pset last, tp, pre_tp;
	int index_p;
	pre_tp = p_var_t; /* 為了不跳warning 隨便初始化，不會實際用到p_var_t的空間*/

	sf_active(sf_vars);
	int pre_cnt = 0;
	/* 紀錄那些var對應到相同(mvcube)的個數*/
	int cnt_same = 1;
	/* 紀錄那些var已經找到配對*/
	int num_unmatch_var = mvcube.num_mv_var;

	foreachi_set(sf_vars, index_p, tp) {
		/* 若已知解已包含tp全部則tp為無解*/
		if (setp_implies(tp, solution)) {
#ifdef mvcube_feas_check_debug
			printf("已知解已包含tp全部則tp為無解 infeasible\n");
			display_pset_eachvar(p, mvcube.mv_var_size, mvcube.num_mv_var);
#endif
			//free_sf_return_false(sf_vars);
			return false;
		}
		if (SIZE(tp) == 1) {
			/* 1. 一個var只對應到一個解時，該解視為已對應到*/
			set_or(solution, solution, tp);
			RESET(tp, ACTIVE);
			num_unmatch_var--;
		}
		else if (SIZE(tp) < mvcube.num_mv_var) {
			/* 將solution 中的解inverse*/
			set_inverse_disjoint(p_var_t, solution);
			//printf("solution:      "); display_pset(solution, mvcube.mv_var_size);
			//printf("original tp:   "); display_pset(tp, mvcube.mv_var_size);
			//printf("solution inve: "); display_pset(p_var_t, mvcube.mv_var_size);

			/* tp去除當前已知解*/
			set_and(tp, p_var_t, tp);
			//printf("removed tp:    "); display_pset(tp, mvcube.mv_var_size);
			if (pre_cnt == SIZE(tp)) {
				if (setp_equal(tp, pre_tp)) {
					cnt_same++;
					/* n個var都有相同的n個解時，解只可能對應到該n個解，故可視為已對應到*/
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
		/* 剩餘皆是full cube 跳過不比對*/
		else {
			break;
		}
		/*2. 尚未配對的var比尚未配對的解還多，則表示無解*/
		if (num_unmatch_var > mvcube.num_mv_var - set_count_ones(solution)) {
#ifdef mvcube_feas_check_debug
			printf("尚未配對的var比尚未配對的解還多 infeasible\n");
			display_pset_eachvar(p, mvcube.mv_var_size, mvcube.num_mv_var);
#endif
			//free_sf_return_false(sf_vars);
			return false;
		}
	}

	/* 還有remove的部分可以做*/

	//sf_free(sf_vars);
	return true;
}