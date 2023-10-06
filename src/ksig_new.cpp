#include "ksig_new.h"

/*計算pset中dc的數量*/
int dc_count(pset a) {
	unsigned int x;
	int tmp, count = 0;

	if ((tmp = cube.inword) != -1) {
		x = (a[tmp] >> 1) & a[tmp];
		x = x & cube.inmask;
		count += count_ones(x);
	}
	for (int i = 1; i < tmp; i++) {
		x = (a[i] >> 1) & a[i];
		x = x & DISJOINT;
		count += count_ones(x);
	}
	return count;
}


/*計算pset中dc的數量，並回傳dc的mask*/
int dc_count(pset r, pset a) {
	int tmp, count = 0;

	if ((tmp = cube.inword) != -1) {
		r[tmp] = (a[tmp] >> 1) & a[tmp] & cube.inmask;
		count += count_ones(r[tmp]);
	}
	for (int i = 1; i < tmp; i++) {
		r[i] = (a[i] >> 1) & a[i] & DISJOINT;
		count += count_ones(r[i]);
	}
	return count;
}

int cdistN(pset a, pset b) {
	unsigned int x;
	int tmp;
	int dist = 0;

	tmp = cube.inword;

	for (int j = 1; j < tmp; j++) {
		x = a[j] ^ b[j];
		x = (x >> 1) & x & DISJOINT;
		dist += count_ones(x);
	}

	if (tmp != -1) {
		x = a[tmp] ^ b[tmp];
		x = (x >> 1) & x & cube.inmask;
		dist += count_ones(x);
	}

	return dist;
}

/*計算兩個pset中變數有dc的數量*/
int set_ordc(pset a, pset b) {
	unsigned int x;
	int tmp, count = 0;

	if ((tmp = cube.inword) != -1) {
		x = ((a[tmp] >> 1) & a[tmp]) | ((b[tmp] >> 1) & b[tmp]);
		x = x & cube.inmask;
		count += count_ones(x);
	}
	for (int i = 1; i < tmp; i++) {
		x = ((a[i] >> 1) & a[i]) | ((b[i] >> 1) & b[i]);
		x = x & DISJOINT;
		count += count_ones(x);
	}
	return count;
}

/*計算兩個pset中變數都為dc的數量*/
int set_anddc(pset a, pset b) {
	unsigned int x;
	int tmp, count = 0;

	if ((tmp = cube.inword) != -1) {
		x = a[tmp] & b[tmp];
		x = (x >> 1) & x;
		x = x & cube.inmask;
		count += count_ones(x);
	}
	for (int i = 1; i < tmp; i++) {
		x = a[i] & b[i];
		x = (x >> 1) & x;
		x = x & DISJOINT;
		count += count_ones(x);
	}
	return count;
}

void parallel_pascal_table(TABLE& pt, int N) {
	pt[0].distx[0] = 1;

	for (int i = 1; i < N; i++) {
		pt[i].distx[0] = 1;
		tbb::parallel_for(1, i + 1, [&](int t) {
			//for (int t = 1; t < i + 1; t++)
			pt[i].distx[t] = pt[i - 1].distx[t] + pt[i - 1].distx[t - 1];
			});
	}
}

void parallel_ksig_table(TABLE& NCkt, TABLE pt, int N) {
	NCkt[0].distx[0] = 0;

	for (int i = 1; i < N; i++) {
		tbb::parallel_for(0, i, [&](int t) {
			//for (int t = 0; t < i; t++)
			NCkt[i].distx[t] = NCkt[i - 1].distx[t] * 2 + pt[i - 1].distx[t] * mp::pow(cpp_int(2), i - 1);
			});
	}
}

void chunk_kso(pset_family pf, ksig_value& kso, TABLE pt, TABLE kt, int tsize) {
	int chunksize = 10; /*資料塊大小*/
	int chunknum = ceil((double)pf->count / chunksize); /*計算on-set的數量為chunksize的幾倍，也為生成執行緒的數目*/
	int cs = ceil((double)pf->count / chunknum); /*計算平均一塊資料塊的積項數量*/

	int TN = chunknum > tsize ? tsize : chunknum;

	TABLE thread_kso(TN, ksig_value(NUMINPUTS));
	atomic<int> Index1 = 0, Index2 = 0;
	int job2_size = (chunknum - 1) * chunknum / 2;
	vector<pair<int, int>> Job2(job2_size, make_pair(-1, -1));

	for (int i = 0, index = 0; i < chunknum; i++) {
		for (int j = i + 1; j < chunknum; j++) {
			Job2[index] = make_pair(i, j);
			index++;
		}
	}

	vector<thread> thread_vec;
	for (int id = 0; id < TN; id++) {
		thread_vec.emplace_back(std::thread([&](int tid) {
			pset one, two;
			int dist, cup_num, cap_num, dc_num;
			int i, j, e1, e2, se1, se2;

			/*KSI*/
			cpp_int powdc2;

			int pid;

			while ((pid = Index1.fetch_add(1)) < chunknum) {
				for (e1 = pid * cs, e2 = 0, one = pf->data + pf->wsize * e1; e1 < pf->count && e2 < cs; e1++, e2++, one += pf->wsize) {
					dc_num = dc_count(one);

					/*KSO*/
					thread_kso[tid] += kt[dc_num];


					for (se1 = e1 + 1, se2 = e2 + 1, two = pf->data + pf->wsize * se1; se1 < pf->count && se2 < cs; se1++, se2++, two += pf->wsize) {
						dist = cdistN(one, two);
						cup_num = set_ordc(one, two);
						cap_num = set_anddc(one, two);
						powdc2 = mp::pow(mp::cpp_int(2), cap_num);

						if (dist > 0) {
							/*KSO*/
							for (i = dist - 1, j = 0; j <= cup_num && i >= 0; i++, j++) {
								thread_kso[tid].distx[i] += pt[cup_num].distx[j] * powdc2;
							}
						}
					}
				}
			}
			while ((pid = Index2.fetch_add(1)) < job2_size) {
				for (e1 = Job2[pid].first * cs, e2 = 0, one = pf->data + pf->wsize * e1; e1 < pf->count && e2 < cs; e1++, e2++, one += pf->wsize) {
					for (se1 = Job2[pid].second * cs, se2 = 0, two = pf->data + pf->wsize * se1; se1 < pf->count && se2 < cs; se1++, se2++, two += pf->wsize) {
						dist = cdistN(one, two);
						cup_num = set_ordc(one, two);
						cap_num = set_anddc(one, two);
						powdc2 = mp::pow(mp::cpp_int(2), cap_num);

						/*kso*/
						if (dist > 0) {
							/*KSO*/
							for (i = dist - 1, j = 0; j <= cup_num && i >= 0; i++, j++) {
								thread_kso[tid].distx[i] += pt[cup_num].distx[j] * powdc2;
							}
						}
					}
				}
			}
			}, id));
	}
	for (auto& thread : thread_vec) {
		thread.join();
	}
	for (int i = 0; i < TN; i++) {
		kso += thread_kso[i];
	}
}

void chunk_ksi(pset_family pf, TABLE& ksi, TABLE pt, TABLE kt, int N, pset sup, int tsize) {

	for (int j = 0; j < NUMINPUTS; j++) {
		if (getvar(sup, j) == '?') {
			for (int i = 0; i < ksi[j].distx.size(); i++) {
				ksi[j].distx[i] = -1;
			}
		}
	}

	int chunksize = 10; /*資料塊大小*/
	int chunknum = ceil((double)pf->count / chunksize); /*計算on-set的數量為chunksize的幾倍，也為生成執行緒的數目*/
	int cs = ceil((double)pf->count / chunknum); /*計算平均一塊資料塊的積項數量*/

	int TN = chunknum > tsize ? tsize : chunknum;

	atomic<int> Index1 = 0, Index2 = 0;
	int job2_size = (chunknum - 1) * chunknum / 2;
	vector<pair<int, int>> Job2(job2_size, make_pair(-1, -1));

	for (int i = 0, index = 0; i < chunknum; i++) {
		for (int j = i + 1; j < chunknum; j++) {
			Job2[index] = make_pair(i, j);
			index++;
		}
	}

	vector<TABLE> thread_ksi(TN, TABLE(NUMINPUTS, ksig_value(N)));
	vector<thread> thread_vec;
	for (int id = 0; id < TN; id++) {
		thread_vec.emplace_back(std::thread([&](int tid) {
			pset one, two, r = set_new(cube.size);
			int dist, cup_num, cap_num, dc_num;
			int i, j, t, k, e1, e2, se1, se2;
			char c1, c2;

			/*KSI*/
			vector<cpp_int> dc_ksig(N), dist_ksig(N);
			cpp_int powdc2;

			int pid;

			while ((pid = Index1.fetch_add(1)) < chunknum) {
				for (e1 = pid * cs, e2 = 0, one = pf->data + pf->wsize * e1; e1 < pf->count && e2 < cs; e1++, e2++, one += pf->wsize) {
					dc_num = dc_count(r, one);
					if (dc_num > 0) {
						/*KSO*/
						powdc2 = mp::pow(mp::cpp_int(2), dc_num - 1);
						for (i = 0; i < dc_num && i < N; i++) {
							dc_ksig[i] = pt[dc_num - 1].distx[i] * powdc2;
						}
						for (i = 1, t = 0; i <= cube.inword && t != dc_num; i++) {
							if (r[i]) {
								for (j = 0; j < varnum; j++) {
									if (is_in_set(r, (i - 1) * BPI + j * 2)) {
										for (k = 0; k < dc_num && k < N; k++) {
											thread_ksi[tid][(i - 1) * varnum + j].distx[k] += dc_ksig[k];
										}
										t++;
									}
								}
							}
						}
					}
					for (se1 = e1 + 1, se2 = e2 + 1, two = pf->data + pf->wsize * se1; se1 < pf->count && se2 < cs; se1++, se2++, two += pf->wsize) {
						dist = cdistN(one, two);
						cup_num = set_ordc(one, two);
						cap_num = set_anddc(one, two);
						powdc2 = mp::pow(mp::cpp_int(2), cap_num);

						if (dist > 0) {
							/*KSI*/
							if (dist <= N) {
								for (i = 0; i < N; i++) {
									if ((i + 1) > dist + cup_num || (i + 1) < dist) {
										dist_ksig[i] = 0;
										dc_ksig[i] = 0;
									}
									else if ((i + 1) - dist == 0) {
										dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * powdc2;
										dc_ksig[i] = 0;
									}
									else {
										dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * powdc2;
										dc_ksig[i] = pt[cup_num - 1].distx[(i + 1) - dist - 1] * powdc2;
									}
								}
								for (i = 0; i < NUMINPUTS; i++) {
									c1 = getvar(one, i);
									c2 = getvar(two, i);

									if (c1 == '2' || c2 == '2') {
										for (j = dist; j <= N; j++) {
											thread_ksi[tid][i].distx[j - 1] += dc_ksig[j - 1];
										}
									}
									else if (c1 != c2) {
										for (j = dist; j <= N; j++) {
											thread_ksi[tid][i].distx[j - 1] += dist_ksig[j - 1];
										}
									}
								}
							}
						}
					}
				}
			}
			FREE(r);
			while ((pid = Index2.fetch_add(1)) < job2_size) {
				for (e1 = Job2[pid].first * cs, e2 = 0, one = pf->data + pf->wsize * e1; e1 < pf->count && e2 < cs; e1++, e2++, one += pf->wsize) {
					for (se1 = Job2[pid].second * cs, se2 = 0, two = pf->data + pf->wsize * se1; se1 < pf->count && se2 < cs; se1++, se2++, two += pf->wsize) {
						dist = cdistN(one, two);
						cup_num = set_ordc(one, two);
						cap_num = set_anddc(one, two);
						powdc2 = mp::pow(mp::cpp_int(2), cap_num);

						if (dist > 0) {
							/*KSI*/
							if (dist <= N) {
								for (i = 0; i < N; i++) {
									if ((i + 1) > dist + cup_num || (i + 1) < dist) {
										dist_ksig[i] = 0;
										dc_ksig[i] = 0;
									}
									else if ((i + 1) - dist == 0) {
										dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * powdc2;
										dc_ksig[i] = 0;
									}
									else {
										dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * powdc2;
										dc_ksig[i] = pt[cup_num - 1].distx[(i + 1) - dist - 1] * powdc2;
									}
								}
								for (i = 0; i < NUMINPUTS; i++) {
									c1 = getvar(one, i);
									c2 = getvar(two, i);

									if (c1 == '2' || c2 == '2') {
										for (j = dist; j <= N; j++) {
											thread_ksi[tid][i].distx[j - 1] += dc_ksig[j - 1];
										}
									}
									else if (c1 != c2) {
										for (j = dist; j <= N; j++) {
											thread_ksi[tid][i].distx[j - 1] += dist_ksig[j - 1];
										}
									}
								}
							}
						}
					}
				}
			}
			}, id));
	}
	for (auto& thread : thread_vec) {
		thread.join();
	}
	for (int i = 0; i < TN; i++) {
		for (int j = 0; j < NUMINPUTS; j++) {
			ksi[j] += thread_ksi[i][j];
		}
	}
}
void chunk_cal_ksignature(pset_family pf, ksig_value& kso, TABLE& ksi, TABLE pt, TABLE kt, int N, pset sup, int tsize) {

	for (int j = 0; j < NUMINPUTS; j++) {
		if (getvar(sup, j) == '?') {
			for (int i = 0; i < ksi[j].distx.size(); i++) {
				ksi[j].distx[i] = -1;
			}
		}
	}

	int chunksize = 10; /*資料塊大小*/
	int chunknum = ceil((double)pf->count / chunksize); /*計算on-set的數量為chunksize的幾倍，也為生成執行緒的數目*/
	int cs = ceil((double)pf->count / chunknum); /*計算平均一塊資料塊的積項數量*/

	int TN = chunknum > tsize ? tsize : chunknum;
	//TN = tsize;
	//cout << TN << endl;
	TABLE thread_kso(TN, ksig_value(NUMINPUTS));

	atomic<int> Index1 = 0, Index2 = 0;
	int job2_size = (chunknum - 1) * chunknum / 2;
	vector<pair<int, int>> Job2(job2_size, make_pair(-1, -1));

	for (int i = 0, index = 0; i < chunknum; i++) {
		for (int j = i + 1; j < chunknum; j++) {
			Job2[index] = make_pair(i, j);
			index++;
		}
	}

	vector<TABLE> thread_ksi(TN, TABLE(NUMINPUTS, ksig_value(N)));
	vector<thread> thread_vec;
	for (int id = 0; id < TN; id++) {
		thread_vec.emplace_back(std::thread([&](int tid) {
			pset one, two, r = set_new(cube.size);
			int dist, cup_num, cap_num, dc_num;
			int i, j, t, k, e1, e2, se1, se2;
			char c1, c2;

			/*KSI*/
			vector<cpp_int> dc_ksig(N), dist_ksig(N);
			cpp_int powdc2;

			int pid;

			while ((pid = Index1.fetch_add(1)) < chunknum) {
				for (e1 = pid * cs, e2 = 0, one = pf->data + pf->wsize * e1; e1 < pf->count && e2 < cs; e1++, e2++, one += pf->wsize) {
					dc_num = dc_count(r, one);
					if (dc_num > 0) {
						/*KSO*/
						thread_kso[tid] += kt[dc_num];
						powdc2 = mp::pow(mp::cpp_int(2), dc_num - 1);
						for (i = 0; i < dc_num && i < N; i++) {
							dc_ksig[i] = pt[dc_num - 1].distx[i] * powdc2;
						}
						for (i = 1, t = 0; i <= cube.inword && t != dc_num; i++) {
							if (r[i]) {
								for (j = 0; j < varnum; j++) {
									if (is_in_set(r, (i - 1) * BPI + j * 2)) {
										for (k = 0; k < dc_num && k < N; k++) {
											thread_ksi[tid][(i - 1) * varnum + j].distx[k] += dc_ksig[k];
										}
										t++;
									}
								}
							}
						}
					}
					for (se1 = e1 + 1, se2 = e2 + 1, two = pf->data + pf->wsize * se1; se1 < pf->count && se2 < cs; se1++, se2++, two += pf->wsize) {
						dist = cdistN(one, two);
						cup_num = set_ordc(one, two);
						cap_num = set_anddc(one, two);
						powdc2 = mp::pow(mp::cpp_int(2), cap_num);

						if (dist > 0) {
							/*KSO*/
							for (i = dist - 1, j = 0; j <= cup_num && i >= 0; i++, j++) {
								thread_kso[tid].distx[i] += pt[cup_num].distx[j] * powdc2;
							}
							/*KSI*/
							if (dist <= N) {
								for (i = 0; i < N; i++) {
									if ((i + 1) > dist + cup_num || (i + 1) < dist) {
										dist_ksig[i] = 0;
										dc_ksig[i] = 0;
									}
									else if ((i + 1) - dist == 0) {
										dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * powdc2;
										dc_ksig[i] = 0;
									}
									else {
										dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * powdc2;
										dc_ksig[i] = pt[cup_num - 1].distx[(i + 1) - dist - 1] * powdc2;
									}
								}
								for (i = 0; i < NUMINPUTS; i++) {
									c1 = getvar(one, i);
									c2 = getvar(two, i);

									if (c1 == '2' || c2 == '2') {
										for (j = dist; j <= N; j++) {
											thread_ksi[tid][i].distx[j - 1] += dc_ksig[j - 1];
										}
									}
									else if (c1 != c2) {
										for (j = dist; j <= N; j++) {
											thread_ksi[tid][i].distx[j - 1] += dist_ksig[j - 1];
										}
									}
								}
							}
						}
					}
				}
			}
			FREE(r);
			while ((pid = Index2.fetch_add(1)) < job2_size) {
				for (e1 = Job2[pid].first * cs, e2 = 0, one = pf->data + pf->wsize * e1; e1 < pf->count && e2 < cs; e1++, e2++, one += pf->wsize) {
					for (se1 = Job2[pid].second * cs, se2 = 0, two = pf->data + pf->wsize * se1; se1 < pf->count && se2 < cs; se1++, se2++, two += pf->wsize) {
						dist = cdistN(one, two);
						cup_num = set_ordc(one, two);
						cap_num = set_anddc(one, two);
						powdc2 = mp::pow(mp::cpp_int(2), cap_num);

						/*kso*/
						if (dist > 0) {
							/*KSO*/
							for (i = dist - 1, j = 0; j <= cup_num && i >= 0; i++, j++) {
								thread_kso[tid].distx[i] += pt[cup_num].distx[j] * powdc2;
							}
							/*KSI*/
							if (dist <= N) {
								for (i = 0; i < N; i++) {
									if ((i + 1) > dist + cup_num || (i + 1) < dist) {
										dist_ksig[i] = 0;
										dc_ksig[i] = 0;
									}
									else if ((i + 1) - dist == 0) {
										dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * powdc2;
										dc_ksig[i] = 0;
									}
									else {
										dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * powdc2;
										dc_ksig[i] = pt[cup_num - 1].distx[(i + 1) - dist - 1] * powdc2;
									}
								}
								for (i = 0; i < NUMINPUTS; i++) {
									c1 = getvar(one, i);
									c2 = getvar(two, i);

									if (c1 == '2' || c2 == '2') {
										for (j = dist; j <= N; j++) {
											thread_ksi[tid][i].distx[j - 1] += dc_ksig[j - 1];
										}
									}
									else if (c1 != c2) {
										for (j = dist; j <= N; j++) {
											thread_ksi[tid][i].distx[j - 1] += dist_ksig[j - 1];
										}
									}
								}
							}
						}
					}
				}
			}
			}, id));
	}
	for (auto& thread : thread_vec) {
		thread.join();
	}
	for (int i = 0; i < TN; i++) {
		kso += thread_kso[i];
		for (int j = 0; j < NUMINPUTS; j++) {
			ksi[j] += thread_ksi[i][j];
		}
	}
}

/* Mode 0: KSO, 1: KSI, 2: KSO+KSI*/
void cal_ksignature(pPLA PLA, TABLE& KSO, vector<TABLE>& KSI, int N, int Mode) {
	TABLE pt(NUMINPUTS, ksig_value(NUMINPUTS)), pt1(NUMINPUTS, ksig_value(NUMINPUTS));
	TABLE kt(NUMINPUTS + 1, ksig_value(NUMINPUTS)), kt1(NUMINPUTS + 1, ksig_value(NUMINPUTS));

	parallel_pascal_table(pt, NUMINPUTS);
	parallel_ksig_table(kt, pt, NUMINPUTS + 1);

	int tsize = thread::hardware_concurrency();

	pset_family f;
	int i;
	pset p, last;
	pset_family supp;
	if (NUMOUTPUTS == 1) {
		supp = sf_new(1, cube.size);
		sf_addset(supp, set_full(cube.size));
	}
	else {
		supp = support_set(PLA->F, PLA->R);
	}

	for (i = 0, p = supp->data; i < NUMOUTPUTS; i++, p += supp->wsize) {
		/* seperate function */
		{
			//fprintf(stderr, "START make_disjoint...\n");
			f = sf_new(PLA->F->count, PLA->F->sf_size);
			f = sep_sup_output(PLA->F, i + cube.first_part[cube.output], p);
			f = sf_contain(f);
			f = make_disjoint(f);
			//fprintf(stderr, "FINISH make_disjoint...\n");
		}

		if (Mode == 0) {
			chunk_kso(f, KSO[i], pt, kt, tsize);
		}
		else if (Mode == 1) {
			chunk_ksi(f, KSI[i], pt, kt, N, p, tsize);
		}
		else {
			chunk_cal_ksignature(f, KSO[i], KSI[i], pt, kt, N, p, tsize);
		}
	}

	FREE(supp);
}