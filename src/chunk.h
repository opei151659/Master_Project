/*
	將onset(offset)切成所需要的區塊大小，切割後的區塊大小可以用不同順序執行，這裡目前只有一種順序(從左至右，從上至下)
*/
#ifndef CHUNK_H
#define CHUNK_H

#include <vector>

/*
	coord (chunk)
	(f_begin, g_begin)___(f_begin, g_end)
	|                                   |
	|                                   |
	|                                   |
	(f_end, g_begin)_______(f_end, g_end)
*/
struct coord {
	int f_begin, f_end, g_begin, g_end;
	void set(int _f_begin, int _f_end, int _g_begin, int _g_end) {
		f_begin = _f_begin;
		f_end = _f_end;
		g_begin = _g_begin;
		g_end = _g_end;
	}
};


/*
	num_length
		|
		|
		|__________num_width
*/
class chunk {
public:
	int chunk_length_size, chunk_width_size, chunk_cnt, f_cnt, g_cnt;
	int num_length, num_width; /*chunks */
	std::vector<coord> chunks;
	chunk() :chunk_length_size(0), chunk_width_size(0), chunk_cnt(0), f_cnt(0), g_cnt(0) {};
	chunk(int _chunk_length_num, int _chunk_width_num, int _f_cnt, int _g_cnt) {
		setup(_chunk_length_num, _chunk_width_num, _f_cnt, _g_cnt);
	}
	void setup(int _chunk_length_num, int _chunk_width_num, int _f_cnt, int _g_cnt) {
		f_cnt = _f_cnt; /* PLA->F 的個數*/
		g_cnt = _g_cnt; /* PLA->R 的個數*/
		num_length = _chunk_length_num;//ceil((float)f_cnt / (float)chunk_length_size); /* 長的個數*/
		num_width = _chunk_width_num;//ceil((float)g_cnt / (float)chunk_width_size); /* 寬的個數*/
		chunk_length_size = floor((float)f_cnt / (float)num_length); /* 每一個chunk的長*/
		chunk_width_size = floor((float)g_cnt / (float)num_width); /* 每一個chunk的寬*/
		chunk_cnt = num_length * num_width;
		chunks.resize(chunk_cnt);
	}
	/* 由上至下由左至右排序
	   _________
	   | 1 2 3 |
	   | 4 5 6 |
	   | 7 8 9 |
	   ————
	
	*/
//#define display_order1
	void order1() {
		int id = 0;
		int f_begin, f_end, g_begin, g_end;
		for (int i = 0; i < num_length; i++) {
			f_begin = i * chunk_length_size;
			if (i == (num_length - 1)) {
				f_end = f_cnt;
			}
			else {
				f_end = f_begin + chunk_length_size;
			}
			
			for (int j = 0; j < num_width; j++) {
				g_begin = j * chunk_width_size;

				if (j == (num_width - 1)) {
					g_end = g_cnt;
				}
				else {
					g_end = g_begin + chunk_width_size;
				}
#ifdef display_order1
				printf("%5d ", id);
#endif
				chunks[id++].set(f_begin, f_end, g_begin, g_end);
			}
#ifdef display_order1
			printf("\n");
#endif
		}
	}
	void display() {
		printf("chunk size: ( %d , %d ) , chunk cnt: %d , num_lenght: %d , num_width: %d\n", chunk_length_size, chunk_width_size, chunk_cnt, num_length, num_width);
		printf("Top left\n");
		int id = 0;
		for (int i = 0; i < num_length; i++) {
			for (int j = 0; j < num_width; j++) {
				printf("(%d, %d) ", chunks[id].f_begin, chunks[id].g_begin);
				id++;
			}printf("\n");
		}
		printf("\n");
		printf("Top right\n");
		id = 0;
		for (int i = 0; i < num_length; i++) {
			for (int j = 0; j < num_width; j++) {
				printf("(%d, %d) ", chunks[id].f_begin, chunks[id].g_end);
				id++;
			}printf("\n");
		}
		printf("\n");
		printf("Bottom left\n");
		id = 0;
		for (int i = 0; i < num_length; i++) {
			for (int j = 0; j < num_width; j++) {
				printf("(%d, %d) ", chunks[id].f_end, chunks[id].g_begin);
				id++;
			}printf("\n");
		}
		printf("\n");
		printf("Bottom right\n");
		id = 0;
		for (int i = 0; i < num_length; i++) {
			for (int j = 0; j < num_width; j++) {
				printf("(%d, %d) ", chunks[id].f_end, chunks[id].g_end);
				id++;
			}printf("\n");
		}
	}
};

#endif