/*
    在平行化計算時常使用到的指令
    parallel_temp_psets: 暫存用的pset與mutex

    執行緒池的實作
*/

#ifndef PARALLEL_TOOL_H
#define PARALLEL_TOOL_H

#include <iostream>
#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include "mainCplusplus.h"

/* parallel */
class parallel_temp_psets {
public:
    std::vector<pset> psets;
    std::mutex& mtx_;

    parallel_temp_psets(int _size, std::mutex& mtx) : mtx_(mtx) {
        psets.resize(_size);
    };
    ~parallel_temp_psets() {
        {
            std::lock_guard<std::mutex> lock(mtx_);
            for (int i = 0; i < psets.size(); i++) {
                if (psets[i] != nullptr) {
                    //printf("parallel release %d\n", i);
                    set_free(psets[i]);
                }
            }
        }
        psets.clear();
        psets.shrink_to_fit();
    };
};

/*
    Barrier 方法較慢 有使用到的部分應該都已移除
*/
class Barrier {
public:
    Barrier() {};
    Barrier(int num_threads) : num_threads_(num_threads), counter_(num_threads), generation_(0) {}

    void setup(int num_threads) {
        num_threads_ = num_threads;
        counter_ = num_threads;
        generation_ = 0;
    }

    void Wait() {
        std::unique_lock<std::mutex> lock(mutex_);
        int gen = generation_;
        if (--counter_ == 0) {
            printf("%d %d\n", generation_, counter_);
            generation_++;
            counter_ = num_threads_;
            cv_.notify_all();
        }
        else {
            printf("%d %d\n", generation_, counter_);
            cv_.wait(lock, [this, gen]() { return generation_ != gen; });
        }
    }

    void finish() {
        std::unique_lock<std::mutex> lock(mutex_);
        --counter_;
        num_threads_--;

    }

private:
    int num_threads_;
    int counter_;
    int generation_;
    std::mutex mutex_;
    std::condition_variable cv_;
};

/*
    執行緒池的實作
    1. 根據需求創建執行緒池
    2. 將所有執行緒池用while保持運轉，如無事則進入idle狀態，不消耗CPU但會吃部分記憶體
    3. 回收所有執行緒
*/
class ThreadPool {
public:
    ThreadPool() : stop(false), threads_num(0) {}

    void setup(size_t numThreads) {
        threads_num = numThreads;
        check_setup.store(0);
        thread_ids.resize(threads_num);
        thread_ids_used_in_one_task.resize(threads_num);
        for (size_t i = 0; i < threads_num; ++i) {
            threads.emplace_back(
                [this, i] {
                    thread_ids[i] = std::this_thread::get_id();
                    check_setup.fetch_add(1);

                    //printf("%d\n", std::this_thread::get_id());
                    while (true) {
                        std::function<void()> task;
                        {
                            std::unique_lock<std::mutex> lock(queue_mutex);
                            condition.wait(lock, [this] {
                                return stop || !tasks.empty();
                                });
                            if (stop && tasks.empty()) {
                                return;
                            }
                            task = std::move(tasks.front());
                            thread_ids_used_in_one_task[i]++;
                            tasks.pop();
                        }
                        task();
                    }
                }
            );
        }
        //check_finish_setup
        while (check_setup < threads_num) {};
       
        //set_thread_ids
        sort(thread_ids.begin(), thread_ids.end(), std::less<std::thread::id>());
       
    }

    int get_thread_ids(std::thread::id id) {
        auto it = std::lower_bound(thread_ids.begin(), thread_ids.end(), id);

        if (it != thread_ids.end() && *it == id) {
            return std::distance(thread_ids.begin(), it);
        }
        else {
            return -1;
        }
    }

    void wait() {
        std::unique_lock<std::mutex> lock(this->counter_mutex);
        this->counter_condition.wait(lock, [this] { return this->counter == 0; });
    }

    template<class F, class... Args>
    void enqueue(F&& f, Args&&... args) {
        auto task = std::bind(std::forward<F>(f), std::forward<Args>(args)...);
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            if (stop) {
                throw std::runtime_error("enqueue on stopped ThreadPool");
            }
            tasks.emplace(task);
        }
        condition.notify_one();
    }



    void task_begin() {
        // increment the counter
        {
            std::unique_lock<std::mutex> lock(counter_mutex);
            ++counter;
        }
        
    }

    void task_end() {
        {
            std::unique_lock<std::mutex> lock(counter_mutex);
            --counter;
            if (this->counter == 0) {
                // all tasks completed, notify the main thread
                this->counter_condition.notify_one();
            }
        }
    }


    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            stop = true;
        }
        condition.notify_all();
        for (std::thread& thread : threads) {
            thread.join();
        }
    }

    void reset_Thread_task_used_cnt() {
        for (int& used : thread_ids_used_in_one_task) {
            used = 0;
        }
    }

    int display_Thread_task_used_cnt() {
        int cnt = 0;
        for (int& used : thread_ids_used_in_one_task) {
            if (used > 0)
                cnt++;
            //cout << used << " ";
        }//cout << endl;
        return cnt;
    }


private:
    int threads_num;
    std::vector<std::thread> threads;
    std::queue<std::function<void()>> tasks;
    std::mutex queue_mutex;
    std::condition_variable condition;
    std::vector<std::thread::id> thread_ids;
    std::atomic<int> check_setup;
    std::vector<int> thread_ids_used_in_one_task;


    std::mutex counter_mutex;
    std::condition_variable counter_condition;
    int counter = 0;
    bool stop;
};



/*
    parallel global vairable
*/
extern ThreadPool pool;



bool parallel_mvcube_feas_check_notation2(pset p, pset_family sf_vars, pset_family sf_sort, pset p_var_t, pset var_supercube, pset solution);
void local_intcpy(unsigned int* d, unsigned int* s, long int n);

#endif