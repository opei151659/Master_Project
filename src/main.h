/*
    原Espresso 中的指令
    在簡化後大部分指令無效 詳細看原版Espresso
*/

enum keys {
    KEY_ESPRESSO, KEY_PLA_verify, KEY_check, KEY_contain, KEY_d1merge,
    KEY_disjoint, KEY_dsharp, KEY_echo, KEY_essen, KEY_exact, KEY_expand,
    KEY_gasp, KEY_intersect, KEY_irred, KEY_lexsort, KEY_make_sparse,
    KEY_map, KEY_mapdc, KEY_minterms, KEY_opo, KEY_opoall,
    KEY_pair, KEY_pairall, KEY_primes, KEY_qm, KEY_reduce, KEY_sharp,
    KEY_simplify, KEY_so, KEY_so_both, KEY_stats, KEY_super_gasp, KEY_taut,
    KEY_test, KEY_equiv, KEY_union, KEY_verify, KEY_MANY_ESPRESSO,
    KEY_separate, KEY_xor, KEY_d1merge_in, KEY_fsm, KEY_signature,
    KEY_unknown
};


/* Lookup table for program options */
struct {
    char *name;
    enum keys key;
    int num_plas;
    bool needs_offset;
    bool needs_dcset;
} option_table [] = {
    /* ways to minimize functions */
    {"ESPRESSO", KEY_ESPRESSO, 1, TRUE, TRUE},    /* must be first */
    {"many", KEY_MANY_ESPRESSO, 1, TRUE, TRUE},
    {"exact", KEY_exact, 1, TRUE, TRUE},
    {"qm", KEY_qm, 1, TRUE, TRUE},
    {"single_output", KEY_so, 1, TRUE, TRUE},
    {"so", KEY_so, 1, TRUE, TRUE},
    {"so_both", KEY_so_both, 1, TRUE, TRUE},
    {"simplify", KEY_simplify, 1, FALSE, FALSE},
    {"echo", KEY_echo, 1, FALSE, FALSE},
    {"signature", KEY_signature, 1, TRUE, TRUE},

    /* output phase assignment and assignment of inputs to two-bit decoders */
    {"opo", KEY_opo, 1, TRUE, TRUE},
    {"opoall", KEY_opoall, 1, TRUE, TRUE},
    {"pair", KEY_pair, 1, TRUE, TRUE},
    {"pairall", KEY_pairall, 1, TRUE, TRUE},

    /* Ways to check covers */
    {"check", KEY_check, 1, TRUE, TRUE},
    {"stats", KEY_stats, 1, FALSE, FALSE},
    {"verify", KEY_verify, 2, FALSE, TRUE},
    {"PLAverify", KEY_PLA_verify, 2, FALSE, TRUE},

    /* hacks */
    {"equiv", KEY_equiv, 1, TRUE, TRUE},
    {"map", KEY_map, 1, FALSE, FALSE},
    {"mapdc", KEY_mapdc, 1, FALSE, FALSE},
    {"fsm", KEY_fsm, 1, FALSE, TRUE},

    /* the basic boolean operations on covers */
    {"contain", KEY_contain, 1, FALSE, FALSE},
    {"d1merge", KEY_d1merge, 1, FALSE, FALSE},
    {"d1merge_in", KEY_d1merge_in, 1, FALSE, FALSE},
    {"disjoint", KEY_disjoint, 1, TRUE, FALSE},
    {"dsharp", KEY_dsharp, 2, FALSE, FALSE},
    {"intersect", KEY_intersect, 2, FALSE, FALSE},
    {"minterms", KEY_minterms, 1, FALSE, FALSE},
    {"primes", KEY_primes, 1, FALSE, TRUE},
    {"separate", KEY_separate, 1, TRUE, TRUE},
    {"sharp", KEY_sharp, 2, FALSE, FALSE},
    {"union", KEY_union, 2, FALSE, FALSE},
    {"xor", KEY_xor, 2, TRUE, TRUE},

    /* debugging only -- call each step of the espresso algorithm */
    {"essen", KEY_essen, 1, FALSE, TRUE},
    {"expand", KEY_expand, 1, TRUE, FALSE},
    {"gasp", KEY_gasp, 1, TRUE, TRUE},
    {"irred", KEY_irred, 1, FALSE, TRUE},
    {"make_sparse", KEY_make_sparse, 1, TRUE, TRUE},
    {"reduce", KEY_reduce, 1, FALSE, TRUE},
    {"taut", KEY_taut, 1, FALSE, FALSE},
    {"super_gasp", KEY_super_gasp, 1, TRUE, TRUE},
    {"lexsort", KEY_lexsort, 1, FALSE, FALSE},
    {"test", KEY_test, 1, TRUE, TRUE},
    {0, KEY_unknown, 0, FALSE, FALSE}             /* must be last */
};


struct {
    char *name;
    int value;
} debug_table[] = {
    {"", EXPAND + ESSEN + IRRED + REDUCE + SPARSE + GASP + SHARP + MINCOV},
    {"compl",   COMPL},  {"essen",       ESSEN},
    {"expand",  EXPAND}, {"expand1",     EXPAND1|EXPAND},
    {"irred",   IRRED},  {"irred1",      IRRED1|IRRED},
    {"reduce",  REDUCE}, {"reduce1",     REDUCE1|REDUCE},
    {"mincov",  MINCOV}, {"mincov1",     MINCOV1|MINCOV},
    {"sparse",  SPARSE}, {"sharp",       SHARP},
    {"taut",    TAUT},   {"gasp",        GASP},
    {"exact",   EXACT},
    {0, 0}
};


struct {
    char *name;
    int *variable;
    int value;
} esp_opt_table[] = {
    {"eat", &echo_comments, FALSE},
    {"eatdots", &echo_unknown_commands, FALSE},
    {"fast", &single_expand, TRUE},
    {"kiss", &kiss, TRUE},
    {"ness", &remove_essential, FALSE},
    {"nirr", &force_irredundant, FALSE},
    {"nunwrap", &unwrap_onset, FALSE},
    {"onset", &recompute_onset, TRUE},
    {"pos", &pos, TRUE},
    {"random", &use_random_order, TRUE},
    {"strong", &use_super_gasp, TRUE},
    {0, 0, 0}
};

// 參數範圍限制
class arg_range {
    unordered_map<string, pair<int, int>> args;
public:
    arg_range() {
        args["-IN"] = make_pair(1, 128);
        args["-RM"] = make_pair(0, 1);
        args["-ROW"] = make_pair(0, INT_MAX);
        args["-COL"] = make_pair(0, INT_MAX);
        args["-T"] = make_pair(1, std::thread::hardware_concurrency()); //小於硬體CPU最大執行緒個數
        args["-CHK"] = make_pair(1, INT_MAX);
        args["-SIG"] = make_pair(0, 15);
        args["-MCQE"] = make_pair(0, 15);
    }
    int check_range_i(string arg, string value) {
        int intValue = 0;
        if (args.find(arg) != args.end()) {
            try {
                intValue = stoi(value); // 將參數轉換為整數
                if (intValue >= args[arg].first && intValue <= args[arg].second) {  // 自訂範圍檢查
                    //cout << arg << " " << intValue << std::endl;
                }
                else {
                    cerr << arg << " out of range (" << args[arg].first << " ~ " << args[arg].second << ")" << endl;
                    exit(EXIT_FAILURE);
                }
            }
            catch (const std::invalid_argument& e) {
                cerr << arg << " wrong format" << e.what() << endl;
                exit(EXIT_FAILURE);
            }
            catch (const std::out_of_range& e) {
                cerr << arg << " out of range (" << args[arg].first << " ~ " << args[arg].second << ")" << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout << "wrong arguments" << endl;
            exit(EXIT_FAILURE);
        }

        return  intValue;
    }


    float check_range_f(string arg, string value) {
        float floatValue = 0.0;
        if (args.find(arg) != args.end()) {
            try {
                floatValue = stof(value); // 將參數轉換為整數
                if (floatValue >= (float)args[arg].first && floatValue <= (float)args[arg].second) {  // 自訂範圍檢查
                    //cout << arg << " " << floatValue << std::endl;
                }
                else {
                    cerr << arg << " out of range (" << (float)args[arg].first << " ~ " << (float)args[arg].second << ")" << endl;
                    exit(EXIT_FAILURE);
                }
            }
            catch (const std::invalid_argument& e) {
                cerr << arg << " wrong format" << e.what() << endl;
                exit(EXIT_FAILURE);
            }
            catch (const std::out_of_range& e) {
                cerr << arg << " out of range (" << args[arg].first << " ~ " << args[arg].second << ")" << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            cout << "wrong arguments" << endl;
            exit(EXIT_FAILURE);
        }

        return floatValue;
    }
};
