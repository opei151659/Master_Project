/*
   原Espresso主程式的C++簡化版
   可以新增功能指令
*/

#include <iomanip>
#include <iostream>
#include "mainCplusplus.h"
extern "C" {
    #include "main.h"       /* table definitions for options */
}
int     opterr = 1,      /* if error message should be printed */
        optind = 1,      /* index into parent argv vector */
        optopt,          /* character checked for validity */
        optreset;        /* reset getopt */
char* optarg;            /* argument associated with option */

#define BADCH   (int)'?'
#define BADARG  (int)':'
#define EMSG    ""

/* ==========user define==========*/
int isINP, isINN, isOUTP, isOUTN; /* is INput/OUTput Permute/iNverse*/
int row_num, col_num;
int ppbm_thread_num, ppbm_chunk_size, ppbm_option;
float remove_sf_rate;
int extand_input_num;
int CAL_SIG_OUTPUT, CAL_SIG_INPUT, CAL_ALL_MCQE;

int SIG_key;
int MCQE_key;
int Total_parallel;
/* ==========user define==========*/
static FILE* last_fp;
/* !!!!!輸入的type類型*/
static int input_type = FD_type; //FR_type; FD_type

void getPLA(int opt, int argc, char** argv, pPLA* PLA, int out_type);
void delete_arg(int* argc,  char** argv, int num);
void init_runtime(void);
void backward_compatibility_hack(int* argc, char** argv, int* option, int* out_type);
void runtime(void);
void usage(void);
bool check_arg(int* argc,  char** argv,  char* s);


using namespace std;

int getopt(int nargc, char* const nargv[], const char* ostr)
{
    static char* place = EMSG;              /* option letter processing */
    const char* oli;                        /* option letter list index */

    if (optreset || !*place) {              /* update scanning pointer */
        optreset = 0;
        if (optind >= nargc || *(place = nargv[optind]) != '-') {
            place = EMSG;
            return (-1);
        }
        if (place[1] && *++place == '-') {      /* found "--" */
            ++optind;
            place = EMSG;
            return (-1);
        }
    }                                       /* option letter okay? */
    if ((optopt = (int)*place++) == (int)':' ||
        !(oli = strchr(ostr, optopt))) {
        /*
        * if the user didn't specify '-' as an option,
        * assume it means -1.
        */
        if (optopt == (int)'-')
            return (-1);
        if (!*place)
            ++optind;
        if (opterr && *ostr != ':')
            (void)printf("illegal option -- %c\n", optopt);
        return (BADCH);
    }
    if (*++oli != ':') {                    /* don't need argument */
        optarg = NULL;
        if (!*place)
            ++optind;
    }
    else {                                  /* need an argument */
        if (*place)                     /* no white space */
            optarg = place;
        else if (nargc <= ++optind) {   /* no arg */
            place = EMSG;
            if (*ostr == ':')
                return (BADARG);
            if (opterr)
                (void)printf("option requires an argument -- %c\n", optopt);
            return (BADCH);
        }
        else                            /* white space */
            optarg = nargv[optind];
        place = EMSG;
        ++optind;
    }
    return (optopt);                        /* dump back option letter */
}

int main(int argc, char** argv) {
    int i, j, first, last, strategy, out_type, option;
    pPLA PLA, PLA1;
    pcover F, Fold, Dold;
    pset last1, p;
    cost_t cost;
    bool error, exact_cover;
    long start;

    start = ptime();

    error = FALSE;
    init_runtime();
#ifdef RANDOM
    srandom(314973);
#endif

    option = 0;         /* default -D: ESPRESSO */
    out_type = F_type;      /* default -o: default is ON-set only */
    debug = 0;          /* default -d: no debugging info */
    verbose_debug = FALSE;  /* default -v: not verbose */
    print_solution = TRUE;  /* default -x: print the solution (!) */
    summary = FALSE;        /* default -s: no summary */
    trace = FALSE;      /* default -t: no trace information */
    strategy = 0;       /* default -S: strategy number */
    first = -1;         /* default -R: select range */
    last = -1;
    remove_essential = TRUE;    /* default -e: */
    force_irredundant = TRUE;
    unwrap_onset = TRUE;
    single_expand = FALSE;
    pos = FALSE;
    recompute_onset = FALSE;
    use_super_gasp = FALSE;
    use_random_order = FALSE;
    kiss = FALSE;
    echo_comments = FALSE; /* printf comments in PLA file (usually is TRUE)*/
    echo_unknown_commands = TRUE;
    exact_cover = FALSE;    /* for -qm option, the default */

    backward_compatibility_hack(&argc, argv, &option, &out_type);

    /* provide version information and summaries */
    if (summary || trace) {
        /* echo command line and arguments */
        printf("#");
        for (i = 0; i < argc; i++) {
            printf(" %s", argv[i]);
        }
        printf("\n");
        printf("# %s\n", VERSION);
    }

    /* READ PLA */
    /* the remaining arguments are argv[optind ... argc-1] */
    optind = 1;
    string str;
    if (argc == 4) {
        str = argv[3];
    }
    PLA = PLA1 = NIL(PLA_t);
    if (argc == 2) {
        getPLA(1, argc, argv, &PLA, out_type);
    }
    /* 輸出特徵值實驗*/
    /* 範例 pla\cm150a.pla pla\cm150a.pla OutputSIG*/
    else if (argc == 4 && str == "OutputSIG") {
        getPLA(optind++, argc, argv, &PLA, out_type);   // PLA 1 (.pla)
        getPLA(optind++, argc, argv, &PLA1, out_type);  // PLA 2 (.pla)
        extand_input_num = 0;
        isINP = 1;
        isINN = 1;
        isOUTP = 0;
        isOUTN = 0;
        remove_sf_rate = 0;
        ppbm_thread_num = 1;
        CAL_SIG_OUTPUT = 1; /*計算輸出特徵值*/
        CAL_SIG_INPUT = 0;
        CAL_ALL_MCQE = 0;
    }
    /* 輸入特徵值實驗*/
    /* 範例 pla\cm150a.pla pla\cm150a.pla InputSIG*/
    else if (argc == 4 && str == "InputSIG") {
        getPLA(optind++, argc, argv, &PLA, out_type);   // PLA 1 (.pla)
        getPLA(optind++, argc, argv, &PLA1, out_type);  // PLA 2 (.pla)
        extand_input_num = 0;
        isINP = 1;
        isINN = 1;
        isOUTP = 0;
        isOUTN = 0;
        remove_sf_rate = 0;
        ppbm_thread_num = 1;
        CAL_SIG_OUTPUT = 0; 
        CAL_SIG_INPUT = 1; /*計算輸入特徵值*/
        CAL_ALL_MCQE = 0;
    }
    else if (argc == 4 && str == "AllMCQE") {
        getPLA(optind++, argc, argv, &PLA, out_type);   // PLA 1 (.pla)
        getPLA(optind++, argc, argv, &PLA1, out_type);  // PLA 2 (.pla)
        extand_input_num = 0;
        isINP = 1;
        isINN = 1;
        isOUTP = 0;
        isOUTN = 0;
        remove_sf_rate = 0;
        ppbm_thread_num = 1;
        CAL_SIG_OUTPUT = 0;
        CAL_SIG_INPUT = 0; /*計算輸入特徵值*/
        CAL_ALL_MCQE = 1;
    }
    /*範例 pla\cm150a.pla pla\cm150a.pla 0 1 1 0 0 0.0 1 0 1 100 2 7 0 0*/
    else if (argc == 17) { // 17 個輸入參數
        getPLA(optind++, argc, argv, &PLA, out_type);   // PLA 1 (.pla)
        getPLA(optind++, argc, argv, &PLA1, out_type);  // PLA 2 (.pla)
        extand_input_num = atoi(argv[optind++]);        // 是否要擴展輸入的input個數，數字小於原本大小包持原樣，最大為128  (0 ~ 128)
        isINP = atoi(argv[optind++]);                   // 是否有input permutation (1/0)
        isINN = atoi(argv[optind++]);                   // 是否有input phase assignment (1/0)
        isOUTP = atoi(argv[optind++]);                  // 是否有output permutation (1/0) 
        isOUTN = atoi(argv[optind++]);                  // 是否有output phase assignment (1/0)
        remove_sf_rate = atof(argv[optind++]);          // 移除on-set 與off-set 的百分比  (0 ~ 1) 0.1表示移除10%
        row_num = atoi(argv[optind++]);                 // row(on-set) 要分割的份數
        col_num = atoi(argv[optind++]);                 // column(off-set) 要分割的份數
        ppbm_thread_num = atoi(argv[optind++]);         // 使用執行緒的最大個數 (最小為1)
        ppbm_chunk_size = atoi(argv[optind++]);         // 執行相同任務多執行緒分配任務的個數
        ppbm_option = atoi(argv[optind++]);             // 使用的平行化函式庫 1:OMP 2:CPP 3:TBB
        SIG_key = atoi(argv[optind++]);                 // 使用的特徵值  0b0001: check_unate, 0b0010 = check_ENE, 0b0100 = check_KSIG, 0b1000 = check_cofactor 可以同時使用多個
        MCQE_key = atoi(argv[optind++]);                // 使用的MCQE  0b0001: Rule1, 0b0010 = Rule2, 0b0100 = Rule3, 0b1000 = Rule4 可以同時使用多個
        Total_parallel = atoi(argv[optind++]);          // 是否執行全部平行化 (1/0)
        CAL_SIG_OUTPUT = 0; /*不計算輸出特徵值*/
    }
    else {
        printf("argv doesn't have pla\n");
    }

    

    if (argc == 17 || argc == 4) {
        /* 限制最大執行緒個數*/
        tbb::global_control MAXTHREADS(tbb::global_control::max_allowed_parallelism, ppbm_thread_num);
        /*tbb::task_arena limited(1); */
        print_solution = FALSE;
        BM(PLA, PLA1, isINP, isINN, isOUTP, isOUTN);
    }

    if (summary || trace) {
        if (PLA != NIL(PLA_t)) PLA_summary(PLA);
        if (PLA1 != NIL(PLA_t)) PLA_summary(PLA1);
    }

    /* Print a runtime summary if trace mode enabled */
    if (trace) {
        runtime();
    }

    /* Print total runtime */
    if (summary || trace) {
        print_trace(PLA->F, option_table[option].name, ptime() - start);
    }

    /* Output the solution */
    if (print_solution) {
        EXECUTE(fprint_pla(stdout, PLA, out_type), WRITE_TIME, PLA->F, cost);
    }

    /* Crash and burn if there was a verify error */
    if (error) {
        fatal("cover verification failed");
    }

    /* cleanup all used memory */
    free_PLA(PLA);
    FREE(cube.part_size);
    setdown_cube();             /* free the cube/cdata structure data */
    sf_cleanup();               /* free unused set structures */
    sm_cleanup();               /* sparse matrix cleanup */

    exit(0);
    return 0;
}

void getPLA(int opt, int argc, char** argv, pPLA* PLA, int out_type)
{
    FILE* fp;
    int needs_dcset, needs_offset;
    char* fname;

    if (opt >= argc) {
        fp = stdin;
        fname = "(stdin)";
    }
    else {
        fname = argv[opt];
        if ((fp = fopen(argv[opt], "r")) == NULL) {
            fprintf(stderr, "%s: Unable to open %s\n", argv[0], fname);
            exit(1);
        }
    }

    needs_dcset = TRUE; //default 
    needs_offset = TRUE;

    if (read_pla(fp, needs_dcset, needs_offset, input_type, PLA) == EOF) {
        fprintf(stderr, "%s: Unable to find PLA on file %s\n", argv[0], fname);
        exit(1);
    }
    (*PLA)->filename = _strdup(fname);
    filename = (*PLA)->filename;
    last_fp = fp;
}

void runtime(void)
{
    int i;
    long total = 1, temp;

    for (i = 0; i < TIME_COUNT; i++) {
        total += total_time[i];
    }
    for (i = 0; i < TIME_COUNT; i++) {
        if (total_calls[i] != 0) {
            temp = 100 * total_time[i];
            printf("# %s\t%2d call(s) for %s (%2ld.%01ld%%)\n",
                total_name[i], total_calls[i], print_time(total_time[i]),
                temp / total, (10 * (temp % total)) / total);
        }
    }
}

void init_runtime(void)
{
    total_name[READ_TIME] = "READ       ";
    total_name[WRITE_TIME] = "WRITE      ";
    total_name[COMPL_TIME] = "COMPL      ";
    total_name[REDUCE_TIME] = "REDUCE     ";
    total_name[EXPAND_TIME] = "EXPAND     ";
    total_name[ESSEN_TIME] = "ESSEN      ";
    total_name[IRRED_TIME] = "IRRED      ";
    total_name[GREDUCE_TIME] = "REDUCE_GASP";
    total_name[GEXPAND_TIME] = "EXPAND_GASP";
    total_name[GIRRED_TIME] = "IRRED_GASP ";
    total_name[MV_REDUCE_TIME] = "MV_REDUCE  ";
    total_name[RAISE_IN_TIME] = "RAISE_IN   ";
    total_name[VERIFY_TIME] = "VERIFY     ";
    total_name[PRIMES_TIME] = "PRIMES     ";
    total_name[MINCOV_TIME] = "MINCOV     ";
}

void subcommands(void)
{
    int i, col;
    printf("                ");
    col = 16;
    for (i = 0; option_table[i].name != 0; i++) {
        if ((col + strlen(option_table[i].name) + 1) > 76) {
            printf(",\n                ");
            col = 16;
        }
        else if (i != 0) {
            printf(", ");
        }
        printf("%s", option_table[i].name);
        col += strlen(option_table[i].name) + 2;
    }
    printf("\n");
}

void usage(void)
{
    printf("%s\n\n", VERSION);
    printf("SYNOPSIS: espresso [options] [file]\n\n");
    printf("  -d        Enable debugging\n");
    printf("  -e[opt]   Select espresso option:\n");
    printf("                fast, ness, nirr, nunwrap, onset, pos, strong,\n");
    printf("                eat, eatdots, kiss, random\n");
    printf("  -o[type]  Select output format:\n");
    printf("                f, fd, fr, fdr, pleasure, eqntott, kiss, cons\n");
    printf("  -rn-m     Select range for subcommands:\n");
    printf("                d1merge: first and last variables (0 ... m-1)\n");
    printf("                minterms: first and last variables (0 ... m-1)\n");
    printf("                opoall: first and last outputs (0 ... m-1)\n");
    printf("  -s        Provide short execution summary\n");
    printf("  -t        Provide longer execution trace\n");
    printf("  -x        Suppress printing of solution\n");
    printf("  -v[type]  Verbose debugging detail (-v '' for all)\n");
    printf("  -D[cmd]   Execute subcommand 'cmd':\n");
    subcommands();
    printf("  -Sn       Select strategy for subcommands:\n");
    printf("                opo: bit2=exact bit1=repeated bit0=skip sparse\n");
    printf("                opoall: 0=minimize, 1=exact\n");
    printf("                pair: 0=algebraic, 1=strongd, 2=espresso, 3=exact\n");
    printf("                pairall: 0=minimize, 1=exact, 2=opo\n");
    printf("                so_espresso: 0=minimize, 1=exact\n");
    printf("                so_both: 0=minimize, 1=exact\n");
}

/*
 *  Hack for backward compatibility (ACK! )
 */

void backward_compatibility_hack(int* argc, char** argv, int* option, int* out_type)
{
    int i, j;

    /* Scan the argument list for something to do (default is ESPRESSO) */
    *option = 0;
    for (i = 1; i < (*argc) - 1; i++) {
        if (strcmp(argv[i], "-do") == 0) {
            for (j = 0; option_table[j].name != 0; j++)
                if (strcmp(argv[i + 1], option_table[j].name) == 0) {
                    *option = j;
                    delete_arg(argc, argv, i + 1);
                    delete_arg(argc, argv, i);
                    break;
                }
            if (option_table[j].name == 0) {
                fprintf(stderr,
                    "espresso: bad keyword \"%s\" following -do\n", argv[i + 1]);
                exit(1);
            }
            break;
        }
    }

    for (i = 1; i < (*argc) - 1; i++) {
        if (strcmp(argv[i], "-out") == 0) {
            for (j = 0; pla_types[j].key != 0; j++)
                if (strcmp(pla_types[j].key + 1, argv[i + 1]) == 0) {
                    *out_type = pla_types[j].value;
                    delete_arg(argc, argv, i + 1);
                    delete_arg(argc, argv, i);
                    break;
                }
            if (pla_types[j].key == 0) {
                fprintf(stderr,
                    "espresso: bad keyword \"%s\" following -out\n", argv[i + 1]);
                exit(1);
            }
            break;
        }
    }

    for (i = 1; i < (*argc); i++) {
        if (argv[i][0] == '-') {
            for (j = 0; esp_opt_table[j].name != 0; j++) {
                if (strcmp(argv[i] + 1, esp_opt_table[j].name) == 0) {
                    delete_arg(argc, argv, i);
                    *(esp_opt_table[j].variable) = esp_opt_table[j].value;
                    break;
                }
            }
        }
    }

    if (check_arg(argc, argv, "-fdr")) input_type = FDR_type;
    if (check_arg(argc, argv, "-fr")) input_type = FR_type;
    if (check_arg(argc, argv, "-f")) input_type = F_type;
}


/* delete_arg -- delete an argument from the argument list */
void delete_arg(int* argc,  char** argv, int num)
{
     int i;
    (*argc)--;
    for (i = num; i < *argc; i++) {
        argv[i] = argv[i + 1];
    }
}


/* check_arg -- scan argv for an argument, and return TRUE if found */
bool check_arg(int* argc,  char** argv,  char* s)
{
     int i;
    for (i = 1; i < *argc; i++) {
        if (strcmp(argv[i], s) == 0) {
            delete_arg(argc, argv, i);
            return TRUE;
        }
    }
    return FALSE;
}