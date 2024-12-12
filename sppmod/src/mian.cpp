// rtklibzj.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>// C++ 中的预处理指令，用于包含标准输入输出库
#include <stdarg.h>//处理可变参数的函数
#include "rtklib.h"



#define PROGNAME    "rnx2rtkp"          /* program name */
#define MAXFILE     16                  /* max number of input files */

static const char* help[] = {
"",
" usage: rnx2rtkp [option]... file file [...]",
"",
" Read RINEX OBS/NAV/GNAV/HNAV/CLK, SP3, SBAS message log files and ccompute ",
" receiver (rover) positions and output position solutions.",
" The first RINEX OBS file shall contain receiver (rover) observations. For the",
" relative mode, the second RINEX OBS file shall contain reference",
" (base station) receiver observations. At least one RINEX NAV/GNAV/HNAV",
" file shall be included in input files. To use SP3 precise ephemeris, specify",
" the path in the files. The extension of the SP3 file shall be .sp3 or .eph.",
" All of the input file paths can include wild-cards (*). To avoid command",
" line deployment of wild-cards, use \"...\" for paths with wild-cards.",
" Command line options are as follows ([]:default). With -k option, the",
" processing options are input from the configuration file. In this case,",
" command line options precede options in the configuration file.",
"",
" -?        print help",
" -k file   input options from configuration file [off]",
" -o file   set output file [stdout]",
" -ts ds ts start day/time (ds=y/m/d ts=h:m:s) [obs start time]",
" -te de te end day/time   (de=y/m/d te=h:m:s) [obs end time]",
" -ti tint  time interval (sec) [all]",
" -p mode   mode (0:single,1:dgps,2:kinematic,3:static,4:moving-base,",
"                 5:fixed,6:ppp-kinematic,7:ppp-static) [2]",
" -m mask   elevation mask angle (deg) [15]",
" -sys s[,s...] nav system(s) (s=G:GPS,R:GLO,E:GAL,J:QZS,C:BDS,I:IRN) [G|R]",
" -f freq   number of frequencies for relative mode (1:L1,2:L1+L2,3:L1+L2+L5) [2]",
" -v thres  validation threshold for integer ambiguity (0.0:no AR) [3.0]",
" -b        backward solutions [off]",
" -c        forward/backward combined solutions [off]",
" -i        instantaneous integer ambiguity resolution [off]",
" -h        fix and hold for integer ambiguity resolution [off]",
" -e        output x/y/z-ecef position [latitude/longitude/height]",
" -a        output e/n/u-baseline [latitude/longitude/height]",
" -n        output NMEA-0183 GGA sentence [off]",
" -g        output latitude/longitude in the form of ddd mm ss.ss' [ddd.ddd]",
" -t        output time in the form of yyyy/mm/dd hh:mm:ss.ss [sssss.ss]",
" -u        output time in utc [gpst]",
" -d col    number of decimals in time [3]",
" -s sep    field separator [' ']",
" -r x y z  reference (base) receiver ecef pos (m) [average of single pos]",
"           rover receiver ecef pos (m) for fixed or ppp-fixed mode",
" -l lat lon hgt reference (base) receiver latitude/longitude/height (deg/m)",
"           rover latitude/longitude/height for fixed or ppp-fixed mode",
" -y level  output soltion status (0:off,1:states,2:residuals) [0]",
" -x level  debug trace level (0:off) [0]"
};

/* show message --------------------------------------------------------------*/
extern int showmsg(const char* format, ...)
{
    va_list arg;
    va_start(arg, format); vfprintf(stderr, format, arg); va_end(arg);
    fprintf(stderr, "\r");
    return 0;
}

extern void settspan(gtime_t ts, gtime_t te) {}
extern void settime(gtime_t time) {}


int main(int argc, char** argv)//参数个,数参数值
{
    //std::cout << "Hello World!\n";
    prcopt_t prcopt = prcopt_default;//处理配置
    solopt_t solopt = solopt_default;//结果配置
    filopt_t filopt = { "" };
    gtime_t ts = { 0 }, te = { 0 };
    double tint = 0.0, es[] = { 2000,1,1,0,0,0 }, ee[] = { 2000,12,31,23,59,59 }, pos[3];
    int i, j, n, ret;
    //char *infile[MAXFILE], char *outfile = "",char *p;
    char* infile[MAXFILE]; // 声明一个指向字符指针的数组
    char* outfile = "";    // 声明一个指向字符的指针，并初始化为空字符串
    char* p;               // 声明一个指向字符的指针

  
    prcopt.mode = PMODE_KINEMA;
    prcopt.navsys = 0;
    prcopt.refpos = 1;
    prcopt.glomodear = 1;
    solopt.timef = 0;
    sprintf_s(solopt.prog, "%s ver.%s %s", PROGNAME, VER_RTKLIB, PATCH_LEVEL);
    sprintf_s(filopt.trace, "%s.trace", PROGNAME);

    /* load options from configuration file */
    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-k") && i + 1 < argc) {
            resetsysopts();
            if (!loadopts(argv[++i], sysopts)) return -1;//
            getsysopts(&prcopt, &solopt, &filopt);
        }
    }
    for (i = 1, n = 0; i < argc; i++) {
        if (!strcmp(argv[i], "-o") && i + 1 < argc) outfile = argv[++i];
        else if (!strcmp(argv[i], "-ts") && i + 2 < argc) {
            sscanf_s(argv[++i], "%lf/%lf/%lf", es, es + 1, es + 2);
            sscanf_s(argv[++i], "%lf:%lf:%lf", es + 3, es + 4, es + 5);
            ts = epoch2time(es);
        }
        else if (!strcmp(argv[i], "-te") && i + 2 < argc) {
            sscanf_s(argv[++i], "%lf/%lf/%lf", ee, ee + 1, ee + 2);
            sscanf_s(argv[++i], "%lf:%lf:%lf", ee + 3, ee + 4, ee + 5);
            te = epoch2time(ee);
        }
        else if (!strcmp(argv[i], "-ti") && i + 1 < argc) tint = atof(argv[++i]);
        else if (!strcmp(argv[i], "-k") && i + 1 < argc) { ++i; continue; }
        else if (!strcmp(argv[i], "-p") && i + 1 < argc) prcopt.mode = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-f") && i + 1 < argc) prcopt.nf = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-sys") && i + 1 < argc) {
            for (p = argv[++i]; *p; p++) {
                switch (*p) {
                case 'G': prcopt.navsys |= SYS_GPS;
                case 'R': prcopt.navsys |= SYS_GLO;
                case 'E': prcopt.navsys |= SYS_GAL;
                case 'J': prcopt.navsys |= SYS_QZS;
                case 'C': prcopt.navsys |= SYS_CMP;
                case 'I': prcopt.navsys |= SYS_IRN;
                }
                if (!(p = strchr(p, ','))) break;
            }
        }
        else if (!strcmp(argv[i], "-m") && i + 1 < argc) prcopt.elmin = atof(argv[++i]) * D2R;
        else if (!strcmp(argv[i], "-v") && i + 1 < argc) prcopt.thresar[0] = atof(argv[++i]);
        else if (!strcmp(argv[i], "-s") && i + 1 < argc) strcpy(solopt.sep, argv[++i]);
        else if (!strcmp(argv[i], "-d") && i + 1 < argc) solopt.timeu = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-b")) prcopt.soltype = 1;
        else if (!strcmp(argv[i], "-c")) prcopt.soltype = 2;
        else if (!strcmp(argv[i], "-i")) prcopt.modear = 2;
        else if (!strcmp(argv[i], "-h")) prcopt.modear = 3;
        else if (!strcmp(argv[i], "-t")) solopt.timef = 1;
        else if (!strcmp(argv[i], "-u")) solopt.times = TIMES_UTC;
        else if (!strcmp(argv[i], "-e")) solopt.posf = SOLF_XYZ;
        else if (!strcmp(argv[i], "-a")) solopt.posf = SOLF_ENU;
        else if (!strcmp(argv[i], "-n")) solopt.posf = SOLF_NMEA;
        else if (!strcmp(argv[i], "-g")) solopt.degf = 1;
        else if (!strcmp(argv[i], "-r") && i + 3 < argc) {
            prcopt.refpos = prcopt.rovpos = 0;
            for (j = 0; j < 3; j++) prcopt.rb[j] = atof(argv[++i]);
            matcpy(prcopt.ru, prcopt.rb, 3, 1);
        }
        else if (!strcmp(argv[i], "-l") && i + 3 < argc) {
            prcopt.refpos = prcopt.rovpos = 0;
            for (j = 0; j < 3; j++) pos[j] = atof(argv[++i]);
            for (j = 0; j < 2; j++) pos[j] *= D2R;
            pos2ecef(pos, prcopt.rb);
            matcpy(prcopt.ru, prcopt.rb, 3, 1);
        }
        else if (!strcmp(argv[i], "-y") && i + 1 < argc) solopt.sstat = atoi(argv[++i]);//strcmp 函数用于比较两个字符串，如果相等则返回 0
        else if (!strcmp(argv[i], "-x") && i + 1 < argc) solopt.trace = atoi(argv[++i]);
        //else if (*argv[i] == '-') printhelp();
        else if (n < MAXFILE) infile[n++] = argv[i];
    }
    if (!prcopt.navsys) {
        prcopt.navsys = SYS_GPS | SYS_GLO;
    }
    if (n <= 0) {
        showmsg("error : no input file");
        return -2;
    }
    ret = postpos(ts, te, tint, 0.0, &prcopt, &solopt, &filopt, infile, n, outfile, "", "");//开始时间结束时间采样频率

    if (!ret) fprintf(stderr, "%40s\r", "");
    return ret;

}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
