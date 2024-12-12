#include <stdarg.h>
#include "rtklib.h"

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define MIN(x,y)    ((x)<=(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

#define VAR_POS     SQR(30.0) /* initial variance of receiver pos (m^2) */
#define VAR_VEL     SQR(10.0) /* initial variance of receiver vel ((m/s)^2) */
#define VAR_ACC     SQR(10.0) /* initial variance of receiver acc ((m/ss)^2) */
#define VAR_HWBIAS  SQR(1.0)  /* initial variance of h/w bias ((m/MHz)^2) */
#define VAR_GRA     SQR(0.001) /* initial variance of gradient (m^2) */
#define INIT_ZWD    0.15     /* initial zwd (m) */

#define PRN_HWBIAS  1E-6     /* process noise of h/w bias (m/MHz/sqrt(s)) */
#define GAP_RESION  120      /* gap to reset ionosphere parameters (epochs) */
#define MAXACC      30.0     /* max accel for doppler slip detection (m/s^2) */

#define VAR_HOLDAMB 0.001    /* constraint to hold ambiguity (cycle^2) */

#define TTOL_MOVEB  (1.0+2*DTTOL)
/* time sync tolerance for moving-baseline (s) */

/* number of parameters (pos,ionos,tropos,hw-bias,phase-bias,real,estimated) */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics==0?3:9)
#define NI(opt)     ((opt)->ionoopt!=IONOOPT_EST?0:MAXSAT)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?2:6))
#define NL(opt)     ((opt)->glomodear!=2?0:NFREQGLO)
#define NB(opt)     ((opt)->mode<=PMODE_DGPS?0:MAXSAT*NF(opt))
#define NR(opt)     (NP(opt)+NI(opt)+NT(opt)+NL(opt))
#define NX(opt)     (NR(opt)+NB(opt))

/* state variable index */
#define II(s,opt)   (NP(opt)+(s)-1)                 /* ionos (s:satellite no) */
#define IT(r,opt)   (NP(opt)+NI(opt)+NT(opt)/2*(r)) /* tropos (r:0=rov,1:ref) */
#define IL(f,opt)   (NP(opt)+NI(opt)+NT(opt)+(f))   /* receiver h/w bias */
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1) /* phase bias (s:satno,f:freq) */

/* global variables ----------------------------------------------------------*/
static int statlevel = 0;          /* rtk status output level (0:off) */
static FILE* fp_stat = NULL;       /* rtk status file pointer */
static char file_stat[1024] = "";  /* rtk status file original path */
static gtime_t time_stat = { 0 };    /* rtk status file time */

/* open solution status file ---------------------------------------------------
* open solution status file and set output level
* args   : char     *file   I   rtk status file
*          int      level   I   rtk status level (0: off)
* return : status (1:ok,0:error)
* notes  : file can constain time keywords (%Y,%y,%m...) defined in reppath().
*          The time to replace keywords is based on UTC of CPU time.
* output : solution status file record format
*
*   $POS,week,tow,stat,posx,posy,posz,posxf,posyf,poszf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          posx/posy/posz    : position x/y/z ecef (m) float
*          posxf/posyf/poszf : position x/y/z ecef (m) fixed
*
*   $VELACC,week,tow,stat,vele,veln,velu,acce,accn,accu,velef,velnf,veluf,accef,accnf,accuf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          vele/veln/velu    : velocity e/n/u (m/s) float
*          acce/accn/accu    : acceleration e/n/u (m/s^2) float
*          velef/velnf/veluf : velocity e/n/u (m/s) fixed
*          accef/accnf/accuf : acceleration e/n/u (m/s^2) fixed
*
*   $CLK,week,tow,stat,clk1,clk2,clk3,clk4
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          clk1     : receiver clock bias GPS (ns)
*          clk2     : receiver clock bias GLO-GPS (ns)
*          clk3     : receiver clock bias GAL-GPS (ns)
*          clk4     : receiver clock bias BDS-GPS (ns)
*
*   $ION,week,tow,stat,sat,az,el,ion,ion-fixed
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          sat      : satellite id
*          az/el    : azimuth/elevation angle(deg)
*          ion      : vertical ionospheric delay L1 (m) float
*          ion-fixed: vertical ionospheric delay L1 (m) fixed
*
*   $TROP,week,tow,stat,rcv,ztd,ztdf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          rcv      : receiver (1:rover,2:base station)
*          ztd      : zenith total delay (m) float
*          ztdf     : zenith total delay (m) fixed
*
*   $HWBIAS,week,tow,stat,frq,bias,biasf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          frq      : frequency (1:L1,2:L2,...)
*          bias     : h/w bias coefficient (m/MHz) float
*          biasf    : h/w bias coefficient (m/MHz) fixed
*
*   $SAT,week,tow,sat,frq,az,el,resp,resc,vsat,snr,fix,slip,lock,outc,slipc,rejc
*          week/tow : gps week no/time of week (s)
*          sat/frq  : satellite id/frequency (1:L1,2:L2,...)
*          az/el    : azimuth/elevation angle (deg)
*          resp     : pseudorange residual (m)
*          resc     : carrier-phase residual (m)
*          vsat     : valid data flag (0:invalid,1:valid)
*          snr      : signal strength (dbHz)
*          fix      : ambiguity flag  (0:no data,1:float,2:fixed,3:hold,4:ppp)
*          slip     : cycle-slip flag (bit1:slip,bit2:parity unknown)
*          lock     : carrier-lock count
*          outc     : data outage count
*          slipc    : cycle-slip count
*          rejc     : data reject (outlier) count
*
*-----------------------------------------------------------------------------*/
extern int rtkopenstat(const char* file, int level)
{
    gtime_t time = utc2gpst(timeget());
    char path[1024];

    trace(3, "rtkopenstat: file=%s level=%d\n", file, level);

    if (level <= 0) return 0;

    reppath(file, path, time, "", "");

    if (!(fp_stat = fopen(path, "w"))) {
        trace(1, "rtkopenstat: file open error path=%s\n", path);
        return 0;
    }
    strcpy(file_stat, file);
    time_stat = time;
    statlevel = level;
    return 1;
}

/* close solution status file --------------------------------------------------
* close solution status file
* args   : none
* return : none
*-----------------------------------------------------------------------------*/
extern void rtkclosestat(void)
{
    trace(3, "rtkclosestat:\n");

    if (fp_stat) fclose(fp_stat);
    fp_stat = NULL;
    file_stat[0] = '\0';
    statlevel = 0;
}

/* write solution status to buffer -------------------------------------------*/
extern int rtkoutstat(rtk_t* rtk, char* buff)
{
    ssat_t* ssat;
    double tow, pos[3], vel[3], acc[3], vela[3] = { 0 }, acca[3] = { 0 }, xa[3];
    int i, j, week, est, nfreq, nf = NF(&rtk->opt);
    char id[32], * p = buff;

    if (rtk->sol.stat <= SOLQ_NONE) {
        return 0;
    }
    /* write ppp solution status to buffer */
    /*if (rtk->opt.mode >= PMODE_PPP_KINEMA) {
        return pppoutstat(rtk, buff);
    }*/
    est = rtk->opt.mode >= PMODE_DGPS;
    nfreq = est ? nf : 1;
    tow = time2gpst(rtk->sol.time, &week);

    /* receiver position */
    if (est) {
        for (i = 0; i < 3; i++) xa[i] = i < rtk->na ? rtk->xa[i] : 0.0;
        p += sprintf(p, "$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", week, tow,
            rtk->sol.stat, rtk->x[0], rtk->x[1], rtk->x[2], xa[0], xa[1],
            xa[2]);
    }
    else {
        p += sprintf(p, "$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", week, tow,
            rtk->sol.stat, rtk->sol.rr[0], rtk->sol.rr[1], rtk->sol.rr[2],
            0.0, 0.0, 0.0);
    }
    /* receiver velocity and acceleration */
    if (est && rtk->opt.dynamics) {
        ecef2pos(rtk->sol.rr, pos);
        ecef2enu(pos, rtk->x + 3, vel);
        ecef2enu(pos, rtk->x + 6, acc);
        if (rtk->na >= 6) ecef2enu(pos, rtk->xa + 3, vela);
        if (rtk->na >= 9) ecef2enu(pos, rtk->xa + 6, acca);
        p += sprintf(p, "$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f\n",
            week, tow, rtk->sol.stat, vel[0], vel[1], vel[2], acc[0], acc[1],
            acc[2], vela[0], vela[1], vela[2], acca[0], acca[1], acca[2]);
    }
    else {
        ecef2pos(rtk->sol.rr, pos);
        ecef2enu(pos, rtk->sol.rr + 3, vel);
        p += sprintf(p, "$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f\n",
            week, tow, rtk->sol.stat, vel[0], vel[1], vel[2],
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }
    /* receiver clocks */
    p += sprintf(p, "$CLK,%d,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f\n",
        week, tow, rtk->sol.stat, 1, rtk->sol.dtr[0] * 1E9, rtk->sol.dtr[1] * 1E9,
        rtk->sol.dtr[2] * 1E9, rtk->sol.dtr[3] * 1E9);

    /* ionospheric parameters */
    if (est && rtk->opt.ionoopt == IONOOPT_EST) {
        for (i = 0; i < MAXSAT; i++) {
            ssat = rtk->ssat + i;
            if (!ssat->vs) continue;
            satno2id(i + 1, id);
            j = II(i + 1, &rtk->opt);
            xa[0] = j < rtk->na ? rtk->xa[j] : 0.0;
            p += sprintf(p, "$ION,%d,%.3f,%d,%s,%.1f,%.1f,%.4f,%.4f\n", week, tow,
                rtk->sol.stat, id, ssat->azel[0] * R2D, ssat->azel[1] * R2D,
                rtk->x[j], xa[0]);
        }
    }
    /* tropospheric parameters */
    if (est && (rtk->opt.tropopt == TROPOPT_EST || rtk->opt.tropopt == TROPOPT_ESTG)) {
        for (i = 0; i < 2; i++) {
            j = IT(i, &rtk->opt);
            xa[0] = j < rtk->na ? rtk->xa[j] : 0.0;
            p += sprintf(p, "$TROP,%d,%.3f,%d,%d,%.4f,%.4f\n", week, tow,
                rtk->sol.stat, i + 1, rtk->x[j], xa[0]);
        }
    }
    /* receiver h/w bias */
    if (est && rtk->opt.glomodear == 2) {
        for (i = 0; i < nfreq; i++) {
            j = IL(i, &rtk->opt);
            xa[0] = j < rtk->na ? rtk->xa[j] : 0.0;
            p += sprintf(p, "$HWBIAS,%d,%.3f,%d,%d,%.4f,%.4f\n", week, tow,
                rtk->sol.stat, i + 1, rtk->x[j], xa[0]);
        }
    }
    return (int)(p - buff);
}

/* save error message --------------------------------------------------------*/
static void errmsg(rtk_t* rtk, const char* format, ...)
{
    char buff[256], tstr[32];
    int n;
    va_list ap;
    time2str(rtk->sol.time, tstr, 2);
    n = sprintf(buff, "%s: ", tstr + 11);
    va_start(ap, format);
    n += vsprintf(buff + n, format, ap);
    va_end(ap);
    n = n < MAXERRMSG - rtk->neb ? n : MAXERRMSG - rtk->neb;
    memcpy(rtk->errbuf + rtk->neb, buff, n);
    rtk->neb += n;
    trace(2, "%s", buff);
}

/* swap solution status file -------------------------------------------------*/
static void swapsolstat(void)
{
    gtime_t time = utc2gpst(timeget());
    char path[1024];

    if ((int)(time2gpst(time, NULL) / INT_SWAP_STAT) ==
        (int)(time2gpst(time_stat, NULL) / INT_SWAP_STAT)) {
        return;
    }
    time_stat = time;

    if (!reppath(file_stat, path, time, "", "")) {
        return;
    }
    if (fp_stat) fclose(fp_stat);

    if (!(fp_stat = fopen(path, "w"))) {
        trace(2, "swapsolstat: file open error path=%s\n", path);
        return;
    }
    trace(3, "swapsolstat: path=%s\n", path);
}

/* output solution status ----------------------------------------------------*/
static void outsolstat(rtk_t* rtk)
{
    ssat_t* ssat;
    double tow;
    char buff[MAXSOLMSG + 1], id[32];
    int i, j, n, week, nfreq, nf = NF(&rtk->opt);

    if (statlevel <= 0 || !fp_stat || !rtk->sol.stat) return;

    trace(3, "outsolstat:\n");

    /* swap solution status file */
    swapsolstat();

    /* write solution status */
    n = rtkoutstat(rtk, buff); buff[n] = '\0';

    fputs(buff, fp_stat);

    if (rtk->sol.stat == SOLQ_NONE || statlevel <= 1) return;

    tow = time2gpst(rtk->sol.time, &week);
    nfreq = rtk->opt.mode >= PMODE_DGPS ? nf : 1;

    /* write residuals and status */
    for (i = 0; i < MAXSAT; i++) {
        ssat = rtk->ssat + i;
        if (!ssat->vs) continue;
        satno2id(i + 1, id);
        for (j = 0; j < nfreq; j++) {
            fprintf(fp_stat, "$SAT,%d,%.3f,%s,%d,%.1f,%.1f,%.4f,%.4f,%d,%.1f,%d,%d,%d,%d,%d,%d\n",
                week, tow, id, j + 1, ssat->azel[0] * R2D, ssat->azel[1] * R2D,
                ssat->resp[j], ssat->resc[j], ssat->vsat[j],
                ssat->snr[j] * SNR_UNIT, ssat->fix[j], ssat->slip[j] & 3,
                ssat->lock[j], ssat->outc[j], ssat->slipc[j], ssat->rejc[j]);
        }
    }
}

/* initialize RTK control ------------------------------------------------------
* initialize RTK control struct
* args   : rtk_t    *rtk    IO  TKk control/result struct
*          prcopt_t *opt    I   positioning options (see rtklib.h)
* return : none
*-----------------------------------------------------------------------------*/
extern void rtkinit(rtk_t* rtk, const prcopt_t* opt)
{
    sol_t sol0 = { {0} };
    ambc_t ambc0 = { {{0}} };
    ssat_t ssat0 = { 0 };
    int i;

    trace(3, "rtkinit :\n");

    rtk->sol = sol0;
    for (i = 0; i < 6; i++) rtk->rb[i] = 0.0;
    //rtk->nx = opt->mode <= PMODE_FIXED ? NX(opt) : pppnx(opt);
    //rtk->na = opt->mode <= PMODE_FIXED ? NR(opt) : pppnx(opt);
    rtk->tt = 0.0;
    rtk->x = zeros(rtk->nx, 1);//
    rtk->P = zeros(rtk->nx, rtk->nx);
    rtk->xa = zeros(rtk->na, 1);
    rtk->Pa = zeros(rtk->na, rtk->na);
    rtk->nfix = rtk->neb = 0;
    for (i = 0; i < MAXSAT; i++) {
        rtk->ambc[i] = ambc0;
        rtk->ssat[i] = ssat0;
    }
    for (i = 0; i < MAXERRMSG; i++) rtk->errbuf[i] = 0;
    rtk->opt = *opt;
}

/* free rtk control ------------------------------------------------------------
* free memory for rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
* return : none
*-----------------------------------------------------------------------------*/
extern void rtkfree(rtk_t* rtk)
{
    trace(3, "rtkfree :\n");

    rtk->nx = rtk->na = 0;
    free(rtk->x); rtk->x = NULL;
    free(rtk->P); rtk->P = NULL;
    free(rtk->xa); rtk->xa = NULL;
    free(rtk->Pa); rtk->Pa = NULL;
}

/* precise positioning ---------------------------------------------------------
* input observation data and navigation message, compute rover position by
* precise positioning
* args   : rtk_t *rtk       IO  RTK control/result struct
*            rtk->sol       IO  solution
*                .time      O   solution time
*                .rr[]      IO  rover position/velocity
*                               (I:fixed mode,O:single mode)
*                .dtr[0]    O   receiver clock bias (s)
*                .dtr[1-5]  O   receiver GLO/GAL/BDS/IRN/QZS-GPS time offset (s)
*                .Qr[]      O   rover position covarinace
*                .stat      O   solution status (SOLQ_???)
*                .ns        O   number of valid satellites
*                .age       O   age of differential (s)
*                .ratio     O   ratio factor for ambiguity validation
*            rtk->rb[]      IO  base station position/velocity
*                               (I:relative mode,O:moving-base mode)
*            rtk->nx        I   number of all states
*            rtk->na        I   number of integer states
*            rtk->ns        O   number of valid satellites in use
*            rtk->tt        O   time difference between current and previous (s)
*            rtk->x[]       IO  float states pre-filter and post-filter
*            rtk->P[]       IO  float covariance pre-filter and post-filter
*            rtk->xa[]      O   fixed states after AR
*            rtk->Pa[]      O   fixed covariance after AR
*            rtk->ssat[s]   IO  satellite {s+1} status
*                .sys       O   system (SYS_???)
*                .az   [r]  O   azimuth angle   (rad) (r=0:rover,1:base)
*                .el   [r]  O   elevation angle (rad) (r=0:rover,1:base)
*                .vs   [r]  O   data valid single     (r=0:rover,1:base)
*                .resp [f]  O   freq(f+1) pseudorange residual (m)
*                .resc [f]  O   freq(f+1) carrier-phase residual (m)
*                .vsat [f]  O   freq(f+1) data vaild (0:invalid,1:valid)
*                .fix  [f]  O   freq(f+1) ambiguity flag
*                               (0:nodata,1:float,2:fix,3:hold)
*                .slip [f]  O   freq(f+1) cycle slip flag
*                               (bit8-7:rcv1 LLI, bit6-5:rcv2 LLI,
*                                bit2:parity unknown, bit1:slip)
*                .lock [f]  IO  freq(f+1) carrier lock count
*                .outc [f]  IO  freq(f+1) carrier outage count
*                .slipc[f]  IO  freq(f+1) cycle slip count
*                .rejc [f]  IO  freq(f+1) data reject count
*                .gf        IO  geometry-free phase (L1-L2) (m)
*                .gf2       IO  geometry-free phase (L1-L5) (m)
*            rtk->nfix      IO  number of continuous fixes of ambiguity
*            rtk->neb       IO  bytes of error message buffer
*            rtk->errbuf    IO  error message buffer
*            rtk->tstr      O   time string for debug
*            rtk->opt       I   processing options
*          obsd_t *obs      I   observation data for an epoch
*                               obs[i].rcv=1:rover,2:reference
*                               sorted by receiver and satellte
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation messages
* return : status (0:no solution,1:valid solution)
* notes  : before calling function, base station position rtk->sol.rb[] should
*          be properly set for relative mode except for moving-baseline
*-----------------------------------------------------------------------------*/
extern int rtkpos(rtk_t* rtk, const obsd_t* obs, int n, const nav_t* nav)
{
    prcopt_t* opt = &rtk->opt;
    sol_t solb = { {0} };
    gtime_t time;
    int i, nu, nr;
    char msg[128] = "";

    trace(3, "rtkpos  : time=%s n=%d\n", time_str(obs[0].time, 3), n);
    trace(4, "obs=\n"); traceobs(4, obs, n);

    /* set base staion position */
    if (opt->refpos <= POSOPT_RINEX && opt->mode != PMODE_SINGLE &&
        opt->mode != PMODE_MOVEB) {
        for (i = 0; i < 6; i++) rtk->rb[i] = i < 3 ? opt->rb[i] : 0.0;
    }
    /* count rover/base station observations */
    for (nu = 0; nu < n && obs[nu].rcv == 1; nu++);//����վ
    for (nr = 0; nu + nr < n && obs[nu + nr].rcv == 2; nr++);//��׼վ

    time = rtk->sol.time; /* previous epoch */

    /* rover position by single point positioning *////���õ��㶨λ��λ��
    if (!pntpos(obs, nu, nav, &rtk->opt, &rtk->sol, NULL, rtk->ssat, msg)) {
        errmsg(rtk, "point pos error (%s)\n", msg);

        if (!rtk->opt.dynamics) {
            outsolstat(rtk);
            return 0;
        }
    }
    if (time.time != 0) rtk->tt = timediff(rtk->sol.time, time);

    /* single point positioning */
    if (opt->mode == PMODE_SINGLE) {//ģʽ���е��㶨λ�Ļ�����������
        outsolstat(rtk);
        return 1;
    }
    /* suppress output of single solution */
    if (!opt->outsingle) {
        rtk->sol.stat = SOLQ_NONE;
    }
    /* precise point positioning */
    /*if (opt->mode >= PMODE_PPP_KINEMA) {
        pppos(rtk, obs, nu, nav);
        outsolstat(rtk);
        return 1;
    }*/
    /* check number of data of base station and age of differential */
    if (nr == 0) {//��ֶ�λ
        errmsg(rtk, "no base station observation data for rtk\n");
        outsolstat(rtk);
        return 1;
    }
    //if (opt->mode == PMODE_MOVEB) { /*  moving baseline */

    //    /* estimate position/velocity of base station */
    //    if (!pntpos(obs + nu, nr, nav, &rtk->opt, &solb, NULL, NULL, msg)) {
    //        errmsg(rtk, "base station position error (%s)\n", msg);
    //        return 0;
    //    }
    //    rtk->sol.age = (float)timediff(rtk->sol.time, solb.time);

    //    if (fabs(rtk->sol.age) > TTOL_MOVEB) {
    //        errmsg(rtk, "time sync error for moving-base (age=%.1f)\n", rtk->sol.age);
    //        return 0;
    //    }
    //    for (i = 0; i < 6; i++) rtk->rb[i] = solb.rr[i];

    //    /* time-synchronized position of base station */
    //    for (i = 0; i < 3; i++) rtk->rb[i] += rtk->rb[i + 3] * rtk->sol.age;
    //}
    else {
        rtk->sol.age = (float)timediff(obs[0].time, obs[nu].time);

        if (fabs(rtk->sol.age) > opt->maxtdiff) {
            errmsg(rtk, "age of differential error (age=%.1f)\n", rtk->sol.age);
            outsolstat(rtk);
            return 1;
        }
    }
    ///* relative potitioning */
    //relpos(rtk, obs, nu, nr, nav);
    //outsolstat(rtk);

    return 1;
}
