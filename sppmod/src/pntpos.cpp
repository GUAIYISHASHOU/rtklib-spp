#include "rtklib.h"

/* constants/macros ----------------------------------------------------------*/

#define SQR(x)      ((x)*(x))

#if 0 /* enable GPS-QZS time offset estimation */
#define NX          (4+5)       /* # of estimated parameters */
#else
#define NX          (4+4)       /* # of estimated parameters */
#endif
#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_ION     5.0         /* ionospheric delay Std (m) */
#define ERR_TROP    3.0         /* tropspheric delay Std (m) */
#define ERR_SAAS    0.3         /* Saastamoinen model error Std (m) */
#define ERR_BRDCI   0.5         /* broadcast ionosphere model error factor */
#define ERR_CBIAS   0.3         /* code bias error Std (m) */
#define REL_HUMI    0.7         /* relative humidity for Saastamoinen model */
#define MIN_EL      (5.0*D2R)   /* min elevation for measurement error (rad) */


///////////////////////////////////////////////////////////////////////////////////////////////
// ephememeris内容
#define SQR(x)   ((x)*(x))

#define RE_GLO   6378136.0        /* radius of earth (m)            ref [2] */
#define MU_GPS   3.9860050E14     /* gravitational constant         ref [1] */
#define MU_GLO   3.9860044E14     /* gravitational constant         ref [2] */
#define MU_GAL   3.986004418E14   /* earth gravitational constant   ref [7] */
#define MU_CMP   3.986004418E14   /* earth gravitational constant   ref [9] */
#define J2_GLO   1.0826257E-3     /* 2nd zonal harmonic of geopot   ref [2] */

#define OMGE_GLO 7.292115E-5      /* earth angular velocity (rad/s) ref [2] */
#define OMGE_GAL 7.2921151467E-5  /* earth angular velocity (rad/s) ref [7] */
#define OMGE_CMP 7.292115E-5      /* earth angular velocity (rad/s) ref [9] */

#define SIN_5 -0.0871557427476582 /* sin(-5.0 deg) */
#define COS_5  0.9961946980917456 /* cos(-5.0 deg) */

#define ERREPH_GLO 5.0            /* error of glonass ephemeris (m) */
#define TSTEP    60.0             /* integration step glonass ephemeris (s) */
#define RTOL_KEPLER 1E-13         /* relative tolerance for Kepler equation */

#define DEFURASSR 0.15            /* default accurary of ssr corr (m) */
#define MAXECORSSR 10.0           /* max orbit correction of ssr (m) */
#define MAXCCORSSR (1E-6*CLIGHT)  /* max clock correction of ssr (m) */
#define MAXAGESSR 90.0            /* max age of ssr orbit and clock (s) */
#define MAXAGESSR_HRCLK 10.0      /* max age of ssr high-rate clock (s) */
#define STD_BRDCCLK 30.0          /* error of broadcast clock (m) */
#define STD_GAL_NAPA 500.0        /* error of galileo ephemeris for NAPA (m) */

#define MAX_ITER_KEPLER 30        /* max number of iteration of Kelpler */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



static int eph_sel[] = { /* GPS,GLO,GAL,QZS,BDS,IRN,SBS */
    0,0,0,0,0,0,0
};//初始化
//extern void setseleph(int sys, int sel)
//{
//    switch (sys) {
//    case SYS_GPS: eph_sel[0] = sel; break;
//    case SYS_GLO: eph_sel[1] = sel; break;
//    case SYS_GAL: eph_sel[2] = sel; break;
//    case SYS_QZS: eph_sel[3] = sel; break;
//    case SYS_CMP: eph_sel[4] = sel; break;
//    case SYS_IRN: eph_sel[5] = sel; break;
//    case SYS_SBS: eph_sel[6] = sel; break;
//    }
//}



/* ionospheric correction ------------------------------------------------------
* compute ionospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          int    sat       I   satellite number
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    ionoopt   I   ionospheric correction option (IONOOPT_???)
*          double *ion      O   ionospheric delay (L1) (m)
*          double *var      O   ionospheric delay (L1) variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int ionocorr(gtime_t time, const nav_t* nav, int sat, const double* pos,
    const double* azel, int ionoopt, double* ion, double* var)//时间星历卫星编号位置
    //高度角方位角 
{
    trace(4, "ionocorr: time=%s opt=%d sat=%2d pos=%.3f %.3f azel=%.3f %.3f\n",
        time_str(time, 3), ionoopt, sat, pos[0] * R2D, pos[1] * R2D, azel[0] * R2D,
        azel[1] * R2D);

    /* GPS broadcast ionosphere model */
    if (ionoopt == IONOOPT_BRDC) {
        *ion = ionmodel(time, nav->ion_gps, pos, azel);
        *var = SQR(*ion * ERR_BRDCI);
        return 1;
    }
    /* SBAS ionosphere model */
 /*   if (ionoopt == IONOOPT_SBAS) {
        return sbsioncorr(time, nav, pos, azel, ion, var);
    }*/
    /* IONEX TEC model */
    if (ionoopt == IONOOPT_TEC) {
        return iontec(time, nav, pos, azel, 1, ion, var);
    }
    /* QZSS broadcast ionosphere model */
    if (ionoopt == IONOOPT_QZS && norm(nav->ion_qzs, 8) > 0.0) {
        *ion = ionmodel(time, nav->ion_qzs, pos, azel);
        *var = SQR(*ion * ERR_BRDCI);
        return 1;
    }
    *ion = 0.0;
    *var = ionoopt == IONOOPT_OFF ? SQR(ERR_ION) : 0.0;
    return 1;
}



/* single-point positioning ----------------------------------------------------
* compute receiver position, velocity, clock bias by single-point positioning
* with pseudorange and doppler observables
* args   : obsd_t *obs      I   observation data
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation data
*          prcopt_t *opt    I   processing options
*          sol_t  *sol      IO  solution
*          double *azel     IO  azimuth/elevation angle (rad) (NULL: no output)//高度角方位角
*          ssat_t *ssat     IO  satellite status              (NULL: no output)//卫星空间信息
*          char   *msg      O   error message for error exit
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
/* pseudorange measurement error variance ------------------------------------*/
static double varerr(const prcopt_t* opt, double el, int sys)
{
    double fact, varr;
    fact = sys == SYS_GLO ? EFACT_GLO : (sys == SYS_SBS ? EFACT_SBS : EFACT_GPS);//不同系统权重不一样
    if (el < MIN_EL) el = MIN_EL;
    varr = SQR(opt->err[0]) * (SQR(opt->err[1]) + SQR(opt->err[2]) / sin(el));
    if (opt->ionoopt == IONOOPT_IFLC) varr *= SQR(3.0); /* iono-free */
    return SQR(fact) * varr;
}
/* get group delay parameter (m) ---------------------------------------------*/
static double gettgd(int sat, const nav_t* nav, int type)
{
    int i, sys = satsys(sat, NULL);

    if (sys == SYS_GLO) {
        for (i = 0; i < nav->ng; i++) {
            if (nav->geph[i].sat == sat) break;
        }
        return (i >= nav->ng) ? 0.0 : -nav->geph[i].dtaun * CLIGHT;
    }
    else {
        for (i = 0; i < nav->n; i++) {
            if (nav->eph[i].sat == sat) break;
        }
        return (i >= nav->n) ? 0.0 : nav->eph[i].tgd[type] * CLIGHT;
    }
}
/* test SNR mask -------------------------------------------------------------*/
static int snrmask(const obsd_t* obs, const double* azel, const prcopt_t* opt)
{
    if (testsnr(0, 0, azel[1], obs->SNR[0] * SNR_UNIT, &opt->snrmask)) {//根据配置文件剔除
        return 0;
    }
    if (opt->ionoopt == IONOOPT_IFLC) {
        if (testsnr(0, 1, azel[1], obs->SNR[1] * SNR_UNIT, &opt->snrmask)) return 0;
    }
    return 1;
}

/* tropospheric correction -----------------------------------------------------
* compute tropospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    tropopt   I   tropospheric correction option (TROPOPT_???)
*          double *trp      O   tropospheric delay (m)
*          double *var      O   tropospheric delay variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int tropcorr(gtime_t time, const nav_t* nav, const double* pos,
    const double* azel, int tropopt, double* trp, double* var)
{
    trace(4, "tropcorr: time=%s opt=%d pos=%.3f %.3f azel=%.3f %.3f\n",
        time_str(time, 3), tropopt, pos[0] * R2D, pos[1] * R2D, azel[0] * R2D,
        azel[1] * R2D);

    /* Saastamoinen model *///对流层改正模型
    if (tropopt == TROPOPT_SAAS || tropopt == TROPOPT_EST || tropopt == TROPOPT_ESTG) {
        *trp = tropmodel(time, pos, azel, REL_HUMI);
        *var = SQR(ERR_SAAS / (sin(azel[1]) + 0.1));
        return 1;
    }
    /* SBAS (MOPS) troposphere model */
  /*  if (tropopt == TROPOPT_SBAS) {
        *trp = sbstropcorr(time, pos, azel, var);
        return 1;
    }*/
    /* no correction */
    *trp = 0.0;
    *var = tropopt == TROPOPT_OFF ? SQR(ERR_TROP) : 0.0;
    return 1;
}

/* psendorange with code bias correction -------------------------------------*/
static double prange(const obsd_t* obs, const nav_t* nav, const prcopt_t* opt,
    double* var)
{
    double P1, P2, gamma, b1, b2;
    int sat, sys;

    sat = obs->sat;
    sys = satsys(sat, NULL);
    P1 = obs->P[0];
    P2 = obs->P[1];
    *var = 0.0;

    if (P1 == 0.0 || (opt->ionoopt == IONOOPT_IFLC && P2 == 0.0)) return 0.0;

    /* P1-C1,P2-C2 DCB correction */
    if (sys == SYS_GPS || sys == SYS_GLO) {
        if (obs->code[0] == CODE_L1C) P1 += nav->cbias[sat - 1][1]; /* C1->P1 */
        if (obs->code[1] == CODE_L2C) P2 += nav->cbias[sat - 1][2]; /* C2->P2 */
    }
    if (opt->ionoopt == IONOOPT_IFLC) { /* dual-frequency */

        if (sys == SYS_GPS || sys == SYS_QZS) { /* L1-L2,G1-G2 */
            gamma = SQR(FREQ1 / FREQ2);
            return (P2 - gamma * P1) / (1.0 - gamma);
        }
        else if (sys == SYS_GLO) { /* G1-G2 */
            gamma = SQR(FREQ1_GLO / FREQ2_GLO);
            return (P2 - gamma * P1) / (1.0 - gamma);
        }
        else if (sys == SYS_GAL) { /* E1-E5b */
            gamma = SQR(FREQ1 / FREQ7);
            if (getseleph(SYS_GAL)) { /* F/NAV */
                P2 -= gettgd(sat, nav, 0) - gettgd(sat, nav, 1); /* BGD_E5aE5b */
            }
            return (P2 - gamma * P1) / (1.0 - gamma);
        }
        else if (sys == SYS_CMP) { /* B1-B2 */
            gamma = SQR(((obs->code[0] == CODE_L2I) ? FREQ1_CMP : FREQ1) / FREQ2_CMP);
            if (obs->code[0] == CODE_L2I) b1 = gettgd(sat, nav, 0); /* TGD_B1I */
            else if (obs->code[0] == CODE_L1P) b1 = gettgd(sat, nav, 2); /* TGD_B1Cp */
            else b1 = gettgd(sat, nav, 2) + gettgd(sat, nav, 4); /* TGD_B1Cp+ISC_B1Cd */
            b2 = gettgd(sat, nav, 1); /* TGD_B2I/B2bI (m) */
            return ((P2 - gamma * P1) - (b2 - gamma * b1)) / (1.0 - gamma);
        }
        else if (sys == SYS_IRN) { /* L5-S */
            gamma = SQR(FREQ5 / FREQ9);
            return (P2 - gamma * P1) / (1.0 - gamma);
        }
    }
    else { /* single-freq (L1/E1/B1) */
        *var = SQR(ERR_CBIAS);

        if (sys == SYS_GPS || sys == SYS_QZS) { /* L1 */
            b1 = gettgd(sat, nav, 0); /* TGD (m) */
            return P1 - b1;
        }
        else if (sys == SYS_GLO) { /* G1 */
            gamma = SQR(FREQ1_GLO / FREQ2_GLO);
            b1 = gettgd(sat, nav, 0); /* -dtaun (m) */
            return P1 - b1 / (gamma - 1.0);
        }
        else if (sys == SYS_GAL) { /* E1 */
            if (getseleph(SYS_GAL)) b1 = gettgd(sat, nav, 0); /* BGD_E1E5a */
            else                    b1 = gettgd(sat, nav, 1); /* BGD_E1E5b */
            return P1 - b1;
        }
        else if (sys == SYS_CMP) { /* B1I/B1Cp/B1Cd */
            if (obs->code[0] == CODE_L2I) b1 = gettgd(sat, nav, 0); /* TGD_B1I */
            else if (obs->code[0] == CODE_L1P) b1 = gettgd(sat, nav, 2); /* TGD_B1Cp */
            else b1 = gettgd(sat, nav, 2) + gettgd(sat, nav, 4); /* TGD_B1Cp+ISC_B1Cd */
            return P1 - b1;
        }
        else if (sys == SYS_IRN) { /* L5 */
            gamma = SQR(FREQ9 / FREQ5);
            b1 = gettgd(sat, nav, 0); /* TGD (m) */
            return P1 - gamma * b1;
        }
    }
    return P1;
}

/* pseudorange residuals -----------------------------------------------------*/
static int rescode(int iter, const obsd_t* obs, int n, const double* rs,
    const double* dts, const double* vare, const int* svh,
    const nav_t* nav, const double* x, const prcopt_t* opt,
    double* v, double* H, double* var, double* azel, int* vsat,
    double* resp, int* ns)//vsat：卫星是否可用
{
    gtime_t time;
    double r, freq, dion = 0.0, dtrp = 0.0, vmeas, vion = 0.0, vtrp = 0.0, rr[3], pos[3], dtr, e[3], P;
    int i, j, nv = 0, sat, sys, mask[NX - 3] = { 0 };

    trace(3, "resprng : n=%d\n", n);

    for (i = 0; i < 3; i++) rr[i] = x[i];//把卫星状态量X赋值给rr
    dtr = x[3];

    ecef2pos(rr, pos);//xyztoBLH

    for (i = *ns = 0; i < n && i < MAXOBS; i++) {//遍历循环所有卫星
        vsat[i] = 0; azel[i * 2] = azel[1 + i * 2] = resp[i] = 0.0;//高度角方位角2*n个
        time = obs[i].time;//接收机钟面时
        sat = obs[i].sat;
        if (!(sys = satsys(sat, NULL))) continue;//检测卫星系统

        /* reject duplicated observation data */
        if (i < n - 1 && i < MAXOBS - 1 && sat == obs[i + 1].sat) {
            trace(2, "duplicated obs data %s sat=%d\n", time_str(time, 3), sat);
            i++;
            continue;
        }
        /* excluded satellite? */
        if (satexclude(sat, vare[i], svh[i], opt)) continue;//删除卫星

        /* geometric distance */
        if ((r = geodist(rs + i * 6, rr, e)) <= 0.0) continue;//计算几何距离

        if (iter > 0) {
            /* test elevation mask */
            if (satazel(pos, e, azel + i * 2) < opt->elmin) continue;//方位角高度角信噪比

            /* test SNR mask *///信噪比
            if (!snrmask(obs + i, azel + i * 2, opt)) continue;//信噪比限制

            /* ionospheric correction */
            if (!ionocorr(time, nav, sat, pos, azel + i * 2, opt->ionoopt, &dion, &vion)) {
                continue;
            }
            if ((freq = sat2freq(sat, obs[i].code[0], nav)) == 0.0) continue;
            dion *= SQR(FREQ1 / freq);//电离层改正
            vion *= SQR(FREQ1 / freq);

            /* tropospheric correction */
            if (!tropcorr(time, nav, pos, azel + i * 2, opt->tropopt, &dtrp, &vtrp)) {
                continue;
            }
        }
        /* psendorange with code bias correction *///具有代码偏差校正功能的 Psendorange
        if ((P = prange(obs + i, nav, opt, &vmeas)) == 0.0) continue;//obs+i第i个卫星观测值

        /* pseudorange residual */
        v[nv] = P - (r + dtr - CLIGHT * dts[i * 2] + dion + dtrp);//计算残差

        /* design matrix */
        for (j = 0; j < NX; j++) {
            H[j + nv * NX] = j < 3 ? -e[j] : (j == 3 ? 1.0 : 0.0);//分成小于三等于三大于三三种情况
        }//NX单点定位中状态量的个数，4+4（xyz+gps钟，其它系统相对于GPS的偏差）
        /* time system offset and receiver bias correction */
        if (sys == SYS_GLO) { v[nv] -= x[4]; H[4 + nv * NX] = 1.0; mask[1] = 1; }//接收机钟系统间的偏差
        else if (sys == SYS_GAL) { v[nv] -= x[5]; H[5 + nv * NX] = 1.0; mask[2] = 1; }
        else if (sys == SYS_CMP) { v[nv] -= x[6]; H[6 + nv * NX] = 1.0; mask[3] = 1; }
        else if (sys == SYS_IRN) { v[nv] -= x[7]; H[7 + nv * NX] = 1.0; mask[4] = 1; }
#if 0 /* enable QZS-GPS time offset estimation */
        else if (sys == SYS_QZS) { v[nv] -= x[8]; H[8 + nv * NX] = 1.0; mask[5] = 1; }
#endif
        else mask[0] = 1;

        vsat[i] = 1; resp[i] = v[nv]; (*ns)++;

        /* variance of pseudorange error */
        var[nv++] = varerr(opt, azel[1 + i * 2], sys) + vare[i] + vmeas + vion + vtrp;//伪距方差计算
        //观测值的方差+卫星位置钟方差+码偏差方差+电离层对流层方差+
        trace(4, "sat=%2d azel=%5.1f %4.1f res=%7.3f sig=%5.3f\n", obs[i].sat,
            azel[i * 2] * R2D, azel[1 + i * 2] * R2D, resp[i], sqrt(var[nv - 1]));
    }
    /* constraint to avoid rank-deficient *///加约束防止秩亏
    for (i = 0; i < NX - 3; i++) {
        if (mask[i]) continue;
        v[nv] = 0.0;
        for (j = 0; j < NX; j++) H[j + nv * NX] = j == i + 3 ? 1.0 : 0.0;
        var[nv++] = 0.01;//假定观测值
    }
    return nv;//
}

/* validate solution ---------------------------------------------------------*/
static int valsol(const double* azel, const int* vsat, int n,
    const prcopt_t* opt, const double* v, int nv, int nx,
    char* msg)
{
    double azels[MAXOBS * 2], dop[4], vv;
    int i, ns;

    trace(3, "valsol  : n=%d nv=%d\n", n, nv);

    /* Chi-square validation of residuals *///卡方分布
    vv = dot(v, v, nv);
    if (nv > nx && vv > chisqr[nv - nx - 1]) {
        sprintf(msg, "chi-square error nv=%d vv=%.1f cs=%.1f", nv, vv, chisqr[nv - nx - 1]);
        return 0;
    }
    /* large GDOP check */
    for (i = ns = 0; i < n; i++) {
        if (!vsat[i]) continue;
        azels[ns * 2] = azel[i * 2];
        azels[1 + ns * 2] = azel[1 + i * 2];
        ns++;
    }
    dops(ns, azels, opt->elmin, dop);//DOP是“Dilution of Precision”“精度衰减因子”
    if (dop[0] <= 0.0 || dop[0] > opt->maxgdop) {
        sprintf(msg, "gdop error nv=%d gdop=%.1f", nv, dop[0]);
        return 0;
    }
    return 1;
}

/* range rate residuals ------------------------------------------------------*/
static int resdop(const obsd_t* obs, int n, const double* rs, const double* dts,
    const nav_t* nav, const double* rr, const double* x,
    const double* azel, const int* vsat, double err, double* v,
    double* H)
{
    double freq, rate, pos[3], E[9], a[3], e[3], vs[3], cosel, sig;
    int i, j, nv = 0;

    trace(3, "resdop  : n=%d\n", n);

    ecef2pos(rr, pos); xyz2enu(pos, E);

    for (i = 0; i < n && i < MAXOBS; i++) {

        freq = sat2freq(obs[i].sat, obs[i].code[0], nav);

        if (obs[i].D[0] == 0.0 || freq == 0.0 || !vsat[i] || norm(rs + 3 + i * 6, 3) <= 0.0) {
            continue;
        }
        /* LOS (line-of-sight) vector in ECEF */
        cosel = cos(azel[1 + i * 2]);
        a[0] = sin(azel[i * 2]) * cosel;
        a[1] = cos(azel[i * 2]) * cosel;
        a[2] = sin(azel[1 + i * 2]);
        matmul("TN", 3, 1, 3, 1.0, E, a, 0.0, e);

        /* satellite velocity relative to receiver in ECEF */
        for (j = 0; j < 3; j++) {
            vs[j] = rs[j + 3 + i * 6] - x[j];
        }
        /* range rate with earth rotation correction */
        rate = dot(vs, e, 3) + OMGE / CLIGHT * (rs[4 + i * 6] * rr[0] + rs[1 + i * 6] * x[0] -
            rs[3 + i * 6] * rr[1] - rs[i * 6] * x[1]);

        /* Std of range rate error (m/s) */
        sig = (err <= 0.0) ? 1.0 : err * CLIGHT / freq;

        /* range rate residual (m/s) */
        v[nv] = (-obs[i].D[0] * CLIGHT / freq - (rate + x[3] - CLIGHT * dts[1 + i * 2])) / sig;

        /* design matrix */
        for (j = 0; j < 4; j++) {
            H[j + nv * 4] = ((j < 3) ? -e[j] : 1.0) / sig;
        }
        nv++;
    }
    return nv;
}

/* estimate receiver position ------------------------------------------------*/
static int estpos(const obsd_t* obs, int n, const double* rs, const double* dts,
    const double* vare, const int* svh, const nav_t* nav,
    const prcopt_t* opt, sol_t* sol, double* azel, int* vsat,
    double* resp, char* msg)
{
    double x[NX] = { 0 }, dx[NX], Q[NX * NX], * v, * H, * var, sig;//x:xyz加gps接收机钟4个值
    int i, j, k, info, stat, nv, ns;

    trace(3, "estpos  : n=%d\n", n);

    v = mat(n + 4, 1); H = mat(NX, n + 4); var = mat(n + 4, 1);//残差+方向余弦矩阵

    for (i = 0; i < 3; i++) x[i] = sol->rr[i];

    for (i = 0; i < MAXITR; i++) {//迭代

        /* pseudorange residuals (m) *///计算所有卫星先验残差
        nv = rescode(i, obs, n, rs, dts, vare, svh, nav, x, opt, v, H, var, azel, vsat, resp,
            &ns);//伪距残差
        trace(4, "spp H=\n");
        tracemat(4, H, NX, n + 4, 8, 3);//打印等级，矩阵，行，列，输出占位个数和精度
        trace(4, "spp v=\n");
        tracemat(4, v, 1, nv, 8, 3);
        if (nv < NX) {//如果残差小于状态量个数
            sprintf(msg, "lack of valid sats ns=%d", nv);
            break;//跳出
        }
        /* weighted by Std *///"std"标准差（Standard Deviation）
        for (j = 0; j < nv; j++) {
            sig = sqrt(var[j]);//方差开根号
            v[j] /= sig;//vj索引的值/sig再赋值给v
            for (k = 0; k < NX; k++) H[k + j * NX] /= sig;//将方向矩阵除std
        }//相当于加权
        /* least square estimation */
        if ((info = lsq(H, v, NX, nv, dx, Q))) {
            sprintf(msg, "lsq error info=%d", info);
            break;
        }
        for (j = 0; j < NX; j++) {
            x[j] += dx[j];//将dx的值赋给x
        }
        if (norm(dx, NX) < 1E-4) {
            sol->type = 0;
            sol->time = timeadd(obs[0].time, -x[3] / CLIGHT);//
            sol->dtr[0] = x[3] / CLIGHT; /* receiver clock bias (s) *///把钟偏放到dtr中
            sol->dtr[1] = x[4] / CLIGHT; /* GLO-GPS time offset (s) */
            sol->dtr[2] = x[5] / CLIGHT; /* GAL-GPS time offset (s) */
            sol->dtr[3] = x[6] / CLIGHT; /* BDS-GPS time offset (s) */
            sol->dtr[4] = x[7] / CLIGHT; /* IRN-GPS time offset (s) */
            for (j = 0; j < 6; j++) sol->rr[j] = j < 3 ? x[j] : 0.0;//rr：位置和速度
            for (j = 0; j < 3; j++) sol->qr[j] = (float)Q[j + j * NX];//存储方差
            sol->qr[3] = (float)Q[1];    /* cov xy */
            sol->qr[4] = (float)Q[2 + NX]; /* cov yz *///协方差
            sol->qr[5] = (float)Q[2];    /* cov zx */
            sol->ns = (uint8_t)ns;//卫星数
            sol->age = sol->ratio = 0.0;

            /* validate solution */
            if ((stat = valsol(azel, vsat, n, opt, v, nv, NX, msg))) {//检查定位结果是否正确
                sol->stat = opt->sateph == EPHOPT_SBAS ? SOLQ_SBAS : SOLQ_SINGLE;//if (opt->sateph==EPHOPT_SBAS?) sol->stat=SOLQ_SBAS:SOLQ_SINGLE
            }
            free(v); free(H); free(var);
            return stat;
        }
    }
    if (i >= MAXITR) sprintf(msg, "iteration divergent i=%d", i);

    free(v); free(H); free(var);
    return 0;
}

/* estimate receiver velocity ------------------------------------------------*/
static void estvel(const obsd_t* obs, int n, const double* rs, const double* dts,
    const nav_t* nav, const prcopt_t* opt, sol_t* sol,
    const double* azel, const int* vsat)
{
    double x[4] = { 0 }, dx[4], Q[16], * v, * H;
    double err = opt->err[4]; /* Doppler error (Hz) */
    int i, j, nv;

    trace(3, "estvel  : n=%d\n", n);

    v = mat(n, 1); H = mat(4, n);

    for (i = 0; i < MAXITR; i++) {

        /* range rate residuals (m/s) */
        if ((nv = resdop(obs, n, rs, dts, nav, sol->rr, x, azel, vsat, err, v, H)) < 4) {
            break;
        }
        /* least square estimation */
        if (lsq(H, v, 4, nv, dx, Q)) break;

        for (j = 0; j < 4; j++) x[j] += dx[j];

        if (norm(dx, 4) < 1E-6) {
            matcpy(sol->rr + 3, x, 3, 1);
            sol->qv[0] = (float)Q[0];  /* xx */
            sol->qv[1] = (float)Q[5];  /* yy */
            sol->qv[2] = (float)Q[10]; /* zz */
            sol->qv[3] = (float)Q[1];  /* xy */
            sol->qv[4] = (float)Q[6];  /* yz */
            sol->qv[5] = (float)Q[2];  /* zx */
            break;
        }
    }
    free(v); free(H);
}

/* RAIM FDE (failure detection and exclution) -------------------------------*/
static int raim_fde(const obsd_t* obs, int n, const double* rs,
    const double* dts, const double* vare, const int* svh,
    const nav_t* nav, const prcopt_t* opt, sol_t* sol,
    double* azel, int* vsat, double* resp, char* msg)
{
    obsd_t* obs_e;
    sol_t sol_e = { {0} };
    char tstr[32], name[16], msg_e[128];
    double* rs_e, * dts_e, * vare_e, * azel_e, * resp_e, rms_e, rms = 100.0;
    int i, j, k, nvsat, stat = 0, * svh_e, * vsat_e, sat = 0;

    trace(3, "raim_fde: %s n=%2d\n", time_str(obs[0].time, 0), n);

    if (!(obs_e = (obsd_t*)malloc(sizeof(obsd_t) * n))) return 0;
    rs_e = mat(6, n); dts_e = mat(2, n); vare_e = mat(1, n); azel_e = zeros(2, n);
    svh_e = imat(1, n); vsat_e = imat(1, n); resp_e = mat(1, n);

    for (i = 0; i < n; i++) {

        /* satellite exclution */
        for (j = k = 0; j < n; j++) {
            if (j == i) continue;
            obs_e[k] = obs[j];
            matcpy(rs_e + 6 * k, rs + 6 * j, 6, 1);
            matcpy(dts_e + 2 * k, dts + 2 * j, 2, 1);
            vare_e[k] = vare[j];
            svh_e[k++] = svh[j];
        }
        /* estimate receiver position without a satellite */
        if (!estpos(obs_e, n - 1, rs_e, dts_e, vare_e, svh_e, nav, opt, &sol_e, azel_e,
            vsat_e, resp_e, msg_e)) {
            trace(3, "raim_fde: exsat=%2d (%s)\n", obs[i].sat, msg);
            continue;
        }//结果合理则不进这里，进行下一步
        for (j = nvsat = 0, rms_e = 0.0; j < n - 1; j++) {
            if (!vsat_e[j]) continue;
            rms_e += SQR(resp_e[j]);
            nvsat++;
        }
        if (nvsat < 5) {
            trace(3, "raim_fde: exsat=%2d lack of satellites nvsat=%2d\n",
                obs[i].sat, nvsat);
            continue;
        }
        rms_e = sqrt(rms_e / nvsat);

        trace(3, "raim_fde: exsat=%2d rms=%8.3f\n", obs[i].sat, rms_e);

        if (rms_e > rms) continue;

        /* save result */
        for (j = k = 0; j < n; j++) {
            if (j == i) continue;
            matcpy(azel + 2 * j, azel_e + 2 * k, 2, 1);
            vsat[j] = vsat_e[k];
            resp[j] = resp_e[k++];
        }
        stat = 1;
        *sol = sol_e;
        sat = obs[i].sat;
        rms = rms_e;
        vsat[i] = 0;
        strcpy(msg, msg_e);
    }
    if (stat) {
        time2str(obs[0].time, tstr, 2); satno2id(sat, name);
        trace(2, "%s: %s excluded by raim\n", tstr + 11, name);
    }
    free(obs_e);
    free(rs_e); free(dts_e); free(vare_e); free(azel_e);
    free(svh_e); free(vsat_e); free(resp_e);
    return stat;
}

extern int pntpos(const obsd_t* obs, int n, const nav_t* nav,
    const prcopt_t* opt, sol_t* sol, double* azel, ssat_t* ssat,
    char* msg)
{
    prcopt_t opt_ = *opt;//复制了下参数
    double* rs, * dts, * var, * azel_, * resp;
    int i, stat, vsat[MAXOBS] = { 0 }, svh[MAXOBS];

    trace(3, "pntpos  : tobs=%s n=%d\n", time_str(obs[0].time, 3), n);

    sol->stat = SOLQ_NONE;//没有结算结果

    if (n <= 0) {
        strcpy(msg, "no observation data");
        return 0;
    }
    sol->time = obs[0].time;
    msg[0] = '\0';

    rs = mat(6, n); dts = mat(2, n); var = mat(1, n); azel_ = zeros(2, n); resp = mat(1, n);//rs:xyz,vxvyvz.

    if (opt_.mode != PMODE_SINGLE) { /* for precise positioning */
        opt_.ionoopt = IONOOPT_BRDC;
        opt_.tropopt = TROPOPT_SAAS;
    }
    /* satellite positons, velocities and clocks *///卫星位置的、速度钟
    satposs(sol->time, obs, n, nav, opt_.sateph, rs, dts, var, svh);

    /* estimate receiver position with pseudorange */
    stat = estpos(obs, n, rs, dts, var, svh, nav, &opt_, sol, azel_, vsat, resp, msg);//dts:钟、钟漂

    /* RAIM FDE *///“Receiver Autonomous Integrity Monitoring”“接收机自主完整性监测”
    if (!stat && n >= 6 && opt->posopt[4]) {
        stat = raim_fde(obs, n, rs, dts, var, svh, nav, &opt_, sol, azel_, vsat, resp, msg);
    }
    /* estimate receiver velocity with Doppler */
    if (stat) {
        estvel(obs, n, rs, dts, nav, &opt_, sol, azel_, vsat);
    }
    if (azel) {
        for (i = 0; i < n * 2; i++) azel[i] = azel_[i];
    }
    if (ssat) {
        for (i = 0; i < MAXSAT; i++) {
            ssat[i].vs = 0;
            ssat[i].azel[0] = ssat[i].azel[1] = 0.0;
            ssat[i].resp[0] = ssat[i].resc[0] = 0.0;
            ssat[i].snr[0] = 0;
        }
        for (i = 0; i < n; i++) {
            ssat[obs[i].sat - 1].azel[0] = azel_[i * 2];
            ssat[obs[i].sat - 1].azel[1] = azel_[1 + i * 2];
            ssat[obs[i].sat - 1].snr[0] = obs[i].SNR[0];
            if (!vsat[i]) continue;
            ssat[obs[i].sat - 1].vs = 1;
            ssat[obs[i].sat - 1].resp[0] = resp[i];
        }
    }
    free(rs); free(dts); free(var); free(azel_); free(resp);
    return stat;
}


////////////////////////////////////////
///////////////////////////////////////
////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
//////////////////////////////// ephememeris内容
/* select ephememeris --------------------------------------------------------*/
//
//static eph_t* seleph(gtime_t time, int sat, int iode, const nav_t* nav)
//{
//    double t, tmax, tmin;
//    int i, j = -1, sys, sel;
//
//    trace(4, "seleph  : time=%s sat=%2d iode=%d\n", time_str(time, 3), sat, iode);
//
//    sys = satsys(sat, NULL);
//    switch (sys) {
//    case SYS_GPS: tmax = MAXDTOE + 1.0; sel = eph_sel[0]; break;
//    case SYS_GAL: tmax = MAXDTOE_GAL; sel = eph_sel[2]; break;
//    case SYS_QZS: tmax = MAXDTOE_QZS + 1.0; sel = eph_sel[3]; break;
//    case SYS_CMP: tmax = MAXDTOE_CMP + 1.0; sel = eph_sel[4]; break;
//    case SYS_IRN: tmax = MAXDTOE_IRN + 1.0; sel = eph_sel[5]; break;
//    default: tmax = MAXDTOE + 1.0; break;
//    }
//    tmin = tmax + 1.0;
//
//    for (i = 0; i < nav->n; i++) {
//        if (nav->eph[i].sat != sat) continue;
//        if (iode >= 0 && nav->eph[i].iode != iode) continue;
//        if (sys == SYS_GAL) {
//            sel = getseleph(SYS_GAL);
//            if (sel == 0 && !(nav->eph[i].code & (1 << 9))) continue; /* I/NAV */
//            if (sel == 1 && !(nav->eph[i].code & (1 << 8))) continue; /* F/NAV */
//            if (timediff(nav->eph[i].toe, time) >= 0.0) continue; /* AOD<=0 */
//        }
//        if ((t = fabs(timediff(nav->eph[i].toe, time))) > tmax) continue;
//        if (iode >= 0) return nav->eph + i;
//        if (t <= tmin) { j = i; tmin = t; } /* toe closest to time */
//    }
//    if (iode >= 0 || j < 0) {
//        trace(3, "no broadcast ephemeris: %s sat=%2d iode=%3d\n",
//            time_str(time, 0), sat, iode);
//        return NULL;
//    }
//    return nav->eph + j;
//}
//
///* select glonass ephememeris ------------------------------------------------*/
//static geph_t* selgeph(gtime_t time, int sat, int iode, const nav_t* nav)
//{
//    double t, tmax = MAXDTOE_GLO, tmin = tmax + 1.0;
//    int i, j = -1;
//
//    trace(4, "selgeph : time=%s sat=%2d iode=%2d\n", time_str(time, 3), sat, iode);
//
//    for (i = 0; i < nav->ng; i++) {
//        if (nav->geph[i].sat != sat) continue;
//        if (iode >= 0 && nav->geph[i].iode != iode) continue;
//        if ((t = fabs(timediff(nav->geph[i].toe, time))) > tmax) continue;
//        if (iode >= 0) return nav->geph + i;
//        if (t <= tmin) { j = i; tmin = t; } /* toe closest to time */
//    }
//    if (iode >= 0 || j < 0) {
//        trace(3, "no glonass ephemeris  : %s sat=%2d iode=%2d\n", time_str(time, 0),
//            sat, iode);
//        return NULL;
//    }
//    return nav->geph + j;
//}
//
///* select sbas ephememeris ---------------------------------------------------*/
//static seph_t* selseph(gtime_t time, int sat, const nav_t* nav)
//{
//    double t, tmax = MAXDTOE_SBS, tmin = tmax + 1.0;
//    int i, j = -1;
//
//    trace(4, "selseph : time=%s sat=%2d\n", time_str(time, 3), sat);
//
//    for (i = 0; i < nav->ns; i++) {
//        if (nav->seph[i].sat != sat) continue;
//        if ((t = fabs(timediff(nav->seph[i].t0, time))) > tmax) continue;
//        if (t <= tmin) { j = i; tmin = t; } /* toe closest to time */
//    }
//    if (j < 0) {
//        trace(3, "no sbas ephemeris     : %s sat=%2d\n", time_str(time, 0), sat);
//        return NULL;
//    }
//    return nav->seph + j;
//}
//
//extern double eph2clk(gtime_t time, const eph_t* eph)
//{
//    double t, ts;
//    int i;
//
//    trace(4, "eph2clk : time=%s sat=%2d\n", time_str(time, 3), eph->sat);
//
//    t = ts = timediff(time, eph->toc);
//
//    for (i = 0; i < 2; i++) {
//        t = ts - (eph->f0 + eph->f1 * t + eph->f2 * t * t);
//    }
//    return eph->f0 + eph->f1 * t + eph->f2 * t * t;
//}
//
///* glonass ephemeris to satellite clock bias -----------------------------------
//* compute satellite clock bias with glonass ephemeris
//* args   : gtime_t time     I   time by satellite clock (gpst)
//*          geph_t *geph     I   glonass ephemeris
//* return : satellite clock bias (s)
//* notes  : see ref [2]
//*-----------------------------------------------------------------------------*/
//extern double geph2clk(gtime_t time, const geph_t* geph)
//{
//    double t, ts;
//    int i;
//
//    trace(4, "geph2clk: time=%s sat=%2d\n", time_str(time, 3), geph->sat);
//
//    t = ts = timediff(time, geph->toe);
//
//    for (i = 0; i < 2; i++) {
//        t = ts - (-geph->taun + geph->gamn * t);
//    }
//    return -geph->taun + geph->gamn * t;
//}
//
///* sbas ephemeris to satellite clock bias --------------------------------------
//* compute satellite clock bias with sbas ephemeris
//* args   : gtime_t time     I   time by satellite clock (gpst)
//*          seph_t *seph     I   sbas ephemeris
//* return : satellite clock bias (s)
//* notes  : see ref [3]
//*-----------------------------------------------------------------------------*/
//extern double seph2clk(gtime_t time, const seph_t* seph)
//{
//    double t;
//    int i;
//
//    trace(4, "seph2clk: time=%s sat=%2d\n", time_str(time, 3), seph->sat);
//
//    t = timediff(time, seph->t0);
//
//    for (i = 0; i < 2; i++) {
//        t -= seph->af0 + seph->af1 * t;
//    }
//    return seph->af0 + seph->af1 * t;
//}
//
///* satellite clock with broadcast ephemeris ----------------------------------*/
//static int ephclk(gtime_t time, gtime_t teph, int sat, const nav_t* nav,
//    double* dts)
//{
//    eph_t* eph;
//    geph_t* geph;
//    seph_t* seph;
//    int sys;
//
//    trace(4, "ephclk  : time=%s sat=%2d\n", time_str(time, 3), sat);
//
//    sys = satsys(sat, NULL);
//
//    if (sys == SYS_GPS || sys == SYS_GAL || sys == SYS_QZS || sys == SYS_CMP || sys == SYS_IRN) {
//        if (!(eph = seleph(teph, sat, -1, nav))) return 0;
//        *dts = eph2clk(time, eph);
//    }
//    else if (sys == SYS_GLO) {
//        if (!(geph = selgeph(teph, sat, -1, nav))) return 0;
//        *dts = geph2clk(time, geph);
//    }
//    else if (sys == SYS_SBS) {
//        if (!(seph = selseph(teph, sat, nav))) return 0;
//        *dts = seph2clk(time, seph);
//    }
//    else return 0;
//
//    return 1;
//}
//
//extern void eph2pos(gtime_t time, const eph_t* eph, double* rs, double* dts,
//    double* var)
//{
//    double tk, M, E, Ek, sinE, cosE, u, r, i, O, sin2u, cos2u, x, y, sinO, cosO, cosi, mu, omge;
//    double xg, yg, zg, sino, coso;
//    int n, sys, prn;
//
//    trace(4, "eph2pos : time=%s sat=%2d\n", time_str(time, 3), eph->sat);
//
//    if (eph->A <= 0.0) {
//        rs[0] = rs[1] = rs[2] = *dts = *var = 0.0;
//        return;
//    }
//    tk = timediff(time, eph->toe);
//
//    switch ((sys = satsys(eph->sat, &prn))) {
//    case SYS_GAL: mu = MU_GAL; omge = OMGE_GAL; break;
//    case SYS_CMP: mu = MU_CMP; omge = OMGE_CMP; break;
//    default:      mu = MU_GPS; omge = OMGE;     break;
//    }
//    M = eph->M0 + (sqrt(mu / (eph->A * eph->A * eph->A)) + eph->deln) * tk;
//
//    for (n = 0, E = M, Ek = 0.0; fabs(E - Ek) > RTOL_KEPLER && n < MAX_ITER_KEPLER; n++) {
//        Ek = E; E -= (E - eph->e * sin(E) - M) / (1.0 - eph->e * cos(E));
//    }
//    if (n >= MAX_ITER_KEPLER) {
//        trace(2, "eph2pos: kepler iteration overflow sat=%2d\n", eph->sat);
//        return;
//    }
//    sinE = sin(E); cosE = cos(E);
//
//    trace(4, "kepler: sat=%2d e=%8.5f n=%2d del=%10.3e\n", eph->sat, eph->e, n, E - Ek);
//
//    u = atan2(sqrt(1.0 - eph->e * eph->e) * sinE, cosE - eph->e) + eph->omg;
//    r = eph->A * (1.0 - eph->e * cosE);
//    i = eph->i0 + eph->idot * tk;
//    sin2u = sin(2.0 * u); cos2u = cos(2.0 * u);
//    u += eph->cus * sin2u + eph->cuc * cos2u;
//    r += eph->crs * sin2u + eph->crc * cos2u;
//    i += eph->cis * sin2u + eph->cic * cos2u;
//    x = r * cos(u); y = r * sin(u); cosi = cos(i);
//
//    /* beidou geo satellite */
//    if (sys == SYS_CMP && (prn <= 5 || prn >= 59)) { /* ref [9] table 4-1 */
//        O = eph->OMG0 + eph->OMGd * tk - omge * eph->toes;
//        sinO = sin(O); cosO = cos(O);
//        xg = x * cosO - y * cosi * sinO;
//        yg = x * sinO + y * cosi * cosO;
//        zg = y * sin(i);
//        sino = sin(omge * tk); coso = cos(omge * tk);
//        rs[0] = xg * coso + yg * sino * COS_5 + zg * sino * SIN_5;
//        rs[1] = -xg * sino + yg * coso * COS_5 + zg * coso * SIN_5;
//        rs[2] = -yg * SIN_5 + zg * COS_5;
//    }
//    else {
//        O = eph->OMG0 + (eph->OMGd - omge) * tk - omge * eph->toes;
//        sinO = sin(O); cosO = cos(O);
//        rs[0] = x * cosO - y * cosi * sinO;
//        rs[1] = x * sinO + y * cosi * cosO;
//        rs[2] = y * sin(i);
//    }
//    tk = timediff(time, eph->toc);
//    *dts = eph->f0 + eph->f1 * tk + eph->f2 * tk * tk;
//
//    /* relativity correction */
//    *dts -= 2.0 * sqrt(mu * eph->A) * eph->e * sinE / SQR(CLIGHT);
//
//    /* position and clock error variance */
//    *var = var_uraeph(sys, eph->sva);
//}
//
///* satellite position and clock by broadcast ephemeris -----------------------*/
//static int ephpos(gtime_t time, gtime_t teph, int sat, const nav_t* nav,
//    int iode, double* rs, double* dts, double* var, int* svh)
//{
//    eph_t* eph;
//    geph_t* geph;
//    seph_t* seph;
//    double rst[3], dtst[1], tt = 1E-3;
//    int i, sys;
//
//    trace(4, "ephpos  : time=%s sat=%2d iode=%d\n", time_str(time, 3), sat, iode);
//
//    sys = satsys(sat, NULL);
//
//    *svh = -1;
//
//    if (sys == SYS_GPS || sys == SYS_GAL || sys == SYS_QZS || sys == SYS_CMP || sys == SYS_IRN) {
//        if (!(eph = seleph(teph, sat, iode, nav))) return 0;
//        eph2pos(time, eph, rs, dts, var);//�������λ��
//        time = timeadd(time, tt);
//        eph2pos(time, eph, rst, dtst, var);//λ��/tt���ٶ�
//        *svh = eph->svh;
//    }
//    else if (sys == SYS_GLO) {
//        if (!(geph = selgeph(teph, sat, iode, nav))) return 0;
//        geph2pos(time, geph, rs, dts, var);
//        time = timeadd(time, tt);
//        geph2pos(time, geph, rst, dtst, var);
//        *svh = geph->svh;
//    }
//    else if (sys == SYS_SBS) {
//        if (!(seph = selseph(teph, sat, nav))) return 0;
//        seph2pos(time, seph, rs, dts, var);
//        time = timeadd(time, tt);
//        seph2pos(time, seph, rst, dtst, var);
//        *svh = seph->svh;
//    }
//    else return 0;
//
//    /* satellite velocity and clock drift by differential approx */
//    for (i = 0; i < 3; i++) rs[i + 3] = (rst[i] - rs[i]) / tt;
//    dts[1] = (dtst[0] - dts[0]) / tt;
//
//    return 1;
//}
///* satellite position and clock with sbas correction -------------------------*/
//static int satpos_sbas(gtime_t time, gtime_t teph, int sat, const nav_t* nav,
//    double* rs, double* dts, double* var, int* svh)
//{
//    const sbssatp_t* sbs = NULL;
//    int i;
//
//    trace(4, "satpos_sbas: time=%s sat=%2d\n", time_str(time, 3), sat);
//
//    /* search sbas satellite correciton */
//    for (i = 0; i < nav->sbssat.nsat; i++) {
//        sbs = nav->sbssat.sat + i;
//        if (sbs->sat == sat) break;
//    }
//    if (i >= nav->sbssat.nsat) {
//        trace(2, "no sbas correction for orbit: %s sat=%2d\n", time_str(time, 0), sat);
//        ephpos(time, teph, sat, nav, -1, rs, dts, var, svh);
//        *svh = -1;
//        return 0;
//    }
//    /* satellite postion and clock by broadcast ephemeris */
//    if (!ephpos(time, teph, sat, nav, sbs->lcorr.iode, rs, dts, var, svh)) return 0;
//
//    /* sbas satellite correction (long term and fast) */
//    if (sbssatcorr(time, sat, nav, rs, dts, var)) return 1;
//    *svh = -1;
//    return 0;
//}
///* variance by ura ssr (ref [10] table 3.3-1 DF389) --------------------------*/
//static double var_urassr(int ura)
//{
//    double std;
//    if (ura <= 0) return SQR(DEFURASSR);
//    if (ura >= 63) return SQR(5.4665);
//    std = (pow(3.0, (ura >> 3) & 7) * (1.0 + (ura & 7) / 4.0) - 1.0) * 1E-3;
//    return SQR(std);
//}
//
///* variance by ura ephemeris -------------------------------------------------*/
//static double var_uraeph(int sys, int ura)
//{
//    const double ura_value[] = {
//        2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
//        3072.0,6144.0
//    };
//    if (sys == SYS_GAL) { /* galileo sisa (ref [7] 5.1.11) */
//        if (ura <= 49) return SQR(ura * 0.01);
//        if (ura <= 74) return SQR(0.5 + (ura - 50) * 0.02);
//        if (ura <= 99) return SQR(1.0 + (ura - 75) * 0.04);
//        if (ura <= 125) return SQR(2.0 + (ura - 100) * 0.16);
//        return SQR(STD_GAL_NAPA);
//    }
//    else { /* gps ura (ref [1] 20.3.3.3.1.1) */
//        return ura < 0 || 14 < ura ? SQR(6144.0) : SQR(ura_value[ura]);
//    }
//}
///* variance by ura ssr (ref [10] table 3.3-1 DF389) --------------------------*/
//static double var_urassr(int ura)
//{
//    double std;
//    if (ura <= 0) return SQR(DEFURASSR);
//    if (ura >= 63) return SQR(5.4665);
//    std = (pow(3.0, (ura >> 3) & 7) * (1.0 + (ura & 7) / 4.0) - 1.0) * 1E-3;
//    return SQR(std);
//}
//
///* satellite position and clock with ssr correction --------------------------*/
//static int satpos_ssr(gtime_t time, gtime_t teph, int sat, const nav_t* nav,
//    int opt, double* rs, double* dts, double* var, int* svh)
//{
//    const ssr_t* ssr;
//    eph_t* eph;
//    double t1, t2, t3, er[3], ea[3], ec[3], rc[3], deph[3], dclk, dant[3] = { 0 }, tk;
//    int i, sys;
//
//    trace(4, "satpos_ssr: time=%s sat=%2d\n", time_str(time, 3), sat);
//
//    ssr = nav->ssr + sat - 1;
//
//    if (!ssr->t0[0].time) {
//        trace(2, "no ssr orbit correction: %s sat=%2d\n", time_str(time, 0), sat);
//        return 0;
//    }
//    if (!ssr->t0[1].time) {
//        trace(2, "no ssr clock correction: %s sat=%2d\n", time_str(time, 0), sat);
//        return 0;
//    }
//    /* inconsistency between orbit and clock correction */
//    if (ssr->iod[0] != ssr->iod[1]) {
//        trace(2, "inconsist ssr correction: %s sat=%2d iod=%d %d\n",
//            time_str(time, 0), sat, ssr->iod[0], ssr->iod[1]);
//        *svh = -1;
//        return 0;
//    }
//    t1 = timediff(time, ssr->t0[0]);
//    t2 = timediff(time, ssr->t0[1]);
//    t3 = timediff(time, ssr->t0[2]);
//
//    /* ssr orbit and clock correction (ref [4]) */
//    if (fabs(t1) > MAXAGESSR || fabs(t2) > MAXAGESSR) {
//        trace(2, "age of ssr error: %s sat=%2d t=%.0f %.0f\n", time_str(time, 0),
//            sat, t1, t2);
//        *svh = -1;
//        return 0;
//    }
//    if (ssr->udi[0] >= 1.0) t1 -= ssr->udi[0] / 2.0;
//    if (ssr->udi[1] >= 1.0) t2 -= ssr->udi[1] / 2.0;
//
//    for (i = 0; i < 3; i++) deph[i] = ssr->deph[i] + ssr->ddeph[i] * t1;
//    dclk = ssr->dclk[0] + ssr->dclk[1] * t2 + ssr->dclk[2] * t2 * t2;
//
//    /* ssr highrate clock correction (ref [4]) */
//    if (ssr->iod[0] == ssr->iod[2] && ssr->t0[2].time && fabs(t3) < MAXAGESSR_HRCLK) {
//        dclk += ssr->hrclk;
//    }
//    if (norm(deph, 3) > MAXECORSSR || fabs(dclk) > MAXCCORSSR) {
//        trace(3, "invalid ssr correction: %s deph=%.1f dclk=%.1f\n",
//            time_str(time, 0), norm(deph, 3), dclk);
//        *svh = -1;
//        return 0;
//    }
//    /* satellite postion and clock by broadcast ephemeris */
//    if (!ephpos(time, teph, sat, nav, ssr->iode, rs, dts, var, svh)) return 0;
//
//    /* satellite clock for gps, galileo and qzss */
//    sys = satsys(sat, NULL);
//    if (sys == SYS_GPS || sys == SYS_GAL || sys == SYS_QZS || sys == SYS_CMP) {
//        if (!(eph = seleph(teph, sat, ssr->iode, nav))) return 0;
//
//        /* satellite clock by clock parameters */
//        tk = timediff(time, eph->toc);
//        dts[0] = eph->f0 + eph->f1 * tk + eph->f2 * tk * tk;
//        dts[1] = eph->f1 + 2.0 * eph->f2 * tk;
//
//        /* relativity correction */
//        dts[0] -= 2.0 * dot(rs, rs + 3, 3) / CLIGHT / CLIGHT;
//    }
//    /* radial-along-cross directions in ecef */
//    if (!normv3(rs + 3, ea)) return 0;
//    cross3(rs, rs + 3, rc);
//    if (!normv3(rc, ec)) {
//        *svh = -1;
//        return 0;
//    }
//    cross3(ea, ec, er);
//
//    /* satellite antenna offset correction */
//    if (opt) {
//        satantoff(time, rs, sat, nav, dant);
//    }
//    for (i = 0; i < 3; i++) {
//        rs[i] += -(er[i] * deph[0] + ea[i] * deph[1] + ec[i] * deph[2]) + dant[i];
//    }
//    /* t_corr = t_sv - (dts(brdc) + dclk(ssr) / CLIGHT) (ref [10] eq.3.12-7) */
//    dts[0] += dclk / CLIGHT;
//
//    /* variance by ssr ura */
//    *var = var_urassr(ssr->ura);
//
//    trace(5, "satpos_ssr: %s sat=%2d deph=%6.3f %6.3f %6.3f er=%6.3f %6.3f %6.3f dclk=%6.3f var=%6.3f\n",
//        time_str(time, 2), sat, deph[0], deph[1], deph[2], er[0], er[1], er[2], dclk, *var);
//
//    return 1;
//}
//
//
//
//
//
//
//
//
//
//extern int satpos(gtime_t time, gtime_t teph, int sat, int ephopt,
//    const nav_t* nav, double* rs, double* dts, double* var,
//    int* svh)
//{
//    trace(4, "satpos  : time=%s sat=%2d ephopt=%d\n", time_str(time, 3), sat, ephopt);
//
//    *svh = 0;
//
//    switch (ephopt) {
//    case EPHOPT_BRDC: return ephpos(time, teph, sat, nav, -1, rs, dts, var, svh);
//    case EPHOPT_SBAS: return satpos_sbas(time, teph, sat, nav, rs, dts, var, svh);
//    case EPHOPT_SSRAPC: return satpos_ssr(time, teph, sat, nav, 0, rs, dts, var, svh);
//    case EPHOPT_SSRCOM: return satpos_ssr(time, teph, sat, nav, 1, rs, dts, var, svh);
//    case EPHOPT_PREC:
//        if (!peph2pos(time, sat, nav, 1, rs, dts, var)) break; else return 1;
//    }
//    *svh = -1;
//    return 0;
//}
//
//extern void satposs(gtime_t teph, const obsd_t* obs, int n, const nav_t* nav,
//    int ephopt, double* rs, double* dts, double* var, int* svh)
//{
//    gtime_t time[2 * MAXOBS] = { {0} };
//    double dt, pr;
//    int i, j;
//
//    trace(3, "satposs : teph=%s n=%d ephopt=%d\n", time_str(teph, 3), n, ephopt);
//
//    for (i = 0; i < n && i < 2 * MAXOBS; i++) {
//        for (j = 0; j < 6; j++) rs[j + i * 6] = 0.0;
//        for (j = 0; j < 2; j++) dts[j + i * 2] = 0.0;
//        var[i] = 0.0; svh[i] = 0;//遍历卫星，初始化
//
//        /* search any pseudorange */
//        for (j = 0, pr = 0.0; j < NFREQ; j++) if ((pr = obs[i].P[j]) != 0.0) break;
//
//        if (j >= NFREQ) {
//            trace(3, "no pseudorange %s sat=%2d\n", time_str(obs[i].time, 3), obs[i].sat);
//            continue;
//        }
//        /* transmission time by satellite clock *///计算卫星1发射时刻
//        time[i] = timeadd(obs[i].time, -pr / CLIGHT);//先减一遍
//
//        /* satellite clock bias by broadcast ephemeris */
//        if (!ephclk(time[i], teph, obs[i].sat, nav, &dt)) {
//            trace(3, "no broadcast clock %s sat=%2d\n", time_str(time[i], 3), obs[i].sat);
//            continue;
//        }
//        time[i] = timeadd(time[i], -dt);//再减一遍
//
//        /* satellite position and clock at transmission time *///计算卫星位置和卫星钟
//        if (!satpos(time[i], teph, obs[i].sat, ephopt, nav, rs + i * 6, dts + i * 2, var + i,
//            svh + i)) {
//            trace(3, "no ephemeris %s sat=%2d\n", time_str(time[i], 3), obs[i].sat);
//            continue;
//        }
//        /* if no precise clock available, use broadcast clock instead */
//        if (dts[i * 2] == 0.0) {
//            if (!ephclk(time[i], teph, obs[i].sat, nav, dts + i * 2)) continue;
//            dts[1 + i * 2] = 0.0;
//            *var = SQR(STD_BRDCCLK);
//        }
//    }
//    for (i = 0; i < n && i < 2 * MAXOBS; i++) {
//        trace(4, "%s sat=%2d rs=%13.3f %13.3f %13.3f dts=%12.3f var=%7.3f svh=%02X\n",
//            time_str(time[i], 6), obs[i].sat, rs[i * 6], rs[1 + i * 6], rs[2 + i * 6],
//            dts[i * 2] * 1E9, var[i], svh[i]);
//    }
//}