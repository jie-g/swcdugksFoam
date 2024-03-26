#include <stdio.h>
#include "slave.h"
#include <math.h>
#include "dma_macros.h"
#include "slave_para.h"
#include <simd.h>

#include "para.h"
// #include <swperf.h>
// #define GHBARPVOL_USE_SIMD
// #define GHBARPVOL_TRAN_FUN
// #define GHBARPVOL_REG_COM

#define REG_PUTR(var, dst) asm volatile ("putr %0,%1\n"::"r"(var),"r"(dst))
#define REG_PUTC(var, dst) asm volatile ("putc %0,%1\n"::"r"(var),"r"(dst))
#define REG_GETR(var) asm volatile ("getr %0\n":"=r"(var))
#define REG_GETC(var) asm volatile ("getc %0\n":"=r"(var))
#define REG_SYNR(mask) \
  asm volatile ("synr %0"::"r"(mask))
#define REG_SYNC(mask) \
  asm volatile ("sync %0"::"r"(mask))

//加这个判断会变慢一些
#define NTHREAD 64
static const double VSMALL = 1.0e-300;
// #define N_slave 456
// typedef long int label//label
typedef long int int64;//label
//scalar  -> double
#define CLOCKRATE 1.45e9
#define CUBE(x) ((x)*(x)*(x))
// #define pi acos(-1)
#define pi 3.1415926535897932384626433832795
static inline unsigned long rpcc(){
	unsigned long addtime;
	asm("rcsr %0,4": "=r" (addtime) : );
	//rcsr %0,4" :  "=r" (time)
	return addtime;
}
// __thread_local volatile int reply_shadow = 0;                        
// __thread_local volatile int count_shadow = 0;                                 
// __thread_local volatile dma_desc pe_get_desc = 0, pe_put_desc = 0;            


typedef union
{
    double   f;
    // uint64_t u;
    struct {int32_t  i0,i1;} s;
}  udi_t;

/*******************************************************************************/
// #include <math.h>
#include <stdint.h>

static inline double pow2n(int n)//计算2的n次方
{
    // if(n<0) return 1/pow(a,-n);
    double res = 1.0;
	double a = 2;
    while(n)
    {
        if(n&1) res *= a;
        a *= a;
        n >>= 1;
    }
    return res;
    // return a*a*a;
}

// static inline double expd(double x){
//     int n = 11;
//     double log2e = 1.44269504088896338700;
//     double xlog = x * log2e;
//     // int64 xint = floor(xlog);//p
// 	unsigned short xint = (int)(xlog);//p
//     double xdeci = xlog - xint;//q
//     // printf("%lf,%d,%lf\n",xlog,xint, xdeci);
//     int i = 0;
//     double exp2p = 1;
//     for(i=0;i<xint;i++){
//         exp2p *= 2;
//     }
// 	// unsigned short exp2p = (unsigned short)(1 << xint);
// 	// double exp2p = (double)(1 << xint);//计算2的x次方
// 	// double exp2p = xint<30?exp2x[xint]:exp2_(xint);//更慢了
// 	// double exp2p = pow(2,xint);
//     double ln2 = 0.69314718055994528623;
//     double qln2 = xdeci * ln2; 
//     // double exp2q = exponential(qln2, 11);
// 	// double exp2q = 1.0;
//     // double term = 1.0;
//     // for (int i = 1; i <= 11; i++) {
//     //     term *= qln2 / i;
//     //     exp2q += term;
//     // }
// 	double exp2q = 1+qln2*(1+qln2*(0.5+qln2*(0.16666666666666666667+qln2*(0.04166666666666666667+qln2*(0.00833333333333333333+ \
//     qln2*(0.00138888888888888889+qln2*(0.00019841269841269841+qln2*(0.00002480158730158730+qln2*(0.00000275573192239859+qln2*0.00000027557319223986)))))))));
    
//     double result = exp2p * exp2q;
//     // printf("%lf,%lf,%.20lf,%lf\n",qln2,exp2p, exp2q, result);
//     // printf("%.20lf\n",result);
//     // double right = exp(x);
//     // printf("right value :%.20lf\n", right);
//     return result;
// }
static inline long double expd(double x){
    int n = 11;
    double log2e = 1.44269504088896338700;
    double xlog = x * log2e;
    // int64 xint = floor(xlog);//p
	unsigned short xint = (int)(xlog);//p
    double xdeci = xlog - xint;//q
    // printf("%lf,%d,%lf\n",xlog,xint, xdeci);
    int i = 0;
    long double exp2p = 1.0;
    for(i=0;i<xint;i++){
        exp2p *= 2;
    }
	// unsigned short exp2p = (unsigned short)(1 << xint);
	// double exp2p = (double)(1 << xint);//计算2的x次方
	// double exp2p = xint<30?exp2x[xint]:exp2_(xint);//更慢了
	// double exp2p = pow(2,xint);
    double ln2 = 0.69314718055994528623;
    long double qln2 = xdeci * ln2; 
    // double exp2q = exponential(qln2, 11);
	// double exp2q = 1.0;
    // double term = 1.0;
    // for (int i = 1; i <= 11; i++) {
    //     term *= qln2 / i;
    //     exp2q += term;
    // }
	// long double exp2q = 1+qln2*(1+qln2*(0.5+qln2*(0.16666666666666666667+qln2*(0.04166666666666666667+qln2*(0.00833333333333333333+ \
    // qln2*(0.00138888888888888889+qln2*(0.00019841269841269841+qln2*(0.00002480158730158730+qln2*(0.00000275573192239859+qln2*0.00000027557319223986)))))))));
    long double exp2q = 1+qln2*(1+qln2*(0.5+qln2*(1.0/6.0+qln2*(1.0/24.0+qln2*(1.0/120.0+ \
    qln2*(1.0/720.0+qln2*(1.0/5040.0+qln2*(1.0/40320.0+ \
	qln2*(1.0/362880.0+qln2*(1.0/3628800.0+qln2*(1.0/39916800.0+qln2*(1.0/479001600.0+ \
	qln2*(1.0/6227020800.0+qln2*(1.0/87178291200.0+qln2*(1.0/1307674368000.0+qln2*0.00000027557319223986)))))))))))))));
    
    long double result = exp2p * exp2q;
    // printf("%lf,%lf,%.20lf,%lf\n",qln2,exp2p, exp2q, result);
    // printf("%.20lf\n",result);
    // double right = exp(x);
    // printf("right value :%.20lf\n", right);
    return result;
}
static inline doublev4 simd_expd(doublev4 x){
    doublev4 log2e = 1.44269504088896338700;
    doublev4 xlog = x * log2e;
	double xlogval[4]__attribute__((aligned(256)));
	long xintval[4]__attribute__((aligned(256)));
	double xintvd[4]__attribute__((aligned(256)));
    // int64 xint = floor(xlog);
	simd_store(xlog,xlogval);
	xintval[0] = floor(xlogval[0]);
	xintval[1] = floor(xlogval[1]);
	xintval[2] = floor(xlogval[2]);
	xintval[3] = floor(xlogval[3]);
	int256 xint ;
  	simd_load(xint,&(xintval[0]));
  	// intv8 a_ = simd_set_intv8 (h,0,2,0,3,0,4,0);
  	// int256 a_i = a_;
	// int256 a_i = simd_set_int256 (xintval[0],xintval[2],xintval[3],xintval[4]);
  	doublev4 exp2p;
	int shfw_mask = 0x67452301;
  	asm ("ldi %0, 1023($31)\n\t"
       "vshff %0, %0, 0, %0\n\t"
       "vaddl %0, %1, %0\n\t"
       "vshfw %0, %0, %2, %0\n\t"
       "vsllw %0, 20, %0\n\t"
      //  "vaddl %0, $31, %1\n\t"
       : "=&r"(exp2p) : "r"(xint), "r"(shfw_mask));
	//IEEE表示中double的偏移为1023，vshff混洗后变为四个1023，再加上xint,vshfw进行字向量混洗，不知道为什么移位20，前面是指数加减？
	xintvd[0] = xintval[0];
	xintvd[1] = xintval[1];
	xintvd[2] = xintval[2];
	xintvd[3] = xintval[3];
	doublev4 xintv4d;
	simd_load(xintv4d,&(xintvd[0]));
	// unsigned short xint = (int)(xlog);//p
    doublev4 xdeci = xlog - xintv4d;//q
    // printf("%lf,%d,%lf\n",xlog,xint, xdeci);
    // int i = 0;
    // double exp2p = 1;
    // for(i=0;i<xint;i++){
    //     exp2p *= 2;
    // }
	// unsigned short exp2p = (unsigned short)(1 << xint);
	// double exp2p = (double)(1 << xint);//计算2的x次方
	// double exp2p = xint<30?exp2x[xint]:exp2_(xint);//更慢了
	// double exp2p = pow(2,xint);
    doublev4 ln2 = 0.69314718055994528623;
    doublev4 qln2 = xdeci * ln2; 
    // double exp2q = exponential(qln2, 11);
	// double exp2q = 1.0;
    // double term = 1.0;
    // for (int i = 1; i <= 11; i++) {
    //     term *= qln2 / i;
    //     exp2q += term;
    // }
	// doublev4 exp2q = 1+qln2*(1+qln2*(0.5+qln2*(0.16666666666666666667+qln2*(0.04166666666666666667+qln2*(0.00833333333333333333+ \
    // qln2*(0.00138888888888888889+qln2*(0.00019841269841269841+qln2*(0.00002480158730158730+qln2*(0.00000275573192239859+qln2*0.00000027557319223986)))))))));
	/*
    doublev4 qln2_2 = qln2*qln2;
    //奇数项
	doublev4 exp2q1 = 1.0+qln2_2*(0.5+qln2_2*(0.04166666666666666667+qln2_2*(0.00138888888888888889 \
	+qln2_2*(0.00002480158730158730+qln2_2*(0.00000027557319223986+qln2_2*0.00000000208767569878)))));
	//偶数项
	doublev4 exp2q2 = qln2*(1.0+qln2_2*(0.16666666666666666667+qln2_2*(0.00833333333333333333+ \
	qln2_2*(0.00019841269841269841+qln2_2*(0.00000275573192239859+qln2_2*(0.00000002505210838544+qln2_2*0.00000000016059043836))))));
    doublev4 result = exp2p * (exp2q1+exp2q2);
	*/
    doublev4 exp2q = 1+qln2*(1+qln2*(0.5+qln2*(0.16666666666666666667+qln2*(0.04166666666666666667+qln2*(0.00833333333333333333+ \
    qln2*(0.00138888888888888889+qln2*(0.00019841269841269841+qln2*(0.00002480158730158730+qln2*(0.00000275573192239859+qln2*0.00000027557319223986)))))))));
    doublev4 result = exp2p * exp2q;

    // printf("%lf,%lf,%.20lf,%lf\n",qln2,exp2p, exp2q, result);
    // printf("%.20lf\n",result);
    // double right = exp(x);
    // printf("right value :%.20lf\n", right);
    return result;
}
// double rsqrtd (double number) //计算1/sqrt
// {
//   long i;
//   double x2, y;
//   const double threehalfs=1.5;//1.5F
//   x2 = number*0.5;//0.5F
//   y = number;
//   i = *(long *) &y;
//   i = 0x5fe6eb50c7b537a9-(i>>1);
//   y = *(double *) &i;
//   y = y*(threehalfs-(x2*y*y));
//   y  = y * ( threehalfs - ( x2 * y * y ) );
//   y  = y * ( threehalfs - ( x2 * y * y ) );
//   y  = y * ( threehalfs - ( x2 * y * y ) );//次数越多，结果越精确
//   return y;
// }
// static inline double powd(double x,int n){
double rsqrtd (double number) //计算sqrt  降速了
{

  long i;
  double x2, y, tmp2, tmp3;
  x2 = number*0.5;//0.5F
  y = number;
  i = *(long *) &y;
  i = 0x5fe6eb50c7b537a9-(i>>1);
  y = *(double *) &i;
  y  = y * ( 1.5 - ( x2 * y * y ) );
  y  = y * ( 1.5 - ( x2 * y * y ) );
  y  = y * ( 1.5 - ( x2 * y * y ) );
  y  = y * ( 1.5 - ( x2 * y * y ) );//次数越多，结果越精确
  tmp2 = number * y;
  tmp3 = number - tmp2*tmp2;
  y = tmp2 + tmp3*y;
  return y;
}

static inline long double powd(long double x,int n){
	
	switch (n) {
        case 1:
            return x;
        case 2:
            return x*x;
        case 3:
            return x*x*x;
    }
	return pow(x,n);
	// double result = x;
	// for(int i=1;i<n;i++)//好像是会损失精度，cavity2迭代会多一轮
	// 	result *= x;
	// return result;
}


/*******************************************************************/
__thread_local volatile unsigned long get_reply,put_reply;
__thread_local volatile unsigned long get_reply_eq[2],put_reply_eq[2];

void Func_ghbarpvol(GHbarPvol *ghbarpvol) {
	dma_init();
	// int N_slave = 640;
#ifndef GHBARPVOL_REG_COM
	int N_slave = 512;//456  448
#else
	int N_slave = 448 ;//456  448 480
#endif
	GHbarPvol slave;
	// unsigned long get_reply,put_reply; 
    int my_id; 
    my_id = athread_get_id(-1);  //获取从核的id号
    get_reply = 0;
    athread_get(PE_MODE , ghbarpvol, &slave, sizeof(GHbarPvol), &get_reply , 0, 0, 0);
    while (get_reply != 1);//等待数据加载完成
	// get_reply = 0;
    // dma_set_size(&pe_get_desc, sizeof(GHbarPvol));
    // dma(pe_get_desc, ghbarpvol, &slave);
    // dma_wait(&get_reply, 1);
	// pe_get(ghbarpvol, &slave, sizeof(GHbarPvol));
	// dma_syn();
	int64 DVsize = ghbarpvol->DVsize;
	int64 size = ghbarpvol->ghBarPvol_size;
	// Load Balance
	int load = size / NTHREAD;
	int rest = size % NTHREAD;
	int stid = (my_id < rest) ? my_id * (load + 1) : my_id * load + rest;
	if (my_id < rest) load++;
	int edid = stid + load;
	int numtemp = edid - stid;
	//Define variable
	//double tauVol_value[N_slave]__attribute__((aligned(256)))  11-9
	double tauVol_value[N_slave]__attribute__ ((aligned(128)));
	double U_value[N_slave*3]__attribute__ ((aligned(128)));
	double q_value[N_slave*3]__attribute__ ((aligned(128)));
	double T_value[N_slave]__attribute__ ((aligned(128)));
	double rho_value[N_slave]__attribute__ ((aligned(128)));
	double gBarPvol_value[N_slave]__attribute__ ((aligned(128)));
	double hBarPvol_value[N_slave]__attribute__ ((aligned(128)));
	double gTildeVol_value[N_slave]__attribute__ ((aligned(128)));
	double hTildeVol_value[N_slave]__attribute__ ((aligned(128)));
#ifdef GHBARPVOL_REG_COM
	double xi_value[N_slave*3]__attribute__ ((aligned(128)));
#endif
	// double relaxFactor[N_slave];
	double dt = slave.dt;
	// double xii_x = ghsurf->xii_x;
	// double xii_y = ghsurf->xii_y;
	// double xii_z = ghsurf->xii_z;
	double xii_x , xii_y, xii_z;
	double R_value = slave.R_value;
	double Pr_value = slave.Pr_value;
	double vUnit_value = slave.vUnit_value;
	int64 D = slave.D;
	int64 K = slave.K;
	int facei = 0;
	int num = 0, xnum = 0;
	int step = 0, xstep = 0;
	int stnum = 0, xstnum = 0;
	int xnumtemp = DVsize;
	double cSqrByRT_value = 0;
	double cqBy5pRT_value = 0;
	double gEqBGK_value = 0;
	double relaxFactor = 0;
	int dvi = 0;
	double temp_1;
	doublev4 xv4;
// #ifdef GHBARPVOL_USE_SIMD
// 	double mantissa[4]__attribute__((aligned(128)));
// 	doublev4 numberv4,mantissav4;
//    	long exponent[4]__attribute__((aligned(128)));//指数
// 	double relaxFactor_[4]__attribute__ ((aligned(128)));
// 	double cSqrByRT_[4]__attribute__ ((aligned(128)));
// 	double cqBy5pRT_[4]__attribute__ ((aligned(128)));
// 	double gEqBGK_[4]__attribute__ ((aligned(128)));
// 	double temp[4]__attribute__ ((aligned(128)));
// 	doublev4 relaxFactorv4=0.0,tauVolv4=0.0,gTildeVolv4=0.0,hTildeVolv4=0.0;
// 	doublev4 gBarPvolv4=0.0,hBarPvolv4=0.0;
// 	doublev4 Uxv4=0.0,Uyv4=0.0,Uzv4=0.0,Tv4=0.0,rhov4=0.0;
// 	doublev4 cSqrByRTv4=0.0,cqBy5pRTv4=0.0,gEqBGKv4=0.0;
// 	doublev4 qxv4=0.0,qyv4=0.0,qzv4=0.0;
// 	doublev4 temp1=0.0,temp2=0.0,temp3=0.0,temp4=0.0,tempv4=0.0;
// #endif
	while(numtemp>0){
		num = (numtemp < N_slave) ? numtemp : N_slave;
		numtemp = numtemp - N_slave;
		stnum = step * N_slave;
		// get data
		pe_get(slave.tauVol_value+stid+stnum, tauVol_value, sizeof(double)*num);
		pe_get(slave.U_value+stid*3+stnum*3, U_value, sizeof(double)*num*3);
		pe_get(slave.q_value+stid*3+stnum*3, q_value, sizeof(double)*num*3);
		pe_get(slave.T_value+stid+stnum, T_value, sizeof(double)*num);
		pe_get(slave.rho_value+stid+stnum, rho_value, sizeof(double)*num);
		dma_syn();
		xstep = 0;
		xstnum = 0;
		xnum = 0;
		xnumtemp = DVsize;
    	for(dvi=0;dvi<DVsize;dvi++){
#ifndef GHBARPVOL_REG_COM
			xii_x = slave.xi_value[dvi*3];
			xii_y = slave.xi_value[dvi*3+1];
			xii_z = slave.xi_value[dvi*3+2];
#else
			//第一列读主存，然后行广播给其他核
			if(my_id % 8 == 0){
				if(dvi==(xstnum+xnum)){
					xnum = (xnumtemp < N_slave) ? xnumtemp : N_slave;
					xnumtemp = xnumtemp - N_slave;
					xstnum = xstep * N_slave;
					//第一列读主存，然后行广播给其他核
					pe_get(slave.xi_value+xstnum*3, xi_value, sizeof(double)*xnum*3);
					dma_syn();
					xstep++;
				}
				((double*)(&xv4))[0] = xi_value[(dvi-xstnum)*3];
				((double*)(&xv4))[1] = xi_value[(dvi-xstnum)*3+1];
				((double*)(&xv4))[2] = xi_value[(dvi-xstnum)*3+2];
				REG_PUTR(xv4,8);
				xii_x = xi_value[(dvi-xstnum)*3];
				xii_y = xi_value[(dvi-xstnum)*3+1];
				xii_z = xi_value[(dvi-xstnum)*3+2];
			}
			else{
			// if(my_id % 8 != 0){
				REG_GETR(xv4);
				xii_x = ((double*)(&xv4))[0];
				xii_y = ((double*)(&xv4))[1];
				xii_z = ((double*)(&xv4))[2];
			}
			// REG_SYNR(0xff);
#endif
			pe_get(slave.gTildeVol_value[dvi]+stid+stnum, gTildeVol_value, sizeof(double)*num);
			pe_get(slave.hTildeVol_value[dvi]+stid+stnum, hTildeVol_value, sizeof(double)*num);
			dma_syn();
#ifndef GHBARPVOL_USE_SIMD

			for(facei=0;facei<num;facei++)
			{
				// relaxFactor = h/(2*tauSurf_value[facei] + h);
				relaxFactor = 1.5*dt/(2.0*tauVol_value[facei] + dt);
				gBarPvol_value[facei] = (1.0 - relaxFactor)*gTildeVol_value[facei] ;
				hBarPvol_value[facei] = (1.0 - relaxFactor)*hTildeVol_value[facei] ;
				cSqrByRT_value = ((U_value[facei*3]-xii_x)*(U_value[facei*3]-xii_x)+(U_value[facei*3+1]-xii_y)*(U_value[facei*3+1]-xii_y)+(U_value[facei*3+2]-xii_z)*(U_value[facei*3+2]-xii_z))/(R_value*T_value[facei]);
				cqBy5pRT_value = ((xii_x - U_value[facei*3])*q_value[facei*3]+(xii_y - U_value[facei*3+1])*q_value[facei*3+1]+(xii_z - U_value[facei*3+2])*q_value[facei*3+2])/(5.0*rho_value[facei]*R_value*T_value[facei]*R_value*T_value[facei]);
#ifndef GHBARPVOL_TRAN_FUN
				gEqBGK_value  = (relaxFactor*rho_value[facei])/pow(sqrt(2.0*pi*R_value*T_value[facei]),D)*exp(-cSqrByRT_value/2.0);
				//替换
#else				
				/******************************以下运算是ok的***************************/
				// gEqBGK_value  = (relaxFactor*rho_value[facei])/pow(sqrt(2.0*pi*R_value*T_value[facei]),D)*exp(-cSqrByRT_value/2.0)/pow(vUnit_value, 3-D);
				// gEqBGK_value  = (relaxFactor*rho_value[facei])/pow(sqrt(2.0*pi*R_value*T_value[facei]),D)*(1/expd( cSqrByRT_value/2.0));
				// gEqBGK_value  = (relaxFactor*rho_value[facei])*(1/expd(cSqrByRT_value/2.0));
				gEqBGK_value  = (relaxFactor*rho_value[facei])*(1/expd(cSqrByRT_value/2.0));
				temp_1 = sqrt(2.0*pi*R_value*T_value[facei]);
				
				// temp = rsqrtd(2.0*pi*R_value*T_value[facei]);
				
				// gEqBGK_value /= pow(temp,D);
				// gEqBGK_value /= temp;
				// gEqBGK_value /= temp;
				// gEqBGK_value /= temp;
				double mantissa, number=1.0;
   				int exponent;
   				// number = 8.0;
   				mantissa = frexp(temp_1, &exponent);
				for(int w = 0;w<exponent;w++){
					number *= 2;
				}
				// number = pow2n(exponent);
				for(int q=0;q<D;q++){
					gEqBGK_value /= (number*mantissa);
				}
				// mantissa = frexp(temp, &exponent);
				// // for(int w = 0;w<exponent;w++){
				// // 	number *= 2;
				// // }
				// udi_t    epart;//和循环计算差距不大
				// epart.s.i0 = 0;
				// epart.s.i1 = ((exponent) + 1023) << 20;
				// for(int q=0;q<D;q++){
				// 	gEqBGK_value /= (epart.f*mantissa);
				// 	// gEqBGK_value /= (number*mantissa);
				// }
				/**********************************************************************/

#endif				
				
				gBarPvol_value[facei] += ( 1.0 + (1.0 - Pr_value)*cqBy5pRT_value*(cSqrByRT_value - D - 2.0) )*gEqBGK_value;
				hBarPvol_value[facei] += ( (K + 3.0 - D) + (1.0 - Pr_value)*cqBy5pRT_value*((cSqrByRT_value - D)*(K + 3.0 - D) - 2*K) )*gEqBGK_value*R_value*T_value[facei];
				
				//原void Foam::discreteVelocity::updateGHtildeVol()函数中的计算
				gTildeVol_value[facei] = -1.0/3*gTildeVol_value[facei] + 4.0/3*gBarPvol_value[facei];
				hTildeVol_value[facei] = -1.0/3*hTildeVol_value[facei] + 4.0/3*hBarPvol_value[facei];
				// // store the gTildePlus in gTilde
    			// gTildeVol_ = -1.0/3*gTildeVol_ + 4.0/3*gBarPvol_;
    			// hTildeVol_ = -1.0/3*hTildeVol_ + 4.0/3*hBarPvol_;

			}
#else
			double relaxFactor_[4]__attribute__ ((aligned(128)));
			double cSqrByRT_[4]__attribute__ ((aligned(128)));
			double cqBy5pRT_[4]__attribute__ ((aligned(128)));
			double gEqBGK_[4]__attribute__ ((aligned(128)));
			double temp[4]__attribute__ ((aligned(128)));
			doublev4 relaxFactorv4=0.0,tauVolv4=0.0,gTildeVolv4=0.0,hTildeVolv4=0.0;
			doublev4 gBarPvolv4=0.0,hBarPvolv4=0.0;
			doublev4 Uxv4=0.0,Uyv4=0.0,Uzv4=0.0,Tv4=0.0,rhov4=0.0;
			doublev4 cSqrByRTv4=0.0,cqBy5pRTv4=0.0,gEqBGKv4=0.0;
			doublev4 qxv4=0.0,qyv4=0.0,qzv4=0.0;
			doublev4 temp1=0.0,temp2=0.0,temp3=0.0,temp4=0.0,tempv4=0.0;
			for(facei=0;facei<(num/4)*4;facei+=4)
			{
				simd_load(tauVolv4, &(tauVol_value[facei]));
				relaxFactorv4 = 1.5*dt/(2.0*tauVolv4 + dt);
				// simd_store(relaxFactorv4, &(relaxFactor_[0]));
				simd_load(gTildeVolv4, &(gTildeVol_value[facei]));
				gBarPvolv4 = (1.0 - relaxFactorv4)*gTildeVolv4;
				simd_load(hTildeVolv4, &(hTildeVol_value[facei]));
				hBarPvolv4 = (1.0 - relaxFactorv4)*hTildeVolv4;
				/**************************************************/
				simd_load(Uxv4, &(U_value[facei*3]));
				simd_load(Uyv4, &(U_value[facei*3+4]));
				simd_load(Uzv4, &(U_value[facei*3+8]));
				temp1 = simd_vshuffle(Uxv4, Uyv4, 0b10010100);//
				temp2 = simd_vshuffle(Uyv4, Uzv4, 0b11101001);//
				Uxv4 = simd_vshuffle(temp2, Uxv4, 0b00101100);//x1-4
				Uyv4 = simd_vshuffle(temp2, temp1, 0b01110010);//y1-4
				Uzv4 = simd_vshuffle(Uzv4, temp1, 0b11000111);//z1-4
				/**************************************************/
				simd_load(Tv4, &(T_value[facei]));

				// simd_store(gBarPvolv4, &(gBarPvol_value[facei]));
				// simd_store(hBarPvolv4, &(hBarPvol_value[facei]));
				
				// Uxv4=simd_set_doublev4(U_value[facei*3],U_value[(facei+1)*3],U_value[(facei+2)*3],U_value[(facei+3)*3]);
				// Uyv4=simd_set_doublev4(U_value[facei*3+1],U_value[(facei+1)*3+1],U_value[(facei+2)*3+1],U_value[(facei+3)*3+1]);
				// Uzv4=simd_set_doublev4(U_value[facei*3+2],U_value[(facei+1)*3+2],U_value[(facei+2)*3+2],U_value[(facei+3)*3+2]);
				
				cSqrByRTv4 = ((Uxv4-xii_x)*(Uxv4-xii_x)+(Uyv4-xii_y)*(Uyv4-xii_y)+(Uzv4-xii_z)*(Uzv4-xii_z))/(R_value*Tv4);
				// simd_store(cSqrByRTv4, &(cSqrByRT_[0]));
				/************************************************************/
				simd_load(qxv4, &(q_value[facei*3]));
				simd_load(qyv4, &(q_value[facei*3+4]));
				simd_load(qzv4, &(q_value[facei*3+8]));
				temp3 = simd_vshuffle(qxv4, qyv4, 0b10010100);//
				temp4 = simd_vshuffle(qyv4, qzv4, 0b11101001);//
				qxv4 = simd_vshuffle(temp4, qxv4, 0b00101100);//x1-4
				qyv4 = simd_vshuffle(temp4, temp3, 0b01110010);//y1-4
				qzv4 = simd_vshuffle(qzv4, temp3, 0b11000111);//z1-4
				/***********************************************************/
				simd_load(rhov4, &(rho_value[facei]));
				// qxv4=simd_set_doublev4(q_value[facei*3],q_value[(facei+1)*3],q_value[(facei+2)*3],q_value[(facei+3)*3]);
				// qyv4=simd_set_doublev4(q_value[facei*3+1],q_value[(facei+1)*3+1],q_value[(facei+2)*3+1],q_value[(facei+3)*3+1]);
				// qzv4=simd_set_doublev4(q_value[facei*3+2],q_value[(facei+1)*3+2],q_value[(facei+2)*3+2],q_value[(facei+3)*3+2]);
				
				cqBy5pRTv4 = ((xii_x - Uxv4)*qxv4+(xii_y - Uyv4)*qyv4+(xii_z - Uzv4)*qzv4)/(5.0*rhov4*R_value*Tv4*R_value*Tv4);
				// simd_store(cqBy5pRTv4, &(cqBy5pRT_[0]));
				
#ifndef GHBARPVOL_TRAN_FUN
				simd_store(cSqrByRTv4, &(cSqrByRT_[0]));
				simd_store(relaxFactorv4, &(relaxFactor_[0]));
				gEqBGK_[0]  = (relaxFactor_[0]*rho_value[facei])/pow(sqrt(2.0*pi*R_value*T_value[facei]),D)*exp(-cSqrByRT_[0]/2.0);
				gEqBGK_[1]  = (relaxFactor_[1]*rho_value[facei+1])/pow(sqrt(2.0*pi*R_value*T_value[facei+1]),D)*exp(-cSqrByRT_[1]/2.0);
				gEqBGK_[2]  = (relaxFactor_[2]*rho_value[facei+2])/pow(sqrt(2.0*pi*R_value*T_value[facei+2]),D)*exp(-cSqrByRT_[2]/2.0);
				gEqBGK_[3]  = (relaxFactor_[3]*rho_value[facei+3])/pow(sqrt(2.0*pi*R_value*T_value[facei+3]),D)*exp(-cSqrByRT_[3]/2.0);
				simd_load(gEqBGKv4, &(gEqBGK_[0]));

#else
				/****************************************************************/
				// gEqBGKv4 = relaxFactorv4 * rhov4 /simd_expd(cSqrByRTv4/2.0);//会报错--问题原因还未知
				/*simd_store(cSqrByRTv4, &(cSqrByRT_[0]));
				simd_store(relaxFactorv4, &(relaxFactor_[0]));
				gEqBGK_[0] = relaxFactor_[0]*rho_value[facei]/expd(cSqrByRT_[0]/2.0);
				gEqBGK_[1] = relaxFactor_[1]*rho_value[facei+1]/expd(cSqrByRT_[1]/2.0);
				gEqBGK_[2] = relaxFactor_[2]*rho_value[facei+2]/expd(cSqrByRT_[2]/2.0);
				gEqBGK_[3] = relaxFactor_[3]*rho_value[facei+3]/expd(cSqrByRT_[3]/2.0);
				simd_load(gEqBGKv4, &(gEqBGK_[0]));
				*/
				simd_store(cSqrByRTv4, &(cSqrByRT_[0]));
				simd_store(relaxFactorv4, &(relaxFactor_[0]));
				gEqBGK_[0]  = (relaxFactor_[0]*rho_value[facei]);
				gEqBGK_[1]  = (relaxFactor_[1]*rho_value[facei+1]);
				gEqBGK_[2]  = (relaxFactor_[2]*rho_value[facei+2]);
				gEqBGK_[3]  = (relaxFactor_[3]*rho_value[facei+3]);
				simd_load(gEqBGKv4, &(gEqBGK_[0]));
				// temp[0] = pow(sqrt(2.0*pi*R_value*T_value[facei]),D);
				// temp[1] = pow(sqrt(2.0*pi*R_value*T_value[facei+1]),D);
				// temp[2] = pow(sqrt(2.0*pi*R_value*T_value[facei+2]),D);
				// temp[3] = pow(sqrt(2.0*pi*R_value*T_value[facei+3]),D);
				// simd_store(tempv4,&(temp[0]));
				// gEqBGKv4 = gEqBGKv4/tempv4;

				// simd_store(cSqrByRTv4, &(cSqrByRT_[0]));
				// gEqBGK_[0] = rho_value[facei]*(exp(-cSqrByRT_[0]/2.0));
				// gEqBGK_[1] = rho_value[facei+1]*(exp(-cSqrByRT_[1]/2.0));
				// gEqBGK_[2] = rho_value[facei+2]*(exp(-cSqrByRT_[2]/2.0));
				// gEqBGK_[3] = rho_value[facei+3]*(exp(-cSqrByRT_[3]/2.0));
				// simd_load(gEqBGKv4, &(gEqBGK_[0]));

				//11-6有问题
				// gEqBGKv4 = rhov4 * (1/simd_expd(cSqrByRTv4/2.0));
			
				temp[0] = sqrt(2.0*pi*R_value*T_value[facei]);
				temp[1] = sqrt(2.0*pi*R_value*T_value[facei+1]);
				temp[2] = sqrt(2.0*pi*R_value*T_value[facei+2]);
				temp[3] = sqrt(2.0*pi*R_value*T_value[facei+3]);
				// simd_load(tempv4,&(temp[0]));

				// temp[0] = 1/rsqrtd(2.0*pi*R_value*T_value[facei]);
				// temp[1] = 1/rsqrtd(2.0*pi*R_value*T_value[facei+1]);
				// temp[2] = 1/rsqrtd(2.0*pi*R_value*T_value[facei+2]);
				// temp[3] = 1/rsqrtd(2.0*pi*R_value*T_value[facei+3]);
				// simd_load(tempv4,&(temp[0]));

				// tempv4 = simd_vsqrtd(2.0*pi*R_value*Tv4);
				// simd_store(tempv4,&(temp[0]));
			
			    // /***************************/
				// // double mantissa[4]__attribute__((aligned(128)));
				// // doublev4 numberv4,mantissav4;
   				// // long exponent[4]__attribute__((aligned(128)));//指数
   				// mantissa[0] = frexp(temp[0], &(exponent[0]));
				// mantissa[1] = frexp(temp[1], &(exponent[1]));
				// mantissa[2] = frexp(temp[2], &(exponent[2]));
				// mantissa[3] = frexp(temp[3], &(exponent[3]));
				// int256 exponentv4;
				// simd_load(exponentv4,&(exponent[0]));
				// int shfw_mask = 0x67452301;
				// asm ("ldi %0, 1023($31)\n\t"
				// 	"vshff %0, %0, 0, %0\n\t"
				// 	"vaddl %0, %1, %0\n\t"
				// 	"vshfw %0, %0, %2, %0\n\t"
				// 	"vsllw %0, 20, %0\n\t"
				// 	//  "vaddl %0, $31, %1\n\t"
				// 	: "=&r"(numberv4) : "r"(exponentv4), "r"(shfw_mask));
				// simd_load(mantissav4,&(mantissa[0]));
				// for(int q=0;q<D;q++){
				// 	gEqBGKv4 /= (numberv4*mantissav4);
				// }
			    /***************************************/
				simd_load(tempv4,&(temp[0]));
				for(int q=0;q<D;q++){
					gEqBGKv4 /= tempv4;
				}
				/*****************************************/
				simd_store(gEqBGKv4, &(gEqBGK_[0]));
				// gEqBGK_[0]  *= exp(-cSqrByRT_[0]/2.0);
				// gEqBGK_[1]  *= exp(-cSqrByRT_[1]/2.0);
				// gEqBGK_[2]  *= exp(-cSqrByRT_[2]/2.0);
				// gEqBGK_[3]  *= exp(-cSqrByRT_[3]/2.0);
				gEqBGK_[0]  *= 1.0/expd(cSqrByRT_[0]/2.0);
				gEqBGK_[1]  *= 1.0/expd(cSqrByRT_[1]/2.0);
				gEqBGK_[2]  *= 1.0/expd(cSqrByRT_[2]/2.0);
				gEqBGK_[3]  *= 1.0/exp (cSqrByRT_[3]/2.0);
				// gEqBGK_[0]  *= 1.0/expd(cSqrByRT_[0]/2.0);
				// gEqBGK_[1]  *= 1.0/expd(cSqrByRT_[1]/2.0);
				// gEqBGK_[2]  *= 1.0/expd(cSqrByRT_[2]/2.0);
				// gEqBGK_[3]  *= 1.0/expd(cSqrByRT_[3]/2.0);
				simd_load(gEqBGKv4, &(gEqBGK_[0]));

				// gEqBGKv4 *= 1.0/simd_expd(cSqrByRTv4/2.0);
				/*******************************************/
				/***************************************************************************/
#endif			

				gBarPvolv4 += ( 1.0 + (1.0 - Pr_value)*cqBy5pRTv4*(cSqrByRTv4 - D - 2.0) )*gEqBGKv4;
				simd_store(gBarPvolv4, &(gBarPvol_value[facei]));
				hBarPvolv4 += ((K + 3.0 - D) + (1.0 - Pr_value)*cqBy5pRTv4*((cSqrByRTv4 - D)*(K + 3.0 - D) - K - K) )*gEqBGKv4*R_value*Tv4;//如果换成2.0*K就会产生误差
				simd_store(hBarPvolv4, &(hBarPvol_value[facei]));				

				
				//原void Foam::discreteVelocity::updateGHtildeVol()函数中的计算
				gTildeVolv4 = -1.0/3*gTildeVolv4 + 4.0/3*gBarPvolv4;
				simd_store(gTildeVolv4, &(gTildeVol_value[facei]));
				hTildeVolv4 = -1.0/3*hTildeVolv4 + 4.0/3*hBarPvolv4;
				simd_store(hTildeVolv4, &(hTildeVol_value[facei]));	
				// // store the gTildePlus in gTilde
    			// gTildeVol_ = -1.0/3*gTildeVol_ + 4.0/3*gBarPvol_;
    			// hTildeVol_ = -1.0/3*hTildeVol_ + 4.0/3*hBarPvol_;
			}
			for(facei=(num/4)*4;facei<num;facei++)
			{
				relaxFactor = 1.5*dt/(2.0*tauVol_value[facei] + dt);
				gBarPvol_value[facei] = (1.0 - relaxFactor)*gTildeVol_value[facei] ;
				hBarPvol_value[facei] = (1.0 - relaxFactor)*hTildeVol_value[facei] ;
				cSqrByRT_value = ((U_value[facei*3]-xii_x)*(U_value[facei*3]-xii_x)+(U_value[facei*3+1]-xii_y)*(U_value[facei*3+1]-xii_y)+(U_value[facei*3+2]-xii_z)*(U_value[facei*3+2]-xii_z))/(R_value*T_value[facei]);
				cqBy5pRT_value = ((xii_x - U_value[facei*3])*q_value[facei*3]+(xii_y - U_value[facei*3+1])*q_value[facei*3+1]+(xii_z - U_value[facei*3+2])*q_value[facei*3+2])/(5.0*rho_value[facei]*R_value*T_value[facei]*R_value*T_value[facei]);
				gEqBGK_value  = (relaxFactor*rho_value[facei])/pow(sqrt(2.0*pi*R_value*T_value[facei]),D)*exp(-cSqrByRT_value/2.0);
				gBarPvol_value[facei] += ( 1.0 + (1.0 - Pr_value)*cqBy5pRT_value*(cSqrByRT_value - D - 2.0) )*gEqBGK_value;
				hBarPvol_value[facei] += ( (K + 3.0 - D) + (1.0 - Pr_value)*cqBy5pRT_value*((cSqrByRT_value - D)*(K + 3.0 - D) - 2*K) )*gEqBGK_value*R_value*T_value[facei];
				
				//原void Foam::discreteVelocity::updateGHtildeVol()函数中的计算
				gTildeVol_value[facei] = -1.0/3*gTildeVol_value[facei] + 4.0/3*gBarPvol_value[facei];
				hTildeVol_value[facei] = -1.0/3*hTildeVol_value[facei] + 4.0/3*hBarPvol_value[facei];
				// // store the gTildePlus in gTilde
    			// gTildeVol_ = -1.0/3*gTildeVol_ + 4.0/3*gBarPvol_;
    			// hTildeVol_ = -1.0/3*hTildeVol_ + 4.0/3*hBarPvol_;
			}
#endif
			pe_put(slave.gBarPvol_value[dvi]+stid+stnum, gBarPvol_value, sizeof(double)*num);
			pe_put(slave.hBarPvol_value[dvi]+stid+stnum, hBarPvol_value, sizeof(double)*num);
			pe_put(slave.gTildeVol_value[dvi]+stid+stnum, gTildeVol_value, sizeof(double)*num);
			pe_put(slave.hTildeVol_value[dvi]+stid+stnum, hTildeVol_value, sizeof(double)*num);
			dma_syn();
		}
		step++;
	}
}
