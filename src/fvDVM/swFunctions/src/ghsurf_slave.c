#include <stdio.h>
#include "slave.h"
#include <math.h>
#include "dma_macros.h"
#include "slave_para.h"
#include <simd.h>

#include "para.h"
// #include <swperf.h>
// #define GHSURF_USE_SIMD//和修改超越函数一起，可以加速，单数误差变大
// #define GHSURF_TRAN_FUN
// #define GHSURF_REG_COM//这里的寄存器通信和双缓冲只能留一个，同时用两个数据上有问题，原因还不知道

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

double rsqrtd (double number) //计算1/sqrt
{
  long i;
  double x2, y;
  const double threehalfs=1.5;//1.5F
  x2 = number*0.5;//0.5F
  y = number;
  i = *(long *) &y;
  i = 0x5fe6eb50c7b537a9-(i>>1);
  y = *(double *) &i;
  y = y*(threehalfs-(x2*y*y));
  y  = y * ( threehalfs - ( x2 * y * y ) );
  y  = y * ( threehalfs - ( x2 * y * y ) );
  y  = y * ( threehalfs - ( x2 * y * y ) );//次数越多，结果越精确
  return y;
}

/*******************************************************************/
__thread_local volatile unsigned long get_reply,put_reply;
__thread_local volatile unsigned long get_reply_eq[2],put_reply_eq[2];


void Func_ghsurf(GHsurf *ghsurf) {
	unsigned long time_sta = 0;
	unsigned long time_end = 0;
	unsigned long time_sta1 = 0;
	unsigned long time_end1 = 0;
	unsigned long time_sta2 = 0;
	unsigned long time_end2 = 0;
	unsigned long time_sta3 = 0;
	unsigned long time_end3 = 0;
	unsigned long time1 = 0;
	unsigned long time2 = 0;
	unsigned long ic=0, ic1 = 0;
	// penv_slave0_cycle_count(&ic);
	penv_slave2_gld_count(&ic);
	time1 = rpcc();
	// int N_slave = 640;
	int N_slave = 456;
	// int N_slave = 384;
	GHsurf slave;
	// unsigned long get_reply,put_reply; 
    int my_id; 
    my_id = athread_get_id(-1);  //获取从核的id号
	dma_init();
    get_reply = 0;
    athread_get(PE_MODE , ghsurf, &slave, sizeof(GHsurf), &get_reply , 0, 0, 0);
    while (get_reply != 1);//等待数据加载完成
	int64 DVsize = ghsurf->DVsize;
	int64 size = ghsurf->gSurf_size;
	// Load Balance
	int load = size / NTHREAD;
	int rest = size % NTHREAD;
	int stid = (my_id < rest) ? my_id * (load + 1) : my_id * load + rest;
	if (my_id < rest) load++;
	int edid = stid + load;
	int numtemp = edid - stid;
	//Define variable
	double tauSurf_value[N_slave]__attribute__((aligned(128)));
	double U_value[N_slave*3+4]__attribute__((aligned(128)));
	double q_value[N_slave*3+4]__attribute__((aligned(128)));
	double T_value[N_slave]__attribute__((aligned(128)));
	double rho_value[N_slave]__attribute__((aligned(128)));
	double gSurf_value[2][N_slave]__attribute__((aligned(128)));
	double hSurf_value[2][N_slave]__attribute__((aligned(128)));
	// double Ux_value[N_slave]__attribute__((aligned(256)));
	// double Uy_value[N_slave]__attribute__((aligned(256)));
	// double Uz_value[N_slave]__attribute__((aligned(256)));
	// double qx_value[N_slave]__attribute__((aligned(256)));
	// double qy_value[N_slave]__attribute__((aligned(256)));
	// double qz_value[N_slave]__attribute__((aligned(256)));
	// double gSurf_value[N_slave];
	// double hSurf_value[N_slave];
	// double relaxFactor[N_slave];
#ifdef GHSURF_REG_COM
	double xi_value[N_slave*3]__attribute__ ((aligned(128)));
#endif

	int xnum = 0;
	int xstep = 0;
	int xstnum = 0;
	int xnumtemp = DVsize;
	doublev4 xv4;

	double h = slave.h;
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
	int num = 0;
	int step = 0;
	int stnum = 0;
	double cSqrByRT_value = 0;
	double cqBy5pRT_value = 0;
	double gEqBGK_value = 0;
	double relaxFactor = 0;
	int dvi = 0;
	int index = 0, next = 0, last = 0;
	while(numtemp>0){
		num = (numtemp < N_slave) ? numtemp : N_slave;
		numtemp = numtemp - N_slave;
		stnum = step * N_slave;
		// time_sta = rpcc();
		//get data
		pe_get(slave.tauSurf_value+stid+stnum, tauSurf_value, sizeof(double)*num);
		pe_get(slave.U_value+stid*3+stnum*3, U_value, sizeof(double)*num*3);
		pe_get(slave.q_value+stid*3+stnum*3, q_value, sizeof(double)*num*3);
		pe_get(slave.T_value+stid+stnum, T_value, sizeof(double)*num);
		pe_get(slave.rho_value+stid+stnum, rho_value, sizeof(double)*num);
		dma_syn();
		U_value[N_slave*3+1] = U_value[N_slave*3+2] = U_value[N_slave*3+3] = U_value[N_slave*3+4]= 0;
		q_value[N_slave*3+1] = q_value[N_slave*3+2] = q_value[N_slave*3+3] = q_value[N_slave*3+4]= 0;
		/*for(int i=0;i<(num/4)*4;i+=4){
			Ux_value[i]=U_value[i*3];Uy_value[i]=U_value[i*3+1];Uz_value[i]=U_value[i*3+2];
			Ux_value[i+1]=U_value[(i+1)*3];Uy_value[i+1]=U_value[(i+1)*3+1];Uz_value[i+1]=U_value[(i+1)*3+2];
			Ux_value[i+2]=U_value[(i+2)*3];Uy_value[i+2]=U_value[(i+2)*3+1];Uz_value[i+2]=U_value[(i+2)*3+2];
			Ux_value[i+3]=U_value[(i+3)*3];Uy_value[i+3]=U_value[(i+3)*3+1];Uz_value[i+3]=U_value[(i+3)*3+2];
		}
		for(int i=(num/4)*4;i<num;i++){
			Ux_value[i]=U_value[i*3];Uy_value[i]=U_value[i*3+1];Uz_value[i]=U_value[i*3+2];
		}
		for(int i=0;i<(num/4)*4;i+=4){
			qx_value[i]=q_value[i*3];qy_value[i]=q_value[i*3+1];qz_value[i]=q_value[i*3+2];
			qx_value[i+1]=q_value[(i+1)*3];qy_value[i+1]=q_value[(i+1)*3+1];qz_value[i+1]=q_value[(i+1)*3+2];
			qx_value[i+2]=q_value[(i+2)*3];qy_value[i+2]=q_value[(i+2)*3+1];qz_value[i+2]=q_value[(i+2)*3+2];
			qx_value[i+3]=q_value[(i+3)*3];qy_value[i+3]=q_value[(i+3)*3+1];qz_value[i+3]=q_value[(i+3)*3+2];
		}
		for(int i=(num/4)*4;i<num;i++){
			qx_value[i]=q_value[i*3];qy_value[i]=q_value[i*3+1];qz_value[i]=q_value[i*3+2];
		}*/
		// for(int i=0;i<num;i++){
		// 	qx_value[i]=q_value[i*3];
		// 	qy_value[i]=q_value[i*3+1];
		// 	qz_value[i]=q_value[i*3+2];
		// }
		// time_end = rpcc();
		// if(my_id==0&&step==0)
		// 	printf("updateGHsurf time slave 1= %.10lf ms",((double)(time_end-time_sta)));
    	xstep = 0;
		xstnum = 0;
		xnum = 0;
		xnumtemp = DVsize;
		for(dvi=0;dvi<DVsize;dvi++){
			index = dvi%2;
			next = (dvi+1)%2;
			last = (dvi-1)%2;
#ifndef GHSURF_REG_COM
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
			// time_sta1 = rpcc();
			//和寄存器通信取一个
			// if(dvi==0){//第一次
				get_reply_eq[index] = 0;
				athread_get(PE_MODE , slave.gSurf_value[dvi]+stid+stnum, gSurf_value[index], sizeof(double)*num, &get_reply_eq[index] , 0, 0, 0);
				athread_get(PE_MODE , slave.hSurf_value[dvi]+stid+stnum, hSurf_value[index], sizeof(double)*num, &get_reply_eq[index] , 0, 0, 0);
			// }
			// if(dvi!=DVsize-1){//除最后一次
			// 	get_reply_eq[next] = 0;
			// 	athread_get(PE_MODE , slave.gSurf_value[dvi+1]+stid+stnum, gSurf_value[next], sizeof(double)*num, &get_reply_eq[next] , 0, 0, 0);
			// 	athread_get(PE_MODE , slave.hSurf_value[dvi+1]+stid+stnum, hSurf_value[next], sizeof(double)*num, &get_reply_eq[next] , 0, 0, 0);
			// }
			while (get_reply_eq[index] != 2);
			// time_end1 = rpcc();
			// if(my_id==0&&step==0&&dvi==0)
			// 	printf("updateGHsurf time slave 2= %.10lf ms",((double)(time_end1-time_sta1)*1000/CLOCKRATE));
			// pe_get(slave.gSurf_value[dvi]+stid+stnum, gSurf_value, sizeof(double)*num);
			// pe_get(slave.hSurf_value[dvi]+stid+stnum, hSurf_value, sizeof(double)*num);
			// dma_syn();
			// time_sta3 = rpcc();
#ifndef GHSURF_USE_SIMD
// #ifdef GHSURF_USE_SIMD
//向量不向量的区别不大
			// for(facei=0;facei<num;facei++)
			// {
			// 	relaxFactor = h/(2*tauSurf_value[facei] + h);
			// 	gSurf_value[index][facei] = (1.0 - relaxFactor)*gSurf_value[index][facei] ;
			// 	hSurf_value[index][facei] = (1.0 - relaxFactor)*hSurf_value[index][facei] ;
			// 	cSqrByRT_value = ((U_value[facei*3]-xii_x)*(U_value[facei*3]-xii_x)+(U_value[facei*3+1]-xii_y)*(U_value[facei*3+1]-xii_y)+(U_value[facei*3+2]-xii_z)*(U_value[facei*3+2]-xii_z))/(R_value*T_value[facei]);
			// 	cqBy5pRT_value = ((xii_x - U_value[facei*3])*q_value[facei*3]+(xii_y - U_value[facei*3+1])*q_value[facei*3+1]+(xii_z - U_value[facei*3+2])*q_value[facei*3+2])/(5.0*rho_value[facei]*R_value*T_value[facei]*R_value*T_value[facei]);
			// 	// gEqBGK_value  = (relaxFactor*rho_value[facei])/pow(sqrt(2.0*pi*R_value*T_value[facei]),D)*exp(-cSqrByRT_value/2.0)/pow(vUnit_value, 3-D);
			// 	// gEqBGK_value  = (relaxFactor*rho_value[facei])/pow(sqrt(2.0*pi*R_value*T_value[facei]),D)*texp(-cSqrByRT_value/2.0);
			// 	gEqBGK_value  = (relaxFactor*rho_value[facei])/pow(sqrt(2.0*pi*R_value*T_value[facei]),D)*(1/expd(cSqrByRT_value/2.0));
			// 	//gEqBGK_value  = (relaxFactor*rho_value[facei])/pow(2.0*pi*R_value*T_value[facei],D/2.0)*exp(-cSqrByRT_value/2.0);
			// 	gSurf_value[index][facei] += ( 1.0 + (1.0 - Pr_value)*cqBy5pRT_value*(cSqrByRT_value - D - 2.0) )*gEqBGK_value;
			// 	hSurf_value[index][facei] += ( (K + 3.0 - D) + (1.0 - Pr_value)*cqBy5pRT_value*((cSqrByRT_value - D)*(K + 3.0 - D) - 2*K) )*gEqBGK_value*R_value*T_value[facei];				
				
			// }
			double temp = 0;
			double gEqBGK1 = 0;
			double gSurf1 = 0;
			double hSurf1 = 0;
			for(facei=0;facei<num;facei++)
			{
				relaxFactor = h/(2*tauSurf_value[facei] + h);
				gSurf_value[index][facei] = (1.0 - relaxFactor)*gSurf_value[index][facei] ;
				hSurf_value[index][facei] = (1.0 - relaxFactor)*hSurf_value[index][facei] ;
				cSqrByRT_value = ((U_value[facei*3]-xii_x)*(U_value[facei*3]-xii_x)+(U_value[facei*3+1]-xii_y)*(U_value[facei*3+1]-xii_y)+(U_value[facei*3+2]-xii_z)*(U_value[facei*3+2]-xii_z))/(R_value*T_value[facei]);
				cqBy5pRT_value = ((xii_x - U_value[facei*3])*q_value[facei*3]+(xii_y - U_value[facei*3+1])*q_value[facei*3+1]+(xii_z - U_value[facei*3+2])*q_value[facei*3+2])/(5.0*rho_value[facei]*R_value*T_value[facei]*R_value*T_value[facei]);
				
#ifndef GHSURF_TRAN_FUN
				// gEqBGK_value  = (relaxFactor*rho_value[facei])/pow(sqrt(2.0*pi*R_value*T_value[facei]),D)*exp(-cSqrByRT_value/2.0)/pow(vUnit_value, 3-D);
				gEqBGK_value  = (relaxFactor*rho_value[facei])/pow(sqrt(2.0*pi*R_value*T_value[facei]),D)*exp(-cSqrByRT_value/2.0);

#else
				gEqBGK_value  = (relaxFactor*rho_value[facei])*(exp(-cSqrByRT_value/2.0));
				// gEqBGK_value  = (relaxFactor*rho_value[facei])*(1/expd(cSqrByRT_value/2.0));//会产生误差
				temp = sqrt(2.0*pi*R_value*T_value[facei]);
				// gEqBGK_value /= pow(temp,D);
				// gEqBGK_value /= temp;
				// gEqBGK_value /= temp;
				// gEqBGK_value /= temp;
				double mantissa, number=1.0;
   				int exponent;//指数
   				// number = 8.0;
   				mantissa = frexp(temp, &exponent);
				for(int w = 0;w<exponent;w++){
					number *= 2;
				}
				// number = pow2n(exponent);
				// udi_t    epart;//和循环计算差距不大
				// epart.s.i0 = 0;
				// epart.s.i1 = ((exponent) + 1023) << 20;
				for(int q=0;q<D;q++){
					// gEqBGK_value /= (epart.f*mantissa);
					gEqBGK_value /= (number*mantissa);
				}
				//gEqBGK_value /= (number*number*number*mantissa*mantissa*mantissa);
				// gEqBGK1 /= powd(temp,D);
				// gEqBGK1 /= (temp*temp*temp);
				// if((gEqBGK_value - gEqBGK1 > 0.1)||(gEqBGK_value - gEqBGK1 < -0.1))
				// if(gEqBGK_value!=gEqBGK1)
					// printf("%d",D);
				// gEqBGK_value = gEqBGK1;
				//gEqBGK_value  = (relaxFactor*rho_value[facei])/pow(2.0*pi*R_value*T_value[facei],D/2.0)*exp(-cSqrByRT_value/2.0);
				// gSurf1 = gSurf_value[index][facei] + ( 1.0 + (1.0 - Pr_value)*cqBy5pRT_value*(cSqrByRT_value - D - 2.0) )*gEqBGK1;
				// hSurf1 = hSurf_value[index][facei] + ( (K + 3.0 - D) + (1.0 - Pr_value)*cqBy5pRT_value*((cSqrByRT_value - D)*(K + 3.0 - D) - 2*K) )*gEqBGK1*R_value*T_value[facei];				
#endif				
				gSurf_value[index][facei] += ( 1.0 + (1.0 - Pr_value)*cqBy5pRT_value*(cSqrByRT_value - D - 2.0) )*gEqBGK_value;
				hSurf_value[index][facei] += ( (K + 3.0 - D) + (1.0 - Pr_value)*cqBy5pRT_value*((cSqrByRT_value - D)*(K + 3.0 - D) - 2*K) )*gEqBGK_value*R_value*T_value[facei];				
				// if((gSurf_value[index][facei] - gSurf1 > 0.000000001)||(gSurf_value[index][facei] - gSurf1 < -0.000000001))
				// 	printf("1");
				// if((hSurf_value[index][facei] - hSurf1 > 0.000000001)||(hSurf_value[index][facei] - hSurf1 < -0.000000001))
				// 	printf("2");

				
			}
#else
			doublev4 Uxv4=0,Uyv4=0,Uzv4=0,Tv4=0,rhov4=0,cSqrByRTv4=0,cqBy5pRTv4=0,gEqBGKv4=0;
			doublev4 qxv4=0,qyv4=0,qzv4=0;
			doublev4 temp1=0,temp2=0,temp3=0,temp4=0,tempv4;
			double cSqrByRT_[4]__attribute__((aligned(128)));
			double cqBy5pRT_[4]__attribute__((aligned(128)));
			double gEqBGK_[4]__attribute__((aligned(128)));
			double relaxFactor_[4]__attribute__((aligned(128)));
			double temp[4]__attribute__((aligned(128)));
			doublev4 relaxFactorv4=0,tauSurfv4=0,gSurfv4=0,hSurfv4=0;
			for(facei = 0;facei<(num/4)*4;facei+=4){
				simd_load(tauSurfv4, &(tauSurf_value[facei]));
				relaxFactorv4 = h/(2.0*tauSurfv4 + h);
				simd_load(gSurfv4, &(gSurf_value[index][facei]));
				gSurfv4 = (1.0 - relaxFactorv4)*gSurfv4;
				simd_load(hSurfv4, &(hSurf_value[index][facei]));
				hSurfv4 = (1.0 - relaxFactorv4)*hSurfv4;
				
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
				/*Uxv4=simd_set_doublev4(U_value[facei*3],U_value[(facei+1)*3],U_value[(facei+2)*3],U_value[(facei+3)*3]);
				Uyv4=simd_set_doublev4(U_value[facei*3+1],U_value[(facei+1)*3+1],U_value[(facei+2)*3+1],U_value[(facei+3)*3+1]);
				Uzv4=simd_set_doublev4(U_value[facei*3+2],U_value[(facei+1)*3+2],U_value[(facei+2)*3+2],U_value[(facei+3)*3+2]);
				*/
				simd_load(Tv4, &(T_value[facei]));
				cSqrByRTv4 = ((Uxv4-xii_x)*(Uxv4-xii_x)+(Uyv4-xii_y)*(Uyv4-xii_y)+(Uzv4-xii_z)*(Uzv4-xii_z))/(R_value*Tv4);
				// qxv4=simd_set_doublev4(q_value[facei*3],q_value[(facei+1)*3],q_value[(facei+2)*3],q_value[(facei+3)*3]);
				// qyv4=simd_set_doublev4(q_value[facei*3+1],q_value[(facei+1)*3+1],q_value[(facei+2)*3+1],q_value[(facei+3)*3+1]);
				// qzv4=simd_set_doublev4(q_value[facei*3+2],q_value[(facei+1)*3+2],q_value[(facei+2)*3+2],q_value[(facei+3)*3+2]);
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
				cqBy5pRTv4 = ((xii_x - Uxv4)*qxv4+(xii_y - Uyv4)*qyv4+(xii_z - Uzv4)*qzv4)/(5.0*rhov4*R_value*Tv4*R_value*Tv4);

				// cqBy5pRT_temp = (5.0*rhov4*R_value*Tv4*R_value*Tv4);
#ifndef GHSURF_TRAN_FUN
				/******************************************************************/
				simd_store(cSqrByRTv4, &(cSqrByRT_[0]));
				simd_store(relaxFactorv4, &(relaxFactor_[0]));
				// gEqBGK_[0]  = (relaxFactor_[0]*rho_value[facei])/pow(sqrt(2.0*pi*R_value*T_value[facei]),D)*(1/expd(cSqrByRT_[0]/2.0));
				// gEqBGK_[1]  = (relaxFactor_[1]*rho_value[facei+1])/pow(sqrt(2.0*pi*R_value*T_value[facei+1]),D)*(1/expd(cSqrByRT_[1]/2.0));
				// gEqBGK_[2]  = (relaxFactor_[2]*rho_value[facei+2])/pow(sqrt(2.0*pi*R_value*T_value[facei+2]),D)*(1/expd(cSqrByRT_[2]/2.0));
				// gEqBGK_[3]  = (relaxFactor_[3]*rho_value[facei+3])/pow(sqrt(2.0*pi*R_value*T_value[facei+3]),D)*(1/expd(cSqrByRT_[3]/2.0));
				gEqBGK_[0]  = (relaxFactor_[0]*rho_value[facei])/pow(sqrt(2.0*pi*R_value*T_value[facei]),D)*(exp(-cSqrByRT_[0]/2.0));
				gEqBGK_[1]  = (relaxFactor_[1]*rho_value[facei+1])/pow(sqrt(2.0*pi*R_value*T_value[facei+1]),D)*(exp(-cSqrByRT_[1]/2.0));
				gEqBGK_[2]  = (relaxFactor_[2]*rho_value[facei+2])/pow(sqrt(2.0*pi*R_value*T_value[facei+2]),D)*(exp(-cSqrByRT_[2]/2.0));
				gEqBGK_[3]  = (relaxFactor_[3]*rho_value[facei+3])/pow(sqrt(2.0*pi*R_value*T_value[facei+3]),D)*(exp(-cSqrByRT_[3]/2.0));
				
				simd_load(gEqBGKv4, &(gEqBGK_[0]));
#else
				simd_store(cSqrByRTv4, &(cSqrByRT_[0]));
				simd_store(relaxFactorv4, &(relaxFactor_[0]));
				gEqBGK_[0]  = (relaxFactor_[0]*rho_value[facei]);
				gEqBGK_[1]  = (relaxFactor_[1]*rho_value[facei+1]);
				gEqBGK_[2]  = (relaxFactor_[2]*rho_value[facei+2]);
				gEqBGK_[3]  = (relaxFactor_[3]*rho_value[facei+3]);
				
				// gEqBGK_[0]  = (relaxFactor_[0]*rho_value[facei])*(1/expd(cSqrByRT_[0]/2.0));
				// gEqBGK_[1]  = (relaxFactor_[1]*rho_value[facei+1])*(1/expd(cSqrByRT_[1]/2.0));
				// gEqBGK_[2]  = (relaxFactor_[2]*rho_value[facei+2])*(1/expd(cSqrByRT_[2]/2.0));
				// gEqBGK_[3]  = (relaxFactor_[3]*rho_value[facei+3])*(1/exp(cSqrByRT_[3]/2.0));
				simd_load(gEqBGKv4, &(gEqBGK_[0]));

				// gEqBGK_[0]  = (relaxFactor_[0]*rho_value[facei])*exp(-cSqrByRT_[0]/2.0);
				// gEqBGK_[1]  = (relaxFactor_[1]*rho_value[facei+1])*exp(-cSqrByRT_[1]/2.0);
				// gEqBGK_[2]  = (relaxFactor_[2]*rho_value[facei+2])*exp(-cSqrByRT_[2]/2.0);
				// gEqBGK_[3]  = (relaxFactor_[3]*rho_value[facei+3])*exp(-cSqrByRT_[3]/2.0);

				// gEqBGKv4 = relaxFactorv4 * rhov4 * (1/simd_expd(cSqrByRTv4/2.0));
				
				temp[0] = sqrt(2.0*pi*R_value*T_value[facei]);
				temp[1] = sqrt(2.0*pi*R_value*T_value[facei+1]);
				temp[2] = sqrt(2.0*pi*R_value*T_value[facei+2]);
				temp[3] = sqrt(2.0*pi*R_value*T_value[facei+3]);
				// temp[0] = 1/rsqrtd(2.0*pi*R_value*T_value[facei]);
				// temp[1] = 1/rsqrtd(2.0*pi*R_value*T_value[facei+1]);
				// temp[2] = 1/rsqrtd(2.0*pi*R_value*T_value[facei+2]);
				// temp[3] = 1/rsqrtd(2.0*pi*R_value*T_value[facei+3]);
				// tempv4 = simd_vsqrtd(2.0*pi*R_value*Tv4);
				/*********************************************************/
				// double mantissa[4]__attribute__((aligned(128)));
				// doublev4 numberv4,mantissav4;
   				// long exponent[4]__attribute__((aligned(128)));//指数
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

#endif
				gSurfv4 += ( 1.0 + (1.0 - Pr_value)*cqBy5pRTv4*(cSqrByRTv4 - D - 2.0) )*gEqBGKv4;
				simd_store(gSurfv4, &(gSurf_value[index][facei]));
				hSurfv4 += ((K + 3.0 - D) + (1.0 - Pr_value)*cqBy5pRTv4*((cSqrByRTv4 - D)*(K + 3.0 - D) - K - K) )*gEqBGKv4*R_value*Tv4;
				simd_store(hSurfv4, &(hSurf_value[index][facei]));
			
			}
			for(facei=(num/4)*4;facei<num;facei++)
			{
				relaxFactor = h/(2*tauSurf_value[facei] + h);
				gSurf_value[index][facei] = (1.0 - relaxFactor)*gSurf_value[index][facei] ;
				hSurf_value[index][facei] = (1.0 - relaxFactor)*hSurf_value[index][facei] ;
				// gSurf_value[facei] = (1.0 - relaxFactor)*gSurf_value[facei] ;
				// hSurf_value[facei] = (1.0 - relaxFactor)*hSurf_value[facei] ;
				cSqrByRT_value = ((U_value[facei*3]-xii_x)*(U_value[facei*3]-xii_x)+(U_value[facei*3+1]-xii_y)*(U_value[facei*3+1]-xii_y)+(U_value[facei*3+2]-xii_z)*(U_value[facei*3+2]-xii_z))/(R_value*T_value[facei]);
				cqBy5pRT_value = ((xii_x - U_value[facei*3])*q_value[facei*3]+(xii_y - U_value[facei*3+1])*q_value[facei*3+1]+(xii_z - U_value[facei*3+2])*q_value[facei*3+2])/(5.0*rho_value[facei]*R_value*T_value[facei]*R_value*T_value[facei]);
				gEqBGK_value  = (relaxFactor*rho_value[facei])/pow(sqrt(2.0*pi*R_value*T_value[facei]),D)*(exp(-cSqrByRT_value/2.0));
				// gEqBGK_value  = (relaxFactor*rho_value[facei])/pow(sqrt(2.0*pi*R_value*T_value[facei]),D)*(1/expd(cSqrByRT_value/2.0));
				
				gSurf_value[index][facei] += ( 1.0 + (1.0 - Pr_value)*cqBy5pRT_value*(cSqrByRT_value - D - 2.0) )*gEqBGK_value;
				hSurf_value[index][facei] += ( (K + 3.0 - D) + (1.0 - Pr_value)*cqBy5pRT_value*((cSqrByRT_value - D)*(K + 3.0 - D) - 2*K) )*gEqBGK_value*R_value*T_value[facei];				
				// gSurf_value[facei] += ( 1.0 + (1.0 - Pr_value)*cqBy5pRT_value*(cSqrByRT_value - D - 2.0) )*gEqBGK_value;
				// hSurf_value[facei] += ( (K + 3.0 - D) + (1.0 - Pr_value)*cqBy5pRT_value*((cSqrByRT_value - D)*(K + 3.0 - D) - 2*K) )*gEqBGK_value*R_value*T_value[facei];				
				
			}
#endif
			
			// time_end3 = rpcc();
			// if(my_id==0&&step==0&&dvi==0)
			// 	printf("updateGHsurf time slave 3= %.10lf ms",((double)(time_end3-time_sta3)));
			// printf("updateGHsurf time slave 3= %.10lf ms",((double)(time_end3-time_sta3)*1000/CLOCKRATE));
			// time_sta2 = rpcc();
			// pe_put(slave.gSurf_value[dvi]+stid+stnum, gSurf_value, sizeof(double)*num);
			// pe_put(slave.hSurf_value[dvi]+stid+stnum, hSurf_value, sizeof(double)*num);
			// dma_syn();
			put_reply_eq[index] = 0; 
			athread_put(PE_MODE,gSurf_value[index],slave.gSurf_value[dvi]+stid+stnum,sizeof(double)*num,&put_reply_eq[index],0,0);                 
			athread_put(PE_MODE,hSurf_value[index],slave.hSurf_value[dvi]+stid+stnum,sizeof(double)*num,&put_reply_eq[index],0,0);                 
			// if(dvi!=0){
			// 	while(put_reply_eq[last]!=2);
			// }
			// if(dvi==DVsize-1){//最后一轮
				while(put_reply_eq[index]!=2);
			// }
			// time_end2 = rpcc();
			// if(my_id==0&&step==0&&dvi==0)
			// 	printf("updateGHsurf time slave 4= %.10lf ms",((double)(time_end2-time_sta2)));
			
		}
		step++;
	}
	time2 = rpcc();
	// penv_slave0_cycle_count(&ic1);
	penv_slave2_gld_count(&ic1);
	// if (my_id == 0){
	// 	// printf("perf time = %lf\n",(double)(ic1 - ic)*1000/CLOCKRATE);
	// 	// printf("rpcc time = %lf\n",(double)(time2 - time1)*1000/CLOCKRATE);
	// 	printf("perf time = %lu\n",(ic1 - ic));
	// }
        

}
