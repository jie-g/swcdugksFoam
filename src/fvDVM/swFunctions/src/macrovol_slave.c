#include <stdio.h>
#include "slave.h"
#include <math.h>
#include "dma_macros.h"
#include "slave_para.h"
#include <simd.h>

#include "para.h"
// #include <swperf.h>
// #define MACROVOL_USE_SIMD
// #define MACROVOL_TRAN_FUN
// #define MACROVOL_REG_COM

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

static inline double expd(double x){
    int n = 11;
    double log2e = 1.44269504088896338700;
    double xlog = x * log2e;
    // int64 xint = floor(xlog);//p
	unsigned short xint = (int)(xlog);//p
    double xdeci = xlog - xint;//q
    // printf("%lf,%d,%lf\n",xlog,xint, xdeci);
    int i = 0;
    double exp2p = 1;
    for(i=0;i<xint;i++){
        exp2p *= 2;
    }
	// unsigned short exp2p = (unsigned short)(1 << xint);
	// double exp2p = (double)(1 << xint);//计算2的x次方
	// double exp2p = xint<30?exp2x[xint]:exp2_(xint);//更慢了
	// double exp2p = pow(2,xint);
    double ln2 = 0.69314718055994528623;
    double qln2 = xdeci * ln2; 
    // double exp2q = exponential(qln2, 11);
	// double exp2q = 1.0;
    // double term = 1.0;
    // for (int i = 1; i <= 11; i++) {
    //     term *= qln2 / i;
    //     exp2q += term;
    // }
	double exp2q = 1+qln2*(1+qln2*(0.5+qln2*(0.16666666666666666667+qln2*(0.04166666666666666667+qln2*(0.00833333333333333333+ \
    qln2*(0.00138888888888888889+qln2*(0.00019841269841269841+qln2*(0.00002480158730158730+qln2*(0.00000275573192239859+qln2*0.00000027557319223986)))))))));
    
    double result = exp2p * exp2q;
    // printf("%lf,%lf,%.20lf,%lf\n",qln2,exp2p, exp2q, result);
    // printf("%.20lf\n",result);
    // double right = exp(x);
    // printf("right value :%.20lf\n", right);
    return result;
}
static inline doublev4 simd_expd(doublev4 x){
    double log2e = 1.44269504088896338700;
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
    double ln2 = 0.69314718055994528623;
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
	doublev4 qln2_2 = qln2*qln2;
	//奇数项
	doublev4 exp2q1 = 1.0+qln2_2*(0.5+qln2_2*(0.04166666666666666667+qln2_2*(0.00138888888888888889 \
	+qln2_2*(0.00002480158730158730+qln2_2*(0.00000027557319223986+qln2_2*0.00000000208767569878)))));
	//偶数项
	doublev4 exp2q2 = qln2*(1.0+qln2_2*(0.16666666666666666667+qln2_2*(0.00833333333333333333+ \
	qln2_2*(0.00019841269841269841+qln2_2*(0.00000275573192239859+qln2_2*(0.00000002505210838544+qln2_2*0.00000000016059043836))))));
    doublev4 result = exp2p * (exp2q1+exp2q2);
	// doublev4 result = exp2p * exp2q;
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

void Func_macrovol(MacroVol *master) {
	int N_slave = 448;
	MacroVol slave;
	// unsigned long get_reply,put_reply; 
    int my_id; 
    my_id = athread_get_id(-1);  //获取从核的id号
    get_reply = 0;
    athread_get(PE_MODE , master, &slave, sizeof(MacroVol), &get_reply , 0, 0, 0);
    while (get_reply != 1);//等待数据加载完成
    // int64 ownersize1 = slave.ownersize;//4*8划分下  1740
	int64 ownersize = master->ownersize;
	int64 DVsize = master->DVsize;
	// Load Balance
	int load = ownersize / NTHREAD;
	int rest = ownersize % NTHREAD;
	int stid = (my_id < rest) ? my_id * (load + 1) : my_id * load + rest;
	if (my_id < rest) load++;
	int edid = stid + load;
	int numtemp = edid - stid;
	//Define variable
	double Sf[N_slave*3];
	double rhoflux_value_own[N_slave];
	double rhoflux_value_nei[N_slave];
	double rhouflux_value_own[N_slave*3];
	double rhouflux_value_nei[N_slave*3];
	double rhoeflux_value_own[N_slave];
	double rhoeflux_value_nei[N_slave];
	double V_own[N_slave];
	double V_nei[N_slave];
	double gSurf_value[N_slave];
	double hSurf_value[N_slave];
	double dt = master->dt;
	// double xii_x = master->xii_x;
	// double xii_y = master->xii_y;
	// double xii_z = master->xii_z;
	double xii_x,xii_y,xii_z;
	// double weight_value = master->weight_value;
	double weight_value;
	int facei = 0;
	double xii_Sf = 0;
	// unsigned long time1 = 0;
	// unsigned long time_sta = 0;
	// unsigned long time_end = 0;
	int num = 0;
	int step = 0;
	int stnum = 0;
	int dvid=0;
	while(numtemp>0){
		num = (numtemp < N_slave) ? numtemp : N_slave;
		numtemp = numtemp - N_slave;
		stnum = step * N_slave;
		get_reply = 0;
		athread_get(PE_MODE , slave.Sf+stid*3+stnum*3, Sf, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.rhoflux_value_own+stid+stnum, rhoflux_value_own, sizeof(double)*num, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.rhoflux_value_nei+stid+stnum, rhoflux_value_nei, sizeof(double)*num, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.rhouflux_value_own+stid*3+stnum*3, rhouflux_value_own, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.rhouflux_value_nei+stid*3+stnum*3, rhouflux_value_nei, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.rhoeflux_value_own+stid+stnum, rhoeflux_value_own, sizeof(double)*num, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.rhoeflux_value_nei+stid+stnum, rhoeflux_value_nei, sizeof(double)*num, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.V_own+stid+stnum, V_own, sizeof(double)*num, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.V_nei+stid+stnum, V_nei, sizeof(double)*num, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.gSurf_value+stid+stnum, gSurf_value, sizeof(double)*num, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.hSurf_value+stid+stnum, hSurf_value, sizeof(double)*num, &get_reply , 0, 0, 0);
		while (get_reply != 3);
		memset(rhoflux_value_own,0,sizeof(double)*N_slave);
		memset(rhoflux_value_nei,0,sizeof(double)*N_slave);
		memset(rhouflux_value_own,0,sizeof(double)*N_slave*3);
		memset(rhouflux_value_nei,0,sizeof(double)*N_slave*3);
		memset(rhoeflux_value_own,0,sizeof(double)*N_slave);
		memset(rhoeflux_value_nei,0,sizeof(double)*N_slave);

		for(dvid=0;dvid<DVsize;dvid++){
			xii_x = slave.xi_value[dvid*3];
			xii_y = slave.xi_value[dvid*3+1];
			xii_z = slave.xi_value[dvid*3+2];
			weight_value = slave.weight_value[dvid];
			get_reply = 0;
			athread_get(PE_MODE , slave.gSurf_value[dvid]+stid+stnum, gSurf_value, sizeof(double)*num, &get_reply , 0, 0, 0);
			athread_get(PE_MODE , slave.hSurf_value[dvid]+stid+stnum, hSurf_value, sizeof(double)*num, &get_reply , 0, 0, 0);
			while (get_reply != 2);

			double rho_value_temp = 0;
			double e_value_temp = 0;
			double u_value_x = 0 , u_value_y = 0 , u_value_z = 0 ;
#ifndef MACROVOL_USE_SIMD
			for(facei=0;facei<num;facei++)
			{
				xii_Sf = xii_x * Sf[facei*3] + xii_y * Sf[facei*3+1] + xii_z * Sf[facei*3+2];
			
				rho_value_temp =  weight_value * gSurf_value[facei];
				rhoflux_value_own[facei] -= ((xii_Sf) * rho_value_temp * dt / V_own[facei]);
				rhoflux_value_nei[facei] += ((xii_Sf) * rho_value_temp * dt / V_nei[facei]);
			
				u_value_x = weight_value * gSurf_value[facei] * xii_x;
				u_value_y = weight_value * gSurf_value[facei] * xii_y;
				u_value_z = weight_value * gSurf_value[facei] * xii_z;
				rhouflux_value_own[facei*3] -= ((xii_Sf) * u_value_x * dt / V_own[facei]);
				rhouflux_value_nei[facei*3] += ((xii_Sf) * u_value_x * dt / V_nei[facei]);
				rhouflux_value_own[facei*3+1] -= ((xii_Sf) * u_value_y * dt / V_own[facei]);
				rhouflux_value_nei[facei*3+1] += ((xii_Sf) * u_value_y * dt / V_nei[facei]);
				rhouflux_value_own[facei*3+2] -= ((xii_Sf) * u_value_z * dt / V_own[facei]);
				rhouflux_value_nei[facei*3+2] += ((xii_Sf) * u_value_z * dt / V_nei[facei]);
				e_value_temp = 0.5 * weight_value * (gSurf_value[facei] * (xii_x * xii_x + xii_y * xii_y + xii_z * xii_z) + hSurf_value[facei]);
				rhoeflux_value_own[facei] -= ((xii_Sf) * e_value_temp * dt / V_own[facei]);
				rhoeflux_value_nei[facei] += ((xii_Sf) * e_value_temp * dt / V_nei[facei]);
			}
#else
			double xii_Sfval[4]__attribute__((aligned(128)));
			double e_value_tempval[4]__attribute__((aligned(128)));
			// for(facei=0;facei<num;facei++)
			for(facei = 0;facei<(num/4)*4;facei+=4)
			{
				xii_Sfval[0] = xii_x * Sf[facei*3] + xii_y * Sf[facei*3+1] + xii_z * Sf[facei*3+2];
				xii_Sfval[1] = xii_x * Sf[(facei+1)*3] + xii_y * Sf[(facei+1)*3+1] + xii_z * Sf[(facei+1)*3+2];
				xii_Sfval[2] = xii_x * Sf[(facei+2)*3] + xii_y * Sf[(facei+2)*3+1] + xii_z * Sf[(facei+2)*3+2];
				xii_Sfval[3] = xii_x * Sf[(facei+3)*3] + xii_y * Sf[(facei+3)*3+1] + xii_z * Sf[(facei+3)*3+2];
			
				// rho_value_temp =  weight_value * gSurf_value[facei];
				rhoflux_value_own[facei] -= ((xii_Sfval[0]) * weight_value * gSurf_value[facei] * dt / V_own[facei]);
				rhoflux_value_own[facei+1] -= ((xii_Sfval[1]) * weight_value * gSurf_value[facei+1] * dt / V_own[facei+1]);
				rhoflux_value_own[facei+2] -= ((xii_Sfval[2]) * weight_value * gSurf_value[facei+2] * dt / V_own[facei+2]);
				rhoflux_value_own[facei+3] -= ((xii_Sfval[3]) * weight_value * gSurf_value[facei+3] * dt / V_own[facei+3]);

				rhoflux_value_nei[facei] += ((xii_Sfval[0]) * weight_value * gSurf_value[facei] * dt / V_nei[facei]);
				rhoflux_value_nei[facei+1] += ((xii_Sfval[1]) * weight_value * gSurf_value[facei+1] * dt / V_nei[facei+1]);
				rhoflux_value_nei[facei+2] += ((xii_Sfval[2]) * weight_value * gSurf_value[facei+2] * dt / V_nei[facei+2]);
				rhoflux_value_nei[facei+3] += ((xii_Sfval[3]) * weight_value * gSurf_value[facei+3] * dt / V_nei[facei+3]);
			
				// u_value_x = weight_value * gSurf_value[facei] * xii_x;
				// u_value_y = weight_value * gSurf_value[facei] * xii_y;
				// u_value_z = weight_value * gSurf_value[facei] * xii_z;
				rhouflux_value_own[facei*3] -= ((xii_Sfval[0]) * weight_value * gSurf_value[facei] * xii_x * dt / V_own[facei]);
				rhouflux_value_nei[facei*3] += ((xii_Sfval[0]) * weight_value * gSurf_value[facei] * xii_x * dt / V_nei[facei]);
				rhouflux_value_own[facei*3+1] -= ((xii_Sfval[0]) * weight_value * gSurf_value[facei] * xii_y * dt / V_own[facei]);
				rhouflux_value_nei[facei*3+1] += ((xii_Sfval[0]) * weight_value * gSurf_value[facei] * xii_y * dt / V_nei[facei]);
				rhouflux_value_own[facei*3+2] -= ((xii_Sfval[0]) * weight_value * gSurf_value[facei] * xii_z * dt / V_own[facei]);
				rhouflux_value_nei[facei*3+2] += ((xii_Sfval[0]) * weight_value * gSurf_value[facei] * xii_z * dt / V_nei[facei]);
				
				rhouflux_value_own[(facei+1)*3] -= ((xii_Sfval[1]) * weight_value * gSurf_value[facei+1] * xii_x * dt / V_own[facei+1]);
				rhouflux_value_nei[(facei+1)*3] += ((xii_Sfval[1]) * weight_value * gSurf_value[facei+1] * xii_x * dt / V_nei[facei+1]);
				rhouflux_value_own[(facei+1)*3+1] -= ((xii_Sfval[1]) * weight_value * gSurf_value[facei+1] * xii_y * dt / V_own[facei+1]);
				rhouflux_value_nei[(facei+1)*3+1] += ((xii_Sfval[1]) * weight_value * gSurf_value[facei+1] * xii_y * dt / V_nei[facei+1]);
				rhouflux_value_own[(facei+1)*3+2] -= ((xii_Sfval[1]) * weight_value * gSurf_value[facei+1] * xii_z * dt / V_own[facei+1]);
				rhouflux_value_nei[(facei+1)*3+2] += ((xii_Sfval[1]) * weight_value * gSurf_value[facei+1] * xii_z * dt / V_nei[facei+1]);
				
				rhouflux_value_own[(facei+2)*3] -= ((xii_Sfval[2]) * weight_value * gSurf_value[facei+2] * xii_x * dt / V_own[facei+2]);
				rhouflux_value_nei[(facei+2)*3] += ((xii_Sfval[2]) * weight_value * gSurf_value[facei+2] * xii_x * dt / V_nei[facei+2]);
				rhouflux_value_own[(facei+2)*3+1] -= ((xii_Sfval[2]) * weight_value * gSurf_value[facei+2] * xii_y * dt / V_own[facei+2]);
				rhouflux_value_nei[(facei+2)*3+1] += ((xii_Sfval[2]) * weight_value * gSurf_value[facei+2] * xii_y * dt / V_nei[facei+2]);
				rhouflux_value_own[(facei+2)*3+2] -= ((xii_Sfval[2]) * weight_value * gSurf_value[facei+2] * xii_z * dt / V_own[facei+2]);
				rhouflux_value_nei[(facei+2)*3+2] += ((xii_Sfval[2]) * weight_value * gSurf_value[facei+2] * xii_z * dt / V_nei[facei+2]);
				
				rhouflux_value_own[(facei+3)*3] -= ((xii_Sfval[3]) * weight_value * gSurf_value[facei+3] * xii_x * dt / V_own[facei+3]);
				rhouflux_value_nei[(facei+3)*3] += ((xii_Sfval[3]) * weight_value * gSurf_value[facei+3] * xii_x * dt / V_nei[facei+3]);
				rhouflux_value_own[(facei+3)*3+1] -= ((xii_Sfval[3]) * weight_value * gSurf_value[facei+3] * xii_y * dt / V_own[facei+3]);
				rhouflux_value_nei[(facei+3)*3+1] += ((xii_Sfval[3]) * weight_value * gSurf_value[facei+3] * xii_y * dt / V_nei[facei+3]);
				rhouflux_value_own[(facei+3)*3+2] -= ((xii_Sfval[3]) * weight_value * gSurf_value[facei+3] * xii_z * dt / V_own[facei+3]);
				rhouflux_value_nei[(facei+3)*3+2] += ((xii_Sfval[3]) * weight_value * gSurf_value[facei+3] * xii_z * dt / V_nei[facei+3]);
				e_value_tempval[0] = 0.5 * weight_value * (gSurf_value[facei] * (xii_x * xii_x + xii_y * xii_y + xii_z * xii_z) + hSurf_value[facei]);
				rhoeflux_value_own[facei] -= ((xii_Sfval[0]) * e_value_tempval[0] * dt / V_own[facei]);
				rhoeflux_value_nei[facei] += ((xii_Sfval[0]) * e_value_tempval[0] * dt / V_nei[facei]);
				e_value_tempval[1] = 0.5 * weight_value * (gSurf_value[facei+1] * (xii_x * xii_x + xii_y * xii_y + xii_z * xii_z) + hSurf_value[facei+1]);
				rhoeflux_value_own[facei+1] -= ((xii_Sfval[1]) * e_value_tempval[1] * dt / V_own[facei+1]);
				rhoeflux_value_nei[facei+1] += ((xii_Sfval[1]) * e_value_tempval[1] * dt / V_nei[facei+1]);
				e_value_tempval[2] = 0.5 * weight_value * (gSurf_value[facei+2] * (xii_x * xii_x + xii_y * xii_y + xii_z * xii_z) + hSurf_value[facei+2]);
				rhoeflux_value_own[facei+2] -= ((xii_Sfval[2]) * e_value_tempval[2] * dt / V_own[facei+2]);
				rhoeflux_value_nei[facei+2] += ((xii_Sfval[2]) * e_value_tempval[2] * dt / V_nei[facei+2]);
				e_value_tempval[3] = 0.5 * weight_value * (gSurf_value[facei+3] * (xii_x * xii_x + xii_y * xii_y + xii_z * xii_z) + hSurf_value[facei+3]);
				rhoeflux_value_own[facei+3] -= ((xii_Sfval[3]) * e_value_tempval[3] * dt / V_own[facei+3]);
				rhoeflux_value_nei[facei+3] += ((xii_Sfval[3]) * e_value_tempval[3] * dt / V_nei[facei+3]);
			}
			for(facei=(num/4)*4;facei<num;facei++)
			{
				xii_Sf = xii_x * Sf[facei*3] + xii_y * Sf[facei*3+1] + xii_z * Sf[facei*3+2];
			
				rho_value_temp =  weight_value * gSurf_value[facei];
				rhoflux_value_own[facei] -= ((xii_Sf) * rho_value_temp * dt / V_own[facei]);
				rhoflux_value_nei[facei] += ((xii_Sf) * rho_value_temp * dt / V_nei[facei]);
			
				u_value_x = weight_value * gSurf_value[facei] * xii_x;
				u_value_y = weight_value * gSurf_value[facei] * xii_y;
				u_value_z = weight_value * gSurf_value[facei] * xii_z;
				rhouflux_value_own[facei*3] -= ((xii_Sf) * u_value_x * dt / V_own[facei]);
				rhouflux_value_nei[facei*3] += ((xii_Sf) * u_value_x * dt / V_nei[facei]);
				rhouflux_value_own[facei*3+1] -= ((xii_Sf) * u_value_y * dt / V_own[facei]);
				rhouflux_value_nei[facei*3+1] += ((xii_Sf) * u_value_y * dt / V_nei[facei]);
				rhouflux_value_own[facei*3+2] -= ((xii_Sf) * u_value_z * dt / V_own[facei]);
				rhouflux_value_nei[facei*3+2] += ((xii_Sf) * u_value_z * dt / V_nei[facei]);
				e_value_temp = 0.5 * weight_value * (gSurf_value[facei] * (xii_x * xii_x + xii_y * xii_y + xii_z * xii_z) + hSurf_value[facei]);
				rhoeflux_value_own[facei] -= ((xii_Sf) * e_value_temp * dt / V_own[facei]);
				rhoeflux_value_nei[facei] += ((xii_Sf) * e_value_temp * dt / V_nei[facei]);
			}
#endif
		}
		put_reply = 0; 
		athread_put(PE_MODE,rhoflux_value_own,slave.rhoflux_value_own+stid+stnum,sizeof(double)*num,&put_reply,0,0);                 
		athread_put(PE_MODE,rhoflux_value_nei,slave.rhoflux_value_nei+stid+stnum,sizeof(double)*num,&put_reply,0,0);                 
		athread_put(PE_MODE,rhouflux_value_own,slave.rhouflux_value_own+stid*3+stnum*3,sizeof(double)*num*3,&put_reply,0,0); 
		athread_put(PE_MODE,rhouflux_value_nei,slave.rhouflux_value_nei+stid*3+stnum*3,sizeof(double)*num*3,&put_reply,0,0);                 
		athread_put(PE_MODE,rhoeflux_value_own,slave.rhoeflux_value_own+stid+stnum,sizeof(double)*num,&put_reply,0,0);                 
		athread_put(PE_MODE,rhoeflux_value_nei,slave.rhoeflux_value_nei+stid+stnum,sizeof(double)*num,&put_reply,0,0);               
		while(put_reply!=6);
		step++;
	}
}

void Func_macroqvol(MacroqVol *macroqvol){
	int N_slave = 448;//int N_slave = 968;
	MacroqVol slave;
	dma_init();
    int my_id; 
    my_id = athread_get_id(-1);  //获取从核的id号
    get_reply = 0;
    athread_get(PE_MODE , macroqvol, &slave, sizeof(MacroqVol), &get_reply , 0, 0, 0);
    while (get_reply != 1);//等待数据加载完成
	int64 qVol_size = macroqvol->qVol_size;	
	int64 DVsize = macroqvol->DVsize;
	// Load Balance
	int load = qVol_size / NTHREAD;
	int rest = qVol_size % NTHREAD;
	int stid = (my_id < rest) ? my_id * (load + 1) : my_id * load + rest;
	if (my_id < rest) load++;
	int edid = stid + load;
	int numtemp = edid - stid;
	// int num = edid - stid;
	//Define variable
	double gTildeVol_value[N_slave];
	double hTildeVol_value[N_slave];
	double rhoVol_value[N_slave];
	// double rhoUvol_value[N_slave*3];
	// double rhoEvol_value[N_slave];
	double qVol_value[N_slave*3];
	// double xi_value[N_slave];
	// double weight_value[N_slave];

	double rhoflux_value[N_slave];
	double rhouflux_value[N_slave*3];
	double rhoeflux_value[N_slave];

	double Tvol_value[N_slave];
	double Uvol_value[N_slave*3];
	double tauVol_value[N_slave];
	double qVol_temp_value[N_slave];
	double R_value = slave.R_value;
	double muRef_value = slave.muRef_value;
	double Tref_value = slave.Tref_value;
	double Pr_value = slave.Pr_value;
	double omega_value = slave.omega_value;
	double dt = slave.dt;
	int64 KInner_value = slave.KInner_value;
	double xii_x,xii_y,xii_z;
	double weight_value;
	// double xii_x = macqsurf->xii_x;
	// double xii_y = macqsurf->xii_y;
	// double xii_z = macqsurf->xii_z;
	// double weight_value = macqsurf->weight_value;
	int facei = 0;
	int num = 0;
	int step = 0;
	int stnum = 0;
	int dvid = 0;
	double magSqr_U = 0;
	double magSqr_c = 0;
	double c_value[3];
	double rhoUvol_x,rhoUvol_y,rhoUvol_z;
	double rhoEvol_value;
	while(numtemp>0){
		num = (numtemp < N_slave) ? numtemp : N_slave;
		numtemp = numtemp - N_slave;
		stnum = step * N_slave;
		//get data
		// get_reply = 0;
		// athread_get(PE_MODE , slave.qVol_value+stid*3+stnum*3, qVol_value, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.rhoUvol_value+stid*3+stnum*3, rhoUvol_value, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.rhoVol_value+stid+stnum, rhoVol_value, sizeof(double)*num, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.rhoEvol_value+stid+stnum, rhoEvol_value, sizeof(double)*num, &get_reply , 0, 0, 0);
		// while (get_reply != 4);
		pe_get(slave.qVol_value+stid*3+stnum*3, qVol_value, sizeof(double)*num*3);
		// pe_get(slave.rhoUvol_value+stid*3+stnum*3, rhoUvol_value, sizeof(double)*num*3);
		pe_get(slave.rhoVol_value+stid+stnum, rhoVol_value, sizeof(double)*num);

		pe_get(slave.rhoflux_value+stid+stnum, rhoflux_value, sizeof(double)*num);
		pe_get(slave.rhouflux_value+stid*3+stnum*3, rhouflux_value, sizeof(double)*num*3);
		pe_get(slave.rhoeflux_value+stid+stnum, rhoeflux_value, sizeof(double)*num);

		pe_get(slave.Uvol_value+stid*3+stnum*3, Uvol_value, sizeof(double)*num*3);
		pe_get(slave.Tvol_value+stid+stnum, Tvol_value, sizeof(double)*num);
		// pe_get(slave.rhoEvol_value+stid+stnum, rhoEvol_value, sizeof(double)*num);
		dma_syn();
		// volVectorField rhoUvol = rhoVol_ * Uvol_;
		// volScalarField rhoEvol = rhoVol_ * (0.5*magSqr(Uvol_) + (KInner_ + 3) / 2.0 * R_ * Tvol_);
        // rhoVol_ += rhoflux_;
		// rhoUvol += rhouflux_;
		// rhoEvol += rhoeflux_;

		// Uvol_ = rhoUvol / rhoVol_;
		// Tvol_ = (rhoEvol - 0.5 * rhoVol_ * magSqr(Uvol_)) / ((KInner_ + 3) / 2.0 * R_ * rhoVol_);
		// tauVol_ = muRef_ * exp(omega_ * log(Tvol_ / Tref_)) / rhoVol_ / Tvol_ / R_;//// updateTau(tauVol_, Tvol_, rhoVol_);
		// volScalarField qVol_temp = 2.0 * tauVol_ / (2.0 * tauVol_ + time_.deltaT() * Pr_);
		
		
		for(facei=0;facei<num;facei++)
		{
			rhoUvol_x = rhoVol_value[facei] * Uvol_value[facei*3];
			rhoUvol_y = rhoVol_value[facei] * Uvol_value[facei*3+1];
			rhoUvol_z = rhoVol_value[facei] * Uvol_value[facei*3+2];
			magSqr_U = Uvol_value[facei*3] * Uvol_value[facei*3] + Uvol_value[facei*3+1] * Uvol_value[facei*3+1] + Uvol_value[facei*3+2] * Uvol_value[facei*3+2];
			rhoEvol_value = rhoVol_value[facei]  * (0.5*magSqr_U + (KInner_value + 3) / 2.0 * R_value * Tvol_value[facei]);
			rhoVol_value[facei] += rhoflux_value[facei];
			rhoUvol_x += rhouflux_value[facei*3];
			rhoUvol_y += rhouflux_value[facei*3+1];
			rhoUvol_z += rhouflux_value[facei*3+2];
			rhoEvol_value += rhoeflux_value[facei];
			Uvol_value[facei*3] = rhoUvol_x / rhoVol_value[facei];
			Uvol_value[facei*3+1] = rhoUvol_y / rhoVol_value[facei];
			Uvol_value[facei*3+2] = rhoUvol_z / rhoVol_value[facei];
			magSqr_U = Uvol_value[facei*3] * Uvol_value[facei*3] + Uvol_value[facei*3+1] * Uvol_value[facei*3+1] + Uvol_value[facei*3+2] * Uvol_value[facei*3+2];
			Tvol_value[facei] = (rhoEvol_value - 0.5 * rhoVol_value[facei] * magSqr_U) / ((KInner_value + 3) / 2.0 * R_value * rhoVol_value[facei]);
			tauVol_value[facei] = muRef_value * exp(omega_value * log(Tvol_value[facei] / Tref_value)) / rhoVol_value[facei] / Tvol_value[facei] / R_value;
			qVol_temp_value[facei] = 2.0 * tauVol_value[facei] / (2.0 * tauVol_value[facei] + dt * Pr_value) ;
		}
		for(dvid=0;dvid<DVsize;dvid++)
		{
			xii_x = slave.xi_value[dvid*3];
			xii_y = slave.xi_value[dvid*3+1];
			xii_z = slave.xi_value[dvid*3+2];
			weight_value = slave.weight_value[dvid];
			// get_reply = 0;
			// athread_get(PE_MODE , slave.gTildeVol_value[dvid]+stid+stnum, gTildeVol_value, sizeof(double)*num, &get_reply , 0, 0, 0);
			// athread_get(PE_MODE , slave.hTildeVol_value[dvid]+stid+stnum, hTildeVol_value, sizeof(double)*num, &get_reply , 0, 0, 0);
			// while (get_reply != 2);
			pe_get(slave.gTildeVol_value[dvid]+stid+stnum, gTildeVol_value, sizeof(double)*num);
			pe_get(slave.hTildeVol_value[dvid]+stid+stnum, hTildeVol_value, sizeof(double)*num);
			dma_syn();
			for(facei=0;facei<num;facei++)
			{
				c_value[0] = xii_x - Uvol_value[facei*3];
				c_value[1] = xii_y - Uvol_value[facei*3+1];
				c_value[2] = xii_z - Uvol_value[facei*3+2];
				// qVol_ += 0.5 * dXiCellSize_ * dv.weight() * c * ( magSqr(c) * dv.gTildeVol() + dv.hTildeVol());
				magSqr_c = c_value[0]*c_value[0]+c_value[1]*c_value[1]+c_value[2]*c_value[2];
				qVol_value[facei*3] += 0.5 * weight_value * c_value[0] * ( magSqr_c * gTildeVol_value[facei] + hTildeVol_value[facei]);
				qVol_value[facei*3+1] += 0.5 * weight_value * c_value[1] * ( magSqr_c * gTildeVol_value[facei] + hTildeVol_value[facei]);
				qVol_value[facei*3+2] += 0.5 * weight_value * c_value[2] * ( magSqr_c * gTildeVol_value[facei] + hTildeVol_value[facei]);
			
			}
		}
		// put_reply = 0; 
		// athread_put(PE_MODE,Tvol_value,slave.Tvol_value+stid+stnum,sizeof(double)*num,&put_reply,0,0);
		// athread_put(PE_MODE,Uvol_value,slave.Uvol_value+stid*3+stnum*3,sizeof(double)*num*3,&put_reply,0,0);
		// athread_put(PE_MODE,tauVol_value,slave.tauVol_value+stid+stnum,sizeof(double)*num,&put_reply,0,0);
		// athread_put(PE_MODE,qVol_temp_value,slave.qVol_temp_value+stid+stnum,sizeof(double)*num,&put_reply,0,0);
    	// athread_put(PE_MODE,qVol_value,slave.qVol_value+stid*3+stnum*3,sizeof(double)*num*3,&put_reply,0,0);                              
  		// while(put_reply!=5);
		pe_put(slave.rhoVol_value+stid+stnum, rhoVol_value, sizeof(double)*num);
		
		pe_put(slave.Tvol_value+stid+stnum, Tvol_value, sizeof(double)*num);
		pe_put(slave.Uvol_value+stid*3+stnum*3, Uvol_value, sizeof(double)*num*3);
		pe_put(slave.tauVol_value+stid+stnum, tauVol_value, sizeof(double)*num);
		pe_put(slave.qVol_temp_value+stid+stnum, qVol_temp_value, sizeof(double)*num);
		pe_put(slave.qVol_value+stid*3+stnum*3, qVol_value, sizeof(double)*num*3);
		dma_syn();
		step++;
	}
}
