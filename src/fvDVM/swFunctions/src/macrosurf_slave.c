#include <stdio.h>
#include "slave.h"
#include <math.h>
#include "dma_macros.h"
#include "slave_para.h"
#include <simd.h>

#include "para.h"
// #include <swperf.h>
// // #define MACROSURF_USE_SIMD
// #define MACROSURF_TRAN_FUN
// #define MACROSURF_REG_COM

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


void Func_macrosurf(MacroSurf *macsurf){
	// int N_mac = 1104;//除8   456
#ifndef MACROSURF_REG_COM
	int N_mac = 768;//1088
#else
	int N_mac = 512;//456  448
#endif

	//
	MacroSurf slave;
    int my_id; 
    my_id = athread_get_id(-1);  //获取从核的id号
    get_reply = 0;
    athread_get(PE_MODE , macsurf, &slave, sizeof(MacroSurf), &get_reply , 0, 0, 0);
    while (get_reply != 1);//等待数据加载完成
    int64 DVsize = macsurf->DVsize;
	int64 ownersize = macsurf->ownersize;	
	// Load Balance
	int load = ownersize / NTHREAD;
	int rest = ownersize % NTHREAD;
	int stid = (my_id < rest) ? my_id * (load + 1) : my_id * load + rest;
	if (my_id < rest) load++;
	int edid = stid + load;
	int numtemp = edid - stid;
	// int num = edid - stid;
	//Define variable
	double rhoSurf_value[N_mac]__attribute__ ((aligned(128)));
	double rhoUsurf_value[N_mac*3]__attribute__ ((aligned(128)));
	// VectorData rhoUsurf_value[N_mac];
	double rhoEsurf_value[N_mac]__attribute__ ((aligned(128)));
	double Usurf_value[N_mac*3]__attribute__ ((aligned(128)));
	double gSurf_value[N_mac]__attribute__ ((aligned(128)));
	double hSurf_value[N_mac]__attribute__ ((aligned(128)));
#ifdef MACROSURF_REG_COM
	double xi_value[N_mac*3]__attribute__ ((aligned(128)));
	double weight_v[N_mac]__attribute__ ((aligned(128)));
#endif
	int xnum = 0;
	int xstep = 0;
	int xstnum = 0;
	int xnumtemp = DVsize;
	doublev4 xv4;
	// double xii_x = macsurf->xii_x;
	// double xii_y = macsurf->xii_y;
	// double xii_z = macsurf->xii_z;
	// double weight_value = macsurf->weight_value;
	// double xi_value[N_mac*3];
	// double weight_[N_mac];
	int facei = 0;
	int dvid = 0;
	int num = 0;
	int step = 0;
	int stnum = 0;
	// int dvnum = 0;
	// int dvstep = 0;
	// int dvnumtemp = 0;
	// int dvstnum = 0;
	double xii_x = 0, xii_y = 0, xii_z = 0, weight_value = 0;
	while(numtemp>0){
		num = (numtemp < N_mac) ? numtemp : N_mac;
		numtemp = numtemp - N_mac;
		stnum = step * N_mac;
		//get data
		get_reply = 0;
		athread_get(PE_MODE , slave.rhoSurf_value+stid+stnum, rhoSurf_value, sizeof(double)*num, &get_reply , 0, 0, 0);
//		athread_get(PE_MODE , slave.rhoUsurf_value+stid*3+stnum*3, rhoUsurf_value, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.Usurf_value+stid*3+stnum*3, Usurf_value, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.rhoUsurf_value+stid+stnum, rhoUsurf_value, sizeof(VectorData)*num, &get_reply , 0, 0, 0);
//		athread_get(PE_MODE , slave.rhoEsurf_value+stid+stnum, rhoEsurf_value, sizeof(double)*num, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.xi_value+stid*3+stnum*3, xi_value, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.weight_value+stid+stnum, weight_value, sizeof(double)*num, &get_reply , 0, 0, 0);
		while (get_reply != 2);
		
		// rhoUsurf = rhoSurf_ * Usurf_;
		// rhoEsurf = rhoSurf_ * magSqr(Usurf_);
		for(facei=0;facei<num;facei++)
		{
			rhoUsurf_value[facei*3] = rhoSurf_value[facei] * Usurf_value[facei*3];
			rhoUsurf_value[facei*3+1] = rhoSurf_value[facei] * Usurf_value[facei*3+1];
			rhoUsurf_value[facei*3+2] = rhoSurf_value[facei] * Usurf_value[facei*3+2];
			rhoEsurf_value[facei] = rhoSurf_value[facei] * (Usurf_value[facei*3]*Usurf_value[facei*3]+Usurf_value[facei*3+1]*Usurf_value[facei*3+1]+Usurf_value[facei*3+2]*Usurf_value[facei*3+2]);
		}


		// dvnumtemp = DVsize;
		// dvstep = 0;
		// while(dvnumtemp>0)
		// {
		// 	dvnum = (dvnumtemp < N_mac) ? dvnumtemp : N_mac;
		// 	dvnumtemp = dvnumtemp - N_mac;
		// 	dvstnum = dvstep * N_mac;
			// get_reply = 0;
			// athread_get(PE_MODE , slave.xi_value+dvstnum*3, xi_value, sizeof(double)*dvnum*3, &get_reply , 0, 0, 0);
			// athread_get(PE_MODE , slave.weight_value+dvstnum, weight_, sizeof(double)*dvnum, &get_reply , 0, 0, 0);
			// while (get_reply != 2);
			// for(dvid=0;dvid<dvnum;dvid++)
			xstep = 0;
			xstnum = 0;
			xnum = 0;
			xnumtemp = DVsize;
			for(dvid=0;dvid<DVsize;dvid++){
#ifndef MACROSURF_REG_COM
			// for(dvid=0;dvid<DVsize;dvid++){
				xii_x = macsurf->xi_value[dvid*3];
				xii_y = macsurf->xi_value[dvid*3+1];
				xii_z = macsurf->xi_value[dvid*3+2];
				weight_value = macsurf->weight_value[dvid];
#else
				//第一列读主存，然后行广播给其他核
				if(my_id % 8 == 0){
					if(dvid==(xstnum+xnum)){
						xnum = (xnumtemp < N_mac) ? xnumtemp : N_mac;
						xnumtemp = xnumtemp - N_mac;
						xstnum = xstep * N_mac;
						//第一列读主存，然后行广播给其他核
						// pe_get(slave.xi_value+xstnum*3, xi_value, sizeof(double)*xnum*3);
						// dma_syn();
						get_reply = 0;
						athread_get(PE_MODE , slave.xi_value+xstnum*3, xi_value, sizeof(double)*xnum*3, &get_reply , 0, 0, 0);
						athread_get(PE_MODE , slave.weight_value+xstnum, weight_v, sizeof(double)*xnum, &get_reply , 0, 0, 0);
						while (get_reply != 2);
						xstep++;
					}
					((double*)(&xv4))[0] = xi_value[(dvid-xstnum)*3];
					((double*)(&xv4))[1] = xi_value[(dvid-xstnum)*3+1];
					((double*)(&xv4))[2] = xi_value[(dvid-xstnum)*3+2];
					((double*)(&xv4))[3] = weight_v[dvid-xstnum];
					REG_PUTR(xv4,8);
					xii_x = xi_value[(dvid-xstnum)*3];
					xii_y = xi_value[(dvid-xstnum)*3+1];
					xii_z = xi_value[(dvid-xstnum)*3+2];
					weight_value = weight_v[dvid-xstnum];
				}
				// else{
				if(my_id % 8 != 0){
					REG_GETR(xv4);
					xii_x = ((double*)(&xv4))[0];
					xii_y = ((double*)(&xv4))[1];
					xii_z = ((double*)(&xv4))[2];
					weight_value = ((double*)(&xv4))[3];
				}
				// REG_SYNR(0xff);
#endif
				
				// xii_x = xi_value[dvid*3];
				// xii_y = xi_value[dvid*3+1];
				// xii_z = xi_value[dvid*3+2];
				// weight_value = weight_[dvid];
				// get_reply = 0;
				// athread_get(PE_MODE , slave.gSurf_value[dvstnum + dvid]+stid+stnum, gSurf_value, sizeof(double)*num, &get_reply , 0, 0, 0);
				// athread_get(PE_MODE , slave.hSurf_value[dvstnum + dvid]+stid+stnum, hSurf_value, sizeof(double)*num, &get_reply , 0, 0, 0);
				// while (get_reply != 2);
				get_reply = 0;
				athread_get(PE_MODE , slave.gSurf_value[dvid]+stid+stnum, gSurf_value, sizeof(double)*num, &get_reply , 0, 0, 0);
				athread_get(PE_MODE , slave.hSurf_value[dvid]+stid+stnum, hSurf_value, sizeof(double)*num, &get_reply , 0, 0, 0);
				while (get_reply != 2);
				for(facei=0;facei<num;facei++)
				{
					// 	rhoSurf_ += dXiCellSize_ * dv.weight() * dv.gSurf();
					// 	rhoUsurf += dXiCellSize_ * dv.weight() * dv.gSurf() * dv.xi();
					// 	rhoEsurf += 0.5 * dXiCellSize_ * dv.weight()* (dv.gSurf() * magSqr(dv.xi())+ dv.hSurf());
					rhoSurf_value[facei] += weight_value * gSurf_value[facei];
					rhoUsurf_value[facei*3] += weight_value * gSurf_value[facei] * xii_x;
					rhoUsurf_value[facei*3+1] += weight_value * gSurf_value[facei] * xii_y;
					rhoUsurf_value[facei*3+2] += weight_value * gSurf_value[facei] * xii_z;
					// rhoUsurf_value[facei].vectordata[0] += weight_value * gSurf_value[facei] * xii_x;
					// rhoUsurf_value[facei].vectordata[1] += weight_value * gSurf_value[facei] * xii_y;
					// rhoUsurf_value[facei].vectordata[2] += weight_value * gSurf_value[facei] * xii_z;
					rhoEsurf_value[facei] += 0.5 * weight_value * (gSurf_value[facei] * (xii_x * xii_x + xii_y * xii_y + xii_z * xii_z) + hSurf_value[facei]);
				}
			}
		// 	dvstep++;
		// }
		put_reply = 0; 
    	athread_put(PE_MODE,rhoSurf_value,slave.rhoSurf_value+stid+stnum,sizeof(double)*num,&put_reply,0,0);                              
 		athread_put(PE_MODE,rhoUsurf_value,slave.rhoUsurf_value+stid*3+stnum*3,sizeof(double)*num*3,&put_reply,0,0);                 
  		// athread_put(PE_MODE,rhoUsurf_value,slave.rhoUsurf_value+stid+stnum,sizeof(VectorData)*num,&put_reply,0,0);                 
		athread_put(PE_MODE,rhoEsurf_value,slave.rhoEsurf_value+stid+stnum,sizeof(double)*num,&put_reply,0,0);                
		while(put_reply!=3);

		step++;
	}
}

void Func_macroqsurf(MacroqSurf *macqsurf){
	int N_slave = 448;//int N_slave = 968;
	MacroqSurf slave;
	dma_init();
    int my_id; 
    my_id = athread_get_id(-1);  //获取从核的id号
    get_reply = 0;
    athread_get(PE_MODE , macqsurf, &slave, sizeof(MacroqSurf), &get_reply , 0, 0, 0);
    while (get_reply != 1);//等待数据加载完成
    int64 ownersize1 = slave.ownersize;
	int64 ownersize = macqsurf->ownersize;	
	int64 DVsize = macqsurf->DVsize;
	// Load Balance
	int load = ownersize / NTHREAD;
	int rest = ownersize % NTHREAD;
	int stid = (my_id < rest) ? my_id * (load + 1) : my_id * load + rest;
	if (my_id < rest) load++;
	int edid = stid + load;
	int numtemp = edid - stid;
	// int num = edid - stid;
	//Define variable
	double gSurf_value[N_slave]__attribute__ ((aligned(128)));
	double hSurf_value[N_slave]__attribute__ ((aligned(128)));
	double rhoSurf_value[N_slave]__attribute__ ((aligned(128)));
	double rhoUsurf_value[N_slave*3]__attribute__ ((aligned(128)));
	double rhoEsurf_value[N_slave]__attribute__ ((aligned(128)));
	double qSurf_value[N_slave*3]__attribute__ ((aligned(128)));
	// double xi_value[N_slave];
	// double weight_value[N_slave];
	double Tsurf_value[N_slave]__attribute__ ((aligned(128)));
	double Usurf_value[N_slave*3]__attribute__ ((aligned(128)));
	double tauSurf_value[N_slave]__attribute__ ((aligned(128)));
	double qSurf_temp_value[N_slave]__attribute__ ((aligned(128)));
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
	// double Tref_omega = pow(Tref_value,omega_value);
	while(numtemp>0){
		num = (numtemp < N_slave) ? numtemp : N_slave;
		numtemp = numtemp - N_slave;
		stnum = step * N_slave;
		//get data
		// get_reply = 0;
		// athread_get(PE_MODE , slave.qSurf_value+stid*3+stnum*3, qSurf_value, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.rhoUsurf_value+stid*3+stnum*3, rhoUsurf_value, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.rhoSurf_value+stid+stnum, rhoSurf_value, sizeof(double)*num, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.rhoEsurf_value+stid+stnum, rhoEsurf_value, sizeof(double)*num, &get_reply , 0, 0, 0);
		// while (get_reply != 4);
		pe_get(slave.qSurf_value+stid*3+stnum*3, qSurf_value, sizeof(double)*num*3);
		pe_get(slave.rhoUsurf_value+stid*3+stnum*3, rhoUsurf_value, sizeof(double)*num*3);
		pe_get(slave.rhoSurf_value+stid+stnum, rhoSurf_value, sizeof(double)*num);
		pe_get(slave.rhoEsurf_value+stid+stnum, rhoEsurf_value, sizeof(double)*num);
		dma_syn();
		// Usurf_ = rhoUsurf / rhoSurf_;
		// Tsurf_ = (rhoEsurf - 0.5 * rhoSurf_ * magSqr(Usurf_)) / ((KInner_ + 3) / 2.0 * R_ * rhoSurf_);
		// tauSurf_ = muRef_ * exp(omega_ * log(Tsurf_ / Tref_)) / rhoSurf_ / Tsurf_ / R_;
		// surfaceScalarField qSurf_temp = 2.0 * tauSurf_ / (2.0 * tauSurf_ + 0.5 * time_.deltaT() * Pr_) ;
		for(facei=0;facei<num;facei++)
		{
			Usurf_value[facei*3] = rhoUsurf_value[facei*3] / rhoSurf_value[facei];
			Usurf_value[facei*3+1] = rhoUsurf_value[facei*3+1] / rhoSurf_value[facei];
			Usurf_value[facei*3+2] = rhoUsurf_value[facei*3+2] / rhoSurf_value[facei];
			magSqr_U = Usurf_value[facei*3] * Usurf_value[facei*3] + Usurf_value[facei*3+1] * Usurf_value[facei*3+1] + Usurf_value[facei*3+2] * Usurf_value[facei*3+2];
			Tsurf_value[facei] = (rhoEsurf_value[facei] - 0.5 * rhoSurf_value[facei] * magSqr_U) / ((KInner_value + 3) / 2.0 * R_value * rhoSurf_value[facei]);
			//tauSurf_value[facei] = muRef_value * exp(omega_value * lnd(Tsurf_value[facei] / Tref_value)) / rhoSurf_value[facei] / Tsurf_value[facei] / R_value;
			// tauSurf_value[facei] = muRef_value * (pow(Tsurf_value[facei],omega_value)*Tref_omega)/ rhoSurf_value[facei] / Tsurf_value[facei] / R_value;
#ifndef MACROSURF_TRAN_FUN
			tauSurf_value[facei] = muRef_value * exp(omega_value * log(Tsurf_value[facei] / Tref_value)) / rhoSurf_value[facei] / Tsurf_value[facei] / R_value;
#else
			tauSurf_value[facei] = muRef_value * expd(omega_value * log(Tsurf_value[facei] / Tref_value)) / rhoSurf_value[facei] / Tsurf_value[facei] / R_value;
#endif
			qSurf_temp_value[facei] = 2.0 * tauSurf_value[facei] / (2.0 * tauSurf_value[facei] + 0.5 * dt * Pr_value) ;
		}
		for(dvid=0;dvid<DVsize;dvid++)
		{
			xii_x = slave.xi_value[dvid*3];
			xii_y = slave.xi_value[dvid*3+1];
			xii_z = slave.xi_value[dvid*3+2];
			weight_value = slave.weight_value[dvid];
			// get_reply = 0;
			// athread_get(PE_MODE , slave.gSurf_value[dvid]+stid+stnum, gSurf_value, sizeof(double)*num, &get_reply , 0, 0, 0);
			// athread_get(PE_MODE , slave.hSurf_value[dvid]+stid+stnum, hSurf_value, sizeof(double)*num, &get_reply , 0, 0, 0);
			// while (get_reply != 2);
			pe_get(slave.gSurf_value[dvid]+stid+stnum, gSurf_value, sizeof(double)*num);
			pe_get(slave.hSurf_value[dvid]+stid+stnum, hSurf_value, sizeof(double)*num);
			dma_syn();
			for(facei=0;facei<num;facei++)
			{
				c_value[0] = xii_x - Usurf_value[facei*3];
				c_value[1] = xii_y - Usurf_value[facei*3+1];
				c_value[2] = xii_z - Usurf_value[facei*3+2];
				magSqr_c = c_value[0]*c_value[0]+c_value[1]*c_value[1]+c_value[2]*c_value[2];
				qSurf_value[facei*3] += 0.5 * weight_value * c_value[0] * (magSqr_c* gSurf_value[facei] + hSurf_value[facei]);
				qSurf_value[facei*3+1] += 0.5 * weight_value * c_value[1] * (magSqr_c* gSurf_value[facei] + hSurf_value[facei]);
				qSurf_value[facei*3+2] += 0.5 * weight_value * c_value[2] * (magSqr_c* gSurf_value[facei] + hSurf_value[facei]);
			
			}
		}
		// put_reply = 0; 
		// athread_put(PE_MODE,Tsurf_value,slave.Tsurf_value+stid+stnum,sizeof(double)*num,&put_reply,0,0);
		// athread_put(PE_MODE,Usurf_value,slave.Usurf_value+stid*3+stnum*3,sizeof(double)*num*3,&put_reply,0,0);
		// athread_put(PE_MODE,tauSurf_value,slave.tauSurf_value+stid+stnum,sizeof(double)*num,&put_reply,0,0);
		// athread_put(PE_MODE,qSurf_temp_value,slave.qSurf_temp_value+stid+stnum,sizeof(double)*num,&put_reply,0,0);
    	// athread_put(PE_MODE,qSurf_value,slave.qSurf_value+stid*3+stnum*3,sizeof(double)*num*3,&put_reply,0,0);                              
  		// while(put_reply!=5);
		pe_put(slave.Tsurf_value+stid+stnum, Tsurf_value, sizeof(double)*num);
		pe_put(slave.Usurf_value+stid*3+stnum*3, Usurf_value, sizeof(double)*num*3);
		pe_put(slave.tauSurf_value+stid+stnum, tauSurf_value, sizeof(double)*num);
		pe_put(slave.qSurf_temp_value+stid+stnum, qSurf_temp_value, sizeof(double)*num);
		pe_put(slave.qSurf_value+stid*3+stnum*3, qSurf_value, sizeof(double)*num*3);
		dma_syn();
		step++;
	}
}