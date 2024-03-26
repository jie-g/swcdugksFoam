#include <stdio.h>
#include "slave.h"
#include <math.h>
#include "dma_macros.h"
#include "slave_para.h"
#include <simd.h>

#include "para.h"

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


void Func_conum(CONUM *conum) {
	int N_slave = 320;
	CONUM slave;
	// unsigned long get_reply,put_reply; 
    int my_id; 
    my_id = athread_get_id(-1);  //获取从核的id号
    get_reply = 0;
    athread_get(PE_MODE , conum, &slave, sizeof(CONUM), &get_reply , 0, 0, 0);
    while (get_reply != 1);//等待数据加载完成
	int64 ownersize = conum->ownersize;
	// Load Balance
	int load = ownersize / NTHREAD;
	int rest = ownersize % NTHREAD;
	int stid = (my_id < rest) ? my_id * (load + 1) : my_id * load + rest;
	if (my_id < rest) load++;
	int edid = stid + load;
	int numtemp = edid - stid;
	//Define variable
	double UbyDx_value[N_slave];
	double Usurf_value[N_slave*3];
	double C_own[N_slave*3];
	double C_nei[N_slave*3];
	double dt = slave.dt;
	double xiMax_value = slave.xiMax_value;
	int64 D = slave.D;
	
	int facei = 0;
	int num = 0;
	int step = 0;
	int stnum = 0;
	double xii_Sf = 0;
	// 	unsigned long time1 = 0;
	// unsigned long time_sta = 0;
	// unsigned long time_end = 0;
	// time_sta=rpcc();
	double DeltaCoeffs=0.0;
	while(numtemp>0){
		num = (numtemp < N_slave) ? numtemp : N_slave;
		numtemp = numtemp - N_slave;
		stnum = step * N_slave;
		//get data
		get_reply = 0;
		athread_get(PE_MODE , slave.Usurf_value+stid*3+stnum*3, Usurf_value, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.C_own+stid*3+stnum*3, C_own, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.C_nei+stid*3+stnum*3, C_nei, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		while (get_reply != 3);

		for(facei=0;facei<num;facei++)
		{
			DeltaCoeffs = 1.0/sqrt((C_nei[facei*3] - C_own[facei*3])*(C_nei[facei*3] - C_own[facei*3])+ \
				(C_nei[facei*3+1] - C_own[facei*3+1])*(C_nei[facei*3+1] - C_own[facei*3+1])+(C_nei[facei*3+2] - C_own[facei*3+2])*(C_nei[facei*3+2] - C_own[facei*3+2]));
			UbyDx_value[facei] = DeltaCoeffs * (sqrt(Usurf_value[facei*3]*Usurf_value[facei*3]+Usurf_value[facei*3+1]* \
				Usurf_value[facei*3+1]+Usurf_value[facei*3+2]*Usurf_value[facei*3+2]) + sqrt((double)(D)) * xiMax_value);
		
    	}
		put_reply = 0; 
    	athread_put(PE_MODE,UbyDx_value,slave.UbyDx_value+stid+stnum,sizeof(double)*num,&put_reply,0,0);                 
  		while(put_reply!=1);
		step++;
	}

}

