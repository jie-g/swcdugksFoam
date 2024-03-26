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


#define NTHREAD 64
static const double VSMALL = 1.0e-300;

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

void Func_gaussgradtemp(GaussGrad1 *gaussgrad) {
	int N_slave = 320;
	GaussGrad1 slave;
	// unsigned long get_reply,put_reply; 
    int my_id; 
    my_id = athread_get_id(-1);  //获取从核的id号
    get_reply = 0;
    athread_get(PE_MODE , gaussgrad, &slave, sizeof(GaussGrad1), &get_reply , 0, 0, 0);
    while (get_reply != 1);//等待数据加载完成
	int64 ownersize = gaussgrad->ownersize;
	int64 gBarPvolsize = gaussgrad->gBarPvolsize;
	// Load Balance
	int load = ownersize / NTHREAD;
	int rest = ownersize % NTHREAD;
	int stid = (my_id < rest) ? my_id * (load + 1) : my_id * load + rest;
	if (my_id < rest) load++;
	int edid = stid + load;
	int numtemp = edid - stid;
	//Define variable
	double Sf[N_slave*3];
	int64 own[N_slave];
	int64 nei[N_slave];
	double gBarPvol_value_own[3][64];
	double gBarPvol_value_nei[3][64];
	double hBarPvol_value_own[3][64];
	double hBarPvol_value_nei[3][64];
	double gBarPgrad_value_own[N_slave*3];
	double gBarPgrad_value_nei[N_slave*3];
	double hBarPgrad_value_own[N_slave*3];
	double hBarPgrad_value_nei[N_slave*3];
	double V_own[N_slave];
	double V_nei[N_slave];
	double interpola_weight[N_slave];
	int subnum[2][3]={1000000,1000000,1000000,1000000,1000000,1000000};//own or nei
	int usecount[2][3]={0};//使用了多少次，找最小的替换掉
	double gBarPvol_own = 0.0;
	double gBarPvol_nei = 0.0;
	double hBarPvol_own = 0.0;
	double hBarPvol_nei = 0.0;
	int ownnum = 0;
	int neinum = 0;
	//own判断
	int iflag = -1;
	int jflag = -1;
	// double ghnum[4]={0};
	int facei = 0;
	int num = 0;
	int step = 0;
	int stnum = 0;
	double xii_Sf = 0;
	// 	unsigned long time1 = 0;
	// unsigned long time_sta = 0;
	// unsigned long time_end = 0;
	// time_sta=rpcc();
	double gsfi=0.0;
	double hsfi=0.0;
	while(numtemp>0){
		num = (numtemp < N_slave) ? numtemp : N_slave;
		numtemp = numtemp - N_slave;
		stnum = step * N_slave;
		//get data
		get_reply = 0;
		athread_get(PE_MODE , slave.Sf_value+stid*3+stnum*3, Sf, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.own+stid+stnum, own, sizeof(int64)*num, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.neighbour+stid+stnum, nei, sizeof(int64)*num, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.hBarPvol_value_own+stid+stnum, hBarPvol_value_own, sizeof(double)*num, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.hBarPvol_value_nei+stid+stnum, hBarPvol_value_nei, sizeof(double)*num, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.V_own+stid+stnum, V_own, sizeof(double)*num, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.V_nei+stid+stnum, V_nei, sizeof(double)*num, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.interpola_weight+stid+stnum, interpola_weight, sizeof(double)*num, &get_reply , 0, 0, 0);
    	while (get_reply != 6);

		memset(gBarPgrad_value_own,0,sizeof(double)*N_slave*3);
		memset(gBarPgrad_value_nei,0,sizeof(double)*N_slave*3);
		memset(hBarPgrad_value_own,0,sizeof(double)*N_slave*3);
		memset(hBarPgrad_value_nei,0,sizeof(double)*N_slave*3);
    	// forAll(owner, facei)
    	// {
        // 	label own = owner[facei];
        // 	label nei = neighbour[facei];
        // 	scalar gsfi = interpola_weight[facei]*(gvfi[own] - gvfi[nei]) + gvfi[nei];
        // 	scalar hsfi = interpola_weight[facei]*(hvfi[own] - hvfi[nei]) + hvfi[nei];
        // 	igGrad[own] += Sf[facei]*gsfi/V[own];
        // 	igGrad[nei] -= Sf[facei]*gsfi/V[nei];
        // 	ihGrad[own] += Sf[facei]*hsfi/V[own];
        // 	ihGrad[nei] -= Sf[facei]*hsfi/V[nei];
    	// }
		for(facei=0;facei<num;facei++)
		{
			gBarPvol_own = 0.0;
			gBarPvol_nei = 0.0;
			hBarPvol_own = 0.0;
			hBarPvol_nei = 0.0;
			ownnum = own[facei];
			neinum = nei[facei];
			//own判断
			iflag = -1;
			jflag = -1;
			for(int i=0;i<3;i++){//每个数三组
				if(subnum[0][i]<=ownnum&&(subnum[0][i]+63)>=ownnum){
					iflag = i;
					jflag = ownnum - subnum[0][i];
				}
			}
			if(iflag != -1){
				gBarPvol_own = gBarPvol_value_own[iflag][jflag];
				hBarPvol_own = hBarPvol_value_own[iflag][jflag];
				usecount[0][iflag]++;
				if(jflag > 61){
					usecount[0][iflag] = 0;
				}
			}
			else{
				//判断哪组用的最少，替换出去
				iflag = 0;
				// if(usecount[0][0] <= usecount[0][1] && usecount[0][0] <= usecount[0][2]){
				// 	iflag = 0;
				// }else 
				if(usecount[0][1] <= usecount[0][0] && usecount[0][1] <= usecount[0][2]){
					iflag = 1;
				}else if(usecount[0][2] <= usecount[0][0] && usecount[0][2] <= usecount[0][1]){
					iflag = 2;
				}
				int com_num = ((gBarPvolsize-ownnum) < 64)?(gBarPvolsize-ownnum) : 64;
				// if(com_num==64){
				// 	jflag = 0;
				// }else{
				// 	jflag = 64-(gBarPvolsize-ownnum);
				// }
				jflag = 64-com_num;
				subnum[0][iflag] = ownnum-jflag;
				usecount[0][iflag] = 0;
				get_reply = 0;
				athread_get(PE_MODE , slave.gBarPvol_value+ownnum-jflag, &(gBarPvol_value_own[iflag][0]), sizeof(double)*64, &get_reply , 0, 0, 0);
				athread_get(PE_MODE , slave.hBarPvol_value+ownnum-jflag, &(hBarPvol_value_own[iflag][0]), sizeof(double)*64, &get_reply , 0, 0, 0);
				while (get_reply != 2);
				gBarPvol_own = gBarPvol_value_own[iflag][jflag];
				hBarPvol_own = hBarPvol_value_own[iflag][jflag];
			}
			//nei判断
			iflag = -1;
			jflag = -1;
			for(int i=0;i<3;i++){//每个数三组
				if(subnum[1][i]<=neinum&&(subnum[1][i]+63)>=neinum){
					iflag = i;
					jflag = neinum - subnum[1][i];
				}
			}
			if(iflag != -1){
				gBarPvol_nei = gBarPvol_value_nei[iflag][jflag];
				hBarPvol_nei = hBarPvol_value_nei[iflag][jflag];
				usecount[1][iflag]++;
				if(jflag > 61){
					usecount[1][iflag] = 0;
				}
			}
			else{
				//判断哪组用的最少，替换出去
				iflag = 0;
				// if(usecount[1][0] <= usecount[1][1] && usecount[1][0] <= usecount[1][2]){
				// 	iflag = 0;
				// }else 
				if(usecount[1][1] <= usecount[1][0] && usecount[1][1] <= usecount[1][2]){
					iflag = 1;
				}else if(usecount[1][2] <= usecount[1][0] && usecount[1][2] <= usecount[1][1]){
					iflag = 2;
				}
				int com_num = ((gBarPvolsize-neinum) < 64)?(gBarPvolsize-neinum) : 64;
				// if(com_num==64){
				// 	jflag = 0;
				// }else{
				// 	jflag = 64-(gBarPvolsize-neinum);
				// }
				jflag = 64-com_num;
				subnum[1][iflag] = neinum-jflag;
				usecount[1][iflag] = 0;
				get_reply = 0;
				athread_get(PE_MODE , slave.gBarPvol_value+neinum-jflag, &(gBarPvol_value_nei[iflag][0]), sizeof(double)*64, &get_reply , 0, 0, 0);
				athread_get(PE_MODE , slave.hBarPvol_value+neinum-jflag, &(hBarPvol_value_nei[iflag][0]), sizeof(double)*64, &get_reply , 0, 0, 0);
				while (get_reply != 2);
				gBarPvol_nei = gBarPvol_value_nei[iflag][jflag];
				hBarPvol_nei = hBarPvol_value_nei[iflag][jflag];
			}
			
			gsfi = interpola_weight[facei]*(gBarPvol_own - gBarPvol_nei) + gBarPvol_nei;
        	hsfi = interpola_weight[facei]*(hBarPvol_own - hBarPvol_nei) + hBarPvol_nei;
			// gsfi = interpola_weight[facei]*(gBarPvol_value_own[facei] - gBarPvol_value_nei[facei]) + gBarPvol_value_nei[facei];
        	// hsfi = interpola_weight[facei]*(hBarPvol_value_own[facei] - hBarPvol_value_nei[facei]) + hBarPvol_value_nei[facei];
        	gBarPgrad_value_own[facei*3]   += Sf[facei*3]  *gsfi/V_own[facei];
			gBarPgrad_value_own[facei*3+1] += Sf[facei*3+1]*gsfi/V_own[facei];
			gBarPgrad_value_own[facei*3+2] += Sf[facei*3+2]*gsfi/V_own[facei];
        	gBarPgrad_value_nei[facei*3]   -= Sf[facei*3]  *gsfi/V_nei[facei];
        	gBarPgrad_value_nei[facei*3+1] -= Sf[facei*3+1]*gsfi/V_nei[facei];
        	gBarPgrad_value_nei[facei*3+2] -= Sf[facei*3+2]*gsfi/V_nei[facei];
        	hBarPgrad_value_own[facei*3]   += Sf[facei*3]  *hsfi/V_own[facei];
			hBarPgrad_value_own[facei*3+1] += Sf[facei*3+1]*hsfi/V_own[facei];
			hBarPgrad_value_own[facei*3+2] += Sf[facei*3+2]*hsfi/V_own[facei];
        	hBarPgrad_value_nei[facei*3]   -= Sf[facei*3]  *hsfi/V_nei[facei];
        	hBarPgrad_value_nei[facei*3+1] -= Sf[facei*3+1]*hsfi/V_nei[facei];
        	hBarPgrad_value_nei[facei*3+2] -= Sf[facei*3+2]*hsfi/V_nei[facei];

    	}
		put_reply = 0; 
    	athread_put(PE_MODE,gBarPgrad_value_own,slave.gBarPgrad_value_own+stid*3+stnum*3,sizeof(double)*num*3,&put_reply,0,0);                 
  		athread_put(PE_MODE,gBarPgrad_value_nei,slave.gBarPgrad_value_nei+stid*3+stnum*3,sizeof(double)*num*3,&put_reply,0,0);                 
  		athread_put(PE_MODE,hBarPgrad_value_own,slave.hBarPgrad_value_own+stid*3+stnum*3,sizeof(double)*num*3,&put_reply,0,0);                 
  		athread_put(PE_MODE,hBarPgrad_value_nei,slave.hBarPgrad_value_nei+stid*3+stnum*3,sizeof(double)*num*3,&put_reply,0,0);               
		while(put_reply!=4);
		step++;
	}

}
void Func_ghbarsurftemp(GHbarSurf1 *ghbarsurf) {
	int N_slave = 256;
	GHbarSurf1 slave;
	// unsigned long get_reply,put_reply; 
    int my_id; 
    my_id = athread_get_id(-1);  //获取从核的id号
    get_reply = 0;
    athread_get(PE_MODE , ghbarsurf, &slave, sizeof(GHbarSurf1), &get_reply , 0, 0, 0);
    while (get_reply != 1);//等待数据加载完成
	int64 ownersize = ghbarsurf->ownersize;
	int64 gBarPvolsize = ghbarsurf->gBarPvolsize;
	// Load Balance
	int load = ownersize / NTHREAD;
	int rest = ownersize % NTHREAD;
	int stid = (my_id < rest) ? my_id * (load + 1) : my_id * load + rest;
	if (my_id < rest) load++;
	int edid = stid + load;
	int numtemp = edid - stid;
	//Define variable
	double Sf[N_slave*3];
	double Cf[N_slave*3];
	// double gBarPvol_value_own[N_slave];
	// double gBarPvol_value_nei[N_slave];
	// double hBarPvol_value_own[N_slave];
	// double hBarPvol_value_nei[N_slave];
	double gBarPgrad_value_own[N_slave*3];
	double gBarPgrad_value_nei[N_slave*3];
	double hBarPgrad_value_own[N_slave*3];
	double hBarPgrad_value_nei[N_slave*3];
	double C_own[N_slave*3];
	double C_nei[N_slave*3];
	double gSurf_value[N_slave];
	double hSurf_value[N_slave];

	int64 own[N_slave];
	int64 nei[N_slave];
	double gBarPvol_value_own[3][64];
	double gBarPvol_value_nei[3][64];
	double hBarPvol_value_own[3][64];
	double hBarPvol_value_nei[3][64];
	int subnum[2][3]={1000000,1000000,1000000,1000000,1000000,1000000};//own or nei
	int usecount[2][3]={0};//使用了多少次，找最小的替换掉
	double gBarPvol_own = 0.0;
	double gBarPvol_nei = 0.0;
	double hBarPvol_own = 0.0;
	double hBarPvol_nei = 0.0;
	int ownnum = 0;
	int neinum = 0;
	//own判断
	int iflag = -1;
	int jflag = -1;
	
	double dt = ghbarsurf->dt;
	double xii_x = ghbarsurf->xii_x;
	double xii_y = ghbarsurf->xii_y;
	double xii_z = ghbarsurf->xii_z;
	int facei = 0;
	int num = 0;
	int step = 0;
	int stnum = 0;
	double xii_Sf = 0;
	// 	unsigned long time1 = 0;
	// unsigned long time_sta = 0;
	// unsigned long time_end = 0;
	// time_sta=rpcc();
	while(numtemp>0){
		num = (numtemp < N_slave) ? numtemp : N_slave;
		numtemp = numtemp - N_slave;
		stnum = step * N_slave;
		//get data
		get_reply = 0;
		athread_get(PE_MODE , slave.Sf+stid*3+stnum*3, Sf, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.Cf+stid*3+stnum*3, Cf, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.gBarPvol_value_own+stid+stnum, gBarPvol_value_own, sizeof(double)*num, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.gBarPvol_value_nei+stid+stnum, gBarPvol_value_nei, sizeof(double)*num, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.hBarPvol_value_own+stid+stnum, hBarPvol_value_own, sizeof(double)*num, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.hBarPvol_value_nei+stid+stnum, hBarPvol_value_nei, sizeof(double)*num, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.own+stid+stnum, own, sizeof(int64)*num, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.neighbour+stid+stnum, nei, sizeof(int64)*num, &get_reply , 0, 0, 0);
		
		athread_get(PE_MODE , slave.gBarPgrad_value_own+stid*3+stnum*3, gBarPgrad_value_own, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.gBarPgrad_value_nei+stid*3+stnum*3, gBarPgrad_value_nei, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.hBarPgrad_value_own+stid*3+stnum*3, hBarPgrad_value_own, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.hBarPgrad_value_nei+stid*3+stnum*3, hBarPgrad_value_nei, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.C_own+stid*3+stnum*3, C_own, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.C_nei+stid*3+stnum*3, C_nei, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.gSurf_value+stid+stnum, gSurf_value, sizeof(double)*num, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.hSurf_value+stid+stnum, hSurf_value, sizeof(double)*num, &get_reply , 0, 0, 0);
    	while (get_reply != 10);
		for(facei=0;facei<num;facei++)
		{
			gBarPvol_own = 0.0;
			gBarPvol_nei = 0.0;
			hBarPvol_own = 0.0;
			hBarPvol_nei = 0.0;
			ownnum = own[facei];
			neinum = nei[facei];
			//own判断
			iflag = -1;
			jflag = -1;
			for(int i=0;i<3;i++){//每个数三组
				if(subnum[0][i]<=ownnum&&(subnum[0][i]+63)>=ownnum){
					iflag = i;
					jflag = ownnum - subnum[0][i];
				}
			}
			if(iflag != -1){
				gBarPvol_own = gBarPvol_value_own[iflag][jflag];
				hBarPvol_own = hBarPvol_value_own[iflag][jflag];
				usecount[0][iflag]++;
				if(jflag > 61){
					usecount[0][iflag] = 0;
				}
			}
			else{
				//判断哪组用的最少，替换出去
				iflag = 0;
				// if(usecount[0][0] <= usecount[0][1] && usecount[0][0] <= usecount[0][2]){
				// 	iflag = 0;
				// }else 
				if(usecount[0][1] <= usecount[0][0] && usecount[0][1] <= usecount[0][2]){
					iflag = 1;
				}else if(usecount[0][2] <= usecount[0][0] && usecount[0][2] <= usecount[0][1]){
					iflag = 2;
				}
				int com_num = ((gBarPvolsize-ownnum) < 64)?(gBarPvolsize-ownnum) : 64;
				// if(com_num==64){
				// 	jflag = 0;
				// }else{
				// 	jflag = 64-(gBarPvolsize-ownnum);
				// }
				jflag = 64-com_num;
				subnum[0][iflag] = ownnum-jflag;
				usecount[0][iflag] = 0;
				get_reply = 0;
				athread_get(PE_MODE , slave.gBarPvol_value+ownnum-jflag, &(gBarPvol_value_own[iflag][0]), sizeof(double)*64, &get_reply , 0, 0, 0);
				athread_get(PE_MODE , slave.hBarPvol_value+ownnum-jflag, &(hBarPvol_value_own[iflag][0]), sizeof(double)*64, &get_reply , 0, 0, 0);
				while (get_reply != 2);
				gBarPvol_own = gBarPvol_value_own[iflag][jflag];
				hBarPvol_own = hBarPvol_value_own[iflag][jflag];
			}
			//nei判断
			iflag = -1;
			jflag = -1;
			for(int i=0;i<3;i++){//每个数三组
				if(subnum[1][i]<=neinum&&(subnum[1][i]+63)>=neinum){
					iflag = i;
					jflag = neinum - subnum[1][i];
				}
			}
			if(iflag != -1){
				gBarPvol_nei = gBarPvol_value_nei[iflag][jflag];
				hBarPvol_nei = hBarPvol_value_nei[iflag][jflag];
				usecount[1][iflag]++;
				if(jflag > 61){
					usecount[1][iflag] = 0;
				}
			}
			else{
				//判断哪组用的最少，替换出去
				iflag = 0;
				// if(usecount[1][0] <= usecount[1][1] && usecount[1][0] <= usecount[1][2]){
				// 	iflag = 0;
				// }else 
				if(usecount[1][1] <= usecount[1][0] && usecount[1][1] <= usecount[1][2]){
					iflag = 1;
				}else if(usecount[1][2] <= usecount[1][0] && usecount[1][2] <= usecount[1][1]){
					iflag = 2;
				}
				int com_num = ((gBarPvolsize-neinum) < 64)?(gBarPvolsize-neinum) : 64;
				// if(com_num==64){
				// 	jflag = 0;
				// }else{
				// 	jflag = 64-(gBarPvolsize-neinum);
				// }
				jflag = 64-com_num;
				subnum[1][iflag] = neinum-jflag;
				usecount[1][iflag] = 0;
				get_reply = 0;
				athread_get(PE_MODE , slave.gBarPvol_value+neinum-jflag, &(gBarPvol_value_nei[iflag][0]), sizeof(double)*64, &get_reply , 0, 0, 0);
				athread_get(PE_MODE , slave.hBarPvol_value+neinum-jflag, &(hBarPvol_value_nei[iflag][0]), sizeof(double)*64, &get_reply , 0, 0, 0);
				while (get_reply != 2);
				gBarPvol_nei = gBarPvol_value_nei[iflag][jflag];
				hBarPvol_nei = hBarPvol_value_nei[iflag][jflag];
			}

			xii_Sf = xii_x * Sf[facei*3] + xii_y * Sf[facei*3+1] + xii_z * Sf[facei*3+2];
			if (xii_Sf >=  VSMALL) // comming from own
			{
				gSurf_value[facei] = gBarPvol_own + \
				(gBarPgrad_value_own[facei*3]*(Cf[facei*3] - C_own[facei*3] - 0.5*xii_x*dt)) \
				+(gBarPgrad_value_own[facei*3+1]*(Cf[facei*3+1] - C_own[facei*3+1] - 0.5*xii_y*dt)) \
				+(gBarPgrad_value_own[facei*3+2]*(Cf[facei*3+2] - C_own[facei*3+2] - 0.5*xii_z*dt));
				hSurf_value[facei] = hBarPvol_own + \
				(hBarPgrad_value_own[facei*3]*(Cf[facei*3] - C_own[facei*3] - 0.5*xii_x*dt)) \
				+(hBarPgrad_value_own[facei*3+1]*(Cf[facei*3+1] - C_own[facei*3+1] - 0.5*xii_y*dt)) \
				+(hBarPgrad_value_own[facei*3+2]*(Cf[facei*3+2] - C_own[facei*3+2] - 0.5*xii_z*dt));
				// iGsurf[facei] = iGbarPvol[own] + (iGbarPgrad[own]&(Cf[facei] - C[own] - 0.5*xii*dt));
            	// iHsurf[facei] = iHbarPvol[own] + (iHbarPgrad[own]&(Cf[facei] - C[own] - 0.5*xii*dt));
			}
			// Debug, no = 0, =0 put to > 0
			else if (xii_Sf < -VSMALL) // comming form nei
			{
				gSurf_value[facei] = gBarPvol_nei + \
				(gBarPgrad_value_nei[facei*3]*(Cf[facei*3] - C_nei[facei*3] - 0.5*xii_x*dt)) \
				+(gBarPgrad_value_nei[facei*3+1]*(Cf[facei*3+1] - C_nei[facei*3+1] - 0.5*xii_y*dt)) \
				+(gBarPgrad_value_nei[facei*3+2]*(Cf[facei*3+2] - C_nei[facei*3+2] - 0.5*xii_z*dt));
				hSurf_value[facei] = hBarPvol_nei + \
				(hBarPgrad_value_nei[facei*3]*(Cf[facei*3] - C_nei[facei*3] - 0.5*xii_x*dt)) \
				+(hBarPgrad_value_nei[facei*3+1]*(Cf[facei*3+1] - C_nei[facei*3+1] - 0.5*xii_y*dt)) \
				+(hBarPgrad_value_nei[facei*3+2]*(Cf[facei*3+2] - C_nei[facei*3+2] - 0.5*xii_z*dt));
				
				// gSurf_value[facei] = iGbarPvol[nei] + (iGbarPgrad[nei]&(Cf[facei] - C[nei] - 0.5*xii*dt));
				// hSurf_value[facei] = iHbarPvol[nei] + (iHbarPgrad[nei]&(Cf[facei] - C[nei] - 0.5*xii*dt));
			}
			else 
			{
				gSurf_value[facei] = 0.5*(gBarPvol_nei + \
				(gBarPgrad_value_nei[facei*3]*(Cf[facei*3] - C_nei[facei*3] - 0.5*xii_x*dt)) \
				+(gBarPgrad_value_nei[facei*3+1]*(Cf[facei*3+1] - C_nei[facei*3+1] - 0.5*xii_y*dt)) \
				+(gBarPgrad_value_nei[facei*3+2]*(Cf[facei*3+2] - C_nei[facei*3+2] - 0.5*xii_z*dt)) \
				+gBarPvol_own + \
				(gBarPgrad_value_own[facei*3]*(Cf[facei*3] - C_own[facei*3] - 0.5*xii_x*dt)) \
				+(gBarPgrad_value_own[facei*3+1]*(Cf[facei*3+1] - C_own[facei*3+1] - 0.5*xii_y*dt)) \
				+(gBarPgrad_value_own[facei*3+2]*(Cf[facei*3+2] - C_own[facei*3+2] - 0.5*xii_z*dt)));
				hSurf_value[facei] = 0.5*(hBarPvol_nei + \
				(hBarPgrad_value_nei[facei*3]*(Cf[facei*3] - C_nei[facei*3] - 0.5*xii_x*dt)) \
				+(hBarPgrad_value_nei[facei*3+1]*(Cf[facei*3+1] - C_nei[facei*3+1] - 0.5*xii_y*dt)) \
				+(hBarPgrad_value_nei[facei*3+2]*(Cf[facei*3+2] - C_nei[facei*3+2] - 0.5*xii_z*dt)) \
				+ hBarPvol_own + \
				(hBarPgrad_value_own[facei*3]*(Cf[facei*3] - C_own[facei*3] - 0.5*xii_x*dt)) \
				+(hBarPgrad_value_own[facei*3+1]*(Cf[facei*3+1] - C_own[facei*3+1] - 0.5*xii_y*dt)) \
				+(hBarPgrad_value_own[facei*3+2]*(Cf[facei*3+2] - C_own[facei*3+2] - 0.5*xii_z*dt)));
				
				// gSurf_value[facei] = 0.5*(iGbarPvol[nei] + ((iGbarPgrad[nei]) &(Cf[facei] - C[nei] - 0.5*xii*dt)) + iGbarPvol[own] + ((iGbarPgrad[own]) &(Cf[facei] - C[own] - 0.5*xii*dt)));
				// hSurf_value[facei] = 0.5*(iHbarPvol[nei] + ((iHbarPgrad[nei]) &(Cf[facei] - C[nei] - 0.5*xii*dt)) + iHbarPvol[own] + ((iHbarPgrad[own]) &(Cf[facei] - C[own] - 0.5*xii*dt)));
			}
		}
		put_reply = 0; 
		athread_put(PE_MODE,gSurf_value,slave.gSurf_value+stid+stnum, sizeof(double)*num, &put_reply , 0, 0);
		athread_put(PE_MODE,hSurf_value,slave.hSurf_value+stid+stnum, sizeof(double)*num, &put_reply , 0, 0);             
		while(put_reply!=2);
		step++;
	}

}
void Func_ghbarsurftemp2(GHbarSurf2 *ghbarsurf) {
	int N_slave = 256;
	GHbarSurf2 slave;
	// unsigned long get_reply,put_reply; 
    int my_id; 
    my_id = athread_get_id(-1);  //获取从核的id号
    get_reply = 0;
    athread_get(PE_MODE , ghbarsurf, &slave, sizeof(GHbarSurf2), &get_reply , 0, 0, 0);
    while (get_reply != 1);//等待数据加载完成
	int64 ownersize = ghbarsurf->ownersize;
	int64 gBarPvolsize = ghbarsurf->gBarPvolsize;
	// Load Balance
	int load = ownersize / NTHREAD;
	int rest = ownersize % NTHREAD;
	int stid = (my_id < rest) ? my_id * (load + 1) : my_id * load + rest;
	if (my_id < rest) load++;
	int edid = stid + load;
	int numtemp = edid - stid;
	//Define variable
	double Sf[N_slave*3];
	double Cf[N_slave*3];
	// double gBarPvol_value_own[N_slave];
	// double gBarPvol_value_nei[N_slave];
	// double hBarPvol_value_own[N_slave];
	// double hBarPvol_value_nei[N_slave];
	// double gBarPgrad_value_own[N_slave*3];
	// double gBarPgrad_value_nei[N_slave*3];
	// double hBarPgrad_value_own[N_slave*3];
	// double hBarPgrad_value_nei[N_slave*3];
	double C_own[N_slave*3];
	double C_nei[N_slave*3];
	double gSurf_value[N_slave];
	double hSurf_value[N_slave];

	int64 own[N_slave];
	int64 nei[N_slave];
	double gBarPvol_value_own[3][64];
	double gBarPvol_value_nei[3][64];
	double hBarPvol_value_own[3][64];
	double hBarPvol_value_nei[3][64];
	double gBarPgrad_value_own[3][192];//64*3
	double gBarPgrad_value_nei[3][192];
	double hBarPgrad_value_own[3][192];
	double hBarPgrad_value_nei[3][192];
	double gBarPgrad_own[3];
	double gBarPgrad_nei[3];
	double hBarPgrad_own[3];
	double hBarPgrad_nei[3];
	int subnum[2][3]={1000000,1000000,1000000,1000000,1000000,1000000};//own or nei
	int usecount[2][3]={0};//使用了多少次，找最小的替换掉
	double gBarPvol_own = 0.0;
	double gBarPvol_nei = 0.0;
	double hBarPvol_own = 0.0;
	double hBarPvol_nei = 0.0;
	int ownnum = 0;
	int neinum = 0;
	//own判断
	int iflag = -1;
	int jflag = -1;
	
	double dt = ghbarsurf->dt;
	double xii_x = ghbarsurf->xii_x;
	double xii_y = ghbarsurf->xii_y;
	double xii_z = ghbarsurf->xii_z;
	int facei = 0;
	int num = 0;
	int step = 0;
	int stnum = 0;
	double xii_Sf = 0;
	// 	unsigned long time1 = 0;
	// unsigned long time_sta = 0;
	// unsigned long time_end = 0;
	// time_sta=rpcc();
	while(numtemp>0){
		num = (numtemp < N_slave) ? numtemp : N_slave;
		numtemp = numtemp - N_slave;
		stnum = step * N_slave;
		//get data
		get_reply = 0;
		athread_get(PE_MODE , slave.Sf+stid*3+stnum*3, Sf, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.Cf+stid*3+stnum*3, Cf, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.gBarPvol_value_own+stid+stnum, gBarPvol_value_own, sizeof(double)*num, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.gBarPvol_value_nei+stid+stnum, gBarPvol_value_nei, sizeof(double)*num, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.hBarPvol_value_own+stid+stnum, hBarPvol_value_own, sizeof(double)*num, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.hBarPvol_value_nei+stid+stnum, hBarPvol_value_nei, sizeof(double)*num, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.own+stid+stnum, own, sizeof(int64)*num, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.neighbour+stid+stnum, nei, sizeof(int64)*num, &get_reply , 0, 0, 0);
		
		// athread_get(PE_MODE , slave.gBarPgrad_value_own+stid*3+stnum*3, gBarPgrad_value_own, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.gBarPgrad_value_nei+stid*3+stnum*3, gBarPgrad_value_nei, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.hBarPgrad_value_own+stid*3+stnum*3, hBarPgrad_value_own, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.hBarPgrad_value_nei+stid*3+stnum*3, hBarPgrad_value_nei, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.C_own+stid*3+stnum*3, C_own, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.C_nei+stid*3+stnum*3, C_nei, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.gSurf_value+stid+stnum, gSurf_value, sizeof(double)*num, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.hSurf_value+stid+stnum, hSurf_value, sizeof(double)*num, &get_reply , 0, 0, 0);
    	while (get_reply != 6);
		for(facei=0;facei<num;facei++)
		{
			gBarPvol_own = 0.0;
			gBarPvol_nei = 0.0;
			hBarPvol_own = 0.0;
			hBarPvol_nei = 0.0;
			ownnum = own[facei];
			neinum = nei[facei];
			//own判断
			iflag = -1;
			jflag = -1;
			for(int i=0;i<3;i++){//每个数三组
				if(subnum[0][i]<=ownnum&&(subnum[0][i]+63)>=ownnum){
					iflag = i;
					jflag = ownnum - subnum[0][i];
				}
			}
			if(iflag != -1){
				gBarPvol_own = gBarPvol_value_own[iflag][jflag];
				hBarPvol_own = hBarPvol_value_own[iflag][jflag];
				gBarPgrad_own[0] =gBarPgrad_value_own[iflag][jflag*3];
				gBarPgrad_own[1] =gBarPgrad_value_own[iflag][jflag*3+1];
				gBarPgrad_own[2] =gBarPgrad_value_own[iflag][jflag*3+2];
				hBarPgrad_own[0] =hBarPgrad_value_own[iflag][jflag*3];
				hBarPgrad_own[1] =hBarPgrad_value_own[iflag][jflag*3+1];
				hBarPgrad_own[2] =hBarPgrad_value_own[iflag][jflag*3+2];
				usecount[0][iflag]++;
				if(jflag > 61){
					usecount[0][iflag] = 0;
				}
			}
			else{
				//判断哪组用的最少，替换出去
				iflag = 0;
				// if(usecount[0][0] <= usecount[0][1] && usecount[0][0] <= usecount[0][2]){
				// 	iflag = 0;
				// }else 
				if(usecount[0][1] <= usecount[0][0] && usecount[0][1] <= usecount[0][2]){
					iflag = 1;
				}else if(usecount[0][2] <= usecount[0][0] && usecount[0][2] <= usecount[0][1]){
					iflag = 2;
				}
				int com_num = ((gBarPvolsize-ownnum) < 64)?(gBarPvolsize-ownnum) : 64;
				// if(com_num==64){
				// 	jflag = 0;
				// }else{
				// 	jflag = 64-(gBarPvolsize-ownnum);
				// }
				jflag = 64-com_num;
				subnum[0][iflag] = ownnum-jflag;
				usecount[0][iflag] = 0;
				get_reply = 0;
				athread_get(PE_MODE , slave.gBarPvol_value+ownnum-jflag, &(gBarPvol_value_own[iflag][0]), sizeof(double)*64, &get_reply , 0, 0, 0);
				athread_get(PE_MODE , slave.hBarPvol_value+ownnum-jflag, &(hBarPvol_value_own[iflag][0]), sizeof(double)*64, &get_reply , 0, 0, 0);
				athread_get(PE_MODE , slave.gBarPgrad_value+ownnum*3-jflag*3, &(gBarPgrad_value_own[iflag][0]), sizeof(double)*64*3, &get_reply , 0, 0, 0);
				athread_get(PE_MODE , slave.hBarPgrad_value+ownnum*3-jflag*3, &(hBarPgrad_value_own[iflag][0]), sizeof(double)*64*3, &get_reply , 0, 0, 0);
				while (get_reply != 4);
				gBarPvol_own = gBarPvol_value_own[iflag][jflag];
				hBarPvol_own = hBarPvol_value_own[iflag][jflag];
				gBarPgrad_own[0] =gBarPgrad_value_own[iflag][jflag*3];
				gBarPgrad_own[1] =gBarPgrad_value_own[iflag][jflag*3+1];
				gBarPgrad_own[2] =gBarPgrad_value_own[iflag][jflag*3+2];
				hBarPgrad_own[0] =hBarPgrad_value_own[iflag][jflag*3];
				hBarPgrad_own[1] =hBarPgrad_value_own[iflag][jflag*3+1];
				hBarPgrad_own[2] =hBarPgrad_value_own[iflag][jflag*3+2];
			}
			//nei判断
			iflag = -1;
			jflag = -1;
			for(int i=0;i<3;i++){//每个数三组
				if(subnum[1][i]<=neinum&&(subnum[1][i]+63)>=neinum){
					iflag = i;
					jflag = neinum - subnum[1][i];
				}
			}
			if(iflag != -1){
				gBarPvol_nei = gBarPvol_value_nei[iflag][jflag];
				hBarPvol_nei = hBarPvol_value_nei[iflag][jflag];
				gBarPgrad_nei[0] =gBarPgrad_value_nei[iflag][jflag*3];
				gBarPgrad_nei[1] =gBarPgrad_value_nei[iflag][jflag*3+1];
				gBarPgrad_nei[2] =gBarPgrad_value_nei[iflag][jflag*3+2];
				hBarPgrad_nei[0] =hBarPgrad_value_nei[iflag][jflag*3];
				hBarPgrad_nei[1] =hBarPgrad_value_nei[iflag][jflag*3+1];
				hBarPgrad_nei[2] =hBarPgrad_value_nei[iflag][jflag*3+2];
				usecount[1][iflag]++;
				if(jflag > 61){
					usecount[1][iflag] = 0;
				}
			}
			else{
				//判断哪组用的最少，替换出去
				iflag = 0;
				// if(usecount[1][0] <= usecount[1][1] && usecount[1][0] <= usecount[1][2]){
				// 	iflag = 0;
				// }else 
				if(usecount[1][1] <= usecount[1][0] && usecount[1][1] <= usecount[1][2]){
					iflag = 1;
				}else if(usecount[1][2] <= usecount[1][0] && usecount[1][2] <= usecount[1][1]){
					iflag = 2;
				}
				int com_num = ((gBarPvolsize-neinum) < 64)?(gBarPvolsize-neinum) : 64;
				// if(com_num==64){
				// 	jflag = 0;
				// }else{
				// 	jflag = 64-(gBarPvolsize-neinum);
				// }
				jflag = 64-com_num;
				subnum[1][iflag] = neinum-jflag;
				usecount[1][iflag] = 0;
				get_reply = 0;
				athread_get(PE_MODE , slave.gBarPvol_value+neinum-jflag, &(gBarPvol_value_nei[iflag][0]), sizeof(double)*64, &get_reply , 0, 0, 0);
				athread_get(PE_MODE , slave.hBarPvol_value+neinum-jflag, &(hBarPvol_value_nei[iflag][0]), sizeof(double)*64, &get_reply , 0, 0, 0);
				athread_get(PE_MODE , slave.gBarPgrad_value+neinum*3-jflag*3, &(gBarPgrad_value_nei[iflag][0]), sizeof(double)*64*3, &get_reply , 0, 0, 0);
				athread_get(PE_MODE , slave.hBarPgrad_value+neinum*3-jflag*3, &(hBarPgrad_value_nei[iflag][0]), sizeof(double)*64*3, &get_reply , 0, 0, 0);
				while (get_reply != 4);
				gBarPvol_nei = gBarPvol_value_nei[iflag][jflag];
				hBarPvol_nei = hBarPvol_value_nei[iflag][jflag];
				gBarPgrad_nei[0] =gBarPgrad_value_nei[iflag][jflag*3];
				gBarPgrad_nei[1] =gBarPgrad_value_nei[iflag][jflag*3+1];
				gBarPgrad_nei[2] =gBarPgrad_value_nei[iflag][jflag*3+2];
				hBarPgrad_nei[0] =hBarPgrad_value_nei[iflag][jflag*3];
				hBarPgrad_nei[1] =hBarPgrad_value_nei[iflag][jflag*3+1];
				hBarPgrad_nei[2] =hBarPgrad_value_nei[iflag][jflag*3+2];
			}

			xii_Sf = xii_x * Sf[facei*3] + xii_y * Sf[facei*3+1] + xii_z * Sf[facei*3+2];
			if (xii_Sf >=  VSMALL) // comming from own
			{
				gSurf_value[facei] = gBarPvol_own + \
				(gBarPgrad_own[0]*(Cf[facei*3] - C_own[facei*3] - 0.5*xii_x*dt)) \
				+(gBarPgrad_own[1]*(Cf[facei*3+1] - C_own[facei*3+1] - 0.5*xii_y*dt)) \
				+(gBarPgrad_own[2]*(Cf[facei*3+2] - C_own[facei*3+2] - 0.5*xii_z*dt));
				hSurf_value[facei] = hBarPvol_own + \
				(hBarPgrad_own[0]*(Cf[facei*3] - C_own[facei*3] - 0.5*xii_x*dt)) \
				+(hBarPgrad_own[1]*(Cf[facei*3+1] - C_own[facei*3+1] - 0.5*xii_y*dt)) \
				+(hBarPgrad_own[2]*(Cf[facei*3+2] - C_own[facei*3+2] - 0.5*xii_z*dt));
				// iGsurf[facei] = iGbarPvol[own] + (iGbarPgrad[own]&(Cf[facei] - C[own] - 0.5*xii*dt));
            	// iHsurf[facei] = iHbarPvol[own] + (iHbarPgrad[own]&(Cf[facei] - C[own] - 0.5*xii*dt));
			}
			// Debug, no = 0, =0 put to > 0
			else if (xii_Sf < -VSMALL) // comming form nei
			{
				gSurf_value[facei] = gBarPvol_nei + \
				(gBarPgrad_nei[0]*(Cf[facei*3] - C_nei[facei*3] - 0.5*xii_x*dt)) \
				+(gBarPgrad_nei[1]*(Cf[facei*3+1] - C_nei[facei*3+1] - 0.5*xii_y*dt)) \
				+(gBarPgrad_nei[2]*(Cf[facei*3+2] - C_nei[facei*3+2] - 0.5*xii_z*dt));
				hSurf_value[facei] = hBarPvol_nei + \
				(hBarPgrad_nei[0]*(Cf[facei*3] - C_nei[facei*3] - 0.5*xii_x*dt)) \
				+(hBarPgrad_nei[1]*(Cf[facei*3+1] - C_nei[facei*3+1] - 0.5*xii_y*dt)) \
				+(hBarPgrad_nei[2]*(Cf[facei*3+2] - C_nei[facei*3+2] - 0.5*xii_z*dt));
				
				// gSurf_value[facei] = iGbarPvol[nei] + (iGbarPgrad[nei]&(Cf[facei] - C[nei] - 0.5*xii*dt));
				// hSurf_value[facei] = iHbarPvol[nei] + (iHbarPgrad[nei]&(Cf[facei] - C[nei] - 0.5*xii*dt));
			}
			else 
			{
				gSurf_value[facei] = 0.5*(gBarPvol_nei + \
				(gBarPgrad_nei[0]*(Cf[facei*3] - C_nei[facei*3] - 0.5*xii_x*dt)) \
				+(gBarPgrad_nei[1]*(Cf[facei*3+1] - C_nei[facei*3+1] - 0.5*xii_y*dt)) \
				+(gBarPgrad_nei[2]*(Cf[facei*3+2] - C_nei[facei*3+2] - 0.5*xii_z*dt)) \
				+gBarPvol_own + \
				(gBarPgrad_own[0]*(Cf[facei*3] - C_own[facei*3] - 0.5*xii_x*dt)) \
				+(gBarPgrad_own[1]*(Cf[facei*3+1] - C_own[facei*3+1] - 0.5*xii_y*dt)) \
				+(gBarPgrad_own[2]*(Cf[facei*3+2] - C_own[facei*3+2] - 0.5*xii_z*dt)));
				hSurf_value[facei] = 0.5*(hBarPvol_nei + \
				(hBarPgrad_nei[0]*(Cf[facei*3] - C_nei[facei*3] - 0.5*xii_x*dt)) \
				+(hBarPgrad_nei[1]*(Cf[facei*3+1] - C_nei[facei*3+1] - 0.5*xii_y*dt)) \
				+(hBarPgrad_nei[2]*(Cf[facei*3+2] - C_nei[facei*3+2] - 0.5*xii_z*dt)) \
				+ hBarPvol_own + \
				(hBarPgrad_own[0]*(Cf[facei*3] - C_own[facei*3] - 0.5*xii_x*dt)) \
				+(hBarPgrad_own[1]*(Cf[facei*3+1] - C_own[facei*3+1] - 0.5*xii_y*dt)) \
				+(hBarPgrad_own[2]*(Cf[facei*3+2] - C_own[facei*3+2] - 0.5*xii_z*dt)));
				
				// gSurf_value[facei] = 0.5*(iGbarPvol[nei] + ((iGbarPgrad[nei]) &(Cf[facei] - C[nei] - 0.5*xii*dt)) + iGbarPvol[own] + ((iGbarPgrad[own]) &(Cf[facei] - C[own] - 0.5*xii*dt)));
				// hSurf_value[facei] = 0.5*(iHbarPvol[nei] + ((iHbarPgrad[nei]) &(Cf[facei] - C[nei] - 0.5*xii*dt)) + iHbarPvol[own] + ((iHbarPgrad[own]) &(Cf[facei] - C[own] - 0.5*xii*dt)));
			}
		}
		put_reply = 0; 
		athread_put(PE_MODE,gSurf_value,slave.gSurf_value+stid+stnum, sizeof(double)*num, &put_reply , 0, 0);
		athread_put(PE_MODE,hSurf_value,slave.hSurf_value+stid+stnum, sizeof(double)*num, &put_reply , 0, 0);             
		while(put_reply!=2);
		step++;
	}

}


void Func_gaussgrad(GaussGrad *gaussgrad) {
	int N_slave = 320;
	GaussGrad slave;
	// unsigned long get_reply,put_reply; 
    int my_id; 
    my_id = athread_get_id(-1);  //获取从核的id号
    get_reply = 0;
    athread_get(PE_MODE , gaussgrad, &slave, sizeof(GaussGrad), &get_reply , 0, 0, 0);
    while (get_reply != 1);//等待数据加载完成
	int64 ownersize = gaussgrad->ownersize;
	// Load Balance
	int load = ownersize / NTHREAD;
	int rest = ownersize % NTHREAD;
	int stid = (my_id < rest) ? my_id * (load + 1) : my_id * load + rest;
	if (my_id < rest) load++;
	int edid = stid + load;
	int numtemp = edid - stid;
	//Define variable
	double Sf[N_slave*3];
	double gBarPvol_value_own[N_slave];
	double gBarPvol_value_nei[N_slave];
	double hBarPvol_value_own[N_slave];
	double hBarPvol_value_nei[N_slave];
	double gBarPgrad_value_own[N_slave*3];
	double gBarPgrad_value_nei[N_slave*3];
	double hBarPgrad_value_own[N_slave*3];
	double hBarPgrad_value_nei[N_slave*3];
	double V_own[N_slave];
	double V_nei[N_slave];
	double interpola_weight[N_slave];
	int facei = 0;
	int num = 0;
	int step = 0;
	int stnum = 0;
	double xii_Sf = 0;
	// 	unsigned long time1 = 0;
	// unsigned long time_sta = 0;
	// unsigned long time_end = 0;
	// time_sta=rpcc();
	double gsfi=0.0;
	double hsfi=0.0;
	while(numtemp>0){
		num = (numtemp < N_slave) ? numtemp : N_slave;
		numtemp = numtemp - N_slave;
		stnum = step * N_slave;
		//get data
		get_reply = 0;
		athread_get(PE_MODE , slave.Sf_value+stid*3+stnum*3, Sf, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.gBarPvol_value_own+stid+stnum, gBarPvol_value_own, sizeof(double)*num, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.gBarPvol_value_nei+stid+stnum, gBarPvol_value_nei, sizeof(double)*num, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.hBarPvol_value_own+stid+stnum, hBarPvol_value_own, sizeof(double)*num, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.hBarPvol_value_nei+stid+stnum, hBarPvol_value_nei, sizeof(double)*num, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.V_own+stid+stnum, V_own, sizeof(double)*num, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.V_nei+stid+stnum, V_nei, sizeof(double)*num, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.interpola_weight+stid+stnum, interpola_weight, sizeof(double)*num, &get_reply , 0, 0, 0);
    	while (get_reply != 8);

		memset(gBarPgrad_value_own,0,sizeof(double)*N_slave*3);
		memset(gBarPgrad_value_nei,0,sizeof(double)*N_slave*3);
		memset(hBarPgrad_value_own,0,sizeof(double)*N_slave*3);
		memset(hBarPgrad_value_nei,0,sizeof(double)*N_slave*3);
    	// forAll(owner, facei)
    	// {
        // 	label own = owner[facei];
        // 	label nei = neighbour[facei];
        // 	scalar gsfi = interpola_weight[facei]*(gvfi[own] - gvfi[nei]) + gvfi[nei];
        // 	scalar hsfi = interpola_weight[facei]*(hvfi[own] - hvfi[nei]) + hvfi[nei];
        // 	igGrad[own] += Sf[facei]*gsfi/V[own];
        // 	igGrad[nei] -= Sf[facei]*gsfi/V[nei];
        // 	ihGrad[own] += Sf[facei]*hsfi/V[own];
        // 	ihGrad[nei] -= Sf[facei]*hsfi/V[nei];
    	// }
		for(facei=0;facei<num;facei++)
		{
			gsfi = interpola_weight[facei]*(gBarPvol_value_own[facei] - gBarPvol_value_nei[facei]) + gBarPvol_value_nei[facei];
        	hsfi = interpola_weight[facei]*(hBarPvol_value_own[facei] - hBarPvol_value_nei[facei]) + hBarPvol_value_nei[facei];
        	gBarPgrad_value_own[facei*3]   += Sf[facei*3]  *gsfi/V_own[facei];
			gBarPgrad_value_own[facei*3+1] += Sf[facei*3+1]*gsfi/V_own[facei];
			gBarPgrad_value_own[facei*3+2] += Sf[facei*3+2]*gsfi/V_own[facei];
        	gBarPgrad_value_nei[facei*3]   -= Sf[facei*3]  *gsfi/V_nei[facei];
        	gBarPgrad_value_nei[facei*3+1] -= Sf[facei*3+1]*gsfi/V_nei[facei];
        	gBarPgrad_value_nei[facei*3+2] -= Sf[facei*3+2]*gsfi/V_nei[facei];
        	hBarPgrad_value_own[facei*3]   += Sf[facei*3]  *hsfi/V_own[facei];
			hBarPgrad_value_own[facei*3+1] += Sf[facei*3+1]*hsfi/V_own[facei];
			hBarPgrad_value_own[facei*3+2] += Sf[facei*3+2]*hsfi/V_own[facei];
        	hBarPgrad_value_nei[facei*3]   -= Sf[facei*3]  *hsfi/V_nei[facei];
        	hBarPgrad_value_nei[facei*3+1] -= Sf[facei*3+1]*hsfi/V_nei[facei];
        	hBarPgrad_value_nei[facei*3+2] -= Sf[facei*3+2]*hsfi/V_nei[facei];

    	}
		put_reply = 0; 
    	athread_put(PE_MODE,gBarPgrad_value_own,slave.gBarPgrad_value_own+stid*3+stnum*3,sizeof(double)*num*3,&put_reply,0,0);                 
  		athread_put(PE_MODE,gBarPgrad_value_nei,slave.gBarPgrad_value_nei+stid*3+stnum*3,sizeof(double)*num*3,&put_reply,0,0);                 
  		athread_put(PE_MODE,hBarPgrad_value_own,slave.hBarPgrad_value_own+stid*3+stnum*3,sizeof(double)*num*3,&put_reply,0,0);                 
  		athread_put(PE_MODE,hBarPgrad_value_nei,slave.hBarPgrad_value_nei+stid*3+stnum*3,sizeof(double)*num*3,&put_reply,0,0);               
		while(put_reply!=4);
		step++;
	}

}

void Func_ghbarsurf(GHbarSurf *ghbarsurf) {
	int N_slave = 256;
	GHbarSurf slave;
	// unsigned long get_reply,put_reply; 
    int my_id; 
    my_id = athread_get_id(-1);  //获取从核的id号
    get_reply = 0;
    athread_get(PE_MODE , ghbarsurf, &slave, sizeof(GHbarSurf), &get_reply , 0, 0, 0);
    while (get_reply != 1);//等待数据加载完成
	int64 ownersize = ghbarsurf->ownersize;
	// Load Balance
	int load = ownersize / NTHREAD;
	int rest = ownersize % NTHREAD;
	int stid = (my_id < rest) ? my_id * (load + 1) : my_id * load + rest;
	if (my_id < rest) load++;
	int edid = stid + load;
	int numtemp = edid - stid;
	//Define variable
	double Sf[N_slave*3];
	double Cf[N_slave*3];
	double gBarPvol_value_own[N_slave];
	double gBarPvol_value_nei[N_slave];
	double hBarPvol_value_own[N_slave];
	double hBarPvol_value_nei[N_slave];
	double gBarPgrad_value_own[N_slave*3];
	double gBarPgrad_value_nei[N_slave*3];
	double hBarPgrad_value_own[N_slave*3];
	double hBarPgrad_value_nei[N_slave*3];
	double C_own[N_slave*3];
	double C_nei[N_slave*3];
	double gSurf_value[N_slave];
	double hSurf_value[N_slave];
	double dt = ghbarsurf->dt;
	double xii_x = ghbarsurf->xii_x;
	double xii_y = ghbarsurf->xii_y;
	double xii_z = ghbarsurf->xii_z;
	int facei = 0;
	int num = 0;
	int step = 0;
	int stnum = 0;
	double xii_Sf = 0;
	// 	unsigned long time1 = 0;
	// unsigned long time_sta = 0;
	// unsigned long time_end = 0;
	// time_sta=rpcc();
	while(numtemp>0){
		num = (numtemp < N_slave) ? numtemp : N_slave;
		numtemp = numtemp - N_slave;
		stnum = step * N_slave;
		//get data
		get_reply = 0;
		athread_get(PE_MODE , slave.Sf+stid*3+stnum*3, Sf, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.Cf+stid*3+stnum*3, Cf, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.gBarPvol_value_own+stid+stnum, gBarPvol_value_own, sizeof(double)*num, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.gBarPvol_value_nei+stid+stnum, gBarPvol_value_nei, sizeof(double)*num, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.hBarPvol_value_own+stid+stnum, hBarPvol_value_own, sizeof(double)*num, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.hBarPvol_value_nei+stid+stnum, hBarPvol_value_nei, sizeof(double)*num, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.gBarPgrad_value_own+stid*3+stnum*3, gBarPgrad_value_own, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.gBarPgrad_value_nei+stid*3+stnum*3, gBarPgrad_value_nei, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.hBarPgrad_value_own+stid*3+stnum*3, hBarPgrad_value_own, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.hBarPgrad_value_nei+stid*3+stnum*3, hBarPgrad_value_nei, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.C_own+stid*3+stnum*3, C_own, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		athread_get(PE_MODE , slave.C_nei+stid*3+stnum*3, C_nei, sizeof(double)*num*3, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.gSurf_value+stid+stnum, gSurf_value, sizeof(double)*num, &get_reply , 0, 0, 0);
		// athread_get(PE_MODE , slave.hSurf_value+stid+stnum, hSurf_value, sizeof(double)*num, &get_reply , 0, 0, 0);
    	while (get_reply != 12);
		for(facei=0;facei<num;facei++)
		{
			xii_Sf = xii_x * Sf[facei*3] + xii_y * Sf[facei*3+1] + xii_z * Sf[facei*3+2];
			if (xii_Sf >=  VSMALL) // comming from own
			{
				gSurf_value[facei] = gBarPvol_value_own[facei] + \
				(gBarPgrad_value_own[facei*3]*(Cf[facei*3] - C_own[facei*3] - 0.5*xii_x*dt)) \
				+(gBarPgrad_value_own[facei*3+1]*(Cf[facei*3+1] - C_own[facei*3+1] - 0.5*xii_y*dt)) \
				+(gBarPgrad_value_own[facei*3+2]*(Cf[facei*3+2] - C_own[facei*3+2] - 0.5*xii_z*dt));
				hSurf_value[facei] = hBarPvol_value_own[facei] + \
				(hBarPgrad_value_own[facei*3]*(Cf[facei*3] - C_own[facei*3] - 0.5*xii_x*dt)) \
				+(hBarPgrad_value_own[facei*3+1]*(Cf[facei*3+1] - C_own[facei*3+1] - 0.5*xii_y*dt)) \
				+(hBarPgrad_value_own[facei*3+2]*(Cf[facei*3+2] - C_own[facei*3+2] - 0.5*xii_z*dt));
				// iGsurf[facei] = iGbarPvol[own] + (iGbarPgrad[own]&(Cf[facei] - C[own] - 0.5*xii*dt));
            	// iHsurf[facei] = iHbarPvol[own] + (iHbarPgrad[own]&(Cf[facei] - C[own] - 0.5*xii*dt));
			}
			// Debug, no = 0, =0 put to > 0
			else if (xii_Sf < -VSMALL) // comming form nei
			{
				gSurf_value[facei] = gBarPvol_value_nei[facei] + \
				(gBarPgrad_value_nei[facei*3]*(Cf[facei*3] - C_nei[facei*3] - 0.5*xii_x*dt)) \
				+(gBarPgrad_value_nei[facei*3+1]*(Cf[facei*3+1] - C_nei[facei*3+1] - 0.5*xii_y*dt)) \
				+(gBarPgrad_value_nei[facei*3+2]*(Cf[facei*3+2] - C_nei[facei*3+2] - 0.5*xii_z*dt));
				hSurf_value[facei] = hBarPvol_value_nei[facei] + \
				(hBarPgrad_value_nei[facei*3]*(Cf[facei*3] - C_nei[facei*3] - 0.5*xii_x*dt)) \
				+(hBarPgrad_value_nei[facei*3+1]*(Cf[facei*3+1] - C_nei[facei*3+1] - 0.5*xii_y*dt)) \
				+(hBarPgrad_value_nei[facei*3+2]*(Cf[facei*3+2] - C_nei[facei*3+2] - 0.5*xii_z*dt));
				
				// gSurf_value[facei] = iGbarPvol[nei] + (iGbarPgrad[nei]&(Cf[facei] - C[nei] - 0.5*xii*dt));
				// hSurf_value[facei] = iHbarPvol[nei] + (iHbarPgrad[nei]&(Cf[facei] - C[nei] - 0.5*xii*dt));
			}
			else 
			{
				gSurf_value[facei] = 0.5*(gBarPvol_value_nei[facei] + \
				(gBarPgrad_value_nei[facei*3]*(Cf[facei*3] - C_nei[facei*3] - 0.5*xii_x*dt)) \
				+(gBarPgrad_value_nei[facei*3+1]*(Cf[facei*3+1] - C_nei[facei*3+1] - 0.5*xii_y*dt)) \
				+(gBarPgrad_value_nei[facei*3+2]*(Cf[facei*3+2] - C_nei[facei*3+2] - 0.5*xii_z*dt)) \
				+gBarPvol_value_own[facei] + \
				(gBarPgrad_value_own[facei*3]*(Cf[facei*3] - C_own[facei*3] - 0.5*xii_x*dt)) \
				+(gBarPgrad_value_own[facei*3+1]*(Cf[facei*3+1] - C_own[facei*3+1] - 0.5*xii_y*dt)) \
				+(gBarPgrad_value_own[facei*3+2]*(Cf[facei*3+2] - C_own[facei*3+2] - 0.5*xii_z*dt)));
				hSurf_value[facei] = 0.5*(hBarPvol_value_nei[facei] + \
				(hBarPgrad_value_nei[facei*3]*(Cf[facei*3] - C_nei[facei*3] - 0.5*xii_x*dt)) \
				+(hBarPgrad_value_nei[facei*3+1]*(Cf[facei*3+1] - C_nei[facei*3+1] - 0.5*xii_y*dt)) \
				+(hBarPgrad_value_nei[facei*3+2]*(Cf[facei*3+2] - C_nei[facei*3+2] - 0.5*xii_z*dt)) \
				+ hBarPvol_value_own[facei] + \
				(hBarPgrad_value_own[facei*3]*(Cf[facei*3] - C_own[facei*3] - 0.5*xii_x*dt)) \
				+(hBarPgrad_value_own[facei*3+1]*(Cf[facei*3+1] - C_own[facei*3+1] - 0.5*xii_y*dt)) \
				+(hBarPgrad_value_own[facei*3+2]*(Cf[facei*3+2] - C_own[facei*3+2] - 0.5*xii_z*dt)));
				
				// gSurf_value[facei] = 0.5*(iGbarPvol[nei] + ((iGbarPgrad[nei]) &(Cf[facei] - C[nei] - 0.5*xii*dt)) + iGbarPvol[own] + ((iGbarPgrad[own]) &(Cf[facei] - C[own] - 0.5*xii*dt)));
				// hSurf_value[facei] = 0.5*(iHbarPvol[nei] + ((iHbarPgrad[nei]) &(Cf[facei] - C[nei] - 0.5*xii*dt)) + iHbarPvol[own] + ((iHbarPgrad[own]) &(Cf[facei] - C[own] - 0.5*xii*dt)));
			}
		}
		put_reply = 0; 
		athread_put(PE_MODE,gSurf_value,slave.gSurf_value+stid+stnum, sizeof(double)*num, &put_reply , 0, 0);
		athread_put(PE_MODE,hSurf_value,slave.hSurf_value+stid+stnum, sizeof(double)*num, &put_reply , 0, 0);             
		while(put_reply!=2);
		step++;
	}

}