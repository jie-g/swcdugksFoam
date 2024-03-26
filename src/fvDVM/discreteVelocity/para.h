#define ATHREAD_TEMP
#define GHBARPVOL//计算单元辅助函数
#define MACROSURF//计算界面宏观量
#define GHSURF//计算界面分布函数，与MACROSURF要一起打开
#define GHTILDEVOL//更新单元分布函数
#define MACROVOL//更新单元宏观量

#define GHBARSURF
// /* ghbarpvol_slave.c */
// // #define GHBARPVOL_USE_SIMD
// #define GHBARPVOL_TRAN_FUN
// #define GHBARPVOL_REG_COM

// /* macrosurf_slave.c */
// // #define MACROSURF_USE_SIMD
// #define MACROSURF_TRAN_FUN
// #define MACROSURF_REG_COM

// /* ghsurf_slave.c */
// // #define GHSURF_USE_SIMD
// #define GHSURF_TRAN_FUN

// /* ghtildevol_slave.c */
// // #define GHTILDEVOL_USE_SIMD
// // #define GHTILDEVOL_TRAN_FUN
// // #define GHTILDEVOL_REG_COM

// /* macrovol_slave.c */
// // #define MACROVOL_USE_SIMD
// // #define MACROVOL_TRAN_FUN
// // #define MACROVOL_REG_COM