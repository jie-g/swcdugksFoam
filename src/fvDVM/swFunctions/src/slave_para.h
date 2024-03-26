typedef long int int64;//label
/*************fvDVM.C*****************/
// typedef struct VectorData{
// 	 double vectordata[3]; 
// }VectorData;
// typedef struct MacroVol{
//     // int64 *owner;
//     // int64 *neighbour;
// 	double *Sf;
// 	double *rhoflux_value_own;
// 	double *rhoflux_value_nei;
// 	double *rhouflux_value_own;
// 	double *rhouflux_value_nei;
// 	double *rhoeflux_value_own;
// 	double *rhoeflux_value_nei;
// 	// double *rhouflux_value;
// //	double *rho_value;
// 	// double *e_value;
// //	double *u_value;
// 	// double *V;
// 	double *V_own;
// 	double *V_nei;
// 	double *gSurf_value;
// 	double *hSurf_value;
// 	// double *xii_value;
// 	double dt;
// 	double xii_x;
// 	double xii_y;
// 	double xii_z;
// 	double weight_value;
// 	int64 ownersize;
// }MacroVol;
typedef struct MacroVol{
	double **gSurf_value;
	double **hSurf_value;
	double *Sf;
	double *rhoflux_value_own;
	double *rhoflux_value_nei;
	double *rhouflux_value_own;
	double *rhouflux_value_nei;
	double *rhoeflux_value_own;
	double *rhoeflux_value_nei;
	double *V_own;
	double *V_nei;
	double *xi_value;
	double *weight_value;
	double dt;
	int64 ownersize;
	int64 DVsize;
}MacroVol;
typedef struct MacroSurf{
	double *rhoSurf_value;
	double *rhoUsurf_value;
	// VectorData *rhoUsurf_value;
	double *rhoEsurf_value;
	double *Usurf_value;
	double **gSurf_value;
	double **hSurf_value;
	double *xi_value;
	double *weight_value;
	int64 ownersize;
	int64 DVsize;
}MacroSurf;
// typedef struct MacroSurf{
// 	double *rhoSurf_value;
// 	double *rhoUsurf_value;
// 	// VectorData *rhoUsurf_value;
// 	double *rhoEsurf_value;
// 	double *gSurf_value;
// 	double *hSurf_value;
// 	double xii_x;
// 	double xii_y;
// 	double xii_z;
// 	int64 ownersize;
// 	double weight_value;
// }MacroSurf;
// typedef struct MacroqSurf{
// 	double *qSurf_value;
// 	double *Usurf_value;
// 	double *gSurf_value;
// 	double *hSurf_value;
// 	double xii_x;
// 	double xii_y;
// 	double xii_z;
// 	int64 ownersize;
// 	double weight_value;
// }MacroqSurf;
typedef struct MacroqSurf{
	double **gSurf_value;
	double **hSurf_value;
	double *rhoSurf_value;
	double *rhoUsurf_value;
	double *rhoEsurf_value;
	double *qSurf_value;
	double *xi_value;
	double *weight_value;
	double *Tsurf_value;
	double *Usurf_value;
	double *tauSurf_value;
	double *qSurf_temp_value;
	double R_value;
	double muRef_value;
	double Tref_value;
	double Pr_value;
	double omega_value;
	double dt;
	int64 KInner_value;
	int64 ownersize;
	int64 DVsize;
}MacroqSurf;
typedef struct MacroqVol{
	double **gTildeVol_value;
	double **hTildeVol_value;
	double *rhoVol_value;
	// double *rhoUvol_value;
	// double *rhoEvol_value;
	double *qVol_value;
	double *xi_value;
	double *weight_value;
	double *Tvol_value;
	double *Uvol_value;
	double *tauVol_value;
	double *qVol_temp_value;

	double *rhoflux_value;
	double *rhouflux_value;
	double *rhoeflux_value;

	double R_value;
	double muRef_value;
	double Tref_value;
	double Pr_value;
	double omega_value;
	double dt;
	int64 KInner_value;
	int64 qVol_size;
	int64 DVsize;
}MacroqVol;
// typedef struct GHtildeVol{
// 	double **gTildeVol_value_own;//ownersize
// 	double **gTildeVol_value_nei;//ownersize
// 	double **hTildeVol_value_own;//ownersize
// 	double **hTildeVol_value_nei;//ownersize
// 	double **gSurf_value;//ownersize
// 	double **hSurf_value;//ownersize
// 	double **xi_value;
// 	double *V_own;//ownersize
// 	double *V_nei;//ownersize
// 	double *Sf;//ownersize*3
// 	double dt;
// 	int64 ownersize;
// 	int64 DVsize;
// }GHtildeVol;//
typedef struct GHsurf{
    double **gSurf_value;
    double **hSurf_value;
    double *U_value;
    double *q_value;
    double *T_value;
    double *rho_value;
    double *tauSurf_value;
	double *xi_value;
	// double *gSurf_value;
    // double *hSurf_value;
    // double *relaxFactor;
    double h;
    // double xii_x;
	// double xii_y;
	// double xii_z;
	double R_value;
	double Pr_value;
	double vUnit_value;
    int64 gSurf_size;
    int64 D;//label
	int64 K;//label
	int64 DVsize;
}GHsurf;
typedef struct GHbarPvol{
    double **gBarPvol_value;
    double **hBarPvol_value;
	double **gTildeVol_value;
    double **hTildeVol_value;
    double *U_value;
    double *q_value;
    double *T_value;
    double *rho_value;
    double *tauVol_value;
	double *xi_value;
	// double *gSurf_value;
    // double *hSurf_value;
    // double *relaxFactor;
    double dt;
    // double xii_x;
	// double xii_y;
	// double xii_z;
	double R_value;
	double Pr_value;
	double vUnit_value;
    int64 ghBarPvol_size;
    int64 D;//label
	int64 K;//label
	int64 DVsize;
}GHbarPvol;
/*************discreteVelocity.C*****************/
typedef struct EquShakhov{
	double *U_value;
	double *q_value;
	double *T_value;
	double *rho_value;
	double *gEq_value;
	double *hEq_value;
	double xii_x;
	double xii_y;
	double xii_z;
	double R_value;
	double Pr_value;
	double vUnit_value;
    int64 D;//label
	int64 K;//label
    int64 gEq_size;
}EquShakhov;
typedef struct GHtildeVol{
	double *Sf;//ownersize*3
	double *gTildeVol_value_own;//ownersize
	double *gTildeVol_value_nei;//ownersize
	double *hTildeVol_value_own;//ownersize
	double *hTildeVol_value_nei;//ownersize
	double *V_own;//ownersize
	double *V_nei;//ownersize
	double *gSurf_value;//ownersize
	double *hSurf_value;//ownersize
	double dt;
	double xii_x;
	double xii_y;
	double xii_z;
	int64 ownersize;
}GHtildeVol;//
// typedef struct GHbarPvol{
//     double *tauVol_value;
//     double *gEq_value;
//     double *hEq_value;
//     double *gTildeVol_value;
//     double *hTildeVol_value;
//     double *gBarPvol_value;
//     double *hBarPvol_value;
//     double *relaxFactor;
//     double dt;
//     int64 size;
// }GHbarPvol;
// typedef struct GHsurf{
//     double *tauSurf_value;
//     double *gEq_value;
//     double *hEq_value;
//     double *gSurf_value;
//     double *hSurf_value;
//     double *relaxFactor;
//     double h;
//     int64 size;
// }GHsurf;
// typedef struct GHsurf{
//     double *U_value;
//     double *q_value;
//     double *T_value;
//     double *rho_value;
//     double *tauSurf_value;
//     double *gSurf_value;
//     double *hSurf_value;
//     // double *relaxFactor;
//     double h;
//     double xii_x;
// 	double xii_y;
// 	double xii_z;
// 	double R_value;
// 	double Pr_value;
// 	double vUnit_value;
//     int64 gSurf_size;
//     int64 D;//label
// 	int64 K;//label
// }GHsurf;
typedef struct equMaxwell{
    double *U_value;
	double *T_value;
    double *rho_value;
	double *geq_value;
	double *heq_value;
	double xii_x;
	double xii_y;
	double xii_z;
	double Ri;//scalar
    int64 rho_size;//lable
    int64 D;//label
	int64 K;//label
}equMaxwell;
// typedef struct VectorData{
// 	 double vectordata[3]; 
// }VectorData;
typedef struct GHbarSurf{
	double *Sf;
	double *Cf;
	double *gBarPvol_value_own;
	double *gBarPvol_value_nei;
	double *hBarPvol_value_own;
	double *hBarPvol_value_nei;
	double *gBarPgrad_value_own;
	double *gBarPgrad_value_nei;
	double *hBarPgrad_value_own;
	double *hBarPgrad_value_nei;
	double *C_own;
	double *C_nei;
	double *gSurf_value;
	double *hSurf_value;

	double dt;
	double xii_x;
	double xii_y;
	double xii_z;
	int64 ownersize;
}GHbarSurf;
typedef struct GaussGrad{
	double *Sf_value;
	double *gBarPvol_value_own;
	double *gBarPvol_value_nei;
	double *hBarPvol_value_own;
	double *hBarPvol_value_nei;
	double *gBarPgrad_value_own;
	double *gBarPgrad_value_nei;
	double *hBarPgrad_value_own;
	double *hBarPgrad_value_nei;
	double *V_own;
	double *V_nei;
	double *interpola_weight;
	// /************1219******************/
	//     double *Cf;
	// 	double *C_own;
	// 	double *C_nei;
	// 	double *gSurf_value;
	// 	double *hSurf_value;
	// 	double dt;
	// 	double xii_x;
	// 	double xii_y;
	// 	double xii_z;
	// /**************************/
	int64 ownersize;
}GaussGrad;

// typedef struct Grad{
// 	double *gBarPgrad_value_own;
// 	double *gBarPgrad_value_nei;
// 	double *hBarPgrad_value_own;
// 	double *hBarPgrad_value_nei;

// 	double *ownLs_value;
// 	double *neiLs_value;
// 	double *gBarPvol_value_own;
// 	double *gBarPvol_value_nei;
// 	double *hBarPvol_value_own;
// 	double *hBarPvol_value_nei;

//     // double *gBarPgrad_value;
// 	// double *hBarPgrad_value;
//     // int64 *owner_value;
// 	// int64 *neighbour_value;
	
// 	int64 ownersize;
// }Grad;
typedef struct CONUM{
	double *UbyDx_value;
	double *Usurf_value;
	double *C_own;
	double *C_nei;
	double dt;
	double xiMax_value;
	int64 D;
	int64 ownersize;
}CONUM;

typedef struct GaussGrad1{
	double *Sf_value;
	double *gBarPvol_value;
	// double *gBarPvol_value_nei;
	double *hBarPvol_value;
	// double *hBarPvol_value_nei;
	double *gBarPgrad_value_own;
	double *gBarPgrad_value_nei;
	double *hBarPgrad_value_own;
	double *hBarPgrad_value_nei;
	double *V_own;
	double *V_nei;
	int64 *own;
	int64 *neighbour;
	double *interpola_weight;
	int64 ownersize;
	int64 gBarPvolsize;
}GaussGrad1;
typedef struct GHbarSurf1{
	double *Sf;
	double *Cf;
	double *gBarPvol_value;
	double *hBarPvol_value;
	// double *gBarPvol_value_own;
	// double *gBarPvol_value_nei;
	// double *hBarPvol_value_own;
	// double *hBarPvol_value_nei;
	double *gBarPgrad_value_own;
	double *gBarPgrad_value_nei;
	double *hBarPgrad_value_own;
	double *hBarPgrad_value_nei;
	
    // double *gBarPgrad_value;
	// double *hBarPgrad_value;
	double *C_own;
	double *C_nei;
	double *gSurf_value;
	double *hSurf_value;
	int64 *own;
	int64 *neighbour;
	double dt;
	double xii_x;
	double xii_y;
	double xii_z;
	int64 ownersize;
	int64 gBarPvolsize;
}GHbarSurf1;
typedef struct GHbarSurf2{
	double *Sf;
	double *Cf;
	double *gBarPvol_value;
	double *hBarPvol_value;
	// double *gBarPvol_value_own;
	// double *gBarPvol_value_nei;
	// double *hBarPvol_value_own;
	// double *hBarPvol_value_nei;
	// double *gBarPgrad_value_own;
	// double *gBarPgrad_value_nei;
	// double *hBarPgrad_value_own;
	// double *hBarPgrad_value_nei;
	
    double *gBarPgrad_value;
	double *hBarPgrad_value;
	double *C_own;
	double *C_nei;
	double *gSurf_value;
	double *hSurf_value;
	int64 *own;
	int64 *neighbour;
	double dt;
	double xii_x;
	double xii_y;
	double xii_z;
	int64 ownersize;
	int64 gBarPvolsize;
}GHbarSurf2;