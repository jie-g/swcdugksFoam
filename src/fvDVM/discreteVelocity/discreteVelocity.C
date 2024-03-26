/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "foam_defs.h"
#include "fvDVM.H"
#include <map>
#include "discreteVelocity.H"
#include "constants.H"
#include "fixedGradientFvPatchField.H"
#include "calculatedMaxwellFvPatchField.H"
#include "symmetryModFvPatchField.H"
#include "DVMsymmetryFvsPatchField.H"

#include "leastSquaresVectors.H"
#include "gaussGrad.H"
#include "surfaceInterpolate.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// #define GHBARSURF


typedef long int int64;//label
//const Foam::word
//Foam::radiation::radiativeIntensityRay::intensityPrefix("ILambda");

#if FOAM_MAJOR <= 3
    #define BOUNDARY_FIELD_REF boundaryField()
#else
    #define BOUNDARY_FIELD_REF boundaryFieldRef()
#endif

static inline unsigned long rpcc(){
	unsigned long addtime;
	asm("rtc %0": "=r" (addtime) : );
	return addtime;
 }
#define CLOCKRATE 1.45e9

#include "para.h"

#ifdef ATHREAD_TEMP
#include "slave_para.h"
extern "C"{
	#include <athread.h>
    void SLAVE_FUN(Func_ghtildevol)(GHtildeVol *ghtildevol);
    void SLAVE_FUN(Func_gaussgrad)(GaussGrad *gaussgrad);
    void SLAVE_FUN(Func_gaussgradtemp)(GaussGrad1 *gaussgrad);
    void SLAVE_FUN(Func_ghbarsurf)(GHbarSurf *ghbarsurf);
    void SLAVE_FUN(Func_ghbarsurftemp)(GHbarSurf1 *ghbarsurf);
    void SLAVE_FUN(Func_ghbarsurftemp2)(GHbarSurf2 *ghbarsurf);
}
#endif
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::discreteVelocity::discreteVelocity
(
    fvDVM& dvm,
    const fvMesh& mesh,
    const Time& time,
    const scalar weight,
    const dimensionedVector xi,
    const label DVid,
    const label symXtargetDVid,
    const label symYtargetDVid,
    const label symZtargetDVid
)
:
    dvm_(dvm),
    mesh_(mesh),
    time_(time),
    weight_(weight),
    xi_(xi),
    myDVid_(DVid),
    symXtargetDVid_(symXtargetDVid),
    symYtargetDVid_(symYtargetDVid),
    symZtargetDVid_(symZtargetDVid),
    gTildeVol_
    (
        IOobject
        (
            "gTildeVol" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "0", dimMass*pow3(dimTime/dimLength)/pow3(dimLength), 0.0
        )
    ),
    hTildeVol_
    (
        IOobject
        (
            "hTildeVol" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "0", dimMass*(dimTime/dimLength)/pow3(dimLength), 0.0
        )
    ),
    gBarPvol_
    (
        IOobject
        (
            "gBarPvol" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "0", dimMass*pow3(dimTime/dimLength)/pow3(dimLength), 0.0
        ),
        "fixedGradient"
    ),
    hBarPvol_
    (
        IOobject
        (
            "hBarPvol" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "0", dimMass*(dimTime/dimLength)/pow3(dimLength), 0.0
        ),
        "fixedGradient"
    ),
    gSurf_
    (
        IOobject
        (
            "gSurf" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "0", dimMass*pow3(dimTime/dimLength)/pow3(dimLength), 0.0
        )
    ),
    hSurf_
    (
        IOobject
        (
            "hSurf" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "0", dimMass*(dimTime/dimLength)/pow3(dimLength), 0.0
        )
    ),
    gBarPgrad_
    (
        IOobject
        (
            "gBarPgrad" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "0", gBarPvol_.dimensions()/dimLength, vector(0,0,0)
        ),
        "zeroGradient"
    ),
    hBarPgrad_
    (
        IOobject
        (
            "hBarPgrad" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "0", hBarPvol_.dimensions()/dimLength, vector(0,0,0)
        ),
        "zeroGradient"
    )
    ,
    tinterpScheme_(new linear<scalar>(mesh_))
{
    initDFtoEq();
    setBCtype();
    initBoundaryField();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::discreteVelocity::~discreteVelocity()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::discreteVelocity::initDFtoEq()
{
    GeometricField<vector, fvPatchField, volMesh> qVolIni
    (
        IOobject
        (
            "qVolIni",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "0", dvm_.qVol().dimensions(), vector(0,0,0)
        )
    );

    equilibriumShakhovVol
    (
        gTildeVol_,
        hTildeVol_,
        dvm_.rhoVol(), 
        dvm_.Uvol(), 
        dvm_.Tvol(), 
        qVolIni
    );

}

void Foam::discreteVelocity::setBCtype()
{
    // Only from rho's B.C can we determine surfaceScalarField f's B.C.
    // for all patchi fo h/g barP, set to zeroGradient
    // May improve for inlet ingoing DF

    // NOTE: symmetryBoundary condition used for VSD (velocity space decompose) parallel run.
    // For normal PSD (physical velocity space decomposition) parallel run, 
    // the B.C. for symmetry Boundary are just symmetryPlane type,  and since symmetryPlane 
    // for the rho Boundary is  a geometry constraint B.C. type, the coressponding B.C. type for 
    // gSurf_ and hSurf_ will be automatically assigned to be symmetryPlane.
    //
    // NOTE: farField for rho BC is not currently used

#if FOAM_MAJOR <= 3
    const GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = dvm_.rhoVol().boundaryField();
#else
    const GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
        rhoBCs = dvm_.rhoVol().BOUNDARY_FIELD_REF;
#endif

    // bondary condition type map from rho BC to g/hSuf BC
    std::map<word, word> bcMap;
    bcMap["fixedValue"] = "mixed"; // For supersonic out boundary only
    bcMap["zeroGradient"] = "zeroGradient"; 
    bcMap["calculatedMaxwell"] = "maxwellWall"; // Maxwell wall
    bcMap["farField"] = "farField"; // Incoming DF are EQ, with undering macro possiblly updated with time
    bcMap["symmetryMod"] = "DVMsymmetry"; // used in VSD
    bcMap["pressureIn"] = "farField"; // Pressure inlet
    bcMap["pressureOut"] = "farField"; // Pressure outlet

    // geometry constraint boundary are not specifid here,
    // they are automatically assigned, including:
    // cyclic, processor, processorCyclic, symmetryPlane

    forAll(rhoBCs, patchi)
    {
        word rhoBCtype = rhoBCs[patchi].type();
        if( bcMap.find(rhoBCtype) != bcMap.end() ) // found
        {
            gSurf_.BOUNDARY_FIELD_REF.set
            (
                patchi, 
                fvsPatchField<scalar>::New
                (
                    bcMap[rhoBCtype], mesh_.boundary()[patchi], gSurf_
                )
            );
            hSurf_.BOUNDARY_FIELD_REF.set
            (
                patchi, 
                fvsPatchField<scalar>::New
                (
                    bcMap[rhoBCtype], mesh_.boundary()[patchi], hSurf_
                )
            );
        }
    }
}

void Foam::discreteVelocity::initBoundaryField()
{
    unsigned long time0 = 0;
	unsigned long time1 = 0;
    time0 = rpcc();
    // fvsPatchField gSurf_ hSurf_  with "mixed type" init field
    // Here we only consider g but not h, since the eq of h is zero, which is 
    // the default value of mixedFvsPatchField
    // Bug: h is not zero
#if FOAM_MAJOR <= 3
    GeometricField<scalar, fvsPatchField, surfaceMesh>::GeometricBoundaryField& 
        gBCs = gSurf_.boundaryField();
    GeometricField<scalar, fvsPatchField, surfaceMesh>::GeometricBoundaryField& 
        hBCs = hSurf_.boundaryField();
#else
    GeometricField<scalar, fvsPatchField, surfaceMesh>::Boundary& 
        gBCs = gSurf_.BOUNDARY_FIELD_REF;
    GeometricField<scalar, fvsPatchField, surfaceMesh>::Boundary& 
        hBCs = hSurf_.BOUNDARY_FIELD_REF;
#endif

    forAll(gBCs, patchi)
    {
        if(gBCs[patchi].type() == "mixed" )
        {
            equilibriumMaxwell
            (
                gBCs[patchi],
                hBCs[patchi],
                dvm_.rhoVol().BOUNDARY_FIELD_REF[patchi],
                dvm_.Uvol().BOUNDARY_FIELD_REF[patchi],
                dvm_.Tvol().BOUNDARY_FIELD_REF[patchi]
            );
        }
    }
    time1 = rpcc();
    // Info << "initBoundaryField time =" <<((double)(time1-time0)*1000/CLOCKRATE) << "ms" << nl << endl;//
    
}

void Foam::discreteVelocity::updateGHbarPvol()
{
    volScalarField gEq
    (
        IOobject
        (
            "gEq",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", gTildeVol_.dimensions(), 0)
    );

    volScalarField hEq
    (
        IOobject
        (
            "hEq",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", hTildeVol_.dimensions(), 0)
    );

    // volScalarField relaxFactor
    // (
    //     IOobject
    //     (
    //         "relaxFactor",
    //         mesh_.time().timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh_,
    //     scalar(0) // dimLess
    // );
    // time1 = rpcc();
    // //- get relaxtion factor 
    // relaxFactor = 1.5*dt/(2.0*dvm_.tauVol() + dt);
    // time2 = rpcc();
    //- get gEq and hEq
    equilibriumShakhovVol
    (
        gEq,
        hEq, 
        dvm_.rhoVol(), 
        dvm_.Uvol(), 
        dvm_.Tvol(), 
        dvm_.qVol() 
    );

    volScalarField relaxFactor
    (
        IOobject
        (
            "relaxFactor",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(0) // dimLess
    );
    // - get delta t
    dimensionedScalar dt = time_.deltaT();
    // - get relaxtion factor 
    relaxFactor = 1.5*dt/(2.0*dvm_.tauVol() + dt);
    gBarPvol_ = (1.0 - relaxFactor)*gTildeVol_ + relaxFactor*gEq;
    hBarPvol_ = (1.0 - relaxFactor)*hTildeVol_ + relaxFactor*hEq;

}

void Foam::discreteVelocity::updateGHbarSurf()
{
    unsigned long time0 = 0;
	unsigned long time1 = 0;
	unsigned long time2 = 0;
	unsigned long time3 = 0;
    unsigned long time4 = 0;
	unsigned long time5 = 0;
    time0 = rpcc();
  // Pre-setup
    // 1. correct the boundary value of gBarPvol_
    //    the gradient at boundary is known
    // 2. get the gradient
    //
    //Info << "Surf first\n" << endl;
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    int64 ownersize = mesh_.owner().size();
    const vectorField& Cf = mesh_.Cf();
    const vectorField& C = mesh_.C();
    const vectorField& Sf = mesh_.Sf();
    const scalarField& V = mesh_.V();

    const Field<scalar>& iGbarPvol = gBarPvol_;
    const Field<scalar>& iHbarPvol = hBarPvol_;
    const Field<vector>& iGbarPgrad = gBarPgrad_;
    const Field<vector>& iHbarPgrad = hBarPgrad_;

    // This is what we want to update in this function
    Field<scalar>& iGsurf = gSurf_;
    Field<scalar>& iHsurf = hSurf_;

    int64 gBarPvolsize = gBarPvol_.size();
    // Info << "gBarPvolsize =" <<(gBarPvolsize) << "ms" << nl << endl;//
    // Info << "ownersize =" <<(ownersize) << "ms" << nl << endl;//

    scalar dt = time_.deltaTValue();
    vector xii = xi_.value();

    scalar *Sf_value = const_cast<scalar*>((mesh_.Sf().begin())->v_);
	scalar *gBarPvol_value = gBarPvol_.begin();
	scalar *hBarPvol_value = hBarPvol_.begin();
    label *owner_value = const_cast<label*>(owner.begin());
    label *neighbour_value = const_cast<label*>(neighbour.begin());

#ifndef GAUSSGRAD
    gBarPvol_.correctBoundaryConditions(); // NOTE: check if the newly defined zeroGradientFvsPatchField 
    hBarPvol_.correctBoundaryConditions();

    gBarPgrad_ = fvc::grad(gBarPvol_); //耗时最长，只有这个函数用gBarPgrad_，并且下次计算的时候清零
    // gBarPgrad_ = calcGrad(gBarPvol_);
    // gBarPgrad_ = calcGrad(gBarPvol_);
    // calcGrad(gBarPvol_);
    hBarPgrad_ = fvc::grad(hBarPvol_);
#else
    gBarPgrad_ = dimensionedVector("0", gBarPgrad_.dimensions(), vector(0, 0, 0));
    hBarPgrad_ = dimensionedVector("0", hBarPgrad_.dimensions(), vector(0, 0, 0));
    const GeometricField<scalar, fvPatchField, volMesh>& gvsf = gBarPvol_;
    const GeometricField<scalar, fvPatchField, volMesh>& hvsf = hBarPvol_;
    Field<vector>& igGrad = gBarPgrad_;
    Field<vector>& ihGrad = hBarPgrad_;
    // const Field<scalar>& gvfi = gvsf.internalField();
    // const Field<scalar>& hvfi = hvsf.internalField();

    scalar *interpola_weight = dvm_.Interpola_weight();
/*
    forAll(owner, facei)
    {

        label own = owner[facei];
        label nei = neighbour[facei];
        // scalar SfdOwn = mag(Sf[facei] & (Cf[facei] - C[own]));
        // scalar SfdNei = mag(Sf[facei] & (C[nei] - Cf[facei]));
        // // sfi[facei] = SfdNei/(SfdOwn + SfdNei)*(vfi[own] - vfi[nei]) + vfi[nei];
        // scalar sfi = SfdNei/(SfdOwn + SfdNei)*(vfi[own] - vfi[nei]) + vfi[nei];
        scalar gsfi = interpola_weight[facei]*(gvfi[own] - gvfi[nei]) + gvfi[nei];
        scalar hsfi = interpola_weight[facei]*(hvfi[own] - hvfi[nei]) + hvfi[nei];
        //  igGrad[own] += Sf[facei]*sfi;
        // igGrad[nei] -= Sf[facei]*sfi;

        igGrad[own] += Sf[facei]*gsfi/V[own];
        igGrad[nei] -= Sf[facei]*gsfi/V[nei];

        ihGrad[own] += Sf[facei]*hsfi/V[own];
        ihGrad[nei] -= Sf[facei]*hsfi/V[nei];

        // igGrad[owner[facei]] += Sf[facei]*issf[facei]/V[owner[facei]];
        // igGrad[neighbour[facei]] -= Sf[facei]*issf[facei]/V[neighbour[facei]];

        // gBarPgrad_[owner[facei]] += Sf[facei]*issf[facei];
        // gBarPgrad_[neighbour[facei]] -= Sf[facei]*issf[facei];
    }
*/    
    /*****************************************************/
    scalar *V_own = dvm_.V_own();
	scalar *V_nei = dvm_.V_nei();
	// scalar *gBarPvol_value_own = (scalar*)malloc(sizeof(scalar) * ownersize);
	// scalar *gBarPvol_value_nei = (scalar*)malloc(sizeof(scalar) * ownersize);
	// scalar *hBarPvol_value_own = (scalar*)malloc(sizeof(scalar) * ownersize);
	// scalar *hBarPvol_value_nei = (scalar*)malloc(sizeof(scalar) * ownersize);
    scalar *gBarPgrad_value_own = (scalar*)malloc(sizeof(scalar) * ownersize*3);
	scalar *gBarPgrad_value_nei = (scalar*)malloc(sizeof(scalar) * ownersize*3);
	scalar *hBarPgrad_value_own = (scalar*)malloc(sizeof(scalar) * ownersize*3);
	scalar *hBarPgrad_value_nei = (scalar*)malloc(sizeof(scalar) * ownersize*3);

    /*for(int facei=0;facei<ownersize;facei++)
	{
		const label own = owner[facei];
		const label nei = neighbour[facei];
		gBarPvol_value_own[facei] = gBarPvol_value[own];
		gBarPvol_value_nei[facei] = gBarPvol_value[nei];
        hBarPvol_value_own[facei] = hBarPvol_value[own];
		hBarPvol_value_nei[facei] = hBarPvol_value[nei];
	}*/
    /*
	//internal faces
    GaussGrad gaussgrad;
	gaussgrad.Sf_value = Sf_value;
	gaussgrad.gBarPvol_value_own = gBarPvol_value_own;
	gaussgrad.gBarPvol_value_nei = gBarPvol_value_nei;
	gaussgrad.hBarPvol_value_own = hBarPvol_value_own;
	gaussgrad.hBarPvol_value_nei = hBarPvol_value_nei;
    gaussgrad.gBarPgrad_value_own = gBarPgrad_value_own;
	gaussgrad.gBarPgrad_value_nei = gBarPgrad_value_nei;
	gaussgrad.hBarPgrad_value_own = hBarPgrad_value_own;
	gaussgrad.hBarPgrad_value_nei = hBarPgrad_value_nei;
	gaussgrad.V_own = V_own;
	gaussgrad.V_nei = V_nei;
	gaussgrad.interpola_weight = interpola_weight;
	gaussgrad.ownersize = ownersize;
	__real_athread_spawn((void *)slave_Func_gaussgrad, &gaussgrad);
    */

    //internal faces
    GaussGrad1 gaussgrad;
	gaussgrad.Sf_value = Sf_value;
	gaussgrad.gBarPvol_value = gBarPvol_value;
	// gaussgrad.gBarPvol_value_nei = gBarPvol_value_nei;
	gaussgrad.hBarPvol_value = hBarPvol_value;
	// gaussgrad.hBarPvol_value_nei = hBarPvol_value_nei;
    gaussgrad.gBarPgrad_value_own = gBarPgrad_value_own;
	gaussgrad.gBarPgrad_value_nei = gBarPgrad_value_nei;
	gaussgrad.hBarPgrad_value_own = hBarPgrad_value_own;
	gaussgrad.hBarPgrad_value_nei = hBarPgrad_value_nei;
	gaussgrad.V_own = V_own;
	gaussgrad.V_nei = V_nei;
    gaussgrad.own = owner_value;
    gaussgrad.neighbour = neighbour_value;
	gaussgrad.interpola_weight = interpola_weight;
	gaussgrad.ownersize = ownersize;
    gaussgrad.gBarPvolsize = gBarPvolsize;
	__real_athread_spawn((void *)slave_Func_gaussgradtemp, &gaussgrad);

    gBarPvol_.correctBoundaryConditions(); // NOTE: check if the newly defined zeroGradientFvsPatchField 
    hBarPvol_.correctBoundaryConditions();
    /********************************************************************************************/
    forAll(mesh_.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh_.boundary()[patchi].faceCells();

        const vectorField& pSf = mesh_.Sf().boundaryField()[patchi];

        // const fvsPatchField<scalar>& pssf = sf.boundaryField()[patchi];
        // sf.boundaryField()[patchi] = vsf.boundaryField()[patchi];
        // const fvsPatchField<scalar>& pssf = gBarPvol_.boundaryField()[patchi];
        forAll(mesh_.boundary()[patchi], facei)
        {
            const label own = pFaceCells[facei];
            // igGrad[pFaceCells[facei]] += pSf[facei]*pssf[facei]/V[own];
            igGrad[pFaceCells[facei]] += pSf[facei]*gvsf.boundaryField()[patchi][facei]/V[own];
            ihGrad[pFaceCells[facei]] += pSf[facei]*hvsf.boundaryField()[patchi][facei]/V[own];
            // igGrad[pFaceCells[facei]] += pSf[facei]*vsf.boundaryField()[patchi][facei];
        }
    }
    athread_join();

    time1 = rpcc();
    for(int facei=0;facei<ownersize;facei++)
	{
		const label own = owner[facei];
		const label nei = neighbour[facei];
		gBarPgrad_[own].v_[0] += gBarPgrad_value_own[facei*3];
		gBarPgrad_[own].v_[1] += gBarPgrad_value_own[facei*3+1];
		gBarPgrad_[own].v_[2] += gBarPgrad_value_own[facei*3+2];
		gBarPgrad_[nei].v_[0] += gBarPgrad_value_nei[facei*3];
		gBarPgrad_[nei].v_[1] += gBarPgrad_value_nei[facei*3+1];
		gBarPgrad_[nei].v_[2] += gBarPgrad_value_nei[facei*3+2];
	}
    for(int facei=0;facei<ownersize;facei++)
	{
		const label own = owner[facei];
		const label nei = neighbour[facei];
		hBarPgrad_[own].v_[0] += hBarPgrad_value_own[facei*3];
		hBarPgrad_[own].v_[1] += hBarPgrad_value_own[facei*3+1];
		hBarPgrad_[own].v_[2] += hBarPgrad_value_own[facei*3+2];
		hBarPgrad_[nei].v_[0] += hBarPgrad_value_nei[facei*3];
		hBarPgrad_[nei].v_[1] += hBarPgrad_value_nei[facei*3+1];
		hBarPgrad_[nei].v_[2] += hBarPgrad_value_nei[facei*3+2];
	}
    time2 = rpcc();
    gBarPgrad_.correctBoundaryConditions();
    hBarPgrad_.correctBoundaryConditions();
    Foam::fv::gaussGrad<scalar>::correctBoundaryConditions(gvsf, gBarPgrad_);
    Foam::fv::gaussGrad<scalar>::correctBoundaryConditions(hvsf, hBarPgrad_);
    // free(gBarPvol_value_own); 
    // free(gBarPvol_value_nei);
    // free(hBarPvol_value_own);
    // free(hBarPvol_value_nei);
    free(gBarPgrad_value_own);
    free(gBarPgrad_value_nei);
    free(hBarPgrad_value_own);
    free(hBarPgrad_value_nei);
#endif
    /****************************************************/

    //Info << "Surf finished\n" << endl;
    // The DVMsymmetry is  processed automatically in fvc::grad operator

    //
    // 3. correct the boundary value of the grad field
    //    to be used at next time
    //gBarPgrad_.correctBoundaryConditions();
    //hBarPgrad_.correctBoundaryConditions();

    // 4. patch the normal component of boundary value of the grad 
    //    to the gradient field of the fixed 
    //    gradient feild of the gBarPvol_ and 
    //    hBarPvol_ ...
    //    NOTE: we need the surfaceNormal Gradient

    forAll(gBarPgrad_.BOUNDARY_FIELD_REF, patchi)
    {
        const vectorField n
        (
            mesh_.Sf().boundaryField()[patchi]/mesh_.magSf().boundaryField()[patchi]
        );
        
        // if(gBarPvol_.BOUNDARY_FIELD_REF[patchi].type() == "processor") {
        //     Info<< gBarPvol_.BOUNDARY_FIELD_REF[patchi].patchInternalField()<<endl;

        //     Info<< gBarPvol_.BOUNDARY_FIELD_REF[patchi].patchInternalField().size()<<endl;
        // }
        
        if ( 
                   gBarPvol_.BOUNDARY_FIELD_REF[patchi].type() != "empty" 
            && gBarPvol_.BOUNDARY_FIELD_REF[patchi].type() != "processor"
            && gBarPvol_.BOUNDARY_FIELD_REF[patchi].type() != "symmetryPlane"
            && gBarPvol_.BOUNDARY_FIELD_REF[patchi].type() != "DVMsymmetry"
            && gBarPvol_.BOUNDARY_FIELD_REF[patchi].type() != "cyclic"
            && gBarPvol_.BOUNDARY_FIELD_REF[patchi].type() != "processorCyclic"
           ) // only for fixed gradient g/hBarPvol
        {
            // normal component of the grad field
            fixedGradientFvPatchField<scalar>& gBarPvolPatch = 
                refCast<fixedGradientFvPatchField<scalar> >
                (gBarPvol_.BOUNDARY_FIELD_REF[patchi]);

            fixedGradientFvPatchField<scalar>& hBarPvolPatch = 
                refCast<fixedGradientFvPatchField<scalar> >
                (hBarPvol_.BOUNDARY_FIELD_REF[patchi]);

            forAll(gBarPvolPatch, pFacei)
            {
                gBarPvolPatch.gradient()[pFacei] =
                    gBarPgrad_.BOUNDARY_FIELD_REF[patchi][pFacei]&n[pFacei];
                hBarPvolPatch.gradient()[pFacei] =
                    hBarPgrad_.BOUNDARY_FIELD_REF[patchi][pFacei]&n[pFacei];
            }
        }
    }
    /****************************************************/

    // const labelUList& owner = mesh_.owner();
    // const labelUList& neighbour = mesh_.neighbour();

    // internal faces first
    // int64 ownersize = mesh_.owner().size();
#ifndef GHBARSURF
    // forAll(owner, facei)
    for(int facei=0;facei<ownersize;facei++)
    {
        label own = owner[facei];
        label nei = neighbour[facei];
        vector temp = Cf[facei] - 0.5*xii*dt;
        if ((xii&Sf[facei]) >=  VSMALL) // comming from own
        {   //printf("1");
            temp -=  C[own];
            iGsurf[facei] = iGbarPvol[own] + (iGbarPgrad[own]&temp);
            iHsurf[facei] = iHbarPvol[own] + (iHbarPgrad[own]&temp);
            // iGsurf[facei] = iGbarPvol[own] + (iGbarPgrad[own]&(Cf[facei] - C[own] - 0.5*xii*dt));
            // iHsurf[facei] = iHbarPvol[own] + (iHbarPgrad[own]&(Cf[facei] - C[own] - 0.5*xii*dt));
        }
        // Debug, no = 0, =0 put to > 0
        else if ((xii&Sf[facei]) < -VSMALL) // comming form nei
        {
            temp -=  C[nei];
            iGsurf[facei] = iGbarPvol[nei] + (iGbarPgrad[nei]&temp);
            iHsurf[facei] = iHbarPvol[nei] + (iHbarPgrad[nei]&temp);
            // iGsurf[facei] = iGbarPvol[nei] + (iGbarPgrad[nei]&(Cf[facei] - C[nei] - 0.5*xii*dt));
            // iHsurf[facei] = iHbarPvol[nei] + (iHbarPgrad[nei]&(Cf[facei] - C[nei] - 0.5*xii*dt));
        }
        else 
        {
            // printf("3");
            iGsurf[facei] = 0.5*(iGbarPvol[nei] + ((iGbarPgrad[nei]) &(Cf[facei] - C[nei] - 0.5*xii*dt)) + iGbarPvol[own] + ((iGbarPgrad[own]) &(Cf[facei] - C[own] - 0.5*xii*dt)));
            iHsurf[facei] = 0.5*(iHbarPvol[nei] + ((iHbarPgrad[nei]) &(Cf[facei] - C[nei] - 0.5*xii*dt)) + iHbarPvol[own] + ((iHbarPgrad[own]) &(Cf[facei] - C[own] - 0.5*xii*dt)));
        }
    }
      // boundary faces
    forAll(gSurf_.BOUNDARY_FIELD_REF, patchi)
    {
        word type = gSurf_.BOUNDARY_FIELD_REF[patchi].type();
        fvsPatchField<scalar>& gSurfPatch = gSurf_.BOUNDARY_FIELD_REF[patchi];
        fvsPatchField<scalar>& hSurfPatch = hSurf_.BOUNDARY_FIELD_REF[patchi];
        const fvsPatchField<vector>& SfPatch =
            mesh_.Sf().boundaryField()[patchi];
        const fvsPatchField<vector>& CfPatch =
            mesh_.Cf().boundaryField()[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();

        const fvPatchScalarField& rhoVolPatch = 
            dvm_.rhoVol().boundaryField()[patchi];
        const fvPatchScalarField& TvolPatch = 
            dvm_.Tvol().boundaryField()[patchi];
        const labelUList& pOwner = mesh_.boundary()[patchi].faceCells();
        
        //- NOTE: outging DF can be treate unifily for all BCs, including processor BC
        if (type == "zeroGradient")
        {
            gSurfPatch == gBarPvol_.BOUNDARY_FIELD_REF[patchi].patchInternalField();
            hSurfPatch == hBarPvol_.BOUNDARY_FIELD_REF[patchi].patchInternalField();
        }
        else if (type == "mixed")
        {
            //check each boundary face in the patch
            forAll(gSurfPatch, facei)
            {
                //out or in ?
                if ((xii&SfPatch[facei]) > 0 ) // outgoing
                {
                    gSurfPatch[facei] = iGbarPvol[faceCells[facei]] 
                      + ((iGbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));
                    hSurfPatch[facei] = iHbarPvol[faceCells[facei]] 
                      + ((iHbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));
                //incoming and parallel to face, not changed.
                }
            }
        }
        else if (type == "farField")
        {
            //check each boundary face in the patch
            forAll(gSurfPatch, facei)
            {
                //out or in ?
                if ((xii&SfPatch[facei]) > 0 ) // outgoing
                {
                    gSurfPatch[facei] = iGbarPvol[faceCells[facei]] 
                      + ((iGbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));
                    hSurfPatch[facei] = iHbarPvol[faceCells[facei]] 
                      + ((iHbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));
                //incoming and parallel to face, not changed.
                }
                else // incomming, set to be equlibrium, give rho and T, extropolate U
                {
                    // set to maxwellian
                    gSurfPatch[facei] = rhoVolPatch[facei]
                       *equilibriumMaxwellByRho
                        (
                            dvm_.Uvol()[pOwner[facei]],
                            TvolPatch[facei]
                        );
                    hSurfPatch[facei] = 
                        gSurfPatch[facei]*(dvm_.R().value()*TvolPatch[facei])
                       *(dvm_.KInner() + 3 - mesh_.nSolutionD());
                }
            }
        }
        else if (type == "maxwellWall")
        {
            calculatedMaxwellFvPatchField<scalar>& rhoPatch = 
                refCast<calculatedMaxwellFvPatchField<scalar> >
                (dvm_.rhoVol().BOUNDARY_FIELD_REF[patchi]); //DEBUG

            forAll(gSurfPatch, facei)
            {
                vector faceSf= SfPatch[facei];
                if ((xii&faceSf) >  0 ) // outgoing
                {
                    gSurfPatch[facei] = iGbarPvol[faceCells[facei]] 
                      + ((iGbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));
                    hSurfPatch[facei] = iHbarPvol[faceCells[facei]] 
                      + ((iHbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));

                    rhoPatch.outGoing()[facei] += //add outgoing normal momentum flux to outGoing container
                        weight_*(xii&faceSf)*gSurfPatch[facei];
                }
            }
        }
        else if (type == "processor"
              || type == "cyclic"
              || type == "processorCyclic"
              ) // parallel
        {
            forAll(gSurfPatch, facei)
            {
                vector faceSf= SfPatch[facei];
                if ((xii&faceSf) >  VSMALL ) // outgoing
                {
                    gSurfPatch[facei] = iGbarPvol[faceCells[facei]] 
                      + ((iGbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));
                    hSurfPatch[facei] = iHbarPvol[faceCells[facei]] 
                      + ((iHbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));
                } 
                else if ((xii&faceSf) <  -VSMALL )//incomming from processor boundaryField
                {
                    gSurfPatch[facei] = gBarPvol_.BOUNDARY_FIELD_REF[patchi][facei]
                      + ((gBarPgrad_.BOUNDARY_FIELD_REF[patchi][facei])
                       &(CfPatch[facei] - mesh_.C().boundaryField()[patchi][facei] - 0.5*xii*dt));
                    hSurfPatch[facei] = hBarPvol_.BOUNDARY_FIELD_REF[patchi][facei]
                      + ((hBarPgrad_.BOUNDARY_FIELD_REF[patchi][facei])
                       &(CfPatch[facei] - mesh_.C().boundaryField()[patchi][facei] - 0.5*xii*dt));
                }
                else 
                {
                    gSurfPatch[facei] = 0.5*(
                            iGbarPvol[faceCells[facei]] + 
                            (   (iGbarPgrad[faceCells[facei]]) & (CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt) )
                            + gBarPvol_.BOUNDARY_FIELD_REF[patchi][facei] + 
                            (   (gBarPgrad_.BOUNDARY_FIELD_REF[patchi][facei]) & (CfPatch[facei] - mesh_.C().boundaryField()[patchi][facei] - 0.5*xii*dt) )
                    );
                    hSurfPatch[facei] = 0.5*(
                            iHbarPvol[faceCells[facei]] + 
                            (   (iHbarPgrad[faceCells[facei]]) & (CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt) )
                            + hBarPvol_.BOUNDARY_FIELD_REF[patchi][facei] + 
                            (   (hBarPgrad_.BOUNDARY_FIELD_REF[patchi][facei]) & (CfPatch[facei] - mesh_.C().boundaryField()[patchi][facei] - 0.5*xii*dt) )
                    );

                }

            }
        }
        /*else if (type == "symmetryPlane" || type == "DVMsymmetry")
        {
            forAll(gSurfPatch, facei)
            {
                vector faceSf= SfPatch[facei];
                if ((xii&faceSf) >  -VSMALL ) // outgoing and **reside(shouldn't be ignored)** DF, 
                                              //incomming shoud be proceed after this function
                {
                    gSurfPatch[facei] = iGbarPvol[faceCells[facei]] 
                      + ((iGbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));
                    hSurfPatch[facei] = iHbarPvol[faceCells[facei]] 
                      + ((iHbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));
                } 
            }
        }*/
    }

#else
    
    /*****************************************************/
    scalar *Cf_value = const_cast<double*>((mesh_.Cf().begin())->v_);
    scalar *C_own = dvm_.C_own();
	scalar *C_nei = dvm_.C_nei();
    scalar *gBarPgrad_value = const_cast<double*>((gBarPgrad_.begin())->v_);
	scalar *hBarPgrad_value = const_cast<double*>((hBarPgrad_.begin())->v_);
    scalar *gSurf_value = gSurf_.begin();
	scalar *hSurf_value = hSurf_.begin();
    time3 = rpcc();
    time4 = rpcc();
	//internal faces
    GHbarSurf2 ghbarsurf;
	ghbarsurf.Sf = Sf_value;
    ghbarsurf.Cf = Cf_value;
	// ghbarsurf.gBarPvol_value_own = gBarPvol_value_own;
	// ghbarsurf.gBarPvol_value_nei = gBarPvol_value_nei;
	// ghbarsurf.hBarPvol_value_own = hBarPvol_value_own;
	// ghbarsurf.hBarPvol_value_nei = hBarPvol_value_nei;
    ghbarsurf.gBarPvol_value = gBarPvol_value;
	ghbarsurf.hBarPvol_value = hBarPvol_value;
    // ghbarsurf.gBarPgrad_value_own = gBarPgrad_value_own;
	// ghbarsurf.gBarPgrad_value_nei = gBarPgrad_value_nei;
	// ghbarsurf.hBarPgrad_value_own = hBarPgrad_value_own;
	// ghbarsurf.hBarPgrad_value_nei = hBarPgrad_value_nei;
    ghbarsurf.gBarPgrad_value = gBarPgrad_value;
	ghbarsurf.hBarPgrad_value = hBarPgrad_value;
    ghbarsurf.own = owner_value;
    ghbarsurf.neighbour = neighbour_value;
	ghbarsurf.C_own = C_own;
	ghbarsurf.C_nei = C_nei;
	ghbarsurf.gSurf_value = gSurf_value;
	ghbarsurf.hSurf_value = hSurf_value;
	ghbarsurf.dt = dt;
	ghbarsurf.xii_x = xii.x();
	ghbarsurf.xii_y = xii.y();
	ghbarsurf.xii_z = xii.z();
	ghbarsurf.ownersize = ownersize;
    ghbarsurf.gBarPvolsize = gBarPvolsize;
	__real_athread_spawn((void *)slave_Func_ghbarsurftemp2, &ghbarsurf);
    
    // boundary faces
    forAll(gSurf_.BOUNDARY_FIELD_REF, patchi)
    {
        word type = gSurf_.BOUNDARY_FIELD_REF[patchi].type();
        fvsPatchField<scalar>& gSurfPatch = gSurf_.BOUNDARY_FIELD_REF[patchi];
        fvsPatchField<scalar>& hSurfPatch = hSurf_.BOUNDARY_FIELD_REF[patchi];
        const fvsPatchField<vector>& SfPatch = mesh_.Sf().boundaryField()[patchi];
        const fvsPatchField<vector>& CfPatch = mesh_.Cf().boundaryField()[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();

        const fvPatchScalarField& rhoVolPatch =  dvm_.rhoVol().boundaryField()[patchi];
        const fvPatchScalarField& TvolPatch = dvm_.Tvol().boundaryField()[patchi];
        const labelUList& pOwner = mesh_.boundary()[patchi].faceCells();
        
        //- NOTE: outging DF can be treate unifily for all BCs, including processor BC
        if (type == "zeroGradient")
        {
            gSurfPatch == gBarPvol_.BOUNDARY_FIELD_REF[patchi].patchInternalField();
            hSurfPatch == hBarPvol_.BOUNDARY_FIELD_REF[patchi].patchInternalField();
        }
        else if (type == "mixed")
        {
            //check each boundary face in the patch
            forAll(gSurfPatch, facei)
            {
                //out or in ?
                if ((xii&SfPatch[facei]) > 0 ) // outgoing
                {
                    gSurfPatch[facei] = iGbarPvol[faceCells[facei]] 
                      + ((iGbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));
                    hSurfPatch[facei] = iHbarPvol[faceCells[facei]] 
                      + ((iHbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));
                //incoming and parallel to face, not changed.
                }
            }
        }
        else if (type == "farField")
        {
            //check each boundary face in the patch
            forAll(gSurfPatch, facei)
            {
                //out or in ?
                if ((xii&SfPatch[facei]) > 0 ) // outgoing
                {
                    gSurfPatch[facei] = iGbarPvol[faceCells[facei]] 
                      + ((iGbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));
                    hSurfPatch[facei] = iHbarPvol[faceCells[facei]] 
                      + ((iHbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));
                //incoming and parallel to face, not changed.
                }
                else // incomming, set to be equlibrium, give rho and T, extropolate U
                {
                    // set to maxwellian
                    gSurfPatch[facei] = rhoVolPatch[facei]
                       *equilibriumMaxwellByRho
                        (
                            dvm_.Uvol()[pOwner[facei]],
                            TvolPatch[facei]
                        );
                    hSurfPatch[facei] = 
                        gSurfPatch[facei]*(dvm_.R().value()*TvolPatch[facei])
                       *(dvm_.KInner() + 3 - mesh_.nSolutionD());
                }
            }
        }
        else if (type == "maxwellWall")
        {
            calculatedMaxwellFvPatchField<scalar>& rhoPatch = 
                refCast<calculatedMaxwellFvPatchField<scalar> >
                (dvm_.rhoVol().BOUNDARY_FIELD_REF[patchi]); //DEBUG

            forAll(gSurfPatch, facei)
            {
                vector faceSf= SfPatch[facei];
                if ((xii&faceSf) >  0 ) // outgoing
                {
                    gSurfPatch[facei] = iGbarPvol[faceCells[facei]] 
                      + ((iGbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));
                    hSurfPatch[facei] = iHbarPvol[faceCells[facei]] 
                      + ((iHbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));

                    rhoPatch.outGoing()[facei] += //add outgoing normal momentum flux to outGoing container
                        weight_*(xii&faceSf)*gSurfPatch[facei];
                }
            }
        }
        else if (type == "processor"
              || type == "cyclic"
              || type == "processorCyclic"
              ) // parallel
        {
            forAll(gSurfPatch, facei)
            {
                vector faceSf= SfPatch[facei];
                if ((xii&faceSf) >  VSMALL ) // outgoing
                {
                    gSurfPatch[facei] = iGbarPvol[faceCells[facei]] 
                      + ((iGbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));
                    hSurfPatch[facei] = iHbarPvol[faceCells[facei]] 
                      + ((iHbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));
                } 
                else if ((xii&faceSf) <  -VSMALL )//incomming from processor boundaryField
                {
                    gSurfPatch[facei] = gBarPvol_.BOUNDARY_FIELD_REF[patchi][facei]
                      + ((gBarPgrad_.BOUNDARY_FIELD_REF[patchi][facei])
                       &(CfPatch[facei] - mesh_.C().boundaryField()[patchi][facei] - 0.5*xii*dt));
                    hSurfPatch[facei] = hBarPvol_.BOUNDARY_FIELD_REF[patchi][facei]
                      + ((hBarPgrad_.BOUNDARY_FIELD_REF[patchi][facei])
                       &(CfPatch[facei] - mesh_.C().boundaryField()[patchi][facei] - 0.5*xii*dt));
                }
                else 
                {
                    gSurfPatch[facei] = 0.5*(
                            iGbarPvol[faceCells[facei]] + 
                            (   (iGbarPgrad[faceCells[facei]]) & (CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt) )
                            + gBarPvol_.BOUNDARY_FIELD_REF[patchi][facei] + 
                            (   (gBarPgrad_.BOUNDARY_FIELD_REF[patchi][facei]) & (CfPatch[facei] - mesh_.C().boundaryField()[patchi][facei] - 0.5*xii*dt) )
                    );
                    hSurfPatch[facei] = 0.5*(
                            iHbarPvol[faceCells[facei]] + 
                            (   (iHbarPgrad[faceCells[facei]]) & (CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt) )
                            + hBarPvol_.BOUNDARY_FIELD_REF[patchi][facei] + 
                            (   (hBarPgrad_.BOUNDARY_FIELD_REF[patchi][facei]) & (CfPatch[facei] - mesh_.C().boundaryField()[patchi][facei] - 0.5*xii*dt) )
                    );

                }

            }
        }
        /*else if (type == "symmetryPlane" || type == "DVMsymmetry")
        {
            forAll(gSurfPatch, facei)
            {
                vector faceSf= SfPatch[facei];
                if ((xii&faceSf) >  -VSMALL ) // outgoing and **reside(shouldn't be ignored)** DF, 
                                              //incomming shoud be proceed after this function
                {
                    gSurfPatch[facei] = iGbarPvol[faceCells[facei]] 
                      + ((iGbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));
                    hSurfPatch[facei] = iHbarPvol[faceCells[facei]] 
                      + ((iHbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));
                } 
            }
        }*/
    }
    athread_join();
    // free(gBarPgrad_value_own);
    // free(gBarPgrad_value_nei);
    // free(hBarPgrad_value_own);
    // free(hBarPgrad_value_nei);
#endif
    time5= rpcc();
}

void Foam::discreteVelocity::updateGHbarSurfMaxwellWallIn()
{
    vector xii = xi_.value();
    forAll(gSurf_.BOUNDARY_FIELD_REF, patchi)
    {
        if (gSurf_.BOUNDARY_FIELD_REF[patchi].type() == "maxwellWall")
        {
            fvsPatchScalarField& gSurfPatch = gSurf_.BOUNDARY_FIELD_REF[patchi];
            fvsPatchScalarField& hSurfPatch = hSurf_.BOUNDARY_FIELD_REF[patchi];
            const fvPatchScalarField& rhoVolPatch = 
                dvm_.rhoVol().BOUNDARY_FIELD_REF[patchi];
            const fvPatchVectorField& UvolPatch = 
                dvm_.Uvol().BOUNDARY_FIELD_REF[patchi];
            const fvPatchScalarField& TvolPatch = 
                dvm_.Tvol().BOUNDARY_FIELD_REF[patchi];
            const fvsPatchVectorField& SfPatch = 
                mesh_.Sf().boundaryField()[patchi];
            forAll(gSurfPatch, facei)
            {
                vector faceSf = SfPatch[facei];
                if ((xii & faceSf) <= 0) // incomming
                {
                    // set to maxwellian
                    gSurfPatch[facei] = rhoVolPatch[facei]
                       *equilibriumMaxwellByRho
                        (
                            UvolPatch[facei],
                            TvolPatch[facei]
                        );

                    //set hSurf at maxwellWall to zero! , WRONG!!!
                    hSurfPatch[facei] = 
                        gSurfPatch[facei]*(dvm_.R().value()*TvolPatch[facei])
                       *(dvm_.KInner() + 3 - mesh_.nSolutionD());
                }
            }
        }
    }
}

void Foam::discreteVelocity::updateGHbarSurfSymmetryIn()
{
    vector xii = xi_.value();

    label nproc = dvm_.mpiReducer().nproc();
    label rank  = dvm_.mpiReducer().rank();
    labelField recvc(nproc);
    labelField displ(nproc);
    label chunck = dvm_.nXi()/nproc;
    label left   = dvm_.nXi()%nproc;

    forAll(gSurf_.BOUNDARY_FIELD_REF, patchi)
    {
        if (gSurf_.BOUNDARY_FIELD_REF[patchi].type() == "DVMsymmetry" 
        &&  gSurf_.BOUNDARY_FIELD_REF[patchi].size() > 0 )
        {
            fvsPatchScalarField& gSurfPatch = gSurf_.BOUNDARY_FIELD_REF[patchi];
            fvsPatchScalarField& hSurfPatch = hSurf_.BOUNDARY_FIELD_REF[patchi];
            
#if FOAM_MAJOR <= 3
            GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
                rhoBCs = dvm_.rhoVol().boundaryField();
#else
            GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
                rhoBCs = dvm_.rhoVol().BOUNDARY_FIELD_REF;
#endif

            symmetryModFvPatchField<scalar>& rhoPatch = 
                refCast<symmetryModFvPatchField<scalar> >(rhoBCs[patchi]);

            const vector faceSf = mesh_.Sf().boundaryField()[patchi][0];
            label ps = gSurfPatch.size();

            forAll(recvc, i)
            {
                recvc[i] = 2*ps*(chunck + (i<left)) ;
                if(i<=left)
                    displ[i] = i*2*ps*(chunck + 1); // (i<=nXi_%nproc)
                else
                    displ[i] = 2*ps*(left*(chunck + 1) + (i-left)*(chunck));
            }

            if ((xii & faceSf) <= 0) // incomming
            {
                vector nomlizedDirec = faceSf/mag(faceSf);
                label targetDVid = round( fabs(nomlizedDirec
                    & vector(symXtargetDVid_, symYtargetDVid_, symZtargetDVid_)));

                label p  = targetDVid%nproc; //locate in p's
                label pi = targetDVid/nproc;
                
                label shift = displ[p] + pi*ps*2;
                memcpy( gSurfPatch.data(),
                        rhoPatch.dfContainer().data() + shift,      ps*sizeof(scalar) );
                memcpy( hSurfPatch.data(),
                        rhoPatch.dfContainer().data() + shift+ps,   ps*sizeof(scalar) );
            }
        }
        else if (gSurf_.BOUNDARY_FIELD_REF[patchi].type() == "symmetryPlane" 
        &&  gSurf_.BOUNDARY_FIELD_REF[patchi].size() > 0 )
        {
            fvsPatchScalarField& gSurfPatch = gSurf_.BOUNDARY_FIELD_REF[patchi];
            fvsPatchScalarField& hSurfPatch = hSurf_.BOUNDARY_FIELD_REF[patchi];

#if FOAM_MAJOR <= 3
            const vector faceSf = mesh_.Sf().boundaryField()[patchi][0];
#else
            const vector faceSf = mesh_.Sf().boundaryField()[patchi][0];
#endif
            if ((xii & faceSf) <= 0) // incomming
            {
                vector nomlizedDirec = faceSf/mag(faceSf);
                label targetDVid = abs(round(nomlizedDirec
                    & vector(symXtargetDVid_, symYtargetDVid_, symZtargetDVid_)));
                forAll(gSurfPatch, facei)
                {
                    gSurfPatch[facei] = dvm_.DVi(targetDVid).gSurf().boundaryField()
                        [patchi][facei];
                    hSurfPatch[facei] = dvm_.DVi(targetDVid).hSurf().boundaryField()
                        [patchi][facei];
                }
            }
        } 
    } 
}

void Foam::discreteVelocity::updateGHsurf()
{
    //- get delta t
    dimensionedScalar h = 0.5*time_.deltaT();
    surfaceScalarField gEq
    (
        IOobject
        (
            "gEq",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( "0", gTildeVol_.dimensions(), 0)
    );
    surfaceScalarField hEq
    (
        IOobject
        (
            "hEq",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", hTildeVol_.dimensions(), 0)
    );
    surfaceScalarField relaxFactor
    (
        IOobject
        (
            "relaxFactor",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(0)
    );
    //- get relaxtion factor 
    relaxFactor = h/(2*dvm_.tauSurf() + h);
    //- get gEq and hEq
    equilibriumShakhovVol
    (
        gEq,
        hEq, 
        dvm_.rhoSurf(), 
        dvm_.Usurf(), 
        dvm_.Tsurf(), 
        dvm_.qSurf() 
    );
    gSurf_ = (1.0 - relaxFactor)*gSurf_ + relaxFactor*gEq;
    hSurf_ = (1.0 - relaxFactor)*hSurf_ + relaxFactor*hEq;

    // forAll(gSurf_.BOUNDARY_FIELD_REF, patchi)
    // {
    //     fvsPatchScalarField& gEqPatch   =    gSurf_.BOUNDARY_FIELD_REF[patchi];
    //     int gEqPatchsize = gEqPatch.size();
    //     for(int i=0;i<gEqPatchsize;i++){
    //         Info << gEqPatch[i] << "  ";
    //     }
    // }

    // NOTE: here the boundar face value are not computed
    // vector xii = xi_.value();
    // We do it mannuly for the outgoing DF!
    // forAll(gSurf_.BOUNDARY_FIELD_REF, patchi)
    // {
    //     fvsPatchScalarField& gSurfPatch = gSurf_.BOUNDARY_FIELD_REF[patchi];
    //     fvsPatchScalarField& hSurfPatch = hSurf_.BOUNDARY_FIELD_REF[patchi];
    //     fvsPatchScalarField& gEqPatch   =    gEq.BOUNDARY_FIELD_REF[patchi];
    //     fvsPatchScalarField& hEqPatch   =    hEq.BOUNDARY_FIELD_REF[patchi];
    //     // printf("here");
    //     const fvsPatchVectorField& SfPatch = mesh_.Sf().boundaryField()[patchi];

    //     fvsPatchScalarField& relaxFactorPatch = 
    //         relaxFactor.BOUNDARY_FIELD_REF[patchi];

    //     //NOTE:  If keep the below block, for the Sod problem, -parallel running would be different from serial run.
    //     //       I also doubt that for cyclic boundary, the same problem will happen.
    //     //       But if I totally deleted it, my experience shows for hypersonic cylinder flow, the free stream BC will have problem.
    //     //       So, I keep it, but exclude the cases of processor/cyclic/processorCyclic.
    //     //       The reason I guess is that, for those 'coupled' type boundary condition, the boundary field will be calculated when 
    //     //       doing GeometricField calculation, such as the relaxiation computation. While for those boundary condition I've defined, 
    //     //       such as mixed or maxwell, I have to explicitly calculate the unchanged boundary values.
      
         
    //     // forAll(gSurfPatch, facei)
    //     // {
    //     //     if ((xii&(SfPatch[facei])) > 0    // Here, if delted , the free stream BC may bo s problem
    //     //     && gBarPvol_.BOUNDARY_FIELD_REF[patchi].type() != "processor"
    //     //     && gBarPvol_.BOUNDARY_FIELD_REF[patchi].type() != "processorCyclic"
    //     //     && gBarPvol_.BOUNDARY_FIELD_REF[patchi].type() != "cyclic")
    //     //     //&& gSurf_.BOUNDARY_FIELD_REF[patchi].type() != "DVMsymmetry")
    //     //     {
    //     //         gSurfPatch[facei] = (1.0 - relaxFactorPatch[facei])
    //     //             *gSurfPatch[facei]
    //     //            + relaxFactorPatch[facei]*gEqPatch[facei];
    //     //         hSurfPatch[facei] = (1.0 - relaxFactorPatch[facei])
    //     //             *hSurfPatch[facei]
    //     //            + relaxFactorPatch[facei]*hEqPatch[facei];
    //     //     }
    //     // }
        

    //     //both in and out DF has been relaxed at DVMsymmetry boundary
    //     if(gSurfPatch.type() == "DVMsymmetry")
    //     {
    //         gSurfPatch = (1.0 - relaxFactorPatch)
    //             *gSurfPatch
    //            + relaxFactorPatch*gEqPatch;
    //         hSurfPatch = (1.0 - relaxFactorPatch)
    //             *hSurfPatch
    //            + relaxFactorPatch*hEqPatch;
    //     }
    // }
    
}

void Foam::discreteVelocity::updateGHtildeVol()
{

    const scalar dt = time_.deltaTValue();
    const vector xii = xi_.value();

	unsigned long time_1st = 0;
    unsigned long time_1en = 0;
    
#ifndef GHTILDEVOL
    // store the gTildePlus in gTilde
    gTildeVol_ = -1.0/3*gTildeVol_ + 4.0/3*gBarPvol_;
    hTildeVol_ = -1.0/3*hTildeVol_ + 4.0/3*hBarPvol_;

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const vectorField Sf = mesh_.Sf();
    const scalarField V = mesh_.V();
    
    // internal faces
    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];
        gTildeVol_[own] -= ((xii&Sf[facei])*gSurf_[facei]*dt/V[own]);
        gTildeVol_[nei] += ((xii&Sf[facei])*gSurf_[facei]*dt/V[nei]);
        hTildeVol_[own] -= ((xii&Sf[facei])*hSurf_[facei]*dt/V[own]);
        hTildeVol_[nei] += ((xii&Sf[facei])*hSurf_[facei]*dt/V[nei]);
    }
    forAll(gSurf_.BOUNDARY_FIELD_REF, patchi)
    {
        const fvsPatchField<scalar>& gSurfPatch =
            gSurf_.BOUNDARY_FIELD_REF[patchi];
        const fvsPatchField<scalar>& hSurfPatch =
            hSurf_.BOUNDARY_FIELD_REF[patchi];
        const fvsPatchField<vector>& SfPatch =
            mesh_.Sf().boundaryField()[patchi];
        const labelUList& pOwner = mesh_.boundary()[patchi].faceCells();
        forAll(pOwner, pFacei)
        {
            const label own = pOwner[pFacei];

            gTildeVol_[own] -= (xii&SfPatch[pFacei]) 
               *gSurfPatch[pFacei]*dt/V[own];
            hTildeVol_[own] -= (xii&SfPatch[pFacei]) 
               *hSurfPatch[pFacei]*dt/V[own];
        }
    }
#else
        const int64 *owner = mesh_.owner().begin();
		int64 ownersize = mesh_.owner().size();
		const int64 *neighbour = mesh_.neighbour().begin();
        scalar *Sf = const_cast<double*>((mesh_.Sf().begin())->v_);

		double *gTildeVol_value = gTildeVol_.begin();
		double *hTildeVol_value = hTildeVol_.begin();
		
        const double *V = mesh_.V().begin();
        
        scalar *V_own = dvm_.V_own();
		scalar *V_nei = dvm_.V_nei();

		scalar *gTildeVol_value_own = (scalar*)malloc(sizeof(scalar) * ownersize);
		scalar *gTildeVol_value_nei = (scalar*)malloc(sizeof(scalar) * ownersize);
		scalar *hTildeVol_value_own = (scalar*)malloc(sizeof(scalar) * ownersize);
		scalar *hTildeVol_value_nei = (scalar*)malloc(sizeof(scalar) * ownersize);
		// memset(gTildeVol_value_own,0,sizeof(scalar)*ownersize);
		// memset(gTildeVol_value_nei,0,sizeof(scalar)*ownersize);
		// memset(hTildeVol_value_own,0,sizeof(scalar)*ownersize);
		// memset(hTildeVol_value_nei,0,sizeof(scalar)*ownersize);
        scalar *gSurf_value = gSurf_.begin();
		scalar *hSurf_value = hSurf_.begin();
		
		//internal faces
    	GHtildeVol ghtildevol;
		ghtildevol.Sf = Sf;
		ghtildevol.gTildeVol_value_own = gTildeVol_value_own;
		ghtildevol.gTildeVol_value_nei = gTildeVol_value_nei;
		ghtildevol.hTildeVol_value_own = hTildeVol_value_own;
		ghtildevol.hTildeVol_value_nei = hTildeVol_value_nei;
		ghtildevol.V_own = V_own;
		ghtildevol.V_nei = V_nei;
		ghtildevol.gSurf_value = gSurf_value;
		ghtildevol.hSurf_value = hSurf_value;
		ghtildevol.dt = dt;
		ghtildevol.xii_x = xii.x();
		ghtildevol.xii_y = xii.y();
		ghtildevol.xii_z = xii.z();
		ghtildevol.ownersize = ownersize;
		__real_athread_spawn((void *)slave_Func_ghtildevol, &ghtildevol);

        forAll(gSurf_.BOUNDARY_FIELD_REF, patchi)
        {
            const fvsPatchField<scalar>& gSurfPatch = gSurf_.BOUNDARY_FIELD_REF[patchi];
            const fvsPatchField<scalar>& hSurfPatch = hSurf_.BOUNDARY_FIELD_REF[patchi];
            const fvsPatchField<vector>& SfPatch = mesh_.Sf().boundaryField()[patchi];
            const labelUList& pOwner = mesh_.boundary()[patchi].faceCells();
            forAll(pOwner, pFacei)
            {
                const label own = pOwner[pFacei];
                
                gTildeVol_[own] -= (xii&SfPatch[pFacei]) *gSurfPatch[pFacei]*dt/V[own];
                hTildeVol_[own] -= (xii&SfPatch[pFacei]) *hSurfPatch[pFacei]*dt/V[own];
            }
        }
		athread_join();
        int facei = 0;
		for(facei=0;facei<ownersize;facei++)
		{
			const label own = owner[facei];
			const label nei = neighbour[facei];
			gTildeVol_value[own] += gTildeVol_value_own[facei];
			gTildeVol_value[nei] += gTildeVol_value_nei[facei];

		}

        for(facei=0;facei<ownersize;facei++)
		{
			const label own = owner[facei];
			const label nei = neighbour[facei];
			hTildeVol_value[own] += hTildeVol_value_own[facei];
			hTildeVol_value[nei] += hTildeVol_value_nei[facei];
		}

		free(gTildeVol_value_own);
		free(gTildeVol_value_nei);
		free(hTildeVol_value_own);
		free(hTildeVol_value_nei);
#endif

}

template <template<class> class PatchType, class GeoMesh> 
Foam::tmp<Foam::GeometricField<Foam::scalar, PatchType, GeoMesh> >
Foam::discreteVelocity::equilibriumMaxwell
(
    const GeometricField<scalar, PatchType, GeoMesh>& rho,
    const GeometricField<vector, PatchType, GeoMesh>& U,
    const GeometricField<scalar, PatchType, GeoMesh>& T
)
{
    const dimensionedScalar R = dvm_.R();
    const label D = mesh_.nSolutionD();
    const dimensionedScalar vUnit("vUnit", dimLength/dimTime, 1);

    tmp<GeometricField<scalar, PatchType, GeoMesh> > tEqu
    (
        new GeometricField<scalar, PatchType, GeoMesh>
        (
            IOobject
            (
                "equ",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar( "0", gTildeVol_.dimensions(), 0)
        )
    );
    GeometricField<scalar, PatchType, GeoMesh>& equ = tEqu();
    //equ = rho/pow(sqrt(2.0*pi*R*T),D)*exp(-magSqr(U - xi_)/(2.0*R*T))
        ///pow(vUnit, 3-D);
    equ = rho/((2.0*pi*R*T))*exp(-magSqr(U - xi_)/(2.0*R*T))/vUnit;
    return tEqu;
}

template <template<class> class PatchType, class GeoMesh> 
void Foam::discreteVelocity::equilibriumShakhovSurf
(
    GeometricField<scalar, PatchType, GeoMesh>& gEq,
    GeometricField<scalar, PatchType, GeoMesh>& hEq,
    const GeometricField<scalar, PatchType, GeoMesh>& rho,
    const GeometricField<vector, PatchType, GeoMesh>& U,
    const GeometricField<scalar, PatchType, GeoMesh>& T,
    const GeometricField<vector, PatchType, GeoMesh>& q
)
{
    label D = mesh_.nSolutionD();
    label K = dvm_.KInner();
    dimensionedScalar vUnit("vUnit", dimLength/dimTime, 1);

    dimensionedScalar R = dvm_.R();
    dimensionedScalar Pr = dvm_.Pr();
    GeometricField<scalar, PatchType, GeoMesh> cSqrByRT = magSqr(U - xi_)/(R*T);
    GeometricField<scalar, PatchType, GeoMesh> cqBy5pRT = ((xi_ - U)&q)/(5.0*rho*R*T*R*T);
    GeometricField<scalar, PatchType, GeoMesh> gEqBGK  = rho/pow(sqrt(2.0*pi*R*T),D)*exp(-cSqrByRT/2.0)/pow(vUnit, 3-D);
    gEq = ( 1.0 + (1.0 - Pr)*cqBy5pRT*(cSqrByRT - D - 2.0) )*gEqBGK;
    hEq = ( (K + 3.0 - D) + (1.0 - Pr)*cqBy5pRT*((cSqrByRT - D)*(K + 3.0 - D) - 2*K) )*gEqBGK*R*T;
}


template <template<class> class PatchType, class GeoMesh> 
void Foam::discreteVelocity::equilibriumShakhovVol
(
    GeometricField<scalar, PatchType, GeoMesh>& gEq,
    GeometricField<scalar, PatchType, GeoMesh>& hEq,
    const GeometricField<scalar, PatchType, GeoMesh>& rho,
    const GeometricField<vector, PatchType, GeoMesh>& U,
    const GeometricField<scalar, PatchType, GeoMesh>& T,
    const GeometricField<vector, PatchType, GeoMesh>& q
)
{
    label D = mesh_.nSolutionD();
    
    label K = dvm_.KInner();
    dimensionedScalar vUnit("vUnit", dimLength/dimTime, 1);
#ifndef ATHREAD
    dimensionedScalar R = dvm_.R();
    dimensionedScalar Pr = dvm_.Pr();
    GeometricField<scalar, PatchType, GeoMesh> cSqrByRT = magSqr(U - xi_)/(R*T);
    GeometricField<scalar, PatchType, GeoMesh> cqBy5pRT = ((xi_ - U)&q)/(5.0*rho*R*T*R*T);
    GeometricField<scalar, PatchType, GeoMesh> gEqBGK  = rho/pow(sqrt(2.0*pi*R*T),D)*exp(-cSqrByRT/2.0)/pow(vUnit, 3-D);
    gEq = ( 1.0 + (1.0 - Pr)*cqBy5pRT*(cSqrByRT - D - 2.0) )*gEqBGK;
    hEq = ( (K + 3.0 - D) + (1.0 - Pr)*cqBy5pRT*((cSqrByRT - D)*(K + 3.0 - D) - 2*K) )*gEqBGK*R*T;
#else
    double R_value = dvm_.R().value();
    double Pr_value = dvm_.Pr();
    double vUnit_value = vUnit.value();
    vector xii = xi_.value();
    int64 gEq_size = gEq.size();
    double *gEq_value = gEq.begin();
    double *hEq_value = hEq.begin();
    ////
    int64 U_size = U.size();
    scalar *U_value = const_cast<double*>((U.begin())->v_ );
    scalar *q_value = const_cast<double*>((q.begin())->v_) ;
    scalar *rho_value = const_cast<double*>(rho.begin());
    scalar *T_value = const_cast<double*>(T.begin());
    ////
    double xii_x = xii.x();
	double xii_y = xii.y();
	double xii_z = xii.z();
    ///
    EquShakhov equshakhov;
	equshakhov.U_value = U_value;
	equshakhov.q_value = q_value;
	equshakhov.T_value = T_value;
    equshakhov.rho_value = rho_value;
	equshakhov.gEq_value = gEq_value;
	equshakhov.hEq_value = hEq_value;
	equshakhov.xii_x = xii.x();
	equshakhov.xii_y = xii.y();
	equshakhov.xii_z = xii.z();
	equshakhov.gEq_size = gEq_size;
	equshakhov.R_value = R_value;
	equshakhov.Pr_value = Pr_value;
	equshakhov.vUnit_value = vUnit_value;
    equshakhov.D = D;//label
	equshakhov.K = K;//label
	__real_athread_spawn((void *)slave_Func_equshakhov, &equshakhov);
	athread_join();
#endif
}
void Foam::discreteVelocity::equilibriumMaxwell
(
    Foam::fvsPatchScalarField& geq,
    Foam::fvsPatchScalarField& heq,
    const Foam::fvPatchScalarField&  rho,
    const Foam::fvPatchVectorField&    U,
    const Foam::fvPatchScalarField&    T
)
{
    scalar Ri = dvm_.R().value(); 
    label D = mesh_.nSolutionD();
    vector xii = xi_.value();

    geq == rho/pow(sqrt(2.0*pi*Ri*T),D)*exp(-magSqr(U - xii)/(2.0*Ri*T));
    heq == (dvm_.KInner() + 3-D)*Ri*T*geq;

}


Foam::scalar Foam::discreteVelocity::equilibriumMaxwellByRho
(
 const  vector U,
 const  scalar T
)
{
    scalar feqByRho;
    scalar Ri = dvm_.R().value(); //but R has dimensionSet
    label D = mesh_.nSolutionD();
    vector xii = xi_.value();
    feqByRho = 1.0/pow(sqrt(2.0*pi*Ri*T),D)*exp(-magSqr(U - xii)/(2.0*Ri*T));
    return feqByRho;
}

// add
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, scalar>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
// Foam::fv::gaussGrad<scalar>::calcGrad
/*Foam::discreteVelocity::calcGrad
(
    const GeometricField<scalar, fvPatchField, volMesh>& vsf//,
    // const word& name
) const
{
    typedef typename outerProduct<vector, scalar>::type GradType;
    tmp<GeometricField<GradType, fvPatchField, volMesh> > tgGrad
    (
        // gradf(tinterpScheme_().interpolate(vsf), name)
        // Foam::fv::gaussGrad<scalar>::gradf(Foam::fv::gaussGrad<scalar>::tinterpScheme_().interpolate(vsf), "interpolate")

        Foam::fv::gaussGrad<scalar>::gradf(tinterpScheme_().Foam::surfaceInterpolationScheme<scalar>::interpolate(vsf), "interpolate") //可以
    );
    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad();

    Foam::fv::gaussGrad<scalar>::correctBoundaryConditions(vsf, gGrad);

    return tgGrad;
}
*/

Foam::discreteVelocity::calcGrad
(
    // GeometricField<scalar, fvPatchField, volMesh>& vsf//,
    const GeometricField<scalar, fvPatchField, volMesh>& vsf//,
    // const word& name
) //const
{
    // typedef typename outerProduct<vector, scalar>::type GradType;
    tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> > tinPo(
        tinterpScheme_().Foam::surfaceInterpolationScheme<scalar>::interpolate(vsf)
    );
    // GeometricField<scalar, fvsPatchField, surfaceMesh>& inPo = tinPo();
    GeometricField<scalar, fvsPatchField, surfaceMesh>& ssf = tinPo();

    // tmp<GeometricField<GradType, fvPatchField, volMesh> > tgGrad
    // (
    //     new GeometricField<GradType, fvPatchField, volMesh>
    //     (
    //         IOobject
    //         (
    //             // name,
    //             "interpolate",
    //             ssf.instance(),
    //             mesh_,
    //             IOobject::NO_READ,
    //             IOobject::NO_WRITE
    //         ),
    //         mesh_,
    //         dimensioned<GradType>
    //         (
    //             "0",
    //             ssf.dimensions()/dimLength,
    //             pTraits<GradType>::zero
    //         ),
    //         zeroGradientFvPatchField<GradType>::typeName
    //     )
    // );
    // GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad();
    GeometricField<vector, fvPatchField, volMesh>& gGrad = gBarPgrad_;
    Field<vector>& igGrad = gGrad;
    // GeometricField<vector, fvPatchField, volMesh>& igGrad = gBarPgrad_;
    // Field<vector>& igGrad = gBarPgrad_;
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const vectorField Sf = mesh_.Sf();
    // const vectorField& V = mesh_.V();
    const scalarField V = mesh_.V();

    
    const Field<scalar>& issf = ssf;

    forAll(owner, facei)
    {
        // vector Sfssf = Sf[facei]*issf[facei];

        // igGrad[owner[facei]] += Sfssf;
        // igGrad[neighbour[facei]] -= Sfssf;
        igGrad[owner[facei]] += Sf[facei]*issf[facei];
        igGrad[neighbour[facei]] -= Sf[facei]*issf[facei];
    }

    forAll(mesh_.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh_.boundary()[patchi].faceCells();

        const vectorField& pSf = mesh_.Sf().boundaryField()[patchi];

        const fvsPatchField<scalar>& pssf = ssf.boundaryField()[patchi];

        forAll(mesh_.boundary()[patchi], facei)
        {
            igGrad[pFaceCells[facei]] += pSf[facei]*pssf[facei];
        }
    }

    // igGrad /= V;
    // gBarPgrad_ /= mesh_.V();
    igGrad /= V;

    gGrad.correctBoundaryConditions();

    // return tgGrad;

    // tmp<GeometricField<GradType, fvPatchField, volMesh> > tgGrad
    // (
    //     Foam::fv::gaussGrad<scalar>::gradf(inPo, "interpolate") //可以
    // );
    // GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad();

    Foam::fv::gaussGrad<scalar>::correctBoundaryConditions(vsf, gGrad);

    // return tgGrad;
    return gBarPgrad_;
}

// ************************************************************************* //
