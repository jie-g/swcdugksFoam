/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
	\\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include <mpi.h>
#include "fvDVM.H"
#include "constants.H"
#include "fvm.H"
#include "calculatedMaxwellFvPatchField.H"
#include "symmetryModFvPatchField.H"
#include "pressureInFvPatchField.H"
#include "pressureOutFvPatchField.H"
#include "scalarIOList.H"
#include "fieldMPIreducer.H"

#include <sys/time.h>
#include <time.h>
#define CLOCKRATE 1.45e9

typedef long int int64;//label

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

#if FOAM_MAJOR <= 3
#define BOUNDARY_FIELD_REF boundaryField()
#else
#define BOUNDARY_FIELD_REF boundaryFieldRef()
#endif

// #define GETCONUM

#include "para.h"

#ifdef ATHREAD_TEMP
#include "slave_para.h"
extern "C"{
	#include <athread.h>
	//parallel computing functions for CPEs
	void SLAVE_FUN(Func_macrovol)(MacroVol *master);
	void SLAVE_FUN(Func_macrosurf)(MacroSurf *macsurf);
	void SLAVE_FUN(Func_macroqsurf)(MacroqSurf *macqsurf);
	void SLAVE_FUN(Func_ghsurf)(GHsurf *ghsurf);
	void SLAVE_FUN(Func_ghsurf1)(GHsurf *ghsurf);
	void SLAVE_FUN(Func_ghbarpvol)(GHbarPvol *ghbarpvol);
	void SLAVE_FUN(Func_macroqvol)(MacroqVol *macroqvol);
	void SLAVE_FUN(Func_conum)(CONUM *conum);
}
#endif
// record time
static inline unsigned long rpcc(){
	unsigned long addtime;
	asm("rtc %0": "=r" (addtime) : );
	return addtime;
 }
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(fvDVM, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvDVM::setDVgrid
(
	scalarField& weights,
	scalarField& Xis,
	scalar xiMin,
	scalar xiMax,
	label nXi
)
{
	// Read from file ./constant/Xis and ./constant/weights
	scalarIOList xiList
	(
		IOobject
		(
			"Xis",
			time_.caseConstant(),
			mesh_,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

	scalarIOList weightList
	(
		IOobject
		(
			"weights",
			time_.caseConstant(),
			mesh_,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

	for (label i = 0; i < nXi; i++)
	{
		weights[i] = weightList[i];
		
		Xis[i] = xiList[i];
	}
}

void Foam::fvDVM::initialiseDV()
{
	const vectorField &C = mesh_.C();

	scalarField weights1D(nXiPerDim_);
	scalarField Xis(nXiPerDim_);
        xiMax_.value()=0.0;

    //we have polyMesh and vMesh at the beginning
	string mesh_name = args_.path() / word("constant") / word("polyMesh");//polyMesh dir
	string vmesh_name = args_.path() / word("constant") / word("vMesh");//vMesh dir
	string temp_name = args_.path() / word("constant") / word("tMesh");//temp name

	//1.change folder name, polyMesh-->tMesh,vMesh-->polyMesh
	//as OF only recognizes certain path in the couse of reading mesh files
	
	if(mpiReducer_.rank() == 0 || args_.optionFound("parallel")) {
          std::rename(mesh_name.c_str(), temp_name.c_str());
          std::rename(vmesh_name.c_str(), mesh_name.c_str());
       }
        if (mpiReducer_.dvParallel())
            MPI_Barrier(MPI_COMM_WORLD);

	//2.read and construct vmesh 
	Foam::fvMesh vmesh
	(
		Foam::IOobject
		(
			Foam::fvMesh::defaultRegion,
			time_.timeName(),
			time_,
			Foam::IOobject::MUST_READ
		)
	);


	//3.change the folder name back, polyMesh-->vMesh,tMesh-->polyMesh,
	//so no one would know we change the name secretly :D
    if (mpiReducer_.dvParallel())
            MPI_Barrier(MPI_COMM_WORLD);
        if(mpiReducer_.rank() == 0 || args_.optionFound("parallel")) {
            std::rename(mesh_name.c_str(), vmesh_name.c_str());
            std::rename(temp_name.c_str(), mesh_name.c_str());
       }


	scalarField weightsGlobal;
	vectorField XisGlobal;
	labelField  symmXtgID;
	labelField  symmYtgID;
	labelField  symmZtgID;

	if (mesh_.nSolutionD() == 3)    //3D(X & Y & Z)
	{
		nXiX_ = nXiY_ = nXiZ_ = vmesh.C().size();
		nXi_ = vmesh.C().size();

		weightsGlobal.setSize(nXi_);
		XisGlobal.setSize(nXi_);
		symmXtgID.setSize(nXi_);
		symmYtgID.setSize(nXi_);
		symmZtgID.setSize(nXi_);

		label i;
		for (i = 0; i < nXi_; i++) {
			vector xi(vmesh.C()[i].x(), vmesh.C()[i].y(), vmesh.C()[i].z());
			scalar weight(vmesh.V()[i]);
			weightsGlobal[i] = weight;
			XisGlobal[i] = xi;
			xiMax_.value()=max(xiMax_.value(),mag(xi));
		}
		/*

				for (label iz = 0; iz < nXiZ_; iz++)
				{
					for (label iy = 0; iy < nXiY_; iy++)
					{
						for (label ix = 0; ix < nXiZ_; ix++)
						{
							scalar weight = weights1D[iz]*weights1D[iy]*weights1D[ix];
							vector xi(Xis[ix], Xis[iy], Xis[iz]);
							weightsGlobal[i] = weight;
							XisGlobal[i] = xi;
							symmXtgID[i] = iz*nXiY_*nXiX_ + iy*nXiX_ + (nXiX_ - ix -1);
							symmYtgID[i] = iz*nXiY_*nXiX_ + (nXiY_ - iy - 1)*nXiX_ + ix;
							symmZtgID[i] = (nXiZ_ - iz -1)*nXiY_*nXiX_ + iy*nXiX_ + ix;
							i++;
						}
					}
				}
		*/
	}
	else
	{
		if (mesh_.nSolutionD() == 2)    //2D (X & Y)
		{
			nXiX_ = nXiY_ = vmesh.C().size();
			nXiZ_ = 1;
			nXi_ = vmesh.C().size();
			weightsGlobal.setSize(nXi_);
			XisGlobal.setSize(nXi_);
			symmXtgID.setSize(nXi_);
			symmYtgID.setSize(nXi_);
			symmZtgID.setSize(nXi_);
			label i;
			for (i = 0; i < nXi_; i++) {
				vector xi(vmesh.C()[i].x(), vmesh.C()[i].y(), 0.0);
				scalar weight(vmesh.V()[i]);
				weightsGlobal[i] = weight;
				XisGlobal[i] = xi;
				xiMax_.value()=max(xiMax_.value(),mag(xi));
			}
			/*
						for (label iy = 0; iy < nXiY_; iy++)
						{
							for (label ix = 0; ix < nXiX_; ix++)
							{
								scalar weight = weights1D[iy]*weights1D[ix]*1;
								vector xi(Xis[ix], Xis[iy], 0.0);
								weightsGlobal[i] = weight;
								XisGlobal[i] = xi;
								symmXtgID[i] = iy*nXiX_ + (nXiX_ - ix -1);
								symmYtgID[i] = (nXiY_ - iy - 1)*nXiX_ + ix;
								symmZtgID[i] = 0;
								i++;
							}
						}
			*/
		}
		else    //1D (X)
		{
			nXiX_ = vmesh.C().size();
			nXiY_ = nXiZ_ = 1;
			nXi_ = vmesh.C().size();
			weightsGlobal.setSize(nXi_);
			XisGlobal.setSize(nXi_);
			symmXtgID.setSize(nXi_);
			symmYtgID.setSize(nXi_);
			symmZtgID.setSize(nXi_);
			label i;
			for (i = 0; i < nXi_; i++) {
				vector xi(vmesh.C()[i].x(), 0.0, 0.0);
				scalar weight(vmesh.V()[i]);
				weightsGlobal[i] = weight;
				XisGlobal[i] = xi;
				xiMax_.value()=max(xiMax_.value(),mag(xi));
			}
			/*
						for (label ix = 0; ix < nXiX_; ix++)
						{
							scalar weight = weights1D[ix]*1*1;
							vector xi(Xis[ix], 0.0, 0.0);
							weightsGlobal[i] = weight;
							XisGlobal[i] = xi;
							symmXtgID[i] = (nXiX_ - ix -1);
							symmYtgID[i] = 0;
							symmZtgID[i] = 0;
							i++;
						}
			*/
		}
	}

	if (mpiReducer_.rank() == 0)
	{
		Info << "fvDVM : Allocated " << XisGlobal.size()
			<< " discrete velocities" << endl;
	}
	label nA = nXi_ / mpiReducer_.csize();//average
	label nB = nXi_ - nA * mpiReducer_.csize();//mod
	label nXiPart = nA + (label)(mpiReducer_.crank() < nB);//cyclic partition
	DV_.setSize(nXiPart);
	if (mpiReducer_.rank() %mpiReducer_.npd()==0) //first column
	{
		std::cout << "myrank    " << mpiReducer_.crank() <<std:: endl;
		std::cout << "nXisPart " << nXiPart << std::endl;
        //std::cout<<mpiReducer_.crank()<<" "<<mpiReducer_.rank()<<endl;
	}

	label chunk = 0;
	label gid = 0;
        

        
	forAll(DV_, i)
	{
		gid = chunk + mpiReducer_.crank();
		DV_.set
		(
			i,
			new discreteVelocity
			(
				*this,
				mesh_,
				time_,
				weightsGlobal[gid],
				dimensionedVector("xi", dimLength / dimTime, XisGlobal[gid]),
				i,
				symmXtgID[gid],
				symmYtgID[gid],
				symmZtgID[gid]
			)
		);
		chunk += mpiReducer_.csize();
               
	}
    /***************************Rearrange the data in advance******************************************/
	//Storing the GeometricField data in an array in advance facilitates data transfer between the MPE and CPEs later on.
	int64 DVsize = DV_.size();
	//scalar *xi_value = (scalar*)malloc(sizeof(scalar) * DVsize*3);
	xii_value = (scalar*)malloc(sizeof(scalar) * DVsize*3);
	for(int DVid=0;DVid<DVsize;DVid++){
		discreteVelocity &dv = DV_[DVid];
		vector xii = dv.xi().value();
		xii_value[DVid*3] = xii.x();
		xii_value[DVid*3+1] = xii.y();
		xii_value[DVid*3+2] = xii.z();
	}
	weight_ = (scalar*)malloc(sizeof(scalar) * DVsize);
	for(int DVid=0;DVid<DVsize;DVid++){
		weight_[DVid] = DV_[DVid].weight();
	}
	gBarPvol_value1 = (scalar**)malloc(sizeof(scalar*) * DVsize );
	hBarPvol_value1 = (scalar**)malloc(sizeof(scalar*) * DVsize );
	gTildeVol_value1 = (scalar**)malloc(sizeof(scalar*) * DVsize );
	hTildeVol_value1 = (scalar**)malloc(sizeof(scalar*) * DVsize );
	for(int DVid=0;DVid<DVsize;DVid++){
		discreteVelocity &dv = DV_[DVid];
		gBarPvol_value1[DVid]=const_cast<scalar*>(dv.gBarPvol_.begin());
		hBarPvol_value1[DVid]=const_cast<scalar*>(dv.hBarPvol_.begin());
		gTildeVol_value1[DVid]=const_cast<scalar*>(dv.gTildeVol_.begin());
		hTildeVol_value1[DVid]=const_cast<scalar*>(dv.hTildeVol_.begin());
	}
	gSurf_value1 = (scalar**)malloc(sizeof(scalar*) * DVsize );
	hSurf_value1 = (scalar**)malloc(sizeof(scalar*) * DVsize );
	for(int DVid=0;DVid<DVsize;DVid++){
		discreteVelocity &dv = DV_[DVid];
		gSurf_value1[DVid]=const_cast<scalar*>(dv.gSurf_.begin());
		hSurf_value1[DVid]=const_cast<scalar*>(dv.hSurf_.begin());
	}
	const scalarField& V = mesh_.V();
	int64 ownersize = mesh_.owner().size();
	V_own1 = (scalar*)malloc(sizeof(scalar) * ownersize);
	V_nei1 = (scalar*)malloc(sizeof(scalar) * ownersize);
	const labelUList& owner = mesh_.owner();
	const labelUList& neighbour = mesh_.neighbour();
	for(int facei=0;facei<ownersize;facei++)
	{
		const label own = owner[facei];
		const label nei = neighbour[facei];
		V_own1[facei] = V[own];
		V_nei1[facei] = V[nei];
	}
	scalar *C_value = const_cast<double*>((mesh_.C().begin())->v_);
	C_own1 = (scalar*)malloc(sizeof(scalar) * ownersize * 3);
	C_nei1 = (scalar*)malloc(sizeof(scalar) * ownersize * 3);
    for(int facei=0;facei<ownersize;facei++)
	{
		const label own = owner[facei];
		const label nei = neighbour[facei];
		C_own1[facei*3] = C_value[own*3];
        C_own1[facei*3+1] = C_value[own*3+1];
        C_own1[facei*3+2] = C_value[own*3+2];

		C_nei1[facei*3] = C_value[nei*3];
        C_nei1[facei*3+1] = C_value[nei*3+1];
        C_nei1[facei*3+2] = C_value[nei*3+2];
	}

	interpola_weight = (scalar*)malloc(sizeof(scalar) * ownersize);
	const vectorField& Cf = mesh_.Cf();
    // const vectorField& C = mesh_.C();
    const vectorField& Sf = mesh_.Sf();
	for(int facei=0;facei<ownersize;facei++)
    {
        label own = owner[facei];
        label nei = neighbour[facei];
        scalar SfdOwn = mag(Sf[facei] & (Cf[facei] - C[own]));
        scalar SfdNei = mag(Sf[facei] & (C[nei] - Cf[facei]));
        interpola_weight[facei] = SfdNei/(SfdOwn + SfdNei);
    }
	/*********************************************************************/
	//clear vmesh data structure to save memory
	vmesh.clearOut();
        
}


void Foam::fvDVM::setCalculatedMaxwellRhoBC()
{
#if FOAM_MAJOR <= 3
	GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField&
		rhoBCs = rhoVol_.boundaryField();
#else
	GeometricField<scalar, fvPatchField, volMesh>::Boundary&
		rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
#endif
	forAll(rhoBCs, patchi)
	{
		if (rhoBCs[patchi].type() == "calculatedMaxwell")
		{
           
			const vectorField& SfPatch = mesh_.Sf().boundaryField()[patchi];
			calculatedMaxwellFvPatchField<scalar>& rhoPatch =
				refCast<calculatedMaxwellFvPatchField<scalar> >(rhoBCs[patchi]);
			const vectorField& Upatch = Uvol_.boundaryField()[patchi];
			const scalarField& Tpatch = Tvol_.boundaryField()[patchi];

			forAll(rhoPatch, facei)
			{   
				//Info<< mesh_.Cf().boundaryField()[patchi][facei]<<endl;
				vector faceSf = SfPatch[facei];
				rhoPatch.inComingByRho()[facei] = 0; // set to zero
                                
				forAll(DV_, dvi) // add one by one
				{
					vector xi = DV_[dvi].xi().value();
					scalar weight = DV_[dvi].weight();
					if ((xi & faceSf) < 0) //inComing
					{
						rhoPatch.inComingByRho()[facei] +=
							-weight * (xi & faceSf)
							* DV_[dvi].equilibriumMaxwellByRho
							(
								Upatch[facei],
								Tpatch[facei]
							);
					}
				}
                               
			}

			if (mpiReducer_.dvParallel()&& mpiReducer_.npd() < mpiReducer_.nproc()) {
					mpiReducer_.reduceField(rhoPatch.inComingByRho());

			}
			 
		}

	}
}

void Foam::fvDVM::setSymmetryModRhoBC()
{
	//prepare the container (set size) to store all DF on the patchi
#if FOAM_MAJOR <= 3
	GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField&
		rhoBCs = rhoVol_.boundaryField();
#else
	GeometricField<scalar, fvPatchField, volMesh>::Boundary&
		rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
#endif
	forAll(rhoBCs, patchi)
	{
		label ps = rhoBCs[patchi].size();
		if (rhoBCs[patchi].type() == "symmetryMod")
		{
			symmetryModFvPatchField<scalar>& rhoPatch =
				refCast<symmetryModFvPatchField<scalar> >(rhoBCs[patchi]);
			rhoPatch.dfContainer().setSize(ps * nXi_ * 2); //*2 means g and h
		}
	}
}


void Foam::fvDVM::updateGHbarPvol()
{
	unsigned long time_sta = 0;
	unsigned long time_end = 0;
	time_sta = rpcc();
#ifndef GHBARPVOL
	forAll(DV_, DVid)
		DV_[DVid].updateGHbarPvol();
#else
	int64 DVsize = DV_.size();
	label D = mesh_.nSolutionD();
	label K = KInner_;
	dimensionedScalar vUnit("vUnit", dimLength/dimTime, 1);
	double R_value = R_.value();
	double Pr_value = Pr_;
	double vUnit_value = vUnit.value();
	scalar *U_value = const_cast<double*>((Uvol_.begin())->v_ );
	scalar *q_value = const_cast<double*>((qVol_.begin())->v_) ;
	scalar *rho_value = const_cast<double*>(rhoVol_.begin());
	scalar *T_value = const_cast<double*>(Tvol_.begin());
	scalar dt = time_.deltaTValue();
	scalar *tauVol_value = const_cast<double*>(tauVol_.begin());

	scalar **gTildeVol_value = (scalar**)malloc(sizeof(scalar*) * DVsize );
	scalar **hTildeVol_value = (scalar**)malloc(sizeof(scalar*) * DVsize );

	for(int DVid=0;DVid<DVsize;DVid++){
		discreteVelocity &dv = DV_[DVid];
		gTildeVol_value[DVid]=const_cast<double*>(dv.gTildeVol_.begin());
		hTildeVol_value[DVid]=const_cast<double*>(dv.hTildeVol_.begin());
	}
	
	int64 size = DV_[0].gBarPvol_.size();
	GHbarPvol ghbarpvol;
	ghbarpvol.U_value = U_value;
	ghbarpvol.q_value = q_value;
	ghbarpvol.T_value = T_value;
	ghbarpvol.rho_value = rho_value;
	ghbarpvol.tauVol_value = tauVol_value;
	ghbarpvol.gBarPvol_value = gBarPvol_value1;
	ghbarpvol.hBarPvol_value = hBarPvol_value1;
	ghbarpvol.gTildeVol_value = gTildeVol_value;
	ghbarpvol.hTildeVol_value = hTildeVol_value;
	ghbarpvol.dt = dt;
	ghbarpvol.ghBarPvol_size = size;
	ghbarpvol.xi_value = xii_value;
	ghbarpvol.R_value = R_value;
	ghbarpvol.Pr_value = Pr_value;
	ghbarpvol.vUnit_value = vUnit_value;
	ghbarpvol.D = D;//label
	ghbarpvol.K = K;//label
	ghbarpvol.DVsize = DVsize;//label
	
	__real_athread_spawn((void *)slave_Func_ghbarpvol, &ghbarpvol);
	
	athread_join();
	
	free(gTildeVol_value);
	free(hTildeVol_value);

    // dimensionedScalar dt = time_.deltaT();
    // relaxFactor = 1.5*dt/(2.0*dvm_.tauVol() + dt);
    // gBarPvol_ = (1.0 - relaxFactor)*gTildeVol_ + relaxFactor*gEq;
    // hBarPvol_ = (1.0 - relaxFactor)*hTildeVol_ + relaxFactor*hEq;
	
/*	forAll(DV_, DVid){
		DV_[DVid].gBarPvol_.correctBoundaryConditions(); // NOTE: check if the newly defined zeroGradientFvsPatchField 
    	DV_[DVid].hBarPvol_.correctBoundaryConditions();
	}
*/   
#endif
	time_end = rpcc();
	if(mpiReducer_.rank() == 0){
		Info << "updateGHbarPvol time =" <<((double)(time_end-time_sta)*1000/CLOCKRATE) << "ms" << nl << endl;
	}
}


void Foam::fvDVM::updateGHbarSurf()
{
	unsigned long time_sta = 0;
	unsigned long time_end = 0;
    time_sta = rpcc();
	forAll(DV_, DVid)
		DV_[DVid].updateGHbarSurf();
	time_end = rpcc();
	if(mpiReducer_.rank() == 0){
		Info << "updateGHbarSurf time =" <<((double)(time_end-time_sta)*1000/CLOCKRATE) << "ms" << nl << endl;
	}
}


void Foam::fvDVM::updateMaxwellWallRho()
{
	unsigned long time_sta = 0;
	unsigned long time_end = 0;
	time_sta = rpcc();
#if FOAM_MAJOR <= 3
	GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField&
		rhoBCs = rhoVol_.boundaryField();
#else
	GeometricField<scalar, fvPatchField, volMesh>::Boundary&
		rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
#endif
	forAll(rhoBCs, patchi)
	{
		if (rhoBCs[patchi].type() == "calculatedMaxwell")
		{
			calculatedMaxwellFvPatchField<scalar>& rhoPatch =
				refCast<calculatedMaxwellFvPatchField<scalar> >(rhoBCs[patchi]);
			if (mpiReducer_.dvParallel()&&mpiReducer_.npd() < mpiReducer_.nproc())
				mpiReducer_.reduceField(rhoPatch.outGoing());
		}
	}
	rhoVol_.correctBoundaryConditions();
	time_end = rpcc(); 

}

void Foam::fvDVM::updateGHbarSurfMaxwellWallIn()
{
	unsigned long time_sta = 0;
	unsigned long time_end = 0;
	time_sta = rpcc();
	forAll(DV_, DVid)
		DV_[DVid].updateGHbarSurfMaxwellWallIn();
	time_end = rpcc();
}

void Foam::fvDVM::updateGHbarSurfSymmetryIn()
{
	unsigned long time_sta = 0;
	unsigned long time_end = 0;
	time_sta = rpcc();
	//1. copy all DV's g/h to rho patch's dfContainer
	//2. MPI_Allgather the rho patch's dfContainer
	//if(args_.optionFound("dvParallel"))
	//{
	label rank = mpiReducer_.rank();
	label nproc = mpiReducer_.nproc();
#if FOAM_MAJOR <= 3
	GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField&
		rhoBCs = rhoVol_.boundaryField();
#else
	GeometricField<scalar, fvPatchField, volMesh>::Boundary&
		rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
#endif
	forAll(rhoBCs, patchi)
	{
		label ps = rhoBCs[patchi].size();
		if (rhoBCs[patchi].type() == "symmetryMod")
		{
			symmetryModFvPatchField<scalar>& rhoPatch =
				refCast<symmetryModFvPatchField<scalar> >(rhoBCs[patchi]);
			//compose the recvcout and displacement array
			//labelField recvc(nproc);
			//labelField displ(nproc);
			Field<int> recvc(nproc);
			Field<int> displ(nproc);
                        label chunck = nXi_ / nproc;
			label left = nXi_ % nproc;
			forAll(recvc, i)
			{
				recvc[i] = 2 * ps * (chunck + (i < left));
				if (i <= left)
					displ[i] = i * 2 * ps * (chunck + 1); // (i<=nXi_%nproc)
				else
					displ[i] = 2 * ps * (left * (chunck + 1) + (i - left) * (chunck));
			}

			// check 12*28+15 dv's g
			label did = 1709;
			label pp = did % nproc;
			label lid = did / nproc;
			//if(rank==pp)
			//{
				//Info << "processing by rank " << rank << endl;
				//Info << "12*28+15 outging g " << DV_[lid].gSurf()[0] << endl;
				//Info << "12*28+15 outging xi " << DV_[lid].xi() <<endl;
				//Info << "12*28+15 outging at boundary " << DV_[lid].gSurf().boundaryField()[patchi][0] << endl;
			//}
			// memcpy each dv's g/h to rho
			forAll(DV_, DVid)
			{
				//label shift = (nXi_ / nproc * rank + DVid)*2*ps;
				label shift = displ[rank] + DVid * 2 * ps;
				memcpy((rhoPatch.dfContainer().data() + shift),
					DV_[DVid].gSurf().boundaryField()[patchi].cdata(), ps * sizeof(scalar));
				memcpy((rhoPatch.dfContainer().data() + shift + ps),
					DV_[DVid].hSurf().boundaryField()[patchi].cdata(), ps * sizeof(scalar));
			}

			// check 
			//if(rank == pp)
				//Info << "dv gid 1709's g = " <<rhoPatch.dfContainer()[displ[pp]+lid*2*ps+32]<< endl;;


			//Allgather
			MPI_Allgatherv(
				//rhoPatch.dfContainer().data() + displ[rank],//2*ps*nXI_/nproc*rank, //send*
				MPI_IN_PLACE,
				2 * ps * DV_.size(), //(how many DV i processed) * 2 * patch size
				MPI_DOUBLE,
				rhoPatch.dfContainer().data(),
				recvc.data(),
				displ.data(),
				MPI_DOUBLE,
				MPI_COMM_WORLD
			);
		}
	}
	forAll(DV_, DVid)
		DV_[DVid].updateGHbarSurfSymmetryIn();
	time_end = rpcc();
}

void Foam::fvDVM::updateMacroSurf()
{
	unsigned long time_sta = 0;
	unsigned long time_end = 0;
	unsigned long updateMacroSurf_sta0 = 0;
	unsigned long updateMacroSurf_sta1 = 0;
	unsigned long updateMacroSurf_end0 = 0;
	unsigned long updateMacroSurf_end1 = 0;
	unsigned long updateMacroSurf_sta2 = 0;
	unsigned long updateMacroSurf_sta3 = 0;
	unsigned long updateMacroSurf_end2 = 0;
	unsigned long updateMacroSurf_end3 = 0;
	time_sta = rpcc();
	// Init to zero before add one DV by one DV
	rhoSurf_ = dimensionedScalar("0", rhoSurf_.dimensions(), 0);
	Usurf_ = dimensionedVector("0", Usurf_.dimensions(), vector(0, 0, 0));
	Tsurf_ = dimensionedScalar("0", Tsurf_.dimensions(), 0);
	qSurf_ = dimensionedVector("0", qSurf_.dimensions(), vector(0, 0, 0));
	stressSurf_ = dimensionedTensor("0",stressSurf_.dimensions(),pTraits<tensor>::zero);
	
#ifndef MACROSURF
	surfaceVectorField rhoUsurf = rhoSurf_ * Usurf_;
	surfaceScalarField rhoEsurf = rhoSurf_ * magSqr(Usurf_);
	// if(mpiReducer_.rank() == 0)
	// 	Info << "Original code "<< nl << endl;
	updateMacroSurf_sta0 = rpcc();
	forAll(DV_, dvi)
	{
		discreteVelocity& dv = DV_[dvi];
		rhoSurf_ += dXiCellSize_ * dv.weight() * dv.gSurf();
		rhoUsurf += dXiCellSize_ * dv.weight() * dv.gSurf() * dv.xi();
		rhoEsurf += 0.5 * dXiCellSize_ * dv.weight()* (dv.gSurf() * magSqr(dv.xi())+ dv.hSurf());
	}
	updateMacroSurf_end0 = rpcc();
#else
	surfaceVectorField rhoUsurf = Usurf_;
	surfaceScalarField rhoEsurf = rhoSurf_ ;

	updateMacroSurf_sta0 = rpcc();
	int64 DVsize = DV_.size();
	int64 ownersize = mesh_.owner().size();
	scalar *rhoSurf_value = rhoSurf_.begin();
	scalar *rhoUsurf_value = (rhoUsurf.begin())->v_ ;
	scalar *rhoEsurf_value = rhoEsurf.begin();
	scalar *Usurf_value = (Usurf_.begin())->v_ ;
	double weight_value = 0;
/*************************************************************************************/
/*		scalar **gSurf_value = (scalar**)malloc(sizeof(scalar*) * DVsize );
		scalar **hSurf_value = (scalar**)malloc(sizeof(scalar*) * DVsize );
		for(int DVid=0;DVid<DVsize;DVid++){
			discreteVelocity &dv = DV_[DVid];
			gSurf_value[DVid]=const_cast<double*>(dv.gSurf_.begin());
		}
		for(int DVid=0;DVid<DVsize;DVid++){
			discreteVelocity &dv = DV_[DVid];
			hSurf_value[DVid]=const_cast<double*>(dv.hSurf_.begin());
		}
*/		

	//internal faces
    MacroSurf macsurf;
	macsurf.rhoSurf_value = rhoSurf_value;
	macsurf.rhoUsurf_value = rhoUsurf_value;
	macsurf.rhoEsurf_value = rhoEsurf_value;
	macsurf.Usurf_value = Usurf_value;
	macsurf.gSurf_value = gSurf_value1;
	macsurf.hSurf_value = hSurf_value1;
	// macsurf.xii_x = xii.x();
	// macsurf.xii_y = xii.y();
	// macsurf.xii_z = xii.z();
	macsurf.xi_value = xii_value;
	macsurf.ownersize = ownersize;
	macsurf.weight_value = weight_;
	macsurf.DVsize = DVsize;
	__real_athread_spawn((void *)slave_Func_macrosurf, &macsurf);
		
	for(int DVi=0;DVi<DVsize;DVi++){
		discreteVelocity &dv = DV_[DVi];
		vector xii = dv.xi().value();
		weight_value = dv.weight();
		forAll(rhoSurf_.BOUNDARY_FIELD_REF, patchi)
		{
			fvsPatchField<scalar>& rhoSurf_Patch = rhoSurf_.BOUNDARY_FIELD_REF[patchi];
			fvsPatchField<vector>& rhoUsurf_Patch = rhoUsurf.BOUNDARY_FIELD_REF[patchi];
			fvsPatchField<scalar>& rhoEsurf_Patch = rhoEsurf.BOUNDARY_FIELD_REF[patchi];
			rhoSurf_Patch += weight_value * dv.gSurf().boundaryField()[patchi] ;
			rhoUsurf_Patch +=  weight_value * dv.gSurf().boundaryField()[patchi] * xii;
			rhoEsurf_Patch  += 0.5 * weight_value * ( dv.gSurf().boundaryField()[patchi] * magSqr(xii) + dv.hSurf().boundaryField()[patchi]);

		}
	}
	athread_join();

	updateMacroSurf_end0 = rpcc();
#endif
	if (mpiReducer_.dvParallel()&& mpiReducer_.npd() < mpiReducer_.nproc())
	{
		mpiReducer_.reduceField(rhoSurf_);
		mpiReducer_.reduceField(rhoUsurf);
		mpiReducer_.reduceField(rhoEsurf);
	}

#ifndef MACROSURF
	//- get Prim. from Consv.
	Usurf_ = rhoUsurf / rhoSurf_;
	Tsurf_ = (rhoEsurf - 0.5 * rhoSurf_ * magSqr(Usurf_)) / ((KInner_ + 3) / 2.0 * R_ * rhoSurf_);
	// updateMacroSurf_end2 = rpcc();
	// updateTau(tauSurf_, Tsurf_, rhoSurf_);
	tauSurf_ = muRef_ * exp(omega_ * log(Tsurf_ / Tref_)) / rhoSurf_ / Tsurf_ / R_;
	surfaceScalarField qSurf_temp = 2.0 * tauSurf_ / (2.0 * tauSurf_ + 0.5 * time_.deltaT() * Pr_) ;
	//- peculiar vel.
	surfaceVectorField c = Usurf_;
	// updateMacroSurf_sta3 = rpcc();
	//-get part heat flux 
	forAll(DV_, dvi)
	{
		discreteVelocity& dv = DV_[dvi];
		c = dv.xi() - Usurf_;
		qSurf_ += 0.5 * dXiCellSize_ * dv.weight() * c * ( magSqr(c) * dv.gSurf() + dv.hSurf());
		//- stressSurf is useless as we never update cell macro by macro flux 
		//- Comment out it as it is expansive
		//stressSurf_ +=  dXiCellSize_*dv.weight()*dv.gSurf()*c*c;
	}
#else
	surfaceScalarField qSurf_temp = tauSurf_;
	// DVsize
	scalar *Tsurf_value = Tsurf_.begin();
	scalar *tauSurf_value = tauSurf_.begin();
	scalar *qSurf_temp_value = qSurf_temp.begin();
	double R_value = R_.value();
	double muRef_value = muRef_.value();
	double Tref_value = Tref_.value();
	double Pr_value = Pr_;
	double omega_value = omega_;
	int64 KInner_value = KInner_;
	scalar dt = time_.deltaTValue();
	scalar *qSurf_value = (qSurf_.begin())->v_ ;
	
	//internal faces
    MacroqSurf macqsurf;
	macqsurf.gSurf_value = gSurf_value1;
	macqsurf.hSurf_value = hSurf_value1;
	macqsurf.rhoSurf_value = rhoSurf_value;
	macqsurf.rhoUsurf_value = rhoUsurf_value;
	macqsurf.rhoEsurf_value = rhoEsurf_value;
	macqsurf.qSurf_value = qSurf_value;
	macqsurf.xi_value = xii_value;
	macqsurf.weight_value = weight_;
	macqsurf.Tsurf_value = Tsurf_value;
	macqsurf.Usurf_value = Usurf_value;
	macqsurf.tauSurf_value = tauSurf_value;
	macqsurf.qSurf_temp_value = qSurf_temp_value;
	macqsurf.R_value = R_value;
	macqsurf.muRef_value = muRef_value;
	macqsurf.Tref_value = Tref_value;
	macqsurf.Pr_value = Pr_value;
	macqsurf.omega_value = omega_value;
	macqsurf.dt = dt;
	macqsurf.KInner_value = KInner_value;
	macqsurf.ownersize = ownersize;
	macqsurf.DVsize = DVsize;
	__real_athread_spawn((void *)slave_Func_macroqsurf, &macqsurf);

	forAll(rhoSurf_.BOUNDARY_FIELD_REF, patchi)
	{
		word type = rhoSurf_.BOUNDARY_FIELD_REF[patchi].type();
		fvsPatchField<scalar>& rhoSurf_Patch = rhoSurf_.BOUNDARY_FIELD_REF[patchi];
		fvsPatchField<vector>& rhoUsurf_Patch = rhoUsurf.BOUNDARY_FIELD_REF[patchi];
		fvsPatchField<scalar>& rhoEsurf_Patch = rhoEsurf.BOUNDARY_FIELD_REF[patchi];
		fvsPatchField<vector>& Usurf_Patch = Usurf_.BOUNDARY_FIELD_REF[patchi];
		fvsPatchField<scalar>& Tsurf_Patch = Tsurf_.BOUNDARY_FIELD_REF[patchi];
		fvsPatchField<scalar>& tauSurf_Patch = tauSurf_.BOUNDARY_FIELD_REF[patchi];
		fvsPatchField<scalar>& qSurf_temp_Patch = qSurf_temp.BOUNDARY_FIELD_REF[patchi];
		Usurf_Patch = rhoUsurf_Patch / rhoSurf_Patch;
		Tsurf_Patch = (rhoEsurf_Patch - 0.5 * rhoSurf_Patch * magSqr(Usurf_Patch)) / ((KInner_value + 3) / 2.0 * R_value * rhoSurf_Patch);
		tauSurf_Patch = muRef_value * exp(omega_value * log(Tsurf_Patch / Tref_value)) / rhoSurf_Patch / Tsurf_Patch / R_value;
		qSurf_temp_Patch *= 2.0 / (2.0 * tauSurf_Patch + 0.5 * time_.deltaTValue() * Pr_value) ;
	}
	for(int DVi=0;DVi<DVsize;DVi++){
		discreteVelocity &dv = DV_[DVi];
		vector xii = dv.xi().value();
		weight_value = dv.weight();
		forAll(qSurf_.BOUNDARY_FIELD_REF, patchi)
		{
			word type = rhoSurf_.BOUNDARY_FIELD_REF[patchi].type();
			fvsPatchField<vector>& Usurf_Patch = Usurf_.BOUNDARY_FIELD_REF[patchi];
			fvsPatchField<vector>& qSurf_Patch = qSurf_.BOUNDARY_FIELD_REF[patchi];
			qSurf_Patch  += 0.5 * weight_value *(xii - Usurf_Patch)* ( magSqr(xii - Usurf_Patch) * dv.gSurf().boundaryField()[patchi]  + dv.hSurf().boundaryField()[patchi]);
	
		}
	}
	athread_join();
	
#endif
	// updateMacroSurf_end3 = rpcc();
	//- Get global heat flux, via MPI_Allreuce
	if (mpiReducer_.dvParallel() && mpiReducer_.npd() < mpiReducer_.nproc()) {
		mpiReducer_.reduceField(qSurf_);
		//mpiReducer_.reduceField(stressSurf_ );
	}
	// updateTau(tauSurf_, Tsurf_, rhoSurf_);

	//- correction for bar to original
	// qSurf_ = 2.0 * tauSurf_ / (2.0 * tauSurf_ + 0.5 * time_.deltaT() * Pr_) * qSurf_;
	qSurf_ *= qSurf_temp ;
	// qSurf_ = qSurf_temp * qSurf_;
	// updateMacroSurf_sta1 = rpcc();

	//- stress at surf is not used, as we dont't update macro in cell by macro flux at surface
	//stressSurf_ =  2.0*tauSurf_/(2.0*tauSurf_ + 0.5*time_.deltaT())*stressSurf_;

	//- heat flux at wall is specially defined. as it ignores the velocity and temperature slip
	//- NOTE: To be changed as it is part macro, but it will not affect the innner fields, so we change it later
/*#if FOAM_MAJOR <= 3
	GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField&
		rhoBCs = rhoVol_.boundaryField();
#else
	GeometricField<scalar, fvPatchField, volMesh>::Boundary&
		rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
#endif
	qWall_ = dimensionedVector("0", qWall_.dimensions(), vector(0, 0, 0));
	stressWall_ = dimensionedTensor
	(
		"0",
		stressWall_.dimensions(),
		pTraits<tensor>::zero
	);
	forAll(rhoBCs, patchi)
	{
		if (rhoBCs[patchi].type() == "calculatedMaxwell")
		{
			fvPatchField<vector>& qPatch = qWall_.BOUNDARY_FIELD_REF[patchi];
			fvPatchField<vector>& Upatch = Uvol_.BOUNDARY_FIELD_REF[patchi];
			fvPatchField<tensor>& stressPatch = stressWall_.BOUNDARY_FIELD_REF[patchi];
			//- tau at surface use the tau at slip temperature as it is.
			fvsPatchField<scalar>& tauPatch = tauSurf_.BOUNDARY_FIELD_REF[patchi];
			forAll(qPatch, facei)
			{
				forAll(DV_, dvi)
				{
					scalar dXiCellSize = dXiCellSize_.value();
					discreteVelocity& dv = DV_[dvi];
					vector xi = dv.xi().value();
					vector c = xi - Upatch[facei];
					qPatch[facei] += 0.5 * dXiCellSize * dv.weight() * c  //sometimes wall moves, then c != \xi
						* (
							magSqr(c) * dv.gSurf().boundaryField()[patchi][facei]
							+ dv.hSurf().boundaryField()[patchi][facei]
							);
					// stressPatch[facei] +=
					// 	dXiCellSize * dv.weight() * dv.gSurf().boundaryField()[patchi][facei] * xi * xi;
				}
				qPatch[facei] = 2.0 * tauPatch[facei] / (2.0 * tauPatch[facei] + 0.5 * time_.deltaT().value() * Pr_) * qPatch[facei];
				//stressPatch[facei] = 2.0 * tauPatch[facei] / (2.0 * tauPatch[facei] + 0.5 * time_.deltaT().value()) * stressPatch[facei];
			}
			if (mpiReducer_.dvParallel() && mpiReducer_.npd() < mpiReducer_.nproc())
			{
				mpiReducer_.reduceField(qPatch);
				// mpiReducer_.reduceField(stressPatch);
			}
		}
	}
*/

	time_end = rpcc();

	if(mpiReducer_.rank() == 0){
		Info << "updateMacroSurf time =" <<((double)(time_end-time_sta)*1000/CLOCKRATE) << "ms" << nl << endl;
	}
}

void Foam::fvDVM::updateGHsurf()
{
	unsigned long time_sta = 0;
	unsigned long time_end = 0;
	unsigned long time_sta1 = 0;
	unsigned long time_end1 = 0;
	time_sta = rpcc();
#ifndef GHSURF
	forAll(DV_, DVid)
		DV_[DVid].updateGHsurf();
#else
	int64 DVsize = DV_.size();
	label D = mesh_.nSolutionD();
	label K = KInner_;
	dimensionedScalar vUnit("vUnit", dimLength/dimTime, 1);
	const surfaceScalarField &rho = rhoSurf_;
	const surfaceVectorField &U = Usurf_;
	const surfaceScalarField &T = Tsurf_;
	const surfaceVectorField &q = qSurf_;
	double R_value = R_.value();
	double Pr_value = Pr_;
	double vUnit_value = vUnit.value();
	scalar *U_value = const_cast<double*>((U.begin())->v_ );
	scalar *q_value = const_cast<double*>((q.begin())->v_) ;
	scalar *rho_value = const_cast<double*>(rho.begin());
	scalar *T_value = const_cast<double*>(T.begin());
	scalar h = 0.5*time_.deltaTValue();
	scalar *tauSurf_value = const_cast<double*>(tauSurf_.begin());

	scalar **gSurf_value = (scalar**)malloc(sizeof(scalar*) * DVsize );
	scalar **hSurf_value = (scalar**)malloc(sizeof(scalar*) * DVsize );
	for(int DVid=0;DVid<DVsize;DVid++){
		discreteVelocity &dv = DV_[DVid];
		gSurf_value[DVid]=const_cast<double*>(dv.gSurf_.begin());
		hSurf_value[DVid]=const_cast<double*>(dv.hSurf_.begin());
	}	
	int64 size = DV_[0].gSurf_.size();
	GHsurf ghsurf;
	ghsurf.U_value = U_value;
	ghsurf.q_value = q_value;
	ghsurf.T_value = T_value;
	ghsurf.rho_value = rho_value;
	ghsurf.tauSurf_value = tauSurf_value;
	ghsurf.gSurf_value = gSurf_value;
	ghsurf.hSurf_value = hSurf_value;
	ghsurf.h = h;
	ghsurf.gSurf_size = size;
	ghsurf.xi_value = xii_value;
	ghsurf.R_value = R_value;
	ghsurf.Pr_value = Pr_value;
	ghsurf.vUnit_value = vUnit_value;
	ghsurf.D = D;//label
	ghsurf.K = K;//label
	ghsurf.DVsize = DVsize;//label
	__real_athread_spawn((void *)slave_Func_ghsurf, &ghsurf);
	surfaceScalarField relaxFactor = tauSurf_;
	
	for(int DVi=0;DVi<DVsize;DVi++){
		// scalar cSqrByRT=0.0, cqBy5pRT=0.0, gEqBGK=0.0, relaxFactor=0.0;
		surfaceScalarField cSqrByRT = rho;
		surfaceScalarField cqBy5pRT = rho;
		surfaceScalarField gEqBGK  = rho;
		
		vector xii = DV_[DVi].xi_.value();
		forAll(DV_[DVi].gSurf_.BOUNDARY_FIELD_REF, patchi)
		{
			word type = DV_[DVi].gSurf_.BOUNDARY_FIELD_REF[patchi].type();
			fvsPatchField<scalar>& gSurfPatch = DV_[DVi].gSurf_.BOUNDARY_FIELD_REF[patchi];
			fvsPatchField<scalar>& hSurfPatch = DV_[DVi].hSurf_.BOUNDARY_FIELD_REF[patchi];
			fvsPatchField<scalar>& relaxFactorPatch = relaxFactor.BOUNDARY_FIELD_REF[patchi];
			fvsPatchField<scalar>& cSqrByRTPatch = cSqrByRT.BOUNDARY_FIELD_REF[patchi];
			fvsPatchField<scalar>& cqBy5pRTPatch = cqBy5pRT.BOUNDARY_FIELD_REF[patchi];
			fvsPatchField<scalar>& gEqBGKPatch = gEqBGK.BOUNDARY_FIELD_REF[patchi];
			const fvsPatchField<scalar>& rhoSurf_Patch = rhoSurf_.BOUNDARY_FIELD_REF[patchi];
			const fvsPatchField<vector>& Usurf_Patch = Usurf_.BOUNDARY_FIELD_REF[patchi];
			const fvsPatchField<scalar>& Tsurf_Patch = Tsurf_.BOUNDARY_FIELD_REF[patchi];
			const fvsPatchField<scalar>& tauSurf_Patch = tauSurf_.BOUNDARY_FIELD_REF[patchi];
			const fvsPatchField<vector>& qSurf_Patch = qSurf_.BOUNDARY_FIELD_REF[patchi];

			if (type == "processor"){
				relaxFactorPatch = h/(2*tauSurf_Patch + h);
				gSurfPatch = (1.0 - relaxFactorPatch)*gSurfPatch;
				hSurfPatch = (1.0 - relaxFactorPatch)*hSurfPatch;
				cSqrByRTPatch = magSqr(Usurf_Patch - xii)/(R_value * Tsurf_Patch);
				cqBy5pRTPatch = ((xii - Usurf_Patch)&qSurf_Patch)/(5.0*rhoSurf_Patch*R_value * Tsurf_Patch*R_value * Tsurf_Patch);		
				// gEqBGK  = (relaxFactor*rhoSurf_Patch[facei])/pow(sqrt(2.0*pi*R_value * Tsurf_Patch[facei]),D)*exp(-cSqrByRT/2.0)/pow(vUnit_value, 3-D);
				gEqBGKPatch  = (relaxFactorPatch*rhoSurf_Patch)/pow(sqrt(2.0*pi*R_value * Tsurf_Patch),D)*exp(-cSqrByRTPatch/2.0);
				gSurfPatch += ( 1.0 + (1.0 - Pr_value)*cqBy5pRTPatch*(cSqrByRTPatch - D - 2.0) )*gEqBGKPatch;
				hSurfPatch += ( (K + 3.0 - D) + (1.0 - Pr_value)*cqBy5pRTPatch*((cSqrByRTPatch - D)*(K + 3.0 - D) - K - K) )*gEqBGKPatch*R_value*Tsurf_Patch;
				
			}
		}
	}
	athread_join();
	free(gSurf_value);
	free(hSurf_value);
#endif
	time_end = rpcc();
	if(mpiReducer_.rank() == 0){
		std::cout << "updateGHsurf time =" <<((double)(time_end-time_sta)*1000/CLOCKRATE) << "ms" << nl << endl;//测整个循环
	}
}

void Foam::fvDVM::updateGHtildeVol()
{
	unsigned long time_sta = 0;
	unsigned long time_end = 0;
	time_sta = rpcc();
	forAll(DV_, DVid){
		DV_[DVid].updateGHtildeVol();
	}
	time_end = rpcc();
	if(mpiReducer_.rank() == 0){
		Info << "updateGHtildeVol time =" <<((double)(time_end-time_sta)*1000/CLOCKRATE) << "ms" << nl << endl;//测整个循环
	}
}

void Foam::fvDVM::updateMacroVol()
{
	unsigned long time_sta = 0;
	unsigned long time_end = 0;
	time_sta = rpcc();
	//- Old macros, used only if we update using macro fluxes.
	// volVectorField rhoUvol = rhoVol_ * Uvol_;
	// volScalarField rhoEvol = rhoVol_ * (0.5*magSqr(Uvol_) + (KInner_ + 3) / 2.0 * R_ * Tvol_);
	qVol_ = dimensionedVector("0", qVol_.dimensions(), vector(0, 0, 0));
	unsigned long time_0 = 0;
	unsigned long time_1 = 0;
	unsigned long time_2 = 0;

	// if (macroFlux_ == "no") // update cell macro by moment from DF
	// {
	// 	//- init to zeros
	// 	rhoVol_ = dimensionedScalar("0", rhoVol_.dimensions(), 0);
	// 	rhoUvol = dimensionedVector("0", rhoUvol.dimensions(), vector(0, 0, 0));
	// 	rhoEvol = dimensionedScalar("0", rhoEvol.dimensions(), 0);

	// 	//- get part macro
	// 	forAll(DV_, dvi)
	// 	{
	// 		discreteVelocity &dv = DV_[dvi];
	// 		rhoVol_ += dXiCellSize_ * dv.weight() * dv.gTildeVol();
	// 		rhoUvol += dXiCellSize_ * dv.weight() * dv.gTildeVol() * dv.xi();
	// 		rhoEvol += 0.5 * dXiCellSize_ * dv.weight() * (magSqr(dv.xi()) * dv.gTildeVol() + dv.hTildeVol());
	// 	}
	// 	//- get global macro via MPI_Allreduce
	// 	if (mpiReducer_.dvParallel() && mpiReducer_.npd() < mpiReducer_.nproc())
	// 	{
	// 		mpiReducer_.reduceField(rhoVol_);
	// 		mpiReducer_.reduceField(rhoUvol);
	// 		mpiReducer_.reduceField(rhoEvol);
	// 	}
		
	// }
	// else // update by macro flux
	// {

#ifndef MACROVOL
	const labelUList &owner = mesh_.owner();
	const labelUList &neighbour = mesh_.neighbour();
	const vectorField Sf = mesh_.Sf();
	const scalarField V = mesh_.V();
	const scalar dt = time_.deltaTValue();

	 //init flux to zero   
	rhoflux_ = dimensionedScalar("0", rhoVol_.dimensions(), 0);
	volVectorField rhouflux_= rhoflux_*Uvol_;  
	volScalarField rhoeflux_ = rhoflux_ * magSqr(Uvol_);
	unsigned long time_sta0 = 0;
	unsigned long time_sta1 = 0;

	forAll(DV_, dvi)
	{   
		//temp var 
     	rho_ = dimensionedScalar("0", rhoSurf_.dimensions(), 0);
		surfaceVectorField u_= rho_*Usurf_;
		surfaceScalarField e_ = rho_ * magSqr(Usurf_);
		discreteVelocity &dv = DV_[dvi];
		vector xii = dv.xi().value();
		rho_ = dXiCellSize_ * dv.weight() * dv.gSurf(); 
		u_ = dXiCellSize_ * dv.weight() * dv.gSurf() * dv.xi();
		e_  = 0.5*dXiCellSize_*dv.weight()*(dv.gSurf()*magSqr(dv.xi()) + dv.hSurf());

		//internal faces
		forAll(owner, facei)
		{
			const label own = owner[facei];
			const label nei = neighbour[facei];

			rhoflux_[own] -= ((xii & Sf[facei]) * rho_[facei] * dt / V[own]);
			rhoflux_[nei] += ((xii & Sf[facei]) * rho_[facei] * dt / V[nei]);
			rhouflux_[own] -= ((xii & Sf[facei]) * u_[facei] * dt / V[own]);
			rhouflux_[nei] += ((xii & Sf[facei]) * u_[facei] * dt / V[nei]);
			rhoeflux_[own] -= ((xii & Sf[facei]) * e_[facei] * dt / V[own]);
			rhoeflux_[nei] += ((xii & Sf[facei]) * e_[facei] * dt / V[nei]);
		}
		forAll(rhoSurf_.boundaryField(), patchi)
		{   
			const fvsPatchField<vector> &SfPatch = mesh_.Sf().boundaryField()[patchi];
			const labelUList &pOwner = mesh_.boundary()[patchi].faceCells();

			forAll(pOwner, pFacei)
			{
				const label own = pOwner[pFacei];
				rhoflux_[own] -= ((xii & SfPatch[pFacei]) * rho_.boundaryField()[patchi][pFacei] * dt / V[own]);
				rhouflux_[own] -= ((xii & SfPatch[pFacei]) * u_.boundaryField()[patchi][pFacei] * dt / V[own]);
				rhoeflux_[own] -= ((xii & SfPatch[pFacei]) * e_.boundaryField()[patchi][pFacei]  * dt / V[own]);
			}
		}

	}

#else
	int64 DVsize = DV_.size();
	const int64 *owner = mesh_.owner().begin();
	int64 ownersize = mesh_.owner().size();
	const int64 *neighbour = mesh_.neighbour().begin();
	scalar *Sf = const_cast<double*>((mesh_.Sf().begin())->v_) ;

	const double *V = mesh_.V().begin();
	const double dt = time_.deltaTValue();
	 //init flux to zero
	rhoflux_ = dimensionedScalar("0", rhoVol_.dimensions(), 0);
	volVectorField rhouflux_= rhoflux_*Uvol_;
	volScalarField rhoeflux_ = rhoflux_ * magSqr(Uvol_);

	double *rhoflux_value = rhoflux_.begin();
	double *rhoeflux_value = rhoeflux_.begin();
	double *rhouflux_value = (rhouflux_.begin())->v_ ;

	double *rhouflux_value_own = (double*)malloc(sizeof(double) * ownersize*3);
	double *rhouflux_value_nei = (double*)malloc(sizeof(double) * ownersize*3);

	scalar *rhoflux_value_own = (scalar*)malloc(sizeof(scalar) * ownersize);
	scalar *rhoflux_value_nei = (scalar*)malloc(sizeof(scalar) * ownersize);
	scalar *rhoeflux_value_own = (scalar*)malloc(sizeof(scalar) * ownersize);
	scalar *rhoeflux_value_nei = (scalar*)malloc(sizeof(scalar) * ownersize);

	rho_ = dimensionedScalar("0", rhoSurf_.dimensions(), 0);
	surfaceVectorField u_= rho_*Usurf_;
	surfaceScalarField e_ = rho_ * magSqr(Usurf_);

	scalar **gSurf_value = (scalar**)malloc(sizeof(scalar*) * DVsize );
	scalar **hSurf_value = (scalar**)malloc(sizeof(scalar*) * DVsize );

	for(int DVid=0;DVid<DVsize;DVid++){
		discreteVelocity &dv = DV_[DVid];
		gSurf_value[DVid]=const_cast<double*>(dv.gSurf_.begin());

	}
	for(int DVid=0;DVid<DVsize;DVid++){
		discreteVelocity &dv = DV_[DVid];
		hSurf_value[DVid]=const_cast<double*>(dv.hSurf_.begin());
	}

	double weight_value = 0;

	//internal faces
    MacroVol master;
	master.Sf = Sf;
	master.rhoflux_value_own = rhoflux_value_own;
	master.rhoflux_value_nei = rhoflux_value_nei;
	master.rhouflux_value_own = rhouflux_value_own;
	master.rhouflux_value_nei = rhouflux_value_nei;
	master.rhoeflux_value_own = rhoeflux_value_own;
	master.rhoeflux_value_nei = rhoeflux_value_nei;
	master.V_own = V_own1;
	master.V_nei = V_nei1;
	master.gSurf_value = gSurf_value;
	master.hSurf_value = hSurf_value;
	master.dt = dt;
	master.xi_value = xii_value;
	master.weight_value = weight_;
	master.ownersize = ownersize;
	master.DVsize = DVsize;

	__real_athread_spawn((void *)slave_Func_macrovol, &master);

			
	forAll(DV_, dvi)
	{

		discreteVelocity &dv = DV_[dvi];
		vector xii = dv.xi().value();
		weight_value = dv.weight();
		forAll(rhoSurf_.boundaryField(), patchi)
		{   	
			const fvsPatchField<vector> &SfPatch = mesh_.Sf().boundaryField()[patchi];
			const labelUList &pOwner = mesh_.boundary()[patchi].faceCells();
			word type = rhoSurf_.BOUNDARY_FIELD_REF[patchi].type();
			const fvsPatchField<scalar>& gSurfPatch = dv.gSurf_.BOUNDARY_FIELD_REF[patchi];
			const fvsPatchField<scalar>& hSurfPatch = dv.hSurf_.BOUNDARY_FIELD_REF[patchi];
			fvsPatchField<scalar>& rho_Patch = rho_.BOUNDARY_FIELD_REF[patchi];
			fvsPatchField<vector>& u_Patch = u_.BOUNDARY_FIELD_REF[patchi];
			fvsPatchField<scalar>& e_Patch = e_.BOUNDARY_FIELD_REF[patchi];
			rho_Patch = weight_value * gSurfPatch;
			u_Patch =  weight_value * gSurfPatch * xii;
			e_Patch  = 0.5 * weight_value * ( gSurfPatch * magSqr(xii) + hSurfPatch);	
			forAll(pOwner, pFacei)
			{
				const label own = pOwner[pFacei];
				rhoflux_[own] -= ((xii & SfPatch[pFacei]) * rho_.boundaryField()[patchi][pFacei] * dt / V[own]);
				rhouflux_[own] -= ((xii & SfPatch[pFacei]) * u_.boundaryField()[patchi][pFacei] * dt / V[own]);
				rhoeflux_[own] -= ((xii & SfPatch[pFacei]) * e_.boundaryField()[patchi][pFacei]  * dt / V[own]);
			}
			
		}
	}

	athread_join();
	for(int facei=0;facei<ownersize;facei++)
	{
		const label own = owner[facei];
		const label nei = neighbour[facei];
		rhoflux_value[own] += rhoflux_value_own[facei];
		rhoflux_value[nei] += rhoflux_value_nei[facei];
	}
	for(int facei=0;facei<ownersize;facei++)
	{
		const label own = owner[facei];
		const label nei = neighbour[facei];
		rhouflux_[own].v_[0] += rhouflux_value_own[facei*3];
		rhouflux_[own].v_[1] += rhouflux_value_own[facei*3+1];
		rhouflux_[own].v_[2] += rhouflux_value_own[facei*3+2];
		rhouflux_[nei].v_[0] += rhouflux_value_nei[facei*3];
		rhouflux_[nei].v_[1] += rhouflux_value_nei[facei*3+1];
		rhouflux_[nei].v_[2] += rhouflux_value_nei[facei*3+2];
	}
    for(int facei=0;facei<ownersize;facei++)
	{
		const label own = owner[facei];
		const label nei = neighbour[facei];
		rhoeflux_value[own] += rhoeflux_value_own[facei];
		rhoeflux_value[nei] += rhoeflux_value_nei[facei];
	}

	free(gSurf_value);
	free(hSurf_value);
	free(rhoflux_value_own);
	free(rhoflux_value_nei);
	free(rhoeflux_value_own);
	free(rhoeflux_value_nei);
	free(rhouflux_value_own);
	free(rhouflux_value_nei);
#endif
	
	if (mpiReducer_.dvParallel()&& mpiReducer_.npd() < mpiReducer_.nproc())
	{
		mpiReducer_.reduceField(rhoflux_);
		mpiReducer_.reduceField(rhouflux_);
		mpiReducer_.reduceField(rhoeflux_);
	}

	// //- get Prim. from Consv.
	// Uvol_ = rhoUvol / rhoVol_;
	// Tvol_ = (rhoEvol - 0.5 * rhoVol_ * magSqr(Uvol_)) / ((KInner_ + 3) / 2.0 * R_ * rhoVol_);
	//不需要
	//- Correct the macro field boundary conition
	// Uvol_.correctBoundaryConditions();
    //     //Info<<Uvol_.boundaryField();
	// Tvol_.correctBoundaryConditions();
	//- Note for maxwell wall, the operation here update 
	//- the boundary rho field but it's meaningless.

	//- The vol macro field's boundary field is meanless!
	//rhoVol_.correctBoundaryConditions(); 

	//- update tau
	// updateTau(tauVol_, Tvol_, rhoVol_);
	//- peculiar vel.
	// volVectorField c = Uvol_;
	// Utime_0 = rpcc();
	// forAll(DV_, dvi)
	// {
	// 	discreteVelocity& dv = DV_[dvi];
	// 	c = dv.xi() - Usurf_;
	// 	qSurf_ += 0.5 * dXiCellSize_ * dv.weight() * c * ( magSqr(c) * dv.gSurf() + dv.hSurf());
	// 	//- stressSurf is useless as we never update cell macro by macro flux 
	// 	//- Comment out it as it is expansive
	// 	//stressSurf_ +=  dXiCellSize_*dv.weight()*dv.gSurf()*c*c;
	// }

#ifndef MACROVOL

	//- Old macros, used only if we update using macro fluxes.
	volVectorField rhoUvol = rhoVol_ * Uvol_;
	volScalarField rhoEvol = rhoVol_ * (0.5*magSqr(Uvol_) + (KInner_ + 3) / 2.0 * R_ * Tvol_);
	
        rhoVol_ += rhoflux_;
		rhoUvol += rhouflux_;
		rhoEvol += rhoeflux_;
	Uvol_ = rhoUvol / rhoVol_;
	Tvol_ = (rhoEvol - 0.5 * rhoVol_ * magSqr(Uvol_)) / ((KInner_ + 3) / 2.0 * R_ * rhoVol_);
	tauVol_ = muRef_ * exp(omega_ * log(Tvol_ / Tref_)) / rhoVol_ / Tvol_ / R_;//// updateTau(tauVol_, Tvol_, rhoVol_);
	volScalarField qVol_temp = 2.0 * tauVol_ / (2.0 * tauVol_ + time_.deltaT() * Pr_);
	//- peculiar vel.
	volVectorField c = Uvol_;
	//-get part heat flux
	forAll(DV_, dvi)
	{
		discreteVelocity& dv = DV_[dvi];
		c = dv.xi() - Uvol_;
		qVol_ += 0.5 * dXiCellSize_ * dv.weight() * c * ( magSqr(c) * dv.gTildeVol() + dv.hTildeVol());
	}
#else
	int64 qVol_size = qVol_.size();
	volScalarField qVol_temp = tauVol_;
	// DVsize
	scalar *Tvol_value = Tvol_.begin();
	scalar *Uvol_value = (Uvol_.begin())->v_ ;
	scalar *tauVol_value = tauVol_.begin();
	scalar *qVol_temp_value = qVol_temp.begin();
	scalar *rhoVol_value = rhoVol_.begin();
	double R_value = R_.value();
	double muRef_value = muRef_.value();
	double Tref_value = Tref_.value();
	double Pr_value = Pr_;
	double omega_value = omega_;
	int64 KInner_value = KInner_;
	scalar **gTildeVol_value = (scalar**)malloc(sizeof(scalar*) * DVsize );
	scalar **hTildeVol_value = (scalar**)malloc(sizeof(scalar*) * DVsize );
	scalar *qVol_value = (qVol_.begin())->v_ ;
		
	for(int DVid=0;DVid<DVsize;DVid++){
		discreteVelocity &dv = DV_[DVid];
		gTildeVol_value[DVid]=const_cast<double*>(dv.gTildeVol().begin());
	}
	for(int DVid=0;DVid<DVsize;DVid++){
		discreteVelocity &dv = DV_[DVid];
		hTildeVol_value[DVid]=const_cast<double*>(dv.hTildeVol().begin());
	}

	//internal faces
    MacroqVol macroqvol;
	macroqvol.gTildeVol_value = gTildeVol_value;
	macroqvol.hTildeVol_value = hTildeVol_value;
	macroqvol.rhoVol_value = rhoVol_value;
	macroqvol.rhoflux_value = rhoflux_value;
	macroqvol.rhoeflux_value = rhoeflux_value;
	macroqvol.rhouflux_value = rhouflux_value;
	macroqvol.qVol_value = qVol_value;
	macroqvol.xi_value = xii_value;
	macroqvol.weight_value = weight_;
	macroqvol.Tvol_value = Tvol_value;
	macroqvol.Uvol_value = Uvol_value;
	macroqvol.tauVol_value = tauVol_value;
	macroqvol.qVol_temp_value = qVol_temp_value;
	macroqvol.R_value = R_value;
	macroqvol.muRef_value = muRef_value;
	macroqvol.Tref_value = Tref_value;
	macroqvol.Pr_value = Pr_value;
	macroqvol.omega_value = omega_value;
	macroqvol.dt = dt;
	macroqvol.KInner_value = KInner_value;
	macroqvol.qVol_size = qVol_size;
	macroqvol.DVsize = DVsize;
	__real_athread_spawn((void *)slave_Func_macroqvol, &macroqvol);
	athread_join();
	//internal faces
	free(gTildeVol_value);
	free(hTildeVol_value);
	
#endif
	// Utime_1 = rpcc();
	//- get global heat flux via MPI_Allreduce
	if (mpiReducer_.dvParallel()&& mpiReducer_.npd() < mpiReducer_.nproc())
		mpiReducer_.reduceField(qVol_);
	// updateTau(tauVol_, Tvol_, rhoVol_);
	//- correction for bar to original
	// qVol_ = 2.0 * tauVol_ / (2.0 * tauVol_ + time_.deltaT() * Pr_) * qVol_;
	qVol_ *= qVol_temp;
	time_end = rpcc();  
	if(mpiReducer_.rank() == 0){
		Info << "updateMacroVol total time =" <<((double)(time_end-time_sta)*1000/CLOCKRATE) << "ms" << nl << endl;
	}
}

void Foam::fvDVM::updatePressureInOutBC()
{
	unsigned long time_sta = 0;
	unsigned long time_end = 0;
	time_sta = rpcc();
	// for pressureIn and pressureOut BC, the boundary value of Uvol(in/out) and Tvol(in/out) should be updated here!
	// boundary faces
#if FOAM_MAJOR <= 3
	GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField&
		rhoBCs = rhoVol_.boundaryField();
#else
	GeometricField<scalar, fvPatchField, volMesh>::Boundary&
		rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
#endif
	forAll(rhoBCs, patchi)
	{
		if (rhoBCs[patchi].type() == "pressureIn")
		{
			const fvsPatchField<vector>& SfPatch = mesh_.Sf().boundaryField()[patchi];
			const fvsPatchField<scalar>& magSfPatch = mesh_.magSf().boundaryField()[patchi];
			pressureInFvPatchField<scalar>& rhoPatch =
				refCast<pressureInFvPatchField<scalar> >(rhoBCs[patchi]);
			fvPatchField<vector>& Upatch = Uvol_.BOUNDARY_FIELD_REF[patchi];
			const fvPatchField<scalar>& Tpatch = Tvol_.boundaryField()[patchi];
			const scalar pressureIn = rhoPatch.pressureIn();
			// now changed rho and U patch
			const labelUList& pOwner = mesh_.boundary()[patchi].faceCells();
			forAll(rhoPatch, facei)
			{
				const scalar  Tin = Tpatch[facei];
				// change density
				rhoPatch[facei] = pressureIn / R_.value() / Tin; // Accturally not changed at all :p

				// inner boundary cell data state data
				label own = pOwner[facei];
				vector Ui = Uvol_[own];
				scalar Ti = Tvol_[own];
				scalar rhoi = rhoVol_[own];
				scalar ai = sqrt(R_.value() * Ti * (KInner_ + 5) / (KInner_ + 3)); // sos

				// change normal velocity component based on the characteristics
				vector norm = SfPatch[facei] / magSfPatch[facei]; // boundary face normal vector
				scalar Un = Ui & norm; // normal component
				scalar UnIn = Un + (pressureIn - rhoi * R_.value() * Ti) / rhoi / ai; // change normal component
				Upatch[facei] = UnIn * norm + (Ui - Un * norm); // tangential component not changed.
			}
		}
		else if (rhoBCs[patchi].type() == "pressureOut")
		{
			const fvsPatchField<vector>& SfPatch = mesh_.Sf().boundaryField()[patchi];
			const fvsPatchField<scalar>& magSfPatch = mesh_.magSf().boundaryField()[patchi];
			pressureOutFvPatchField<scalar>& rhoPatch =
				refCast<pressureOutFvPatchField<scalar> >(rhoBCs[patchi]);
			fvPatchField<vector>& Upatch = Uvol_.BOUNDARY_FIELD_REF[patchi];
			fvPatchField<scalar>& Tpatch = Tvol_.BOUNDARY_FIELD_REF[patchi];
			const scalar pressureOut = rhoPatch.pressureOut();
			// now changed rho and U patch
			const labelUList& pOwner = mesh_.boundary()[patchi].faceCells();
			forAll(rhoPatch, facei)
			{
				// inner cell data state data
				label own = pOwner[facei];
				vector Ui = Uvol_[own];
				scalar Ti = Tvol_[own];
				scalar rhoi = rhoVol_[own];
				scalar ai = sqrt(R_.value() * Ti * (KInner_ + 5) / (KInner_ + 3)); // sos

				// change outlet density
				rhoPatch[facei] = rhoi + (pressureOut - rhoi * R_.value() * Ti) / ai / ai; // Accturally not changed at all :p
				Tpatch[facei] = pressureOut / (R_.value() * rhoi);

				// change normal velocity component based on the characteristics
				vector norm = SfPatch[facei] / magSfPatch[facei]; // boundary face normal vector
				scalar Un = Ui & norm; // normal component
				scalar UnIn = Un + (rhoi * R_.value() * Ti - pressureOut) / rhoi / ai; // change normal component
				Upatch[facei] = UnIn * norm + (Ui - Un * norm); // tangential component not changed.
			}
		}
	}
	time_end = rpcc();        
}

template<template<class> class PatchType, class GeoMesh>
void Foam::fvDVM::updateTau
(
	GeometricField<scalar, PatchType, GeoMesh>& tau,
	const GeometricField<scalar, PatchType, GeoMesh>& T,
	const GeometricField<scalar, PatchType, GeoMesh>& rho
)
{
	tau = muRef_ * exp(omega_ * log(T / Tref_)) / rho / T / R_;
}


void Foam::fvDVM::writeDFonCell(label cellId)
{
	std::ostringstream convert;
	convert << cellId;
	scalarIOList df
	(
		IOobject
		(
			"DF" + convert.str(),
			"0",
			mesh_,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		)
	);
	//set size of df
	df.setSize(nXi_);

	scalarList dfPart(DV_.size());
	//put cellId's DF to dfPart
	forAll(dfPart, dfi)
		dfPart[dfi] = DV_[dfi].gTildeVol()[cellId];

	label nproc = mpiReducer_.nproc();
	//gather
	//tmp list for recv
	scalarList dfRcv(nXi_);

	//Compose displc and recvc
	//labelField recvc(nproc);
	//labelField displ(nproc);
	Field<int> recvc(nproc);
        Field<int> displ(nproc);
	label chunck = nXi_ / nproc;
	label left = nXi_ % nproc;
	forAll(recvc, i)
	{
		recvc[i] = chunck + (i < left);
		if (i <= left)
			displ[i] = i * (chunck + 1); // (i<=nXi_%nproc)
		else
			displ[i] = left * (chunck + 1) + (i - left) * (chunck);
	}
	MPI_Gatherv(dfPart.data(), dfPart.size(), MPI_DOUBLE,
		dfRcv.data(), recvc.data(), displ.data(),
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//reposition
	if (mpiReducer_.rank() == 0)
	{
		forAll(df, i)
		{
			label p = i % nproc;
			label ldi = i / nproc;
			df[i] = dfRcv[displ[p] + ldi];
		}
		df.write();
	}
}

void Foam::fvDVM::writeDFonCells()
{
	if (time_.outputTime())
		forAll(DFwriteCellList_, i)
		writeDFonCell(i);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvDVM::fvDVM
(
	volScalarField& rho,
	volVectorField& U,
	volScalarField& T,
	int* argc,
	char*** argv,
	Foam::argList& args
)
	:
	IOdictionary
	(
		IOobject
		(
			"DVMProperties",
			T.time().constant(),
			T.mesh(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	),
	mesh_(rho.mesh()),
	time_(rho.time()),
	rhoVol_(rho),
	Uvol_(U),
	Tvol_(T),
	args_(args),
	fvDVMparas_(subOrEmptyDict("fvDVMparas")),
	gasProperties_(subOrEmptyDict("gasProperties")),
	nXiPerDim_(readLabel(fvDVMparas_.lookup("nDV"))),
	xiMax_(fvDVMparas_.lookup("xiMax")),
	xiMin_(fvDVMparas_.lookup("xiMin")),
	dXi_((xiMax_ - xiMin_) / (nXiPerDim_ - 1)),
	dXiCellSize_
	(
		"dXiCellSize",
		pow(dimLength / dimTime, 3),
		scalar(1.0)
	),
	macroFlux_(fvDVMparas_.lookupOrDefault("macroFlux", word("no"))),
    //res_(fvDVMparas_.lookupOrDefault("res", 1.0e-20)),
	//checkSteps_(fvDVMparas_.lookupOrDefault("checkSteps", 100)),
	R_(gasProperties_.lookup("R")),
	omega_(readScalar(gasProperties_.lookup("omega"))),
	Tref_(gasProperties_.lookup("Tref")),
	muRef_(gasProperties_.lookup("muRef")),
	Pr_(readScalar(gasProperties_.lookup("Pr"))),
	KInner_((gasProperties_.lookupOrDefault("KInner", 0))),
	mpiReducer_(args, argc, argv), // args comes from setRootCase.H in dugksFoam.C;
	DV_(0),
	rhoSurf_
	(
		IOobject
		(
			"rhoSurf",
			mesh_.time().timeName(),
			mesh_,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh_,
		dimensionedScalar("0", rho.dimensions(), 0)
	),
	rhoflux_
	(
		IOobject
		(
			"rhoflux",
			mesh_.time().timeName(),
			mesh_,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh_,
		dimensionedScalar("0", rho.dimensions(), 0)
	),
    rho_
	(
		IOobject
		(
			"rhotemp",
			mesh_.time().timeName(),
			mesh_,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh_,
		dimensionedScalar("0", rho.dimensions(), 0)
	),
	Tsurf_
	(
		IOobject
		(
			"Tsurf",
			mesh_.time().timeName(),
			mesh_,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh_,
		dimensionedScalar("0", T.dimensions(), 0)
	),
	Usurf_
	(
		IOobject
		(
			"Usurf",
			mesh_.time().timeName(),
			mesh_,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh_,
		dimensionedVector("0", U.dimensions(), vector(0, 0, 0))
	),
	qSurf_
	(
		IOobject
		(
			"qSurf",
			mesh_.time().timeName(),
			mesh_,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh_,
		dimensionedVector("0", dimMass / pow(dimTime, 3), vector(0, 0, 0))
	),
	stressSurf_
	(
		IOobject
		(
			"stressSurf",
			mesh_.time().timeName(),
			mesh_,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh_,
		dimensionedTensor("0", dimensionSet(1, -1, -2, 0, 0, 0, 0), pTraits<tensor>::zero)
	),
	qVol_
	(
		IOobject
		(
			"q",
			mesh_.time().timeName(),
			mesh_,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh_,
		dimensionedVector("0", dimMass / pow(dimTime, 3), vector(0, 0, 0))
	),
	tauVol_
	(
		IOobject
		(
			"tauVol",
			mesh_.time().timeName(),
			mesh_,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh_,
		dimensionedScalar("0", dimTime, 0)
	),
	tauSurf_
	(
		IOobject
		(
			"tauSurf",
			mesh_.time().timeName(),
			mesh_,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh_,
		dimensionedScalar("0", dimTime, 0)
	),
	qWall_
	(
		IOobject
		(
			"qWall",
			mesh_.time().timeName(),
			mesh_,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh_,
		dimensionedVector("0", dimMass / pow(dimTime, 3), vector(0, 0, 0))
	),
	stressWall_
	(
		IOobject
		(
			"stressWall",
			mesh_.time().timeName(),
			mesh_,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh_,
		dimensionedTensor("0", dimensionSet(1, -1, -2, 0, 0, 0, 0), pTraits<tensor>::zero)
	)
{
	DFwriteCellList_ = lookupOrDefault<labelList>("DFwriteCellList", labelList()),
	initialiseDV();
	setCalculatedMaxwellRhoBC();
	setSymmetryModRhoBC();
	// set initial rho in pressureIn/Out BC
	updatePressureInOutBC();
	updateTau(tauVol_, Tvol_, rhoVol_); //calculate the tau at cell when init
	Usurf_ = fvc::interpolate(Uvol_, "linear"); // for first time Dt calculation.
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvDVM::~fvDVM()
{
	free(xii_value);
	free(weight_);

	free(gBarPvol_value1);
	free(hBarPvol_value1);
	free(gTildeVol_value1);
	free(hTildeVol_value1);

	free(gSurf_value1);
	free(hSurf_value1);

	free(V_own1);
	free(V_nei1);
	free(C_own1);
	free(C_nei1);

	free(interpola_weight);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::fvDVM::evolution(label It )
{
	
    unsigned long time_sta1 = 0;
	unsigned long time_end1 = 0;
	time_sta1 = rpcc();    
	updateGHbarPvol();
	updateGHbarSurf();
	updateMaxwellWallRho();
	updateGHbarSurfMaxwellWallIn();
/*	updateGHbarSurfSymmetryIn();*/
	updateMacroSurf();
	updateGHsurf();
	updateGHtildeVol();
	updateMacroVol();
/*	updatePressureInOutBC();*/
	time_end1 = rpcc();
	if(mpiReducer_.rank() == 0){
        Info << "Time2 = " <<(double)(time_end1-time_sta1)*1000/CLOCKRATE<< "ms" << nl << endl;//测整个循环
    }
}

void Foam::fvDVM::getCoNum(scalar& maxCoNum, scalar& meanCoNum)
{
	scalar dt = time_.deltaTValue();
#ifndef GETCONUM
	
	scalarField UbyDx =
		mesh_.surfaceInterpolation::deltaCoeffs()
		* (mag(Usurf_) + sqrt(scalar(mesh_.nSolutionD())) * xiMax_);
	
#else
	tmp<surfaceScalarField> deltaCoeffs_
	(
		new surfaceScalarField
		(
			IOobject
			(
				"deltaCoeffs",
				mesh_.pointsInstance(),
				mesh_,
				IOobject::NO_READ,
				IOobject::NO_WRITE,
				false // Do not register
			),
			mesh_,
			dimless/dimLength
		)
    );
	scalarField UbyDx = deltaCoeffs_();
	int64 ownersize = mesh_.owner().size();
	label D = mesh_.nSolutionD();
	scalar xiMax_value = xiMax_.value();
	scalar *Usurf_value = const_cast<scalar*>((Usurf_.begin())->v_ );
	scalar *UbyDx_value = UbyDx.begin();
	// scalar *C_own = dvm_.C_own();
	// scalar *C_nei = dvm_.C_nei();
	CONUM conum;
	conum.UbyDx_value = UbyDx_value;
	conum.Usurf_value = Usurf_value;
	conum.C_own = C_own1;
	conum.C_nei = C_nei1;
	conum.dt = dt;
	conum.xiMax_value = xiMax_value;
	conum.D = D;
	conum.ownersize = ownersize;
	__real_athread_spawn((void *)slave_Func_conum, &conum);
	athread_join();
#endif
	
	maxCoNum = gMax(UbyDx) * dt;//全局最大值，需要进行进程间通信
	meanCoNum = gSum(UbyDx) / UbyDx.size() * dt;
}


const fieldMPIreducer& Foam::fvDVM::mpiReducer() const
{
	return mpiReducer_;
}


// ************************************************************************* //
