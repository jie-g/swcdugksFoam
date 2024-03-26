/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Application
    dugksFoam

Description
    Discrete Unified Gas Kinetic Scheme(Zhaoli Guo, Kun Xu) solver.
    Lianhua Zhu

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvDVM.H"
#include "swlu.h"
#define ATHREAD

#ifdef ATHREAD
extern "C"{
	#include <athread.h>
    void penv_slave0_cycle_init();
    void penv_host0_cycle_init();
    void penv_host0_cycle_count(unsigned long *ic);
    void penv_slave2_gld_init();
}
#endif
static inline unsigned long rpcc(){
	unsigned long addtime;
	asm("rtc %0": "=r" (addtime) : );
	return addtime;
}
#define CLOCKRATE 1.45e9
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    unsigned long time_sta = 0;
	unsigned long time_end = 0;
    unsigned long time_sta1 = 0;
	unsigned long time_end1 = 0;
    double time = 0;
    time_sta = rpcc();
#ifdef ATHREAD
    unsigned long sta = 0,end = 0;
    penv_slave0_cycle_init();
    penv_slave2_gld_init();
#endif
    
    /************************************************************/
    //Add global valid option
    Foam::argList::addBoolOption("dvParallel", "Use discrete velocity domain decomposition\n");
    Foam::argList::addOption(
        "pd",
        "label",
        "num of phy domain");

    //convergence monitor
    scalar TemperatureChange = 1.0;
    scalar rhoChange = 1.0;
    scalar Uchange = 1.0;

#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"
#include "readTimeControlsExplicit.H"

    fvDVM dvm(rho, U, T, &argc, &argv, args);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef ATHREAD
	athread_init();
#endif

    if (dvm.mpiReducer().rank() == 0)
        Info << "\nStarting time loop\n"
             << endl;
    label It = 0;
    //while (runTime.run() && (TemperatureChange > convergeTol))
    while (runTime.run() && (Uchange > convergeTol))
    {

#include "CourantNo.H" // calculate the Co num
#include "readTimeControlsExplicit.H"
#include "setDeltaTvar.H"

        runTime++;
        It++;

        if (dvm.mpiReducer().rank() == 0)
            Info << "Time = " << runTime.timeName() << nl << endl;

        time_sta1 = rpcc();
        dvm.evolution(It);
        time_end1 = rpcc();
        time += (double)(time_end1-time_sta1)*1000/CLOCKRATE;
        //write DF
        ///dvm.writeDFonCells();
        //Pout<<" "<<T<<endl;

        /*
        if (dvm.mpiReducer().rank() < dvm.mpiReducer().nproc())
        {
            runTime.write();
             //std::cout<<dvm.mpiReducer().rank()<<std::endl;
            Info<< "Step =" << It << "  ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;

            if(It%convergeCheckSteps == 0 && It >= convergeCheckSteps)
            {
                tmp<Foam::GeometricField<scalar, Foam::fvPatchField, Foam::volMesh> > 
                    deltaTem = mag(T-Told);
                tmp<Foam::GeometricField<scalar, Foam::fvPatchField, Foam::volMesh> > 
		      	deltaRho = mag(rho-rhoOld);
                tmp<Foam::GeometricField<scalar, Foam::fvPatchField, Foam::volMesh> > 
                    deltaU = mag(U-Uold);

                TemperatureChange = gSum(deltaTem())/gSum(T);
                rhoChange         = gSum(deltaRho())/gSum(rho);
                Uchange           = gSum(deltaU())/gSum(mag(U)());
                
               //if(dvm.mpiReducer().rank()==0)
                
                Info << "Temperature changes = " << TemperatureChange << endl;
                Info << "Density     changes = " << rhoChange         << endl;
                Info << "Velocity    changes = " << Uchange << nl     << endl;
                Told = T;
                rhoOld = rho;
                Uold = U;
            }
        }
        */

        //Info Pout printf std::cout
        runTime.write();
        if (dvm.mpiReducer().rank() == 0){
            // printf("Step = %d ,ExecutionTime = %lf s,ClockTime = %lf s\n",It,runTime.elapsedCpuTime(),runTime.elapsedClockTime());
            Pout << "Step =" << It << "  ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                 << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                 << nl << endl;
        }
        if (It % convergeCheckSteps == 0 && It >= convergeCheckSteps)
        {
            tmp<Foam::GeometricField<scalar, Foam::fvPatchField, Foam::volMesh> >
                deltaTem = mag(T - Told);
            tmp<Foam::GeometricField<scalar, Foam::fvPatchField, Foam::volMesh> >
                deltaRho = mag(rho - rhoOld);
            tmp<Foam::GeometricField<scalar, Foam::fvPatchField, Foam::volMesh> >
                deltaU = mag(U - Uold);

            TemperatureChange = gSum(deltaTem()) / gSum(T);
            rhoChange = gSum(deltaRho()) / gSum(rho);
            Uchange = gSum(deltaU()) / gSum(mag(U)());

            Told = T;
            rhoOld = rho;
            Uold = U;
        }
        

        if (dvm.mpiReducer().rank() == 0)
        {
            Info << "Temperature changes = " << TemperatureChange << endl;
            Info << "Density     changes = " << rhoChange << endl;
            Info << "Velocity    changes = " << Uchange << nl << endl;
        }



    }
    
#ifdef ATHREAD
	printf("ATHREAD OFF\n");
	athread_halt();
#endif


    if (dvm.mpiReducer().rank() == 0)
        Info << "End\n"
             << endl;
    
    time_end = rpcc();
    if (dvm.mpiReducer().rank() == 0){
        Info << "Time1 = " << ((double)(time_end-time_sta)*1000/CLOCKRATE) << "ms" << nl << endl;
        Info << "Time2 = " << time << "ms" << nl << endl;
    }
            

    return 0;
}

// ************************************************************************* //
