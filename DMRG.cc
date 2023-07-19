#include "math.h"
#include <complex.h>
#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <vector>
#include <string>
#include <climits>

#include "itensor/all.h"

using namespace itensor;

int main(int argc, char const *argv[])
{
	std::clock_t begin = std::clock();
    int N = atoi(argv[1]); // then bigger sizes
    double h = atof(argv[2]);
    double J1 = atof(argv[3]);
    double J2 = atof(argv[4]); // change from 0 to 10

    std::cout << "N = " << N << ", h = " << h << ", J1 = " << J1 << ", J2 = " << J2 << "\n";

    std::vector<itensor::MPO> obsMagVec;
    std::vector<itensor::MPO> NNNeighborVec;
    std::vector<itensor::MPO> NNeighborVec;
    std::vector<itensor::MPO> xObsVec;

    auto sites = SpinHalf(N);
    auto psi = MPS(sites);

    /*****************************
    *  Hamiltonian construction  *

    j1 = 0.1 and j2 = 0.1sd
	*****************************/


    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N; ++j)
        {
            ampo += -J1*4.0,"Sz",j,"Sz",j+1;
        }

    for(int j = 1; j < N-1; ++j)
        {
            ampo += -J2*4.0,"Sz",j,"Sz",j+2;
        }


    for (int j = 1; j <= N; ++j)
        {
            ampo += h*2.0,"Sx",j;
        }


    //  Magnetization measurement loop start
    for(int i=1; i<=N; i++)
    {
        auto aobs  = AutoMPO(sites);
        aobs += 2.0,"Sz",i;
        obsMagVec.push_back(MPO(aobs));
    }   // Magnetization measurement loop end

    //  Correlation NNN loop start
    for (int i = 1; i <= (N-2); i++)
    {
        auto aobs = AutoMPO(sites);
        aobs += 4.0, "Sz", i, "Sz", i + 2;
        NNNeighborVec.push_back(MPO(aobs));
    }   // Correlation NNN loop end

    //  Correlation NN loop start
    for (int i = 1; i <= (N-1); i++)
    {
        auto aobs = AutoMPO(sites);
        aobs += 4.0, "Sz", i, "Sz", i + 1;
        NNeighborVec.push_back(MPO(aobs));
    }   // Correlation NNN loop end

    //  x loop start
    for (int i = 1; i <= (N); i++)
    {
        auto aobs = AutoMPO(sites);
        aobs += 2.0, "Sx", i;
        xObsVec.push_back(MPO(aobs));
    }   // x loop end
                
    auto H = MPO(ampo);

    auto sweeps = Sweeps(1);
    sweeps.cutoff() = 1E-8;

    double bd = 100;
    int MaxBond = 75000;
    
    double energyFin = 100000;
    double energyIni;

    double sweepBDChangeThresholdValue = 0.005;

    for (int cont = 1; bd <= MaxBond ; ++cont)
        {

            energyIni = energyFin;
            energyFin = dmrg(psi,H,sweeps,{"Quiet",true});
            printfln("GroundStateEnergy = %.20f",energyFin);

            /* Checking for BD increase */
            if ((energyIni-energyFin)/energyIni < sweepBDChangeThresholdValue)
            {
                bd += 100;
                sweeps.maxm() = bd;
                printfln("BD = %.20i",bd);
            }
        }// LOOP FOR SWEEPS - END

    // for (int cont = 1; cont < 5000 ; ++cont)
    //     {
    //         energyFin = dmrg(psi,H,sweeps,{"Quiet",true});
    //     }// LOOP FOR SWEEPS - END
    //     printfln("GroundStateEnergy = %.20f",energyFin);

    
        // Sz measurement (Magnon Density) loop for each site start
        for(int j=1;j<=N;j++)
        {
            auto ev = overlap(psi,obsMagVec[j-1],psi);

            std::cout << "\nMagnonDensityAtSite " << j << " is " << ev;
        }   // Sz measurement (Magnon Density) loop for each site end
        std::cout << "\n\n";


    double NNNeighborCor = 0;
    // NNNeighbor loop for each site start
    for (int j = 0; j < NNNeighborVec.size(); j++)
    {
        auto ev = overlap(psi, NNNeighborVec[j], psi);

        NNNeighborCor += ev;
        std::cout << "\nNNNeighbor " << j+1 << " is " << ev;
    }   // NNNeighbor loop for each site end
    std::cout << "\nNNNeighborCorr*L " << NNNeighborCor;
    std::cout << "\nNNNeighborCorr " << -NNNeighborCor/(N-2) << "\n\n";

    double NNeighborCor = 0;

    // NNeighbor loop for each site start
    for (int j = 0; j < NNeighborVec.size(); j++)
    {
        auto ev = overlap(psi, NNeighborVec[j], psi);

        NNeighborCor += ev;
        std::cout << "\nNNeighbor " << j+1<< " is " << ev;
    }   // NNeighbor loop for each site end
    std::cout << "\nNNeighborCorr*L " << NNeighborCor;
    std::cout << "\nNNeighborCorr " << -NNeighborCor/(N-1) << "\n\n";
    std::cout << "\n\n";

    double xObsVal = 0;
    // xObsVec loop for each site start
    for (int j = 0; j < xObsVec.size(); j++)
    {
        auto ev = overlap(psi, xObsVec[j], psi);

        xObsVal += ev;
        std::cout << "\nxObs " << j+1 << " is " << ev;
    }   // NNeighbor loop for each site end
    std::cout << "\nxObsCorr*L " << xObsVal;
	std::cout << "\nxObsCorr " << -xObsVal/(N);
    std::cout << "\n\n";

    Print(psi);

    std::stringstream psiString; 
    std::stringstream sitesString; 
    psiString   << "./testFolder/" << "PsiGSfile_N-" << N << "_h-" << h << "_J1-" << J1 << "_J2-" << J2;
    sitesString << "./testFolder/" << "PsiSitesfile_N-" << N << "_h-" << h << "_J1-" << J1 << "_J2-" << J2;
    const std::string sitesTmp = sitesString.str(); const char* sitesName = sitesTmp.c_str();
    const std::string psiTmp = psiString.str(); const char* psiName = psiTmp.c_str();

    writeToFile(sitesName, sites);
    writeToFile(psiName,psi);

    std::clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "\nCODE FINISHED - Elapsed time: " << elapsed_secs << " s";

    std::cout << "\n\n"; 
    return 0;
}
