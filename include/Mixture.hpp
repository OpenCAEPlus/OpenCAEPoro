#pragma once
#include "OpenCAEPoro_consts.hpp"
#include <iostream>
#include <vector>

using namespace std;

// total fluid, contains all phase
class Mixture
{
    friend class Bulk;
    friend class Well;

public:
    Mixture() = default;
    virtual ~Mixture(){};

    // return type
    int getType() { return MixtureType; }
    // black oil
    virtual bool empty_PVDG() = 0;

    virtual void Flash_Sj(const OCP_DBL Pin, const OCP_DBL Pbbin, const OCP_DBL Tin,
                          const OCP_DBL* Sjin, OCP_DBL Vpore, const OCP_DBL* Ziin)   = 0;
    virtual void Flash_Ni(const OCP_DBL Pin, const OCP_DBL Tin, const OCP_DBL* Niin) = 0;

    virtual void getProp(){};

    // return xi
    virtual OCP_DBL xiPhase(OCP_DBL Pin, OCP_DBL T, OCP_DBL* Ziin) = 0;

    // return rho
    virtual OCP_DBL rhoPhase(OCP_DBL Pin, OCP_DBL T, OCP_DBL* Ziin) = 0;

    // return gamma
    virtual OCP_DBL gammaPhaseO(OCP_DBL Pin, OCP_DBL Pbbin)              = 0;
    virtual OCP_DBL gammaPhaseW(OCP_DBL Pin)                            = 0;
    virtual OCP_DBL gammaPhaseG(OCP_DBL Pin)                            = 0;
    virtual OCP_DBL gammaPhaseOG(OCP_DBL Pin, OCP_DBL Tin, OCP_DBL* Ziin) = 0;

    // check
    void checkNi(const OCP_DBL* Ni)
    {
        bool flag = false;
        for (int i = 0; i < Nc; i++) {
            if (Ni[i] < 0) {
                cout << "###WARNING:  ";
                ERRORcheck("Ni < 0 ");
            }

            if (Ni[i] > 0) flag = true;
        }
        if (!flag) {
            cout << "###ERROR:  ";
            ERRORcheck("All Ni <= 0 ");
            exit(0);
        }
    }

protected:
    int MixtureType;

    int    Np; // num of phase
    int    Nc; // num of component
    OCP_DBL P;  // Pressure
    OCP_DBL T;  // Temperature

    std::vector<OCP_DBL> Ni;         // molar of component : Nc
    std::vector<bool>   PhaseExist; // existence of phase : Np
    std::vector<OCP_DBL> S;          // saturation of phase : Np
    std::vector<OCP_DBL> Rho;        // mass density of phase : Np
    std::vector<OCP_DBL> Xi;         // molar density of phase: Np
    std::vector<OCP_DBL> Cij;        // Nij / Nj : Np*Nc
    std::vector<OCP_DBL> Mu;         // viscosity of phase: Np
    std::vector<OCP_DBL> V;          // volume of phase

    OCP_DBL              Vf;  // volume of fluids
    OCP_DBL              Vfp; //
    std::vector<OCP_DBL> Vfi; // dVf / dNi   : Nc
};
