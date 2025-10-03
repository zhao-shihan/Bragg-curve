#include "TFile.h"
#include "TGraph.h"

#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

constexpr double um{1};
constexpr double mm{1000 * um};
constexpr double MeV{1};

class ParticleSliceX {
public:
    ParticleSliceX(double deltaX, const TGraph& dEdXGraph) :
        fDeltaX{deltaX},
        fDEdXGraph{&dEdXGraph},
        fRange{},
        fEnergyPerNucleon{},
        fDEDX{} {}

    auto Initialize(double energyPerNucleon) -> void {
        fRange = 0;
        fEnergyPerNucleon = energyPerNucleon;
        fDEDX = fDEdXGraph->Eval(fEnergyPerNucleon);
    }

    auto Range() -> auto {
        return fRange;
    }

    auto EnergyPerNucleon() -> auto {
        return fEnergyPerNucleon;
    }

    auto DEDX() -> auto {
        return fDEDX;
    }

    auto Stopped() -> auto {
        return fEnergyPerNucleon <= 0;
    }

    auto Step() -> void {
        fRange += fDeltaX;
        fEnergyPerNucleon -= fDEDX * fDeltaX;
        if (Stopped()) {
            fEnergyPerNucleon = 0;
            fDEDX = 0;
        } else {
            fDEDX = fDEdXGraph->Eval(fEnergyPerNucleon);
        }
    }

private:
    double fDeltaX;
    const TGraph* fDEdXGraph;

    double fRange;
    double fEnergyPerNucleon;
    double fDEDX;
};

class ReverseParticleSliceX {
public:
    ReverseParticleSliceX(double deltaX, const TGraph& dEdXGraph) :
        fDeltaX{deltaX},
        fDEdXGraph{&dEdXGraph},
        fDEDX{},
        fEnergyPerNucleon{},
        fRange{} {}

    auto Initialize(double range) -> void {
        fDEDX = 0;
        fEnergyPerNucleon = 0;
        fRange = range;
    }

    auto DEDX() -> auto {
        return fDEDX;
    }

    auto EnergyPerNucleon() -> auto {
        return fEnergyPerNucleon;
    }

    auto Range() -> auto {
        return fRange;
    }

    auto Returned() -> auto {
        return fRange <= 0;
    }

    auto ReverseStep() -> void {
        fDEDX = fDEdXGraph->Eval(fEnergyPerNucleon);
        fEnergyPerNucleon += fDEDX * fDeltaX;
        fRange -= fDeltaX;
        if (Returned()) {
            fRange = 0;
        }
    }

private:
    double fDeltaX;
    const TGraph* fDEdXGraph;

    double fDEDX;
    double fEnergyPerNucleon;
    double fRange;
};

std::unique_ptr<TGraph> ParseDEDXData(const char filePath[]) {
    std::ifstream dataFile{filePath};
    if (not dataFile.good()) {
        std::cerr << "Cannot open input file:" << filePath << std::endl;
        return {};
    }
    auto dEdXGraph{std::make_unique<TGraph>()};
    std::string str;
    for (int i{}; i < 3; ++i) {
        std::getline(dataFile, str); // skip the first three lines
    }

    // energy: MeV/u, Range: micron - do not depends on A.

    double energy, r[10], dummy;
    while (std::getline(dataFile, str)) {
        std::stringstream ss(str);
        ss >> energy >> r[0] >> dummy >> r[1] >> dummy >> r[2] >> dummy >> r[3] >> dummy >> r[4] /* used for calculation */
            >> dummy >> r[5] >> dummy >> r[6] >> dummy >> r[7] >> dummy >> r[8] >> dummy >> r[9];
        dEdXGraph->AddPoint(energy, r[4]);
    }
    return dEdXGraph;
}

auto dEdxSliceX(const char dataFilePath[], double targetRange, double deltaX) {
    const auto dEdXGraph{ParseDEDXData(dataFilePath)};
    dEdXGraph->DrawClone();

    // Find energy for targetRange (mm) range
    ReverseParticleSliceX reverseParticleSliceX{deltaX, *dEdXGraph};
    reverseParticleSliceX.Initialize(targetRange * mm);
    while (not reverseParticleSliceX.Returned()) {
        reverseParticleSliceX.ReverseStep();
    }
    const auto energyPerNucleon{reverseParticleSliceX.EnergyPerNucleon()};
    std::cout << "Particle with " << energyPerNucleon / MeV << " MeV/u will has " << targetRange << " mm range.\n";

    // Calculate Bragg curve
    TGraph braggCurve;
    braggCurve.SetName("Bragg_curve");
    ParticleSliceX particleSliceX{deltaX, *dEdXGraph};
    particleSliceX.Initialize(energyPerNucleon);
    braggCurve.AddPoint(particleSliceX.Range() / mm, particleSliceX.DEDX() / (MeV / mm));
    while (not particleSliceX.Stopped()) {
        particleSliceX.Step();
        braggCurve.AddPoint(particleSliceX.Range() / mm, particleSliceX.DEDX() / (MeV / mm));
    }
    braggCurve.DrawClone();
}
