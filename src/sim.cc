// ============================================================================
// CERN BL4S 2027 — Bismuth Scintillator Test-Beam Simulation (DESY-II / ELSA)
// Single-file serial Geant4 application controlled via CLI flags.
// Physics: FTFP_BERT + G4EmStandardPhysics_option4 (high-precision thin-slab EM)
// ============================================================================

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "FTFP_BERT.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4Tubs.hh"
#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VisAttributes.hh"
#include "Randomize.hh"

#ifdef BL4S_VIS
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#endif

// ----------------------------------------------------------------------------
// 1. Simulation Configuration (Populated from CLI Arguments)
// ----------------------------------------------------------------------------
struct SimConfig {
    std::string material = "BiNanoPDMS";  // EJ200 | BiNanoPDMS | BiPh3PVK | LYSO
    double loading_wt = 0.25;             // Inorganic / BiPh3 weight fraction
    double thickness_mm = 2.0;            // DUT disc thickness [mm]
    double diameter_mm = 25.0;            // DUT disc diameter [mm]
    double energy_GeV = 3.0;              // Beam mean energy [GeV]
    double espread_GeV = 0.15;            // DESY-II momentum spread [GeV]
    double beam_diam_mm = 20.0;           // DESY-II beam spot diameter [mm]
    std::string particle = "e-";          // e- or e+
    int events = 10000;
    long seed = 2027;
    std::string out_csv = "results/sim_output.csv";
    bool vis = false;
};

static SimConfig ParseArgs(int argc, char** argv) {
    SimConfig cfg;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto nextVal = [&]() -> std::string {
            if (i + 1 < argc) return argv[++i];
            std::cerr << "Missing value for argument: " << arg << "\n";
            std::exit(1);
        };
        if (arg == "--material") cfg.material = nextVal();
        else if (arg == "--loading") cfg.loading_wt = std::stod(nextVal());
        else if (arg == "--thickness") cfg.thickness_mm = std::stod(nextVal());
        else if (arg == "--diameter") cfg.diameter_mm = std::stod(nextVal());
        else if (arg == "--energy") cfg.energy_GeV = std::stod(nextVal());
        else if (arg == "--espread") cfg.espread_GeV = std::stod(nextVal());
        else if (arg == "--beam-diam") cfg.beam_diam_mm = std::stod(nextVal());
        else if (arg == "--particle") cfg.particle = nextVal();
        else if (arg == "--events") cfg.events = std::stoi(nextVal());
        else if (arg == "--seed") cfg.seed = std::stol(nextVal());
        else if (arg == "--out") cfg.out_csv = nextVal();
        else if (arg == "--vis") cfg.vis = true;
        else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: ./bl4s_sim [--material EJ200|BiNanoPDMS|BiPh3PVK|LYSO] "
                      << "[--loading 0.25] [--thickness 2.0] [--energy 3.0] [--espread 0.15] "
                      << "[--particle e-] [--events 25000] [--seed 2027] [--out file.csv] [--vis]\n";
            std::exit(0);
        }
    }
    return cfg;
}

// ----------------------------------------------------------------------------
// 2. Detector Construction (DESY-II / ELSA Beamline Telescope + DUT + Cal)
// ----------------------------------------------------------------------------

/// Logical volumes the stepping action attributes energy deposits to.
struct DetectorVolumes {
    G4LogicalVolume* s1 = nullptr;
    G4LogicalVolume* s2 = nullptr;
    G4LogicalVolume* dut = nullptr;
    G4LogicalVolume* cal = nullptr;
};

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    explicit DetectorConstruction(const SimConfig& cfg) : fCfg(cfg) {}

    G4VPhysicalVolume* Construct() override {
        auto* nist = G4NistManager::Instance();
        auto* air = nist->FindOrBuildMaterial("G4_AIR");
        auto* silicon = nist->FindOrBuildMaterial("G4_Si");
        auto* ej200 = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
        auto* pbGlass = nist->FindOrBuildMaterial("G4_GLASS_LEAD");

        fDutMaterial = BuildDutMaterial(nist, ej200);

        // World Volume (Air, 40 x 40 x 120 cm)
        auto* worldSolid = new G4Box("World", 20 * cm, 20 * cm, 60 * cm);
        auto* worldLV = new G4LogicalVolume(worldSolid, air, "WorldLV");
        worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
        auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World", nullptr, false, 0);

        // Trigger Paddles S1 (z = -60 mm) and S2 (z = +60 mm): 10 x 10 x 2 mm EJ-200
        auto* trigSolid = new G4Box("TriggerSolid", 5.0 * mm, 5.0 * mm, 1.0 * mm);
        fVol.s1 = new G4LogicalVolume(trigSolid, ej200, "S1_LV");
        fVol.s2 = new G4LogicalVolume(trigSolid, ej200, "S2_LV");
        auto* trigVis = new G4VisAttributes(G4Colour(0.15, 0.45, 0.95, 0.6));
        trigVis->SetForceSolid(true);
        fVol.s1->SetVisAttributes(trigVis);
        fVol.s2->SetVisAttributes(trigVis);

        new G4PVPlacement(nullptr, {0, 0, -60.0 * mm}, fVol.s1, "S1", worldLV, false, 0, true);
        new G4PVPlacement(nullptr, {0, 0, +60.0 * mm}, fVol.s2, "S2", worldLV, false, 0, true);

        // 6-Plane Silicon Beam Telescope (3 upstream, 3 downstream): 20 x 20 x 0.05 mm
        auto* telSolid = new G4Box("TelSolid", 10.0 * mm, 10.0 * mm, 0.025 * mm);
        auto* telLV = new G4LogicalVolume(telSolid, silicon, "TelLV");
        auto* telVis = new G4VisAttributes(G4Colour(0.55, 0.55, 0.65, 0.4));
        telVis->SetForceSolid(true);
        telLV->SetVisAttributes(telVis);

        const std::vector<double> telZ = {-42.0, -30.0, -18.0, +18.0, +30.0, +42.0};
        for (size_t i = 0; i < telZ.size(); ++i) {
            new G4PVPlacement(nullptr, {0, 0, telZ[i] * mm}, telLV, "TelescopePlane", worldLV, false, i, true);
        }

        // Device Under Test (DUT): Cylindrical Scintillator Disc at z = 0
        auto* dutSolid = new G4Tubs("DUTSolid", 0.0, 0.5 * fCfg.diameter_mm * mm,
                                    0.5 * fCfg.thickness_mm * mm, 0.0, 360.0 * deg);
        fVol.dut = new G4LogicalVolume(dutSolid, fDutMaterial, "DUT_LV");
        auto* dutVis = new G4VisAttributes(G4Colour(0.95, 0.45, 0.10, 0.85));
        dutVis->SetForceSolid(true);
        fVol.dut->SetVisAttributes(dutVis);
        new G4PVPlacement(nullptr, {0, 0, 0}, fVol.dut, "DUT", worldLV, false, 0, true);

        // SiPM Silicon Sensor Pad (3 x 3 x 0.5 mm) flush on back of DUT
        auto* sipmSolid = new G4Box("SiPMSolid", 1.5 * mm, 1.5 * mm, 0.25 * mm);
        auto* sipmLV = new G4LogicalVolume(sipmSolid, silicon, "SiPM_LV");
        auto* sipmVis = new G4VisAttributes(G4Colour(0.10, 0.80, 0.30, 0.9));
        sipmVis->SetForceSolid(true);
        sipmLV->SetVisAttributes(sipmVis);
        const double sipmZ = (0.5 * fCfg.thickness_mm + 0.25) * mm;
        new G4PVPlacement(nullptr, {0, 0, sipmZ}, sipmLV, "SiPM", worldLV, false, 0, true);

        // Lead-Glass Electromagnetic Calorimeter (10 x 10 x 37 cm) at z = +285 mm
        auto* calSolid = new G4Box("CalSolid", 5.0 * cm, 5.0 * cm, 18.5 * cm);
        fVol.cal = new G4LogicalVolume(calSolid, pbGlass, "Cal_LV");
        auto* calVis = new G4VisAttributes(G4Colour(0.75, 0.70, 0.20, 0.35));
        calVis->SetForceSolid(true);
        fVol.cal->SetVisAttributes(calVis);
        new G4PVPlacement(nullptr, {0, 0, 285.0 * mm}, fVol.cal, "Calorimeter", worldLV, false, 0, true);

        return worldPV;
    }

    const DetectorVolumes& Volumes() const { return fVol; }
    G4Material* DutMaterial() const { return fDutMaterial; }

private:
    G4Material* BuildDutMaterial(G4NistManager* nist, G4Material* ej200) {
        if (fCfg.material == "EJ200") return ej200;

        auto* elH  = nist->FindOrBuildElement("H");
        auto* elC  = nist->FindOrBuildElement("C");
        auto* elN  = nist->FindOrBuildElement("N");
        auto* elO  = nist->FindOrBuildElement("O");
        auto* elSi = nist->FindOrBuildElement("Si");
        auto* elBr = nist->FindOrBuildElement("Br");
        auto* elAg = nist->FindOrBuildElement("Ag");
        auto* elSb = nist->FindOrBuildElement("Sb");
        auto* elCs = nist->FindOrBuildElement("Cs");
        auto* elBi = nist->FindOrBuildElement("Bi");

        if (fCfg.material == "BiNanoPDMS") {
            // 3D Double Perovskite Cs2AgBi0.9Sb0.1Br6 (scaled x10 for integer stoichiometry)
            const double rhoCrystal = 4.80;  // g/cm3
            auto* crystal = new G4Material("Cs2AgBi09Sb01Br6", rhoCrystal * g / cm3, 5);
            crystal->AddElement(elCs, 20);
            crystal->AddElement(elAg, 10);
            crystal->AddElement(elBi, 9);
            crystal->AddElement(elSb, 1);
            crystal->AddElement(elBr, 60);

            // PDMS Sylgard 184 repeat unit: C2H6OSi
            const double rhoPDMS = 1.02;  // g/cm3
            auto* pdms = new G4Material("PDMS", rhoPDMS * g / cm3, 4);
            pdms->AddElement(elC, 2);
            pdms->AddElement(elH, 6);
            pdms->AddElement(elO, 1);
            pdms->AddElement(elSi, 1);

            // Exact Rule-of-Mixtures composite density
            const double w = std::clamp(fCfg.loading_wt, 0.01, 0.99);
            const double rhoComp = 1.0 / (w / rhoCrystal + (1.0 - w) / rhoPDMS);
            auto* comp = new G4Material("BiNanoPDMS", rhoComp * g / cm3, 2);
            comp->AddMaterial(crystal, w);
            comp->AddMaterial(pdms, 1.0 - w);
            return comp;
        }

        if (fCfg.material == "BiPh3PVK") {
            // Triphenylbismuth: Bi(C6H5)3 -> C18H15Bi
            const double rhoBiPh3 = 1.585;  // g/cm3
            auto* biph3 = new G4Material("BiPh3", rhoBiPh3 * g / cm3, 3);
            biph3->AddElement(elC, 18);
            biph3->AddElement(elH, 15);
            biph3->AddElement(elBi, 1);

            // Poly(N-vinylcarbazole) PVK: C14H11N
            const double rhoPVK = 1.20;  // g/cm3
            auto* pvk = new G4Material("PVK", rhoPVK * g / cm3, 3);
            pvk->AddElement(elC, 14);
            pvk->AddElement(elH, 11);
            pvk->AddElement(elN, 1);

            const double w = std::clamp(fCfg.loading_wt, 0.01, 0.50);
            const double rhoComp = 1.0 / (w / rhoBiPh3 + (1.0 - w) / rhoPVK);
            auto* comp = new G4Material("BiPh3PVK", rhoComp * g / cm3, 2);
            comp->AddMaterial(biph3, w);
            comp->AddMaterial(pvk, 1.0 - w);
            return comp;
        }

        if (fCfg.material == "LYSO") {
            auto* elLu = nist->FindOrBuildElement("Lu");
            auto* elY  = nist->FindOrBuildElement("Y");
            auto* lyso = new G4Material("LYSO", 7.10 * g / cm3, 4);
            lyso->AddElement(elLu, 18);
            lyso->AddElement(elY, 2);
            lyso->AddElement(elSi, 10);
            lyso->AddElement(elO, 50);
            return lyso;
        }

        std::cerr << "Unknown material: " << fCfg.material << "\n";
        std::exit(1);
    }

    SimConfig fCfg;
    DetectorVolumes fVol;
    G4Material* fDutMaterial = nullptr;
};

// ----------------------------------------------------------------------------
// 3. Primary Generator Action (DESY-II / ELSA Electron/Positron Beam)
// ----------------------------------------------------------------------------
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    explicit PrimaryGeneratorAction(const SimConfig& cfg) : fCfg(cfg) {
        fGun = new G4ParticleGun(1);
        auto* particleDef = G4ParticleTable::GetParticleTable()->FindParticle(cfg.particle);
        if (!particleDef) {
            std::cerr << "Unknown particle: " << cfg.particle << "\n";
            std::exit(1);
        }
        fGun->SetParticleDefinition(particleDef);
        fGun->SetParticleMomentumDirection({0.0, 0.0, 1.0});
    }

    ~PrimaryGeneratorAction() override { delete fGun; }

    void GeneratePrimaries(G4Event* event) override {
        // Uniform circular beam spot of diameter beam_diam_mm (centered on trigger)
        const double rMax = 0.5 * fCfg.beam_diam_mm * mm;
        const double r = rMax * std::sqrt(G4UniformRand());
        const double phi = 2.0 * pi * G4UniformRand();
        fGun->SetParticlePosition({r * std::cos(phi), r * std::sin(phi), -100.0 * mm});

        // Gaussian momentum/energy spread (DESY-II typical sigma ~ 0.15 GeV)
        double eGeV = fCfg.energy_GeV;
        if (fCfg.espread_GeV > 0.0) {
            eGeV = std::max(0.1, G4RandGauss::shoot(fCfg.energy_GeV, fCfg.espread_GeV));
        }
        fGun->SetParticleEnergy(eGeV * GeV);
        fGun->GeneratePrimaryVertex(event);
    }

private:
    SimConfig fCfg;
    G4ParticleGun* fGun = nullptr;
};

// ----------------------------------------------------------------------------
// 4. Event Record, Run and Stepping Actions (clean CSV + metadata writer)
// ----------------------------------------------------------------------------

/// One simulated event. SteppingAction fills the raw quantities, EventAction
/// derives the trigger and scattering angles, RunAction writes the row — so no
/// per-quantity accessor boilerplate is needed.
struct EventRecord {
    double edepS1_MeV = 0.0;
    double edepS2_MeV = 0.0;
    double edepDut_MeV = 0.0;
    double edepCal_GeV = 0.0;
    int bremCount = 0;
    double bremEnergy_MeV = 0.0;
    bool enteredDut = false;
    bool exitedDut = false;
    G4ThreeVector dirIn{0, 0, 1};
    G4ThreeVector dirOut{0, 0, 1};
    int triggered = 0;
    double scatterX_mrad = 0.0;
    double scatterY_mrad = 0.0;
    double scatterTotal_mrad = 0.0;
};

class RunAction : public G4UserRunAction {
public:
    RunAction(const SimConfig& cfg, const DetectorConstruction* det) : fCfg(cfg), fDet(det) {}

    void BeginOfRunAction(const G4Run*) override {
        if (fCfg.vis) return;
        fCsv.open(fCfg.out_csv);
        if (!fCsv.is_open()) {
            std::cerr << "Failed to open CSV output: " << fCfg.out_csv << "\n";
            std::exit(1);
        }
        fCsv << "event_id,edep_S1_MeV,edep_S2_MeV,edep_DUT_MeV,edep_Cal_GeV,"
             << "scatter_x_mrad,scatter_y_mrad,scatter_total_mrad,"
             << "brem_photons,brem_energy_MeV,triggered\n";
    }

    void EndOfRunAction(const G4Run*) override {
        if (fCfg.vis) return;
        fCsv.close();

        // Companion metadata with Geant4's own material properties
        const auto* mat = fDet->DutMaterial();
        const double radlen_mm = mat->GetRadlen() / mm;
        std::ofstream meta(fCfg.out_csv + ".meta.json");
        meta << std::fixed << std::setprecision(6)
             << "{\n"
             << "  \"material\": \"" << fCfg.material << "\",\n"
             << "  \"loading_wt\": " << fCfg.loading_wt << ",\n"
             << "  \"thickness_mm\": " << fCfg.thickness_mm << ",\n"
             << "  \"energy_GeV\": " << fCfg.energy_GeV << ",\n"
             << "  \"events\": " << fCfg.events << ",\n"
             << "  \"density_g_cm3\": " << mat->GetDensity() / (g / cm3) << ",\n"
             << "  \"radlen_mm\": " << radlen_mm << ",\n"
             << "  \"x_over_X0\": " << (fCfg.thickness_mm / radlen_mm) << "\n"
             << "}\n";
    }

    void WriteRow(int eventId, const EventRecord& r) {
        if (!fCsv.is_open()) return;
        fCsv << eventId << "," << std::fixed << std::setprecision(5)
             << r.edepS1_MeV << "," << r.edepS2_MeV << "," << r.edepDut_MeV << "," << r.edepCal_GeV << ","
             << r.scatterX_mrad << "," << r.scatterY_mrad << "," << r.scatterTotal_mrad << ","
             << r.bremCount << "," << r.bremEnergy_MeV << "," << r.triggered << "\n";
    }

private:
    SimConfig fCfg;
    const DetectorConstruction* fDet = nullptr;
    std::ofstream fCsv;
};

class EventAction : public G4UserEventAction {
public:
    explicit EventAction(RunAction* runAction) : fRunAction(runAction) {}

    void BeginOfEventAction(const G4Event*) override { fRec = EventRecord{}; }
    void EndOfEventAction(const G4Event* event) override;

    EventRecord& Record() { return fRec; }

private:
    RunAction* fRunAction = nullptr;
    EventRecord fRec;
};

void EventAction::EndOfEventAction(const G4Event* event) {
    // Hardware coincidence trigger: S1 AND S2 (> 50 keV each)
    fRec.triggered = (fRec.edepS1_MeV > 0.05 && fRec.edepS2_MeV > 0.05) ? 1 : 0;

    if (fRec.enteredDut && fRec.exitedDut) {
        const double cosTheta = std::clamp(fRec.dirIn.dot(fRec.dirOut), -1.0, 1.0);
        fRec.scatterTotal_mrad = std::acos(cosTheta) * 1000.0;
        fRec.scatterX_mrad = (std::atan2(fRec.dirOut.x(), fRec.dirOut.z()) -
                              std::atan2(fRec.dirIn.x(), fRec.dirIn.z())) * 1000.0;
        fRec.scatterY_mrad = (std::atan2(fRec.dirOut.y(), fRec.dirOut.z()) -
                              std::atan2(fRec.dirIn.y(), fRec.dirIn.z())) * 1000.0;
    }

    fRunAction->WriteRow(event->GetEventID(), fRec);
}

class SteppingAction : public G4UserSteppingAction {
public:
    SteppingAction(const DetectorConstruction* detector, EventAction* eventAction)
        : fDet(detector), fEvt(eventAction) {}

    void UserSteppingAction(const G4Step* step) override {
        const auto* preVol = step->GetPreStepPoint()->GetPhysicalVolume();
        if (!preVol) return;
        const auto* postVol = step->GetPostStepPoint()->GetPhysicalVolume();

        // Volumes() is only populated once Construct() has run, so resolve the
        // handles per step rather than caching them at registration time.
        const auto& vol = fDet->Volumes();
        auto* preLV = preVol->GetLogicalVolume();
        auto* postLV = postVol ? postVol->GetLogicalVolume() : nullptr;
        auto& rec = fEvt->Record();

        // Lead-glass calorimeter: its response is not plotted, and tracking a
        // multi-GeV EM shower through 18.5 cm of lead glass dominates runtime.
        // Record the entering kinetic energy as the calorimeter signal and stop
        // the shower at the front face.
        if (postLV == vol.cal && preLV != vol.cal) {
            rec.edepCal_GeV += step->GetPostStepPoint()->GetKineticEnergy() / GeV;
            step->GetTrack()->SetTrackStatus(fStopAndKill);
            return;
        }

        const double edepMeV = step->GetTotalEnergyDeposit() / MeV;
        if (edepMeV > 0.0) {
            if (preLV == vol.s1) rec.edepS1_MeV += edepMeV;
            else if (preLV == vol.s2) rec.edepS2_MeV += edepMeV;
            else if (preLV == vol.dut) rec.edepDut_MeV += edepMeV;
        }

        if (preLV != vol.dut) return;

        // Primary-beam scattering across the DUT: entry and exit directions
        if (step->GetTrack()->GetTrackID() == 1) {
            if (!rec.enteredDut) {
                rec.dirIn = step->GetPreStepPoint()->GetMomentumDirection().unit();
                rec.enteredDut = true;
            }
            if (postLV != vol.dut) {
                rec.dirOut = step->GetPostStepPoint()->GetMomentumDirection().unit();
                rec.exitedDut = true;
            }
        }

        // Bremsstrahlung gammas (> 1 MeV) created inside the DUT
        const auto* secondaries = step->GetSecondaryInCurrentStep();
        if (!secondaries) return;
        for (const auto* sec : *secondaries) {
            if (sec->GetDefinition() != G4Gamma::Definition()) continue;
            const double keMeV = sec->GetKineticEnergy() / MeV;
            if (keMeV > 1.0) {
                rec.bremCount++;
                rec.bremEnergy_MeV += keMeV;
            }
        }
    }

private:
    const DetectorConstruction* fDet = nullptr;
    EventAction* fEvt = nullptr;
};

// ----------------------------------------------------------------------------
// 5. Optional Interactive Viewer
// ----------------------------------------------------------------------------
#ifdef BL4S_VIS
static void RunViewer(int argc, char** argv) {
    auto* visManager = new G4VisExecutive();
    visManager->Initialize();
    auto* ui = new G4UIExecutive(argc, argv);

    auto* uiManager = G4UImanager::GetUIpointer();
    uiManager->ApplyCommand("/vis/open OGL 900x700-0+0");
    uiManager->ApplyCommand("/vis/drawVolume");
    uiManager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 65 25");
    uiManager->ApplyCommand("/vis/scene/add/trajectories smooth");
    uiManager->ApplyCommand("/vis/scene/endOfEventAction accumulate");
    uiManager->ApplyCommand("/run/beamOn 20");

    ui->SessionStart();
    delete ui;
    delete visManager;
}
#endif

// ----------------------------------------------------------------------------
// 6. Main Entry Point
// ----------------------------------------------------------------------------
int main(int argc, char** argv) {
    const SimConfig cfg = ParseArgs(argc, argv);

#ifndef BL4S_VIS
    if (cfg.vis) {
        std::cerr << "This build has no visualization support (BL4S_VIS=OFF).\n"
                  << "Reconfigure with: cmake -DBL4S_VIS=ON ..\n";
        return 1;
    }
#endif

    CLHEP::HepRandom::setTheSeed(cfg.seed);

    auto* runManager = new G4RunManager();
    auto* detector = new DetectorConstruction(cfg);
    runManager->SetUserInitialization(detector);

    auto* physics = new FTFP_BERT(0);
    physics->ReplacePhysics(new G4EmStandardPhysics_option4());
    runManager->SetUserInitialization(physics);

    auto* runAction = new RunAction(cfg, detector);
    auto* eventAction = new EventAction(runAction);
    runManager->SetUserAction(new PrimaryGeneratorAction(cfg));
    runManager->SetUserAction(runAction);
    runManager->SetUserAction(eventAction);
    runManager->SetUserAction(new SteppingAction(detector, eventAction));
    runManager->Initialize();

    if (cfg.vis) {
#ifdef BL4S_VIS
        RunViewer(argc, argv);
#endif
    } else {
        runManager->BeamOn(cfg.events);
    }

    delete runManager;
    return 0;
}
