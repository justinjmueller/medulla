import os
import ROOT as r


r.gInterpreter.Declare('''
#include <vector>

std::vector<double> checkPOT(const char* filename){

    // Return: [pot_sum, gates_on_sum, gates_off_sum]
    double pot_sum = 0.0;
    double gates_on_sum = 0.0;
    double gates_off_sum = 0.0;

    TFile* myfile = TFile::Open(filename, "READ");
    if (!myfile || myfile->IsZombie()) {
        if (myfile) myfile->Close();
        return std::vector<double>{pot_sum, gates_on_sum, gates_off_sum};
    }

    // Try the tree name used in the SpinePlot configs first.
    TTree* t = (TTree*)myfile->Get("recTree");
    if (!t) {
        // Fallback: some files may store the default tree as "selected;1" but Get("selected")
        // should work; if not, try a common alternative.
        t = (TTree*)myfile->Get("Selected");
    }
    if (!t) {
        myfile->Close();
        return std::vector<double>{pot_sum, gates_on_sum, gates_off_sum};
    }

    // Branch name requested by user.
    // If this branch is a float in the file, ROOT will convert to double here.
    Float_t pot = 0.0;
    ULong64_t gates_on = 0.0;
    Double_t gates_off = 0.0;

    if (t->SetBranchAddress("rec.hdr.pot", &pot) < 0) {
        myfile->Close();
        return std::vector<double>{pot_sum, gates_on_sum, gates_off_sum};
    }
    if (t->SetBranchAddress("rec.hdr.nnumiinfo", &gates_on) < 0) {
        myfile->Close();
        return std::vector<double>{pot_sum, gates_on_sum, gates_off_sum};
    }
    if (t->SetBranchAddress("rec.hdr.noffbeamnumi", &gates_off) < 0) {
        myfile->Close();
        return std::vector<double>{pot_sum, gates_on_sum, gates_off_sum};
    }

    const Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        pot_sum += pot;
        gates_on_sum += gates_on;
        gates_off_sum += gates_off;
    }

    myfile->Close();
    return std::vector<double>{pot_sum, gates_on_sum, gates_off_sum};
}
'''
)


def getPOT(dataDir):
    total_pot = 0.0
    total_gates_on = 0.0
    total_gates_off = 0.0

    for filename in sorted(os.listdir(dataDir)):
        if not filename.endswith('.root'):
            continue

        fullpath = os.path.join(dataDir, filename)
        pot, gates_on, gates_off = r.checkPOT(fullpath)

        total_pot += float(pot)
        total_gates_on += float(gates_on)
        total_gates_off += float(gates_off)
        print(f"{filename}: POT={float(pot):.6g}, Gates on: {int(gates_on)}, Gates off: {int(gates_off)}")

    print(f"TOTAL: POT={total_pot:.6g}, Gates On: {total_gates_on}, Gates off: {total_gates_off}")


if __name__ == "__main__":
    
    #dataDir = '/pnfs/icarus/persistent/users/dcarber/spine/combined_files/NuMI_data_flat_cafs_2/' #beam on
    dataDir = "/pnfs/icarus/persistent/users/dcarber/spine/NuMI_CV_ext/v09_89_01_02p02_2/" #beam on, CV ext
    #dataDir = '/pnfs/icarus/persistent/users/dcarber/spine/combined_files/NuMI_offbeam_data_flat_cafs/' #beam off
    getPOT(dataDir)
