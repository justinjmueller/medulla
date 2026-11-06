// compute_corrected_pot.C
//
// Scans combined flat CAFs and identifies "failed subruns" where
// SPINE produced zero truth interactions (ndlp_true == 0 for all events)
// despite GENIE neutrinos being present (mc.nnu > 0).
//
// Subtracts the failed-subrun POT from the total to get the corrected POT.
//
// Usage:
//   root -l -b -q 'compute_corrected_pot_list.C("<dir or list.txt>")' 2>&1 | tee corrected_pot.log
//
// List variant of compute_corrected_pot.C (same logic). If the argument ends in .txt it is
// read as one file path or root:// URL per line (for nodes without /pnfs mounted, and for
// samples whose files are not named *.flat.root); otherwise it lists <dir>/*.flat.root as
// the original does. A final "RESULT key=value ..." line is printed for scripts to parse.

#include <iostream>
#include <vector>

void compute_corrected_pot_list(std::string dir) {

    //const char* dir = "/pnfs/icarus/persistent/users/dcarber/spine/NuMI_CV_ext/v09_89_01_02p02_2/combined_files";
    const int N_FILES = 0;  // 0 = all

    const double POT_DOCUMENTED = 7.78521e20;
    const double RATE_OLD = 23452.0 / 1.47114935e20;
 
    std::cout << "Listing files from: " << dir << std::endl;
    std::cout << "  (this may take a minute over /pnfs...)" << std::endl;

    TString file_list;
    if (TString(dir.c_str()).EndsWith(".txt")) {
        TString cmd = TString::Format("grep -v '^#' %s | grep -v '^ *$'", dir.c_str());
        file_list = gSystem->GetFromPipe(cmd.Data());
    } else {
        TString cmd = TString::Format("ls %s/*.flat.root 2>/dev/null", dir.c_str());
        file_list = gSystem->GetFromPipe(cmd.Data());
    }
    TObjArray* files = file_list.Tokenize("\n");
    int n_total_files = files->GetEntries();
    std::cout << "File listing complete." << std::endl;
    std::cout << "Found " << n_total_files << " files" << std::endl;

    int n_to_process = (N_FILES > 0 && N_FILES < n_total_files) ? N_FILES : n_total_files;
    std::cout << "Processing " << n_to_process << " files" << std::endl;

    double total_pot = 0.0;
    double good_pot = 0.0;
    double failed_pot = 0.0;
    double total_pot_th1d = 0.0;
    long long total_events = 0;
    long long total_failed_events = 0;
    long long total_subruns = 0;
    long long total_failed_subruns = 0;
    int n_files_ok = 0;
    int n_files_fail = 0;

    for (int i = 0; i < n_to_process; ++i) {
        TString fname = ((TObjString*)files->At(i))->GetString();
        TFile* f = TFile::Open(fname, "READ");
        if (!f || f->IsZombie()) {
            std::cout << "FAILED_FILE " << fname << std::endl;
            n_files_fail++;
            if (f) delete f;
            continue;
        }

        // Grab TotalPOT from histogram
        TH1D* h_pot = (TH1D*)f->Get("TotalPOT");
        double file_pot_th1d = h_pot ? h_pot->Integral() : 0.0;
        total_pot_th1d += file_pot_th1d;

        TTree* tree = (TTree*)f->Get("recTree");
        if (!tree) {
            std::cout << "FAILED_FILE " << fname << std::endl;
            n_files_fail++;
            f->Close();
            delete f;
            continue;
        }

        int n_entries = tree->GetEntries();
        total_events += n_entries;

        // KEY OPTIMIZATION: disable all branches, enable only the 4 we need
        tree->SetBranchStatus("*", 0);
        tree->SetBranchStatus("rec.hdr.pot", 1);
        tree->SetBranchStatus("rec.hdr.first_in_subrun", 1);
        tree->SetBranchStatus("rec.mc.nnu", 1);
        tree->SetBranchStatus("rec.ndlp_true", 1);

        float pot_val;
        Char_t first_in_subrun;
        unsigned long long mc_nnu;
        int ndlp_true;

        tree->SetBranchAddress("rec.hdr.pot", &pot_val);
        tree->SetBranchAddress("rec.hdr.first_in_subrun", &first_in_subrun);
        tree->SetBranchAddress("rec.mc.nnu", &mc_nnu);
        tree->SetBranchAddress("rec.ndlp_true", &ndlp_true);

        // Read only the 4 branches into vectors
        std::vector<float> v_pot(n_entries);
        std::vector<Char_t> v_fis(n_entries);
        std::vector<unsigned long long> v_nnu(n_entries);
        std::vector<int> v_ndlp(n_entries);

        for (int j = 0; j < n_entries; ++j) {
            if (j % 100000 == 0) std::cout << "    file " << i+1 << " entry " << j << "/" << n_entries << std::endl;
            tree->GetEntry(j);
            v_pot[j] = pot_val;
            v_fis[j] = first_in_subrun;
            v_nnu[j] = mc_nnu;
            v_ndlp[j] = ndlp_true;
        }

        // Find subrun start indices
        std::vector<int> subrun_starts;
        for (int j = 0; j < n_entries; ++j) {
            if (v_fis[j]) subrun_starts.push_back(j);
        }

        int n_subruns_file = subrun_starts.size();
        total_subruns += n_subruns_file;

        if (n_subruns_file == 0) {
            double fp = 0;
            for (int j = 0; j < n_entries; ++j) fp += v_pot[j];
            total_pot += fp;

            bool all_zero = true;
            bool has_nu = false;
            for (int j = 0; j < n_entries; ++j) {
                if (v_ndlp[j] > 0) all_zero = false;
                if (v_nnu[j] > 0) has_nu = true;
            }
            if (all_zero && has_nu) {
                failed_pot += fp;
                total_failed_subruns++;
                total_failed_events += n_entries;
                total_subruns++;
            } else {
                good_pot += fp;
                total_subruns++;
            }
        } else {
            for (int s = 0; s < n_subruns_file; ++s) {
                int start = subrun_starts[s];
                int end = (s + 1 < n_subruns_file) ? subrun_starts[s+1] : n_entries;
                int subrun_size = end - start;

                double subrun_pot = (double)v_pot[start];
                total_pot += subrun_pot;

                bool all_zero = true;
                bool has_nu = false;
                for (int j = start; j < end; ++j) {
                    if (v_ndlp[j] > 0) all_zero = false;
                    if (v_nnu[j] > 0) has_nu = true;
                    if (!all_zero) break;
                }

                if (all_zero && has_nu) {
                    failed_pot += subrun_pot;
                    total_failed_subruns++;
                    total_failed_events += subrun_size;
                } else {
                    good_pot += subrun_pot;
                }
            }
        }

        n_files_ok++;
        f->Close();
        delete f;

        double frac = (total_pot > 0) ? 100.0 * failed_pot / total_pot : 0;
        std::cout << Form("  [%2d/%d] %s  entries=%d  subruns=%d  failed=%lld  cumul_failed_POT=%.2f%%",
                          i+1, n_to_process,
                          gSystem->BaseName(fname.Data()),
                          n_entries, n_subruns_file,
                          total_failed_subruns, frac) << std::endl;
    }

    // Results
    double corrected_pot = good_pot;

    std::cout << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << "CORRECTED POT RESULTS" << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << Form("Files processed:            %d", n_files_ok) << std::endl;
    std::cout << Form("Files failed:               %d", n_files_fail) << std::endl;
    std::cout << std::endl;
    std::cout << Form("Total events:               %lld", total_events) << std::endl;
    std::cout << Form("Failed events:              %lld  (%.2f%%)", total_failed_events,
                      100.0 * total_failed_events / total_events) << std::endl;
    std::cout << std::endl;
    std::cout << Form("Total subruns:              %lld", total_subruns) << std::endl;
    std::cout << Form("Failed subruns:             %lld  (%.2f%%)", total_failed_subruns,
                      100.0 * total_failed_subruns / total_subruns) << std::endl;
    std::cout << std::endl;
    std::cout << Form("POT (per-event sum):        %.6e", total_pot) << std::endl;
    std::cout << Form("POT (TH1D sum):             %.6e", total_pot_th1d) << std::endl;
    if (total_pot_th1d > 0) {
        std::cout << Form("Ratio per-event/TH1D:       %.4f", total_pot / total_pot_th1d) << std::endl;
    }
    std::cout << std::endl;
    std::cout << Form("Failed subrun POT:          %.6e  (%.2f%%)", failed_pot,
                      100.0 * failed_pot / total_pot) << std::endl;
    std::cout << Form("Good subrun POT:            %.6e  (%.2f%%)", good_pot,
                      100.0 * good_pot / total_pot) << std::endl;
    std::cout << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << Form("CORRECTED POT = %.15e", corrected_pot) << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << std::endl;
    std::cout << Form("Documented POT:             %.6e", POT_DOCUMENTED) << std::endl;
    std::cout << Form("Corrected POT:              %.6e", corrected_pot) << std::endl;
    std::cout << Form("Ratio corrected/documented: %.4f", corrected_pot / POT_DOCUMENTED) << std::endl;
    std::cout << std::endl;

    double rate_new_doc = 113768.0 / POT_DOCUMENTED;
    double rate_new_corr = 113768.0 / corrected_pot;

    std::cout << "Cross-check against signal-tree deficit:" << std::endl;
    std::cout << Form("  Old rate:              %.4e events/POT", RATE_OLD) << std::endl;
    std::cout << Form("  New rate (doc POT):    %.4e events/POT  (ratio = 0.9167)", rate_new_doc) << std::endl;
    std::cout << Form("  New rate (corr POT):   %.4e events/POT  (ratio = %.4f)",
                      rate_new_corr, rate_new_corr / RATE_OLD) << std::endl;
    std::cout << std::endl;
    std::cout << "  If ratio ~ 1.000 -> corrected POT fully explains the deficit" << std::endl;
    std::cout << "  If ratio < 1.000 -> residual deficit remains" << std::endl;
    std::cout << std::endl;
    std::cout << Form("RESULT files_listed=%d files_ok=%d files_fail=%d events=%lld failed_events=%lld "
                      "subruns=%lld failed_subruns=%lld pot_event=%.15e pot_th1d=%.15e "
                      "failed_pot=%.15e corrected_pot=%.15e",
                      n_total_files, n_files_ok, n_files_fail, total_events, total_failed_events,
                      total_subruns, total_failed_subruns, total_pot, total_pot_th1d,
                      failed_pot, corrected_pot) << std::endl;
}
