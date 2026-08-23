// Batched ROOT-file integrity check used by `campaign.py scan`.
//
// Reads a list of file paths (one per line) from filelist_path, opens each
// with TFile::Open, and writes the full path of any file that fails to
// open or reports IsZombie() to report_path (one per line). Running this
// as a single ROOT process over many files avoids the per-file ROOT
// startup cost (~1-2s) that a one-process-per-file approach would incur.
void scan_check(const char* filelist_path, const char* report_path)
{
    std::ifstream in(filelist_path);
    std::ofstream out(report_path);
    std::string line;
    while(std::getline(in, line))
    {
        if(line.empty()) continue;
        TFile *f = TFile::Open(line.c_str());
        if(!f || f->IsZombie())
            out << line << "\n";
        if(f) f->Close();
    }
}
