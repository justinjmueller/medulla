void to_gundam()
{ 

  /////////////////////////////////////////////////////////////
  /// Input
  /////////////////////////////////////////////////////////////

  // Selected
  string infile_sel_string = "/exp/icarus/app/users/lkashur/medulla_dev/srcs/medulla/build/output_sel.root";
  string selnutree_string = "/events/cv/selected_nu";
  string selcostree_string = "/events/cv/selected_cos";
  string potmc_string = "/events/cv/POT"; // Add to to_gundam.cc
  string livetimemc_string = "events/cv/Livetime";
  string seloffbeamtree_string = "/events/offbeam/selected_cos";
  string livetimeoffbeam_string = "events/offbeam/Livetime";
  string selonbeamtree_string = "/events/onbeam/selected_nu";
  string potonbeam_string = "/events/onbeam/POT";
  string livetimeonbeam_string = "/events/onbeam/Livetime";

  // Signal
  string infile_sig_string = "/exp/icarus/app/users/lkashur/medulla_dev/srcs/medulla/build/output_sig.root";
  string signaltree_string = "events/cv/signal";

  /////////////////////////////////////////////////////////////
  /// Output
  /////////////////////////////////////////////////////////////
  string outfile_mc_string = "mc_offbeam_syst_gundaminput.root"; // mc = CV + off-beam
  string outfile_data_string = "onbeam_syst_gundaminput.root"; // data = on-beam
  string outfile_signal_string = "signal_syst_gundaminput.root";

  /////////////////////////////////////////////////////////////
  /// Store input in ROOT RDataFrames for further processing
  /////////////////////////////////////////////////////////////
  TChain ch("selected");
  ch.Add((infile_sel_string + selnutree_string).c_str());
  ch.Add((infile_sel_string + selcostree_string).c_str());
  ch.Add((infile_sel_string + seloffbeamtree_string).c_str());
  ROOT::RDataFrame rdf_mc(ch);
  ROOT::RDataFrame rdf_data(selonbeamtree_string, infile_sel_string);
  ROOT::RDataFrame rdf_signal(signaltree_string, infile_sig_string);

  /////////////////////////////////////////////////////////////
  /// Make any necessary modifications and save output
  /////////////////////////////////////////////////////////////
  auto rdf_mc_save = rdf_mc
    .Snapshot("selected", outfile_mc_string.c_str());
  
  auto rdf_data_save = rdf_data
    .Snapshot("selected", outfile_data_string.c_str());

  auto rdf_signal_save = rdf_signal
    .Snapshot("signal", outfile_signal_string.c_str());

  // Write POT and Livetime of original samples to output file
  std::unique_ptr<TFile> infile( TFile::Open(infile_sel_string.c_str(), "READ") );
  std::unique_ptr<TFile> outfile_mc( TFile::Open(outfile_mc_string.c_str(), "UPDATE") ); // "UPDATE" because we created these files above
  std::unique_ptr<TFile> outfile_data( TFile::Open(outfile_data_string.c_str(), "UPDATE") );

  TH1D *POT_mc = (TH1D*)infile->Get(potmc_string.c_str());
  TH1D *POT_mc_clone = (TH1D*)POT_mc->Clone("POT_mc");
  TH1D *POT_onbeam = (TH1D*)infile->Get(potonbeam_string.c_str());
  TH1D *POT_onbeam_clone = (TH1D*)POT_onbeam->Clone("POT_onbeam");
  TH1D *Livetime_mc = (TH1D*)infile->Get(livetimemc_string.c_str());
  TH1D *Livetime_mc_clone = (TH1D*)Livetime_mc->Clone("Livetime_mc");
  TH1D *Livetime_offbeam = (TH1D*)infile->Get(livetimeoffbeam_string.c_str());
  TH1D *Livetime_offbeam_clone = (TH1D*)Livetime_offbeam->Clone("Livetime_offbeam");
  TH1D *Livetime_onbeam = (TH1D*)infile->Get(livetimeonbeam_string.c_str());
  TH1D *Livetime_onbeam_clone = (TH1D*)Livetime_onbeam->Clone("Livetime_onbeam");

  outfile_mc->cd();
  POT_mc_clone->Write();
  Livetime_mc_clone->Write();
  Livetime_offbeam_clone->Write();

  outfile_data->cd();
  POT_onbeam_clone->Write();
  Livetime_onbeam_clone->Write();
  
}
