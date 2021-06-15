// The Hist2D amplitude generates events based on a 2D histogram, which would typically be created from data to simulate a realistic background distribution.  This is a simple macro to write an example formatted histogram for the Hist2D amplitude.

void exampleHist2D() {

	// simple model of mass spectrum with gaussian centered at 1 GeV
	TF1 *fMass = new TF1("mass","gaus", 0.0, 2.0);
	fMass->SetParameters(100, 1.0, 0.1);

	// simple model of linear cross section energy dependence 
	TF1 *fEgamma = new TF1("brem","[0] - [1]*x", 3.0, 12.0);
	fEgamma->SetParameters(1.0, 0.01);

	// fill histogram with random distribution of Mass vs E_gamma
	TH2F *h2 = new TH2F("MVsE", "; E_{#gamma} (GeV); Mass (GeV)", 120, 0.0, 12.0, 200, 0.0, 2.0);	
	for(int i=0; i<100000; i++) 
		h2->Fill(fEgamma->GetRandom(), fMass->GetRandom());
	h2->Draw("colz");

	TFile *fout = TFile::Open("exampleHist2D.root","recreate");
	h2->Write();
	fout->Close();

	return;
}
