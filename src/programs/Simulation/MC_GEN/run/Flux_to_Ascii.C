#include<fstream>

void Flux_to_Ascii(TString flux_path)// = "flux_SP17_30274_31057.root")
{
	TString flux_name = "tagged_flux";
	TFile* f = new TFile(flux_path);
	histo_flux = (TH1F*) f->Get(flux_name);

	TString outName = flux_path;
	outName.Replace(outName.Length()-5, outName.Length()-1, ".ascii");

	ofstream outfile;
	outfile.open(outName, ios::out);

	for(int i = 1; i< histo_flux->GetNbinsX()+1;i++)
	{
		outfile << histo_flux->GetBinCenter(i) << " " << histo_flux->GetBinContent(i) << endl;
	}
	outfile.close();
}
