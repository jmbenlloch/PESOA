void histos_nFixed_noCher(){

	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);

    TFile *fIn1 = new TFile("/home/jmbenlloch/next/petalo/work/histo/lXe_refl97_VUV3mm_xy2.4cm_z5cm_n1.7_noCher_QE_1_SPTR_0_ASIC_0_DT300_histos.root", "read");
    TFile *fIn2 = new TFile("/home/jmbenlloch/next/petalo/work/histo/lyso_refl97_VUV3mm_xy2.4cm_z5cm_n1.8_noCher_QE_1_SPTR_0_ASIC_0_DT300_histos.root", "read");

    TH1F *h1 = (TH1F*) fIn1.Get("DTOF.DTOF3");
    TH1F *h2 = (TH1F*) fIn2.Get("DTOF.DTOF3");

	TF1* gauF1 = new TF1("gauF1","gaus",-100,100);
	gauF1->SetLineColor(kBlue);
	gauF1->SetLineWidth(1);
	TF1* gauF2 = new TF1("gauF2","gaus",-100,100);
	gauF2->SetLineColor(kRed);
	gauF2->SetLineWidth(1);

	TCanvas *c1 = new TCanvas("c1","multipads",900,700);

	h1->Scale(1/h1->Integral(), "width");
	h2->Scale(1/h2->Integral(), "width");

	h1->Fit("gauF1","","e",-10,10);
	h2->Fit("gauF2","","e",-25,25);

	h1->SetTitle("DTOF (n fixed, Scint)");
    h1->Draw();
    h2->SetLineColor(kRed);
    h2->Draw("same");

	std::cout << "sigma: "  << h1->GetFunction("gauF1")->GetParameter(2) << std::endl;

	TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
	leg->SetFillColor(0);
	//leg->SetHeader("test legend");
	leg->AddEntry(h1, "LXe", "lp");
	leg->AddEntry(h2, "LYSO", "lp");
	leg->AddEntry(h1, Form("Sigma LXe %.3g", h1->GetFunction("gauF1")->GetParameter(2)), "lp");
	leg->AddEntry(h2, Form("Sigma LYSO %.3g", h2->GetFunction("gauF2")->GetParameter(2)), "lp");
	leg->Draw("same");

	c1->Print("/home/jmbenlloch/next/petalo/work/histo/crt_nFixed_noCher.pdf");
	c1->Print("/home/jmbenlloch/next/petalo/work/histo/crt_nFixed_noCher.png");
}
