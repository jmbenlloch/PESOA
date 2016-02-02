void histos_nfirst(){

	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);

    TFile *fIn1 = new TFile("/home/jmbenlloch/next/petalo/work/histo/test_nfirst_histos.root", "read");
    TH1F *h1 = (TH1F*) fIn1.Get("NFirst.Time_1_PEBox1");
    TH1F *h2 = (TH1F*) fIn1.Get("NFirst.Time_2_PEBox1");
    TH1F *h3 = (TH1F*) fIn1.Get("NFirst.Time_3_PEBox1");
    TH1F *h4 = (TH1F*) fIn1.Get("NFirst.Time_4_PEBox1");
    TH1F *h5 = (TH1F*) fIn1.Get("NFirst.Time_5_PEBox1");
    TH1F *h6 = (TH1F*) fIn1.Get("NFirst.Time_6_PEBox1");
    TH1F *h7 = (TH1F*) fIn1.Get("NFirst.Time_7_PEBox1");
    TH1F *h8 = (TH1F*) fIn1.Get("NFirst.Time_8_PEBox1");

	TCanvas *c1 = new TCanvas("c1","multipads",900,700);

/*	h1->Scale(1/h1->Integral(), "width");
	h2->Scale(1/h2->Integral(), "width");
	h3->Scale(1/h3->Integral(), "width");
	h4->Scale(1/h4->Integral(), "width");
	h5->Scale(1/h5->Integral(), "width");
	h6->Scale(1/h6->Integral(), "width");
	h7->Scale(1/h7->Integral(), "width");
	h8->Scale(1/h8->Integral(), "width");*/

    h1->SetLineColor(kRed);
    h2->SetLineColor(kBlue);
    h3->SetLineColor(kGreen);
    h4->SetLineColor(kYellow);
    h5->SetLineColor(kPink);
    h6->SetLineColor(kOrange);
    h7->SetLineColor(kCyan);
    h8->SetLineColor(kViolet);

    h1->Draw();
    h2->Draw("same");
    h3->Draw("same");
    h4->Draw("same");
//    h5->Draw("same");
//    h6->Draw("same");
//    h7->Draw("same");
//    h8->Draw("same");

	c1->Print("/home/jmbenlloch/next/petalo/work/histo/nfirst.pdf");
	c1->Print("/home/jmbenlloch/next/petalo/work/histo/nfirst.png");
}
