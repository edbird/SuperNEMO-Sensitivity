

void graphdraw()
{

    TFile *fgr = new TFile("fgraph.root");
    TGraph *g_chi2_ch0 = (TGraph*)fgr->Get("g_ch0");
    TGraph *g_chi2_ch1 = (TGraph*)fgr->Get("g_ch1");
    TGraph *g_chi2_ch10 = (TGraph*)fgr->Get("g_ch10");

    TCanvas *canvas = new TCanvas("c", "c");
    canvas->SetLogy();
    g_chi2_ch0->GetXaxis()->SetTitle("Number of bins");
    g_chi2_ch0->GetYaxis()->SetTitle("#chi^{2}(SSD, HSD)");
    g_chi2_ch0->GetXaxis()->SetTitleFont(43);
    g_chi2_ch0->GetYaxis()->SetTitleFont(43);
    g_chi2_ch0->GetXaxis()->SetTitleSize(20);
    g_chi2_ch0->GetYaxis()->SetTitleSize(20);
    g_chi2_ch0->GetXaxis()->SetLabelFont(43);
    g_chi2_ch0->GetYaxis()->SetLabelFont(43);
    g_chi2_ch0->GetXaxis()->SetLabelSize(20);
    g_chi2_ch0->GetYaxis()->SetLabelSize(20);
    g_chi2_ch0->SetMinimum(1.0e-1);
    g_chi2_ch0->SetMaximum(1.0e+4);
    g_chi2_ch0->SetTitle(0);
    g_chi2_ch0->Draw("al");
    g_chi2_ch1->Draw("lsame");
    //g_chi2_ch10->Draw("lsame");

}
