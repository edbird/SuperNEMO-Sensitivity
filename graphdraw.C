

void graphdraw()
{

    TFile *fgr = new TFile("fgraph.root");
    TGraph *g_chi2_ch0 = (TGraph*)fgr->Get("g_ch0");
    TGraph *g_chi2_ch1 = (TGraph*)fgr->Get("g_ch1");
    TGraph *g_chi2_ch10 = (TGraph*)fgr->Get("g_ch10");

    TGraph *g_chi2_ch0_odz = (TGraph*)fgr->Get("g_ch0_odz");
    TGraph *g_chi2_ch1_odz = (TGraph*)fgr->Get("g_ch1_odz");
    TGraph *g_chi2_ch10_odz = (TGraph*)fgr->Get("g_ch10_odz");

    TGraph *g_chi2_ch0_sz = (TGraph*)fgr->Get("g_ch0_sz");
    TGraph *g_chi2_ch1_sz = (TGraph*)fgr->Get("g_ch1_sz");
    TGraph *g_chi2_ch10_sz = (TGraph*)fgr->Get("g_ch10_sz");

    g_chi2_ch0->SetLineColor(kRed);
    g_chi2_ch1->SetLineColor(kMagenta);
    g_chi2_ch10->SetLineColor(kBlue);

    g_chi2_ch0_odz->SetLineColor(kRed);
    g_chi2_ch1_odz->SetLineColor(kMagenta);
    g_chi2_ch10_odz->SetLineColor(kBlue);

    g_chi2_ch0_sz->SetLineColor(kRed);
    g_chi2_ch1_sz->SetLineColor(kMagenta);
    g_chi2_ch10_sz->SetLineColor(kBlue);

    g_chi2_ch0->SetLineStyle(1);
    g_chi2_ch1->SetLineStyle(1);
    g_chi2_ch10->SetLineStyle(1);

    g_chi2_ch0_odz->SetLineStyle(2);
    g_chi2_ch1_odz->SetLineStyle(2);
    g_chi2_ch10_odz->SetLineStyle(2);

    g_chi2_ch0_sz->SetLineStyle(3);
    g_chi2_ch1_sz->SetLineStyle(3);
    g_chi2_ch10_sz->SetLineStyle(3);
    



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
    g_chi2_ch10->Draw("lsame");

    g_chi2_ch0_odz->Draw("lsame");
    g_chi2_ch1_odz->Draw("lsame");
    g_chi2_ch10_odz->Draw("lsame");

    g_chi2_ch0_sz->Draw("lsame");
    g_chi2_ch1_sz->Draw("lsame");
    g_chi2_ch10_sz->Draw("lsame");


    double lposx = 0.65;
    double lposy = 0.15;
    double lsizex = 0.30;
    double lsizey = 0.50;
    TLegend *leg = new TLegend(lposx, lposy, lposx + lsizex, lposy + lsizey);
    leg->AddEntry(g_chi2_ch0, "#chi^{2} CH0", "L");
    leg->AddEntry(g_chi2_ch1, "#chi^{2} CH1", "L");
    leg->AddEntry(g_chi2_ch10, "#chi^{2} CH10", "L");
    leg->AddEntry(g_chi2_ch0_sz, "#chi^{2} CH0 Stat", "L");
    leg->AddEntry(g_chi2_ch1_sz, "#chi^{2} CH1 Stat", "L");
    leg->AddEntry(g_chi2_ch10_sz, "#chi^{2} CH10 Stat", "L");
    leg->AddEntry(g_chi2_ch0_odz, "#chi^{2} CH0 dz", "L");
    leg->AddEntry(g_chi2_ch1_odz, "#chi^{2} CH1 dz", "L");
    leg->AddEntry(g_chi2_ch10_odz, "#chi^{2} CH10 dz", "L");
    leg->SetBorderSize(5);
    leg->SetShadowColor(kGray + 2);
    leg->SetTextFont(63);
    leg->SetTextSize(20);
    leg->SetMargin(0.15);
    leg->Draw("BR");

}
