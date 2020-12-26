#ifndef DRAW_H
#define DRAW_H


void draw_hist_ratio_pull(
    const int channel,
    TH1D* h_data,
    TH1D* h_model,
    std::vector<double> *chi2_sys,
    const TString &filename,
    const TString &filedir,
    const double chi2)
{

    const Int_t kGrey = kGray;

    // set error, data
//    for(Int_t i = 1; i <= h_data->GetNbinsX(); ++ i)
//    {
//        const double content = h_data->GetBinContent(i);
//        if(content >= 0.0)
//        {
//            h_data->SetBinError(i, std::sqrt(content));
//        }
//        else
//        {
//            h_data->SetBinError(i, 0.0);
//        }
//    }

    // set error, model
    for(Int_t i = 1; i <= h_model->GetNbinsX(); ++ i)
    {
        const double content = h_model->GetBinContent(i);
        const double content_data = h_data->GetBinContent(i);
        double total_error = 0.0;
        double total_data_error = 0.0;
        if(content >= 0.0)
        {
            //total_error = std::sqrt(content);
            total_error = content;
        }
        if(content_data >= 0.0)
        {
            //total_data_error = std::sqrt(content_data);
            total_data_error = content_data;
        }
        if(chi2_sys != nullptr)
        {
            //total_error += std::sqrt(chi2_sys->at(i - 1));
            //total_error += std::sqrt(chi2_sys->at(i - 1));
            total_error += chi2_sys->at(i - 1);
            //total_data_error += std::sqrt(chi2_sys->at(i - 1));
            //total_data_error += std::sqrt(chi2_sys->at(i - 1));
            total_data_error += chi2_sys->at(i - 1);
        }
        total_error = std::sqrt(total_error);
        total_data_error = std::sqrt(total_data_error);
        h_model->SetBinError(i, total_error);
        h_data->SetBinError(i, total_data_error);
    }

    double chi2_other = 0.0;
    for(Int_t i = 1; i <= h_model->GetNbinsX(); ++ i)
    {
        if(h_model->GetBinContent(i) <= 0.0) continue;
        double e = h_model->GetBinError(i);
        double data = h_data->GetBinContent(i);
        double model = h_model->GetBinContent(i);
        double d = data - model;
        double chi = std::pow(d / e, 2.0);
        chi2_other += chi;
        /*
        std::cout << "bin: i=" << i 
                  << " chi2=" << chi
                  << " = (( " << data
                  << " - " << model
                  << " ) / " << e
                  << " )^2 =>> error composition: "
                  << std::sqrt(model) * std::sqrt(model) << " + "
                  << chi2_sys->at(i - 1)
                  << " | model="  << model
                  << " delta=" << d
                  << std::endl;
        */
    }
    std::cout << __func__ << " chi2_other=" << chi2_other << std::endl;

    const int n_bins_x = h_data->GetNbinsX();
    const double low = h_data->GetXaxis()->GetBinLowEdge(1);
    const double high = h_data->GetXaxis()->GetBinUpEdge(n_bins_x);
    TH1D *h_ratio = new TH1D("h_ratio", "h_ratio", n_bins_x, low, high);
    TH1D *h_pull = new TH1D("h_pull", "h_pull", n_bins_x, low, high);

    for(Int_t i = 1; i <= h_data->GetNbinsX(); ++ i)
    {
        const double model = h_model->GetBinContent(i);
        const double model_error = h_model->GetBinError(i);
        const double data = h_data->GetBinContent(i);
        double ratio = 1.0;
        double ratio_error = 0.0;
        double pull = 0.0;
        if(model != 0.0)
        {
            ratio = data / model;
            ratio_error = model_error * data / (model * model);
        }
        h_ratio->SetBinContent(i, ratio);
        h_ratio->SetBinError(i, ratio_error);
        if(model > 0.0)
        {
            pull = (data - model) / model_error;
        }
        h_pull->SetBinContent(i, pull);
    }
    

    ///////////////////////////////////////////////////////////////////////////
    // pad min max
    ///////////////////////////////////////////////////////////////////////////

    std::map<int, double> axismaxmap;
    for(int i = 1; i <= 100; ++ i)
    {
        if(channel == 0)
        {
            if(i == 1)
            {
                axismaxmap[i] = 1.2e+06;
            }
            else if(i == 2)
            {
                axismaxmap[i] = 1.0e+06;
            }
            else if(i == 3)
            {
                axismaxmap[i] = 800.0e+03;
            }
            else if(i == 4)
            {
                axismaxmap[i] = 700.0e+03;
            }
            else if(i <= 6)
            {
                axismaxmap[i] = 500.0e+03;
            }
            else if(i == 7)
            {
                axismaxmap[i] = 400.0e+03;
            }
            else if(i == 8)
            {
                axismaxmap[i] = 350.0e+03;
            }
            else if(i <= 10)
            {
                axismaxmap[i] = 300.0e+03;
            }
            else if(i <= 13)
            {
                axismaxmap[i] = 250.0e+03;
            }
            else if(i <= 17)
            {
                axismaxmap[i] = 200.0e+03;
            }
            else if(i <= 24)
            {
                axismaxmap[i] = 150.0e+03;
            }
            else if(i <= 28)
            {
                axismaxmap[i] = 120.0e+03;
            }
            else if(i <= 60)
            {
                axismaxmap[i] = 100.0e+03;
            }
            else
            {
                axismaxmap[i] = 50.0e+03;
            }
        }
        else if(channel == 1)
        {
            if(i == 1)
            {
                axismaxmap[i] = 2.5e+06;
            }
            else if(i == 2)
            {
                axismaxmap[i] = 2.5e+06;
            }
            else if(i == 3)
            {
                axismaxmap[i] = 2.0e+06;
            }
            else if(i == 4)
            {
                axismaxmap[i] = 1800.0e+03;
            }
            else if(i == 5)
            {
                axismaxmap[i] = 1500.0e+03;
            }
            else if(i == 6)
            {
                axismaxmap[i] = 1200.0e+03;
            }
            else if(i == 7)
            {
                axismaxmap[i] = 1200.0e+03;
            }
            else if(i <= 9)
            {
                axismaxmap[i] = 1000.0e+03;
            }
            else if(i <= 14)
            {
                axismaxmap[i] = 800.0e+03;
            }
            else if(i <= 19)
            {
                axismaxmap[i] = 500.0e+03;
            }
            else if(i <= 25)
            {
                axismaxmap[i] = 400.0e+03;
            }
            else if(i <= 30)
            {
                axismaxmap[i] = 300.0e+03;
            }
            else
            {
                axismaxmap[i] = 200.0e+03;
            }
        }
        else
        {
            axismaxmap[i] = 1.0e+06;
        }
    }

    double PAD_L_Y_MAX = 1.0;
    double PAD_L_Y_MIN = -1.0;
    if(channel == 0)
    {
        PAD_L_Y_MAX = 1.5;
        PAD_L_Y_MIN = -1.5;
    }
    else
    {
        PAD_L_Y_MAX = 5.0;
        PAD_L_Y_MIN = -5.0;
    }
    double PAD_M_Y_MAX = 1.5;
    double PAD_M_Y_MIN = 0.5;
    if(channel == 1)
    {
        PAD_M_Y_MAX = 4.0;
        PAD_M_Y_MIN = 0.0;
    }
    const double PAD_U_Y_MAX = axismaxmap[h_model->GetNbinsX()];
    const double PAD_U_Y_MIN = 0.0;


    ///////////////////////////////////////////////////////////////////////////
    // set drawing properties
    ///////////////////////////////////////////////////////////////////////////

    h_model->SetTitle(0);
    h_ratio->SetTitle(0);
    h_pull->SetTitle(0);

    h_model->SetStats(0);
    h_ratio->SetStats(0);
    h_pull->SetStats(0);


    // tick size

    h_model->GetXaxis()->SetTickSize(0.04);
    h_ratio->GetXaxis()->SetTickSize(0.10);
    h_pull->GetXaxis()->SetTickSize(0.08);

    // x label / x title

    h_model->GetXaxis()->SetLabelSize(0.0);
    h_model->GetXaxis()->SetTitleSize(0.0);
    h_ratio->GetXaxis()->SetLabelSize(0.0);
    h_ratio->GetXaxis()->SetTitleSize(0.0);

    // y ticks / label

    h_model->GetYaxis()->SetTickSize(0.0);
    h_model->GetYaxis()->SetLabelSize(0.0);

    h_ratio->GetYaxis()->SetTickSize(0.0);
    h_ratio->GetYaxis()->SetLabelSize(0.0);

    h_pull->GetYaxis()->SetTickSize(0.0);
    h_pull->GetYaxis()->SetLabelSize(0.0);

    // y title

    h_model->GetYaxis()->SetTitleFont(43);
    h_model->GetYaxis()->SetTitleSize(20);
    TString ystring;
    int binw = (int)std::round(1000.0 * h_model->GetXaxis()->GetBinWidth(1));
    ystring.Form("Events / %d keV", binw);
    h_model->GetYaxis()->SetTitle(ystring);
    h_model->GetYaxis()->SetTitleOffset(1.0);

    h_ratio->GetYaxis()->SetTitleFont(43);
    h_ratio->GetYaxis()->SetTitleSize(20);
    h_ratio->GetYaxis()->SetTitle("SSD / HSD");
    h_ratio->GetYaxis()->SetTitleOffset(1.0);

    h_pull->GetYaxis()->SetTitleFont(43);
    h_pull->GetYaxis()->SetTitleSize(20);
    h_pull->GetYaxis()->SetTitle("Pull (#sigma)     ");
    h_pull->GetYaxis()->SetTitleOffset(1.0);

    // x title
    h_pull->GetXaxis()->SetTitleFont(43);
    h_pull->GetXaxis()->SetTitleSize(20);
    h_pull->GetXaxis()->SetTitle("T_{e} [MeV]");
    h_pull->GetXaxis()->SetTitleOffset(3.0);

    h_data->SetLineColor(kMagenta);
    h_model->SetLineColor(kBlue);
    h_data->SetLineWidth(2);
    h_model->SetLineWidth(2);
//    h_model->SetFillColor(kBlue);
//    h_model->SetFillStyle(3002);



    h_ratio->SetLineWidth(2);
    h_ratio->SetLineColor(kBlack);


    h_pull->SetLineColor(kGray + 2);
    h_pull->SetLineWidth(2);
    h_pull->SetFillColor(kGrey + 2);
    h_pull->SetFillStyle(3001);


    // maximum / minimum

    h_model->SetMaximum(PAD_U_Y_MAX);
    h_model->SetMinimum(PAD_U_Y_MIN);

    h_ratio->SetMaximum(PAD_M_Y_MAX);
    h_ratio->SetMinimum(PAD_M_Y_MIN);
    
    h_pull->SetMinimum(PAD_L_Y_MIN);
    h_pull->SetMaximum(PAD_L_Y_MAX);



    ///////////////////////////////////////////////////////////////////////////
    // create pads
    ///////////////////////////////////////////////////////////////////////////

    TCanvas *canvas = new TCanvas(filename, filename);
    canvas->SetFillColor(kWhite);
    
    TPad *pad0 = new TPad("pad0", "pad0", 0.0, 0.5, 1.0, 1.0);
    pad0->SetBottomMargin(0.0);
    pad0->SetTicks(2, 2);
    pad0->SetRightMargin(0.05);
    pad0->Draw();

    canvas->cd();
    TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 0.5);
    pad1->SetTopMargin(0.0);
    pad1->SetBottomMargin(0.0);
    pad1->SetTicks(2, 2);
    pad1->SetRightMargin(0.05);
    pad1->Draw();

    canvas->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.3);
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(1.0 / 3.0);
    pad2->SetTicks(2, 2);
    pad2->SetRightMargin(0.05);
    pad2->Draw();

    ///////////////////////////////////////////////////////////////////////////
    // clone model
    ///////////////////////////////////////////////////////////////////////////
    TH1D *h_model_u = (TH1D*)h_model->Clone("h_model_u");
    TH1D *h_model_m = (TH1D*)h_model->Clone("h_model_m");
    TH1D *h_model_l = (TH1D*)h_model->Clone("h_model_l");
    for(Int_t i = 1; i <= h_model->GetNbinsX(); ++ i)
    {
        double content = h_model->GetBinContent(i);
        double error = h_model->GetBinError(i);
        double content_u = content + error;
        double content_l = content - error;
        h_model_u->SetBinContent(i, content_u);
        h_model_l->SetBinContent(i, content_l);
    }

    h_model_u->SetLineWidth(1);
    h_model_l->SetLineWidth(1);
//    h_model_l->SetFillColor(kWhite);
//    h_model_l->SetFillStyle(1001);
    h_model->SetFillStyle(3002);
    h_model->SetFillColor(kBlue);

    pad0->cd();
    h_model->Draw("axis");
    h_model->Draw("E2same");
    h_model_u->Draw("histsame");
    h_model_m->Draw("histPsame");
    h_model_l->Draw("histsame");
    h_data->Draw("Esame");

    pad1->cd();
    TLine *l_ratio_1 = new TLine(low, 1.0, high, 1.0);
    l_ratio_1->SetLineWidth(2);
    l_ratio_1->SetLineColor(kGrey + 2);
    h_ratio->Draw("axis");
    l_ratio_1->Draw();
    h_ratio->Draw("Esame");

    pad2->cd();
    TLine *l_pull_0 = new TLine(low, 0.0, high, 0.0);
    l_pull_0->SetLineWidth(2);
    l_pull_0->SetLineColor(kBlack);
    h_pull->Draw("axis");
    l_pull_0->Draw();
    h_pull->Draw("histsame");



    
    ///////////////////////////////////////////////////////////////////////////
    // axis
    ///////////////////////////////////////////////////////////////////////////

    TGaxis *axis_u_y_1 = new TGaxis(0.0, PAD_U_Y_MIN, 0.0, PAD_U_Y_MAX, PAD_U_Y_MIN, PAD_U_Y_MAX, 310, "S-");
    TGaxis *axis_u_y_2 = new TGaxis(H_E_MAX, PAD_U_Y_MIN, H_E_MAX, PAD_U_Y_MAX, PAD_U_Y_MIN, PAD_U_Y_MAX, 310, "SL+");

    TGaxis *axis_m_y_1 = new TGaxis(0.0, PAD_M_Y_MIN, 0.0, PAD_M_Y_MAX, PAD_M_Y_MIN, PAD_M_Y_MAX, 205, "S-");
    TGaxis *axis_m_y_2 = new TGaxis(H_E_MAX, PAD_M_Y_MIN, H_E_MAX, PAD_M_Y_MAX, PAD_M_Y_MIN, PAD_M_Y_MAX, 205, "SL+");
    
    TGaxis *axis_l_y_1 = new TGaxis(0.0, PAD_L_Y_MIN, 0.0, PAD_L_Y_MAX, PAD_L_Y_MIN, PAD_L_Y_MAX, 205, "S-");
    TGaxis *axis_l_y_2 = new TGaxis(H_E_MAX, PAD_L_Y_MIN, H_E_MAX, PAD_L_Y_MAX, PAD_L_Y_MIN, PAD_L_Y_MAX, 205, "SL+");

    axis_u_y_1->SetLabelFont(43);
    axis_u_y_2->SetLabelFont(43);
    axis_m_y_1->SetLabelFont(43);
    axis_m_y_2->SetLabelFont(43);
    axis_l_y_1->SetLabelFont(43);
    axis_l_y_2->SetLabelFont(43);

    axis_u_y_1->SetLabelSize(15);
    axis_u_y_2->SetLabelSize(15);
    axis_m_y_1->SetLabelSize(15);
    axis_m_y_2->SetLabelSize(15);
    axis_l_y_1->SetLabelSize(15);
    axis_l_y_2->SetLabelSize(15);

/*
    axis_u_y_1->SetTickSize(0.04);
    axis_u_y_2->SetTickSize(0.04);
    axis_m_y_1->SetTickSize(0.10);
    axis_m_y_2->SetTickSize(0.10);
    axis_l_y_1->SetTickSize(0.08);
    axis_l_y_2->SetTickSize(0.08);
*/

    axis_u_y_1->ChangeLabel(1, 0, 0);
    axis_u_y_2->ChangeLabel(1, 0, 0);
    axis_u_y_2->ChangeLabel(-1, 0, 0);
    if(channel == 0)
    {
        axis_m_y_1->ChangeLabel(-1, 0, 0);
        axis_m_y_2->ChangeLabel(-1, 0, 0);
        axis_m_y_1->ChangeLabel(1, 0, 0);
        axis_m_y_2->ChangeLabel(1, 0, 0);
    }
//    axis_l_y_1->ChangeLabel(-1, 0, 0);
//    axis_l_y_2->ChangeLabel(-1, 0, 0);

    pad0->cd();
    axis_u_y_1->Draw();
    axis_u_y_2->Draw();

    pad1->cd();
    axis_m_y_1->Draw();
    axis_m_y_2->Draw();

    pad2->cd();
    axis_l_y_1->Draw();
    axis_l_y_2->Draw();


    ///////////////////////////////////////////////////////////////////////////
    // legend
    ///////////////////////////////////////////////////////////////////////////

    pad0->cd();
    
    double lposx = 0.625;
    double lposy = 0.65;
    double lsizex = 0.30;
    double lsizey = 0.30;
    TLegend *leg = new TLegend(lposx, lposy, lposx + lsizex, lposy + lsizey);
    TString text_data;
    TString text_model;
    text_data.Form("^{82}Se (HSD) [#xi_{31}=0?]");
    text_model.Form("^{82}Se (SSD) [#xi_{31}=???]");
    leg->AddEntry(h_data, text_data, "F");
    leg->AddEntry(h_model, text_model, "F");
    leg->SetBorderSize(5);
    leg->SetShadowColor(kGray + 2);
    leg->SetTextFont(63);
    leg->SetTextSize(18);
    leg->SetMargin(0.20);
    leg->Draw("BL");


    ///////////////////////////////////////////////////////////////////////////
    // latex labels
    ///////////////////////////////////////////////////////////////////////////

    pad0->cd();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(63);
    latex.SetTextSize(18);

    TString s0;
    s0.Form("SuperNEMO");
    latex.DrawLatex(0.15, 0.80, s0);

    TString s1;
    s1.Form("#chi^{2} = %.1f", chi2);
    latex.DrawLatex(0.15, 0.70, s1);

    
    bool saveas_pdf = true;
    bool saveas_png = true;
    if(saveas_pdf)
    {
        TString fname = filedir + "/" + filename + ".pdf";
        canvas->SaveAs(fname);
    }
    if(saveas_png)
    {
        TString fname = filedir + "/" + filename + ".png";
        canvas->SaveAs(fname);
    }

}


#endif
