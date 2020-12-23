



const double G0_ps_integral_MeV = 0.331909E-46;
const double G0_ps_integral_yrinv = 0.159132E-17;
const double G2_ps_integral_MeV = 0.146777E-46;
const double G2_ps_integral_yrinv = 0.703714E-18;

const double T12_year = 9.39e+19;
double N0;
const double N_A = 6.022e+23;
const double mass = 6500.0; // gram
const double N_iso = 82.0;
const double Q_bb = 2.99512; // MeV
const double E_min = 0.12; // keV
const std::string isotope = "082Se";
const double efficiency = 15.0e-2; // 10 % - 20 % ~ manu
const double runtime_year = 10.0;
const double N2e = 1000000.0; // number of observed 2e events (after efficiency)

const double low_energy_threshold = 0.3; // low energy cut threshold


TRandom3 gRNG;


TH2D *h_psf_G0 = nullptr;
TH2D *h_psf_G2 = nullptr;



///////////////////////////////////////////////////////////////////////////////
//

TH1D *h_ch0_HSD = nullptr;
TH1D *h_ch0_SSD = nullptr;

TH1D *h_ch1_HSD = nullptr;
TH1D *h_ch1_SSD = nullptr;

TH2D *h_ch10_HSD = nullptr;
TH2D *h_ch10_SSD = nullptr;

const double H_E_MAX = 4.0;
Int_t N_BINS_CH0 = 40;
Int_t N_BINS_CH1 = 40;
Int_t N_BINS_CH10 = 40;

std::vector<double> V_PHYS_SYS_CH0;
std::vector<double> V_PHYS_SYS_CH1;
std::vector<double> V_PHYS_SYS_CH10;





void load_psf_data()
{

    const Int_t N_BINS_XY = (Q_bb - 1.0e-3 * E_min) / 1.0e-3;
    //const Int_t N_BINS_XY = H_E_MAX / 1.0e-3;
    std::cout << "N_BINS_XY=" << N_BINS_XY << std::endl;

    h_psf_G0 = new TH2D("h_psf_G0", "h_psf_G0",
                        N_BINS_XY, 0.0, H_E_MAX,
                        N_BINS_XY, 0.0, H_E_MAX);

    h_psf_G2 = new TH2D("h_psf_G2", "h_psf_G2",
                        N_BINS_XY, 0.0, H_E_MAX,
                        N_BINS_XY, 0.0, H_E_MAX);

    
    ///////////////////////////////////////////////////////////////////////////
    // open file
    ///////////////////////////////////////////////////////////////////////////

    std::string file_base = "/home/ecb/100Mo-150Nd/gA_theoretical_files/psf-nuclei";

    std::string fname_isotope_0 = file_base + "/"
                              + isotope + "/"
                              + "0-N0" + "/"
                              + "nEqNull.dat";
    std::string fname_isotope_1 = file_base + "/"
                              + isotope + "/"
                              + "1-N2" + "/"
                              + "nEqTwo.dat";

    ifstream if_isotope_0(fname_isotope_0);
    ifstream if_isotope_1(fname_isotope_1);

    if(if_isotope_0.is_open())
        cout << "Files Opened Succefully! " << endl;
    else
        cerr << " Problem Openening Files!" << endl;

    
    ///////////////////////////////////////////////////////////////////////////
    // read
    ///////////////////////////////////////////////////////////////////////////

    double e1, e2, w;

    int counts = 0;
    while(if_isotope_0.good())
    {
        if_isotope_0 >> e1 >> e2 >> w;
        h_psf_G0->Fill(e1, e2, w);
        ++ counts;
    }
    std::cout << fname_isotope_0 << " count: " << counts << std::endl;

    counts = 0;
    while(if_isotope_1.good())
    {
        if_isotope_1 >> e1 >> e2 >> w;
        h_psf_G2->Fill(e1, e2, w);
        ++ counts;
    }
    std::cout << fname_isotope_1 << " count: " << counts << std::endl;

    std::cout << "Integral" << std::endl;

    std::cout << h_psf_G0->Integral() << std::endl;
    std::cout << h_psf_G2->Integral() << std::endl;

    h_psf_G0->Scale(1.0 / h_psf_G0->Integral());
    h_psf_G2->Scale(1.0 / h_psf_G2->Integral());

    std::cout << h_psf_G0->Integral() << std::endl;
    std::cout << h_psf_G2->Integral() << std::endl;

}



bool flag_rebuild_ch0_HSD = true;
bool flag_rebuild_ch0_SSD = true;
bool flag_rebuild_ch1_HSD = true;
bool flag_rebuild_ch1_SSD = true;
bool flag_rebuild_ch10_HSD = true;
bool flag_rebuild_ch10_SSD = true;

void rebuild_histograms()
{

    ///////////////////////////////////////////////////////////////////////////
    // ch0

    if(flag_rebuild_ch0_HSD)
    {
        h_ch0_HSD = new TH1D("h_ch0_HSD", "h_ch0_HSD",
                             N_BINS_CH0, 0.0, H_E_MAX);
    }

    if(flag_rebuild_ch0_HSD)
    {
        h_ch0_SSD = new TH1D("h_ch0_SSD", "h_ch0_SSD",
                             N_BINS_CH0, 0.0, H_E_MAX);
    }

    ///////////////////////////////////////////////////////////////////////////
    // ch 1

    if(flag_rebuild_ch1_HSD)
    {
        h_ch1_HSD = new TH1D("h_ch1_HSD", "h_ch1_HSD",
                             N_BINS_CH1, 0.0, H_E_MAX);
    }

    if(flag_rebuild_ch1_HSD)
    {
        h_ch1_SSD = new TH1D("h_ch1_SSD", "h_ch1_SSD",
                             N_BINS_CH1, 0.0, H_E_MAX);
    }

    ///////////////////////////////////////////////////////////////////////////
    // ch10

    if(flag_rebuild_ch10_HSD)
    {
        h_ch10_HSD = new TH2D("h_ch10_HSD", "h_ch10_HSD",
                             N_BINS_CH10, 0.0, H_E_MAX,
                             N_BINS_CH10, 0.0, H_E_MAX);
    }

    if(flag_rebuild_ch10_HSD)
    {
        h_ch10_SSD = new TH2D("h_ch10_SSD", "h_ch10_SSD",
                             N_BINS_CH10, 0.0, H_E_MAX,
                             N_BINS_CH10, 0.0, H_E_MAX);
    }


    for(Int_t j = 1; j <= h_psf_G0->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_psf_G0->GetNbinsX(); ++ i)
        {
            const double content_G0 = h_psf_G0->GetBinContent(i, j);
            const double content_G2 = h_psf_G2->GetBinContent(i, j);

            const double e1 = h_psf_G0->GetXaxis()->GetBinCenter(i);
            const double e2 = h_psf_G0->GetYaxis()->GetBinCenter(j);
    
            const double xi_31_HSD = 0.0; // TODO check
            const double xi_31_SSD = 0.3; // TODO: this is a guess, get correct value from Fedor

            const double G0 = G0_ps_integral_MeV;
            const double G2 = G2_ps_integral_MeV;

            const double weight_HSD = 1.0 / (G0 + xi_31_HSD * G2) * (content_G0 + xi_31_HSD * content_G2);
            const double weight_SSD = 1.0 / (G0 + xi_31_SSD * G2) * (content_G0 + xi_31_SSD * content_G2);

            // channel 0
            if(flag_rebuild_ch0_HSD)
            {
                h_ch0_HSD->Fill(e1 + e2, weight_HSD);
            }
            if(flag_rebuild_ch0_SSD)
            {
                h_ch0_SSD->Fill(e1 + e2, weight_SSD);
            }

            // channel 1
            if(flag_rebuild_ch1_HSD)
            {
                h_ch1_HSD->Fill(e1, weight_HSD);
                h_ch1_HSD->Fill(e2, weight_HSD);
            }
            if(flag_rebuild_ch1_SSD)
            {
                h_ch1_SSD->Fill(e1, weight_SSD);
                h_ch1_SSD->Fill(e2, weight_SSD);
            }

            // channel 10
            if(flag_rebuild_ch10_HSD)
            {
                if(e1 < e2)
                {
                    h_ch10_HSD->Fill(e1, e2, weight_HSD);
                }
                else
                {
                    h_ch10_HSD->Fill(e2, e1, weight_HSD);
                }
            }
            if(flag_rebuild_ch10_SSD)
            {
                if(e1 < e2)
                {
                    h_ch10_SSD->Fill(e1, e2, weight_SSD);
                }
                else
                {
                    h_ch10_SSD->Fill(e2, e1, weight_SSD);
                }
            }
        }
    }

    h_ch0_HSD->Scale(1.0 / h_ch0_HSD->Integral());
    h_ch0_SSD->Scale(1.0 / h_ch0_SSD->Integral());
    h_ch1_HSD->Scale(2.0 / h_ch1_HSD->Integral());
    h_ch1_SSD->Scale(2.0 / h_ch1_SSD->Integral());
    h_ch10_HSD->Scale(1.0 / h_ch10_HSD->Integral());
    h_ch10_SSD->Scale(1.0 / h_ch10_SSD->Integral());


}



///////////////////////////////////////////////////////////////////////////////
// chi2 1D
///////////////////////////////////////////////////////////////////////////////

double calc_chi2(TH1D* h_data, TH1D* h_model)
{
    double chi2 = 0.0;
    for(Int_t i = 1; i <= h_model->GetNbinsY(); ++ i)
    {
        double content_model = h_model->GetBinContent(i);
        double content_data = h_data->GetBinContent(i);
        double error_model = std::sqrt(content_model);
        if((content_model == 0.0) && (content_data == 0.0))
        {
            continue;
        }
        double chi = std::pow((content_data - content_model) / error_model, 2.0);
        chi2 += chi;
    }
    return chi2;
}


///////////////////////////////////////////////////////////////////////////////
// chi2 2D
///////////////////////////////////////////////////////////////////////////////

double calc_chi2(TH2D* h_data, TH2D* h_model)
{
    double chi2 = 0.0;
    for(Int_t j = 1; j <= h_model->GetNbinsX(); ++ j)
    {
        for(Int_t i = 1; i <= h_model->GetNbinsY(); ++ i)
        {
            double content_model = h_model->GetBinContent(i, j);
            double content_data = h_data->GetBinContent(i, j);
            double error_model = std::sqrt(content_model);
            if((content_model == 0.0) && (content_data == 0.0))
            {
                continue;
            }
            double chi = std::pow((content_data - content_model) / error_model, 2.0);
            chi2 += chi;
        }
    }
    return chi2;
}


///////////////////////////////////////////////////////////////////////////////
// chi2 2D with systematics
///////////////////////////////////////////////////////////////////////////////

double calc_chi2(TH2D* h_data, TH2D* h_model, std::vector<double> sys)
{
    double chi2 = 0.0;
    for(Int_t j = 1; j <= h_model->GetNbinsX(); ++ j)
    {
        for(Int_t i = 1; i <= h_model->GetNbinsY(); ++ i)
        {
            double content_model = h_model->GetBinContent(i, j);
            double content_data = h_data->GetBinContent(i, j);
            double error_model = std::sqrt(content_model);
            if((content_model == 0.0) && (content_data == 0.0))
            {
                continue;
            }
            double chi = std::pow((content_data - content_model) / error_model, 2.0);
            chi2 += chi;
        }
    }
    return chi2;
}



int main()
{

    gStyle->SetLabelFont(43);
    gStyle->SetLabelSize(15);
    gStyle->SetTitleFont(43);
    gStyle->SetTitleSize(20);

    load_psf_data();
    
    rebuild_histograms();

    TCanvas *c = new TCanvas("c", "c");
    h_ch10_HSD->Draw("colz");
    //h_ch10_SSD->Draw("colz");


    // draw total E for 3, 4 bins
    N_BINS_CH0 = 3;
    rebuild_histograms();
    const int N_events = N2e;
    h_ch0_HSD->Scale(N_events / h_ch0_HSD->Integral());
    h_ch0_SSD->Scale(N_events / h_ch0_SSD->Integral());
    h_ch1_HSD->Scale(2.0 * N_events / h_ch1_HSD->Integral());
    h_ch1_SSD->Scale(2.0 * N_events / h_ch1_SSD->Integral());
    h_ch10_HSD->Scale(N_events / h_ch10_HSD->Integral());
    h_ch10_SSD->Scale(N_events / h_ch10_SSD->Integral());
    new TCanvas;
    h_ch0_HSD->Draw("hist");
    h_ch0_SSD->Draw("histsame");

    N_BINS_CH0 = 4;
    rebuild_histograms();
    h_ch0_HSD->Scale(N_events / h_ch0_HSD->Integral());
    h_ch0_SSD->Scale(N_events / h_ch0_SSD->Integral());
    h_ch1_HSD->Scale(2.0 * N_events / h_ch1_HSD->Integral());
    h_ch1_SSD->Scale(2.0 * N_events / h_ch1_SSD->Integral());
    h_ch10_HSD->Scale(N_events / h_ch10_HSD->Integral());
    h_ch10_SSD->Scale(N_events / h_ch10_SSD->Integral());
    new TCanvas;
    h_ch0_HSD->Draw("hist");
    h_ch0_SSD->Draw("histsame");




    if(0)
    {
        int N_BINS_MIN = 1;
        int N_BINS_MAX = 100;
        TGraph *g_chi2_ch0 = new TGraph(N_BINS_MAX + 1 - N_BINS_MIN);
        TGraph *g_chi2_ch1 = new TGraph(N_BINS_MAX + 1 - N_BINS_MIN);
        TGraph *g_chi2_ch10 = new TGraph(N_BINS_MAX + 1 - N_BINS_MIN);
        int i = 0;
        for(int n_bins = N_BINS_MIN; n_bins <= N_BINS_MAX; ++ n_bins)
        {
            N_BINS_CH0 = n_bins;    
            N_BINS_CH1 = n_bins;    
            N_BINS_CH10 = n_bins;    
            
            rebuild_histograms();
            const int N_events = N2e;
            h_ch0_HSD->Scale(N_events / h_ch0_HSD->Integral());
            h_ch0_SSD->Scale(N_events / h_ch0_SSD->Integral());
            h_ch1_HSD->Scale(2.0 * N_events / h_ch1_HSD->Integral());
            h_ch1_SSD->Scale(2.0 * N_events / h_ch1_SSD->Integral());
            h_ch10_HSD->Scale(N_events / h_ch10_HSD->Integral());
            h_ch10_SSD->Scale(N_events / h_ch10_SSD->Integral());

            double chi2_ch0 = calc_chi2(h_ch0_SSD, h_ch0_HSD);
            double chi2_ch1 = calc_chi2(h_ch1_SSD, h_ch1_HSD);
            double chi2_ch10 = calc_chi2(h_ch10_SSD, h_ch10_HSD);
            
            g_chi2_ch0->SetPoint(i, n_bins, chi2_ch0);
            g_chi2_ch1->SetPoint(i, n_bins, chi2_ch1);
            g_chi2_ch10->SetPoint(i, n_bins, chi2_ch10);

            ++ i;
        }
        g_chi2_ch0->SetLineColor(kRed);
        g_chi2_ch1->SetLineColor(kGreen);
        g_chi2_ch10->SetLineColor(kBlue);
        g_chi2_ch0->SetLineWidth(2.0);
        g_chi2_ch1->SetLineWidth(2.0);
        g_chi2_ch10->SetLineWidth(2.0);
        
        TCanvas *c_chi2 = new TCanvas("c_chi2", "c_chi2");
        c_chi2->SetLogy();
        //g_chi2_ch10->SetMinimum(1.0e-3);
        //g_chi2_ch10->SetMaximum(1.0e+4);
        g_chi2_ch10->SetTitle(0);
        g_chi2_ch10->Draw("al");
        g_chi2_ch0->Draw("l");
        g_chi2_ch1->Draw("l");
    }


    
    Int_t N_BINS_CH0 = 40;
    Int_t N_BINS_CH1 = 40;
    Int_t N_BINS_CH10 = 40;
    rebuild_histograms();


    // draw chi2 as a function of number of events
    if(0)
    {   
        double N_events_min = 100;
        double N_events_max = 1e+6;
        int N_steps = 100;
        double exp_min = std::log(N_events_min) / std::log(10.0);
        double exp_max = std::log(N_events_max) / std::log(10.0);
        double step = (exp_max - exp_min) / (double)N_steps;

        int i_max = 0;
        for(double exp = exp_min; exp <= exp_max; exp += step) ++ i_max;

        TGraph *g_chi2_ch0 = new TGraph(i_max);
        TGraph *g_chi2_ch1 = new TGraph(i_max);
        TGraph *g_chi2_ch10 = new TGraph(i_max);

        int i = 0;
        for(double exp = exp_min; exp <= exp_max; exp += step)
        {
            double N_events = std::pow(10.0, exp);
            
            h_ch0_HSD->Scale(N_events / h_ch0_HSD->Integral());
            h_ch0_SSD->Scale(N_events / h_ch0_SSD->Integral());
            h_ch1_HSD->Scale(2.0 * N_events / h_ch1_HSD->Integral());
            h_ch1_SSD->Scale(2.0 * N_events / h_ch1_SSD->Integral());
            h_ch10_HSD->Scale(N_events / h_ch10_HSD->Integral());
            h_ch10_SSD->Scale(N_events / h_ch10_SSD->Integral());

            double chi2_ch0 = calc_chi2(h_ch0_SSD, h_ch0_HSD);
            double chi2_ch1 = calc_chi2(h_ch1_SSD, h_ch1_HSD);
            double chi2_ch10 = calc_chi2(h_ch10_SSD, h_ch10_HSD);

            g_chi2_ch0->SetPoint(i, N_events, chi2_ch0);
            g_chi2_ch1->SetPoint(i, N_events, chi2_ch1);
            g_chi2_ch10->SetPoint(i, N_events, chi2_ch10);

            ++ i;
        }

        g_chi2_ch0->SetLineColor(kRed);
        g_chi2_ch1->SetLineColor(kGreen);
        g_chi2_ch10->SetLineColor(kBlue);
        g_chi2_ch0->SetLineWidth(2.0);
        g_chi2_ch1->SetLineWidth(2.0);
        g_chi2_ch10->SetLineWidth(2.0);
        
        TCanvas *c_chi2 = new TCanvas("c_chi2", "c_chi2");
        c_chi2->SetLogx();
        c_chi2->SetLogy();
        g_chi2_ch10->SetMinimum(1.0e-3);
        g_chi2_ch10->SetMaximum(1.0e+4);
        g_chi2_ch10->SetTitle(0);
        g_chi2_ch10->Draw("al");
        g_chi2_ch0->Draw("l");
        g_chi2_ch1->Draw("l");
    }


    return 0;

}
