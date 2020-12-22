



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

TH1D *h_psf_G0_px = nullptr;
TH1D *h_psf_G0_py = nullptr;
TH1D *h_psf_G0_Cx = nullptr;
TH1D *h_psf_G0_Cy = nullptr;
TH1D *h_psf_G0_Cinvx = nullptr;
TH1D *h_psf_G0_Cinvy = nullptr;

TH1D *h_psf_G2_px = nullptr;
TH1D *h_psf_G2_py = nullptr;
TH1D *h_psf_G2_Cx = nullptr;
TH1D *h_psf_G2_Cy = nullptr;
TH1D *h_psf_G2_Cinvx = nullptr;
TH1D *h_psf_G2_Cinvy = nullptr;




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
Int_t N_BINS_CH10 = 1000;




void marginalize(
    const TH2D* const h_input,
    TH1D* h_output_x,
    TH1D* h_output_y)
{
    // x
    for(Int_t i = 1; i <= h_input->GetNbinsX(); ++ i)
    {
        double sum = 0.0;
        for(Int_t j = 1; j <= h_input->GetNbinsY(); ++ j)
        {
            sum += h_input->GetBinContent(i, j);
        }
        h_output_x->SetBinContent(i, sum);
    }
    h_output_x->Scale(1.0 / h_output_x->Integral());

    // y
    for(Int_t j = 1; j <= h_input->GetNbinsY(); ++ j)
    {
        double sum = 0.0;
        for(Int_t i = 1; i <= h_input->GetNbinsX(); ++ i)
        {
            sum += h_input->GetBinContent(i, j);
        }
        h_output_y->SetBinContent(j, sum);
    }
    h_output_y->Scale(1.0 / h_output_y->Integral());
}



void convert_probability_to_CDF_inverse(
    const TH1D* const h_prob,
    TH1D* h_CDF,
    TH1D* h_CDF_inverse)
{

    // CDF
    double sum = 0.0;
    for(Int_t i = 1; i <= h_prob->GetNbinsX(); ++ i)
    {
        double content = h_prob->GetBinContent(i);
        double sum1 = sum + 0.5 * content;
        h_CDF->SetBinContent(i, sum1);
        sum += content;
    }
    std::cout << "CDF: " << h_CDF->GetBinContent(1) << " " << h_CDF->GetBinContent(h_CDF->GetNbinsX()) << std::endl;

    // graph CDF inv
    TGraph *g_CDF_inverse = new TGraph(h_CDF->GetNbinsX());
    g_CDF_inverse->SetPoint(0, 0.0, 0.0);

    for(Int_t i = 1; i <= h_CDF->GetNbinsX(); ++ i)
    {
        const double content = h_CDF->GetBinContent(i);
        const double bin_center = h_CDF->GetXaxis()->GetBinCenter(i);
        g_CDF_inverse->SetPoint(i, content, bin_center);
    }

    // histogram CDF inv
    for(Int_t i = 1; i <= h_CDF_inverse->GetNbinsX(); ++ i)
    {
        const double x = h_CDF_inverse->GetXaxis()->GetBinCenter(i);
        const double content = g_CDF_inverse->Eval(x);
        h_CDF_inverse->SetBinContent(i, content);
    }

}










void load_psf_data()
{

    //const Int_t N_BINS_XY = (Q_bb - 1.0e-3 * E_min) / 1.0e-3;
    const Int_t N_BINS_XY = H_E_MAX / 1.0e-3;
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

    
    ///////////////////////////////////////////////////////////////////////////
    // build probability, CDF and CDF inverse histograms

    h_psf_G0_px = new TH1D("h_psf_G0_px", "h_psf_G0_px", N_BINS_XY, 0.0, H_E_MAX);
    h_psf_G0_py = new TH1D("h_psf_G0_py", "h_psf_G0_py", N_BINS_XY, 0.0, H_E_MAX);

    h_psf_G2_px = new TH1D("h_psf_G2_px", "h_psf_G2_px", N_BINS_XY, 0.0, H_E_MAX);
    h_psf_G2_py = new TH1D("h_psf_G2_py", "h_psf_G2_py", N_BINS_XY, 0.0, H_E_MAX);

    marginalize(h_psf_G0, h_psf_G0_px, h_psf_G0_py);
    marginalize(h_psf_G2, h_psf_G2_px, h_psf_G2_py);

    /*
    new TCanvas;
    h_psf_G0_px->Draw("hist");
    new TCanvas;
    h_psf_G0_py->Draw("hist");
    */

    h_psf_G0_Cx = new TH1D("h_psf_G0_Cx", "h_psf_G0_Cx", N_BINS_XY, 0.0, H_E_MAX);
    h_psf_G0_Cy = new TH1D("h_psf_G0_Cy", "h_psf_G0_Cy", N_BINS_XY, 0.0, H_E_MAX);

    h_psf_G0_Cinvx = new TH1D("h_psf_G0_Cinvx", "h_psf_G0_Cinvx", N_BINS_XY, 0.0, 1.0);
    h_psf_G0_Cinvy = new TH1D("h_psf_G0_Cinvy", "h_psf_G0_Cinvy", N_BINS_XY, 0.0, 1.0);

    convert_probability_to_CDF_inverse(h_psf_G0_px, h_psf_G0_Cx, h_psf_G0_Cinvx);
    convert_probability_to_CDF_inverse(h_psf_G0_py, h_psf_G0_Cy, h_psf_G0_Cinvy);

    /*
    new TCanvas;
    h_psf_G0_Cinvx->Draw("hist");
    new TCanvas;
    h_psf_G0_Cinvy->Draw("hist");
    */

    h_psf_G2_Cx = new TH1D("h_psf_G2_Cx", "h_psf_G2_Cx", N_BINS_XY, 0.0, H_E_MAX);
    h_psf_G2_Cy = new TH1D("h_psf_G2_Cy", "h_psf_G2_Cy", N_BINS_XY, 0.0, H_E_MAX);

    h_psf_G2_Cinvx = new TH1D("h_psf_G2_Cinvx", "h_psf_G2_Cinvx", N_BINS_XY, 0.0, 1.0);
    h_psf_G2_Cinvy = new TH1D("h_psf_G2_Cinvy", "h_psf_G2_Cinvy", N_BINS_XY, 0.0, 1.0);

    convert_probability_to_CDF_inverse(h_psf_G2_px, h_psf_G2_Cx, h_psf_G2_Cinvx);
    convert_probability_to_CDF_inverse(h_psf_G2_py, h_psf_G2_Cy, h_psf_G2_Cinvy);


    /*
    TCanvas *c = new TCanvas("c", "c");
    h_psf_G0_px->Draw("hist");
    TCanvas *c2 = new TCanvas("c2", "c2");
    h_psf_G0_Cx->Draw("hist");
    TCanvas *c3 = new TCanvas("c3", "c3");
    h_psf_G0_Cinvx->Draw("hist");
    */

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



    for(unsigned long long event = 0; event < N2e; ++ event)
    {
        double ux = gRNG.Uniform();
        double uy = gRNG.Uniform();

        double e1_G0 = h_psf_G0_Cinvx->Interpolate(ux);
        double e2_G0 = h_psf_G0_Cinvy->Interpolate(uy);

//        double e1_G2 = h_psf_G2_Cinvx->Interpolate(ux);
//        double e2_G2 = h_psf_G2_Cinvy->Interpolate(uy);

        /*
        std::cout << e1_G0 << " " << e2_G0 << std::endl;
        std::cout << e1_G2 << " " << e2_G2 << std::endl;
        std::cin.get();
        */

        const double xi_31_HSD = 0.0; // TODO assume this is correct
        const double xi_31_SSD = 0.3; // TODO no idea what this is

        const double G0 = G0_ps_integral_MeV;
        const double G2 = G2_ps_integral_MeV;

        double weight_HSD = 1.0;
        double weight_SSD = weight_HSD * ((G0 + xi_31_SSD * G2) / (G0));

        /*
        double weight_G0_HSD = 1.0 / (G0 + xi_31_HSD * G2);
        double weight_G2_HSD = xi_31_HSD * weight_G0_HSD;

        double weight_G0_SSD = 1.0 / (G0 + xi_31_SSD * G2);
        double weight_G2_SSD = xi_31_SSD * weight_G0_SSD;
        */

        ///////////////////////////////////////////////////////////////////////
        // fill histograms - total electron energy (ch0)
        
        if(flag_rebuild_ch0_HSD)
        {
            //h_ch0_HSD->Fill(e1_G0 + e2_G0, weight_G0_HSD);
            //h_ch0_HSD->Fill(e1_G2 + e2_G2, weight_G2_HSD);
            h_ch0_HSD->Fill(e1_G0 + e2_G0, weight_HSD);
        }

        if(flag_rebuild_ch0_SSD)
        {
            //h_ch0_SSD->Fill(e1_G0 + e2_G0, weight_G0_SSD);
            //h_ch0_SSD->Fill(e1_G2 + e2_G2, weight_G2_SSD);
            h_ch0_SSD->Fill(e1_G0 + e2_G0, weight_SSD);
        }


        ///////////////////////////////////////////////////////////////////////
        // fill histograms - single electron energy (ch1)
        
        if(flag_rebuild_ch1_HSD)
        {
            //h_ch1_HSD->Fill(e1_G0, weight_G0_HSD);
            //h_ch1_HSD->Fill(e2_G0, weight_G0_HSD);
            //h_ch1_HSD->Fill(e1_G2, weight_G2_HSD);
            //h_ch1_HSD->Fill(e2_G2, weight_G2_HSD);
            h_ch1_HSD->Fill(e1_G0, weight_HSD);
            h_ch1_HSD->Fill(e2_G0, weight_HSD);
        }

        if(flag_rebuild_ch1_SSD)
        {
            //h_ch1_SSD->Fill(e1_G0, weight_G0_SSD);
            //h_ch1_SSD->Fill(e2_G0, weight_G0_SSD);
            //h_ch1_SSD->Fill(e1_G2, weight_G2_SSD);
            //h_ch1_SSD->Fill(e2_G2, weight_G2_SSD);
            h_ch1_SSD->Fill(e1_G0, weight_SSD);
            h_ch1_SSD->Fill(e2_G0, weight_SSD);
        }


        ///////////////////////////////////////////////////////////////////////
        // fill histograms - high / low electron energy 2d (ch10)

        if(flag_rebuild_ch10_HSD)
        {
            //h_ch10_HSD->Fill(e1_G0, e2_G0, weight_G0_HSD);
            //h_ch10_HSD->Fill(e1_G2, e2_G2, weight_G2_HSD);
            if(e1_G0 < e2_G0)
            {
                h_ch10_HSD->Fill(e1_G0, e2_G0, weight_HSD);
            }
            else
            {
                h_ch10_HSD->Fill(e2_G0, e1_G0, weight_HSD);
            }
        }

        if(flag_rebuild_ch10_SSD)
        {
            //h_ch10_SSD->Fill(e1_G0, e2_G0, weight_G0_SSD);
            //h_ch10_SSD->Fill(e1_G2, e2_G2, weight_G2_SSD);
            if(e1_G0 < e2_G0)
            {
                h_ch10_SSD->Fill(e1_G0, e2_G0, weight_SSD);
            }
            else
            {
                h_ch10_SSD->Fill(e2_G0, e1_G0, weight_SSD);
            }
        }

    }
    


}



int main()
{


    load_psf_data();
    
    rebuild_histograms();

    TCanvas *c = new TCanvas("c", "c");
    h_ch10_HSD->Draw("colz");
    //h_ch10_SSD->Draw("colz");

    return 0;

}
