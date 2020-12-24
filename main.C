



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

TH1D *h_ch0_HSD_sys = nullptr;
TH1D *h_ch1_HSD_sys = nullptr;
TH2D *h_ch10_HSD_sys = nullptr;

const double H_E_MAX = 3.2;
Int_t N_BINS_CH0 = 32;
Int_t N_BINS_CH1 = 32;
Int_t N_BINS_CH10 = 32;

//std::vector<double> V_PHYS_SYS_CH0;
//std::vector<double> V_PHYS_SYS_CH1;
//std::vector<double> V_PHYS_SYS_CH10;

std::vector<double> alpha_coeff_ch0_sys1; // 5.55 % efficiency
std::vector<double> alpha_coeff_ch0_sys2; // 0.2 % energy shift
std::vector<double> alpha_coeff_ch1_sys1; // 5.55 % efficiency
std::vector<double> alpha_coeff_ch1_sys2; // 0.2 % energy shift
std::vector<double> alpha_coeff_ch10_sys1; // 5.55 % efficiency
std::vector<double> alpha_coeff_ch10_sys2; // 0.2 % energy shift

std::vector<double> chi2_SYS_CH0;
std::vector<double> chi2_SYS_CH1;
std::vector<double> chi2_SYS_CH10;




#include "draw.h"



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

//    h_psf_G0->Scale(1.0 / h_psf_G0->Integral());
//    h_psf_G2->Scale(1.0 / h_psf_G2->Integral());

//    std::cout << h_psf_G0->Integral() << std::endl;
//    std::cout << h_psf_G2->Integral() << std::endl;

}



bool flag_rebuild_ch0_HSD = true;
bool flag_rebuild_ch0_SSD = true;
bool flag_rebuild_ch1_HSD = true;
bool flag_rebuild_ch1_SSD = true;
bool flag_rebuild_ch10_HSD = true;
bool flag_rebuild_ch10_SSD = true;

bool flag_rebuild_ch0_HSD_sys = false;
bool flag_rebuild_ch1_HSD_sys = false;
bool flag_rebuild_ch10_HSD_sys = false;


double integral_ch0_HSD = 1.0;
double integral_ch1_HSD = 1.0;
double integral_ch10_HSD = 1.0;


double gSYS_efficiency = 0.0;
double gSYS_coherent_energy_shift = 0.0;


void rebuild_histograms(const double N_events)
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

    if(flag_rebuild_ch0_HSD_sys)
    {
        h_ch0_HSD_sys = new TH1D("h_ch0_HSD_sys", "h_ch0_HSD_sys",
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

    if(flag_rebuild_ch1_HSD_sys)
    {
        h_ch1_HSD_sys = new TH1D("h_ch1_HSD_sys", "h_ch1_HSD_sys",
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

    if(flag_rebuild_ch10_HSD_sys)
    {
        h_ch10_HSD_sys = new TH2D("h_ch10_HSD_sys", "h_ch10_HSD_sys",
                             N_BINS_CH10, 0.0, H_E_MAX,
                             N_BINS_CH10, 0.0, H_E_MAX);
    }

    
    ///////////////////////////////////////////////////////////////////////////
    // read data
    ///////////////////////////////////////////////////////////////////////////

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

            const double weight_HSD = (1.0 / (G0 + xi_31_HSD * G2)) * (content_G0 + xi_31_HSD * content_G2);
            const double weight_SSD = (1.0 / (G0 + xi_31_SSD * G2)) * (content_G0 + xi_31_SSD * content_G2);

            // channel 0
            if(flag_rebuild_ch0_HSD)
            {
                h_ch0_HSD->Fill(e1 + e2, weight_HSD);
            }
            if(flag_rebuild_ch0_SSD)
            {
                h_ch0_SSD->Fill(e1 + e2, weight_SSD);
            }
            if(flag_rebuild_ch0_HSD_sys)
            {
                double e1_t = e1 * (1.0 + gSYS_coherent_energy_shift);
                double e2_t = e2 * (1.0 + gSYS_coherent_energy_shift);
                double weight_HSD_t = weight_HSD * (1.0 + gSYS_efficiency);
                h_ch0_HSD_sys->Fill(e1_t + e2_t, weight_HSD_t);
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
            if(flag_rebuild_ch1_HSD_sys)
            {
                double e1_t = e1 * (1.0 + gSYS_coherent_energy_shift);
                double e2_t = e2 * (1.0 + gSYS_coherent_energy_shift);
                double weight_HSD_t = weight_HSD * (1.0 + gSYS_efficiency);
                h_ch1_HSD_sys->Fill(e1_t, weight_HSD_t);
                h_ch1_HSD_sys->Fill(e2_t, weight_HSD_t);
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
            if(flag_rebuild_ch10_HSD_sys)
            {
                if(e1 < e2)
                {
                    double e1_t = e1 * (1.0 + gSYS_coherent_energy_shift);
                    double e2_t = e2 * (1.0 + gSYS_coherent_energy_shift);
                    double weight_HSD_t = weight_HSD * (1.0 + gSYS_efficiency);
                    h_ch10_HSD_sys->Fill(e1_t, e2_t, weight_HSD_t);
                }
                else
                {
                    double e1_t = e1 * (1.0 + gSYS_coherent_energy_shift);
                    double e2_t = e2 * (1.0 + gSYS_coherent_energy_shift);
                    double weight_HSD_t = weight_HSD * (1.0 + gSYS_efficiency);
                    h_ch10_HSD_sys->Fill(e2_t, e1_t, weight_HSD_t);
                }
            }
        }
    }


    if(flag_rebuild_ch0_HSD)
    {
        integral_ch0_HSD = h_ch0_HSD->Integral();
    }
    if(flag_rebuild_ch1_HSD)
    {
        integral_ch1_HSD = h_ch1_HSD->Integral();
    }
    if(flag_rebuild_ch10_HSD)
    {
        integral_ch10_HSD = h_ch10_HSD->Integral();
    }

    
    if(flag_rebuild_ch0_HSD)
    {
        h_ch0_HSD->Scale(N_events / integral_ch0_HSD);
    }
    if(flag_rebuild_ch0_SSD)
    {
        h_ch0_SSD->Scale(N_events / integral_ch0_HSD);
    }
    if(flag_rebuild_ch0_HSD_sys)
    {
        h_ch0_HSD_sys->Scale(N_events / integral_ch0_HSD);
    }

    if(flag_rebuild_ch1_HSD)
    {
        h_ch1_HSD->Scale(2.0 * N_events / integral_ch1_HSD);
    }
    if(flag_rebuild_ch1_SSD)
    {
        h_ch1_SSD->Scale(2.0 * N_events / integral_ch1_HSD);
    }
    if(flag_rebuild_ch1_HSD_sys)
    {
        h_ch1_HSD_sys->Scale(2.0 * N_events / integral_ch1_HSD);
    }

    if(flag_rebuild_ch10_HSD)
    {
        h_ch10_HSD->Scale(N_events / integral_ch10_HSD);
    }
    if(flag_rebuild_ch10_SSD)
    {
        h_ch10_SSD->Scale(N_events / integral_ch10_HSD);
    }
    if(flag_rebuild_ch10_HSD_sys)
    {
        h_ch10_HSD_sys->Scale(N_events / integral_ch10_HSD);
    }


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
// chi2 1D with systematics
///////////////////////////////////////////////////////////////////////////////

double calc_chi2(
    TH1D* h_data,
    TH1D* h_model,
    const std::vector<double> &sys1,
    const std::vector<double> &sys2)
{
    
    std::vector<double> delta; // d - m
    std::vector<double> sigma; // sigma i

    for(Int_t i = 1; i <= h_model->GetNbinsX(); ++ i)
    {
        double content_model = h_model->GetBinContent(i);
        double content_data = h_data->GetBinContent(i);
        double error_model = std::sqrt(content_model);
        if(content_model == 0.0) continue;

        delta.push_back(content_data - content_model);
        sigma.push_back(error_model);
    }
    TMatrixD *matrix = new TMatrixD(delta.size(), delta.size());

    {
        int j = 0;
        for(int jj = 1; jj <= h_model->GetNbinsX(); ++ jj)
        //for(int j = 0; j < matrix->GetNrows(); ++ j)
        {
            if(h_model->GetBinContent(jj) <= 0.0) continue; 

            int i = 0;
            for(int ii = 1; ii <= h_model->GetNbinsX(); ++ ii)
            //for(int i = 0; i < matrix->GetNcols(); ++ i)
            {
                if(h_model->GetBinContent(ii) <= 0.0) continue;

                double c = 0.0;
                if(i == j)
                {
                    double sigma_i = sigma.at(i);
                    double sigma_j = sigma.at(j);
                    c += sigma_i * sigma_j;
                }
                double alpha_i = sys1.at(ii - 1);
                double alpha_j = sys1.at(jj - 1);
                c += alpha_i * alpha_j;
                alpha_i = sys2.at(ii - 1);
                alpha_j = sys2.at(jj - 1);
                c += alpha_i * alpha_j;

                //std::cout << "j=" << j << " i=" << i << " " << c << std::endl;
                matrix->operator[](j).operator[](i) = c;

                ++ i;
            }

            ++ j;
        }

    }
    matrix->Invert();

    double chi2 = 0.0;
    for(int j = 0; j < matrix->GetNrows(); ++ j)
    {
        for(int i = 0; i < matrix->GetNcols(); ++ i)
        {
            double v_i = delta.at(i);
            double V_ij = matrix->operator[](j).operator[](i);
            double v_j = delta.at(j);
            chi2 += v_i * V_ij * v_j;
        }
    }

    return chi2;
}


///////////////////////////////////////////////////////////////////////////////
// chi2 2D with systematics
///////////////////////////////////////////////////////////////////////////////

double calc_chi2(TH2D* h_data,
    TH2D* h_model,
    const std::vector<double> &sys1,
    const std::vector<double> &sys2)
{
    
    std::vector<double> delta; // d - m
    std::vector<double> sigma; // sigma i

    int count = 0;
    int count_sqrt = 0;
    for(Int_t j = 1; j <= h_model->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_model->GetNbinsX(); ++ i)
        {
            double content_model = h_model->GetBinContent(i);
            double content_data = h_data->GetBinContent(i);
            double error_model = std::sqrt(content_model);
            
            sigma.push_back(error_model);

            if(content_model == 0.0) continue;
            ++ count;

            delta.push_back(content_data - content_model);
        }
    }
    TMatrixD *matrix = new TMatrixD(count, count);

    for(Int_t i = 1; i <= h_model->GetNbinsX(); ++ i)
    {
        double content_model = h_model->GetBinContent(i);
        if(content_model > 0.0) ++ count_sqrt;
    }
    

    {
        int l = 0;
        for(int ll = 1; ll <= h_model->GetNbinsY(); ++ ll)
        {
            int k = 0;
            for(int kk = 1; kk <= h_model->GetNbinsX(); ++ kk)
            {
                if(h_model->GetBinContent(kk, ll) <= 0.0) continue; 

                int j = 0;
                for(int jj = 1; jj <= h_model->GetNbinsY(); ++ jj)
                //for(int j = 0; j < matrix->GetNrows(); ++ j)
                {
                    int i = 0;
                    for(int ii = 1; ii <= h_model->GetNbinsX(); ++ ii)
                    //for(int i = 0; i < matrix->GetNcols(); ++ i)
                    {
                        if(h_model->GetBinContent(ii, jj) <= 0.0) continue;

                        double c = 0.0;
                        if((i == k) && (j == l))
                        {
                            double sigma_i = sigma.at((ii - 1) + h_model->GetNbinsX() * (jj - 1));
                            double sigma_j = sigma.at((kk - 1) + h_model->GetNbinsX() * (ll - 1));
                            c += sigma_i * sigma_j;
                        }
                        int index_i = (ii - 1) + h_model->GetNbinsX() * (jj - 1);
                        int index_j = (kk - 1) + h_model->GetNbinsX() * (kk - 1);
                        double alpha_i = sys1.at(index_i);
                        double alpha_j = sys1.at(index_j);
                        c += alpha_i * alpha_j;
                        alpha_i = sys2.at(index_i);
                        alpha_j = sys2.at(index_j);
                        c += alpha_i * alpha_j;

                        //std::cout << "j=" << j << " i=" << i << " " << c << std::endl;
                        matrix->operator[](k + count_sqrt * l).operator[](i + count_sqrt * j) = c;

                        ++ i;
                    }

                    ++ j;
                }

                ++ k;
            }

            ++ l;
        }

    }
    matrix->Invert();

    double chi2 = 0.0;
    for(int j = 0; j < matrix->GetNrows(); ++ j)
    {
        for(int i = 0; i < matrix->GetNcols(); ++ i)
        {
            double v_i = delta.at(i);
            double V_ij = matrix->operator[](j).operator[](i);
            double v_j = delta.at(j);
            chi2 += v_i * V_ij * v_j;
        }
    }

    return chi2;
}


///////////////////////////////////////////////////////////////////////////////
// chi2 2D
///////////////////////////////////////////////////////////////////////////////

double calc_chi2(TH2D* h_data, TH2D* h_model)
{
    throw "todo";
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


void reset_systematic_arrays()
{

    //V_PHYS_SYS_CH0.clear();
    //V_PHYS_SYS_CH1.clear();
    //V_PHYS_SYS_CH10.clear();

    alpha_coeff_ch0_sys1.clear();
    alpha_coeff_ch0_sys2.clear();
    alpha_coeff_ch1_sys1.clear();
    alpha_coeff_ch1_sys2.clear();
    alpha_coeff_ch10_sys1.clear();
    alpha_coeff_ch10_sys2.clear();

    chi2_SYS_CH0.clear();
    chi2_SYS_CH1.clear();
    chi2_SYS_CH10.clear();

    /*
    for(int c = 0; c < N_BINS_CH0 * N_BINS_CH0; ++ c)
    {
        V_PHYS_SYS_CH0.push_back(0.0);
    }
    for(int c = 0; c < N_BINS_CH1 * N_BINS_CH1; ++ c)
    {
        V_PHYS_SYS_CH1.push_back(0.0);
    }
    for(int c = 0; c < N_BINS_CH10 * N_BINS_CH10 * N_BINS_CH10 * N_BINS_CH10; ++ c)
    {
        V_PHYS_SYS_CH10.push_back(0.0);
    }
    */
    
    for(int c = 0; c < N_BINS_CH0; ++ c)
    {
        alpha_coeff_ch0_sys1.push_back(0.0);
        alpha_coeff_ch0_sys2.push_back(0.0);
    }
    for(int c = 0; c < N_BINS_CH1; ++ c)
    {
        alpha_coeff_ch1_sys1.push_back(0.0);
        alpha_coeff_ch1_sys2.push_back(0.0);
    }
    for(int c = 0; c < N_BINS_CH10 * N_BINS_CH10; ++ c)
    {
        alpha_coeff_ch10_sys1.push_back(0.0);
        alpha_coeff_ch10_sys2.push_back(0.0);
    }


    for(int c = 0; c < N_BINS_CH0; ++ c)
    {
        chi2_SYS_CH0.push_back(0.0);
    }
    for(int c = 0; c < N_BINS_CH1; ++ c)
    {
        chi2_SYS_CH1.push_back(0.0);
    }
    for(int c = 0; c < N_BINS_CH10 * N_BINS_CH10; ++ c)
    {
        chi2_SYS_CH10.push_back(0.0);
    }
}

void systematic_arrays_add(std::vector<double> &alpha_coeff_ch0,
                           std::vector<double> &alpha_coeff_ch1,
                           std::vector<double> &alpha_coeff_ch10)
{

    ///////////////////////////////////////////////////////////////////////////
    // ch0
    ///////////////////////////////////////////////////////////////////////////

    /*
    for(Int_t j = 1; j <= N_BINS_CH0; ++ j)
    {
        for(Int_t i = 1; i <= N_BINS_CH0; ++ i)
        {
            double alpha_i = h_ch0_HSD_sys->GetBinContent(i) - h_ch0_HSD->GetBinContent(i);
            double alpha_j = h_ch0_HSD_sys->GetBinContent(j) - h_ch0_HSD->GetBinContent(j);
            V_PHYS_SYS_CH0.at((i - 1) + (j - 1) * N_BINS_CH0) += alpha_i * alpha_j;
        }
    }
    */

    for(Int_t i = 1; i <= N_BINS_CH0; ++ i)
    {
        double alpha_i = h_ch0_HSD_sys->GetBinContent(i) - h_ch0_HSD->GetBinContent(i);
        alpha_coeff_ch0.at(i - 1) = alpha_i;
    }

    for(Int_t i = 1; i <= N_BINS_CH0; ++ i)
    {
        double delta = h_ch0_HSD_sys->GetBinContent(i) - h_ch0_HSD->GetBinContent(i);
        //std::cout << "delta=" << delta << std::endl;
        chi2_SYS_CH0.at(i - 1) += std::abs(delta);
    }


    ///////////////////////////////////////////////////////////////////////////
    // ch1
    ///////////////////////////////////////////////////////////////////////////

    /*
    for(Int_t j = 1; j <= N_BINS_CH1; ++ j)
    {
        for(Int_t i = 1; i <= N_BINS_CH1; ++ i)
        {
            double alpha_i = h_ch1_HSD_sys->GetBinContent(i) - h_ch1_HSD->GetBinContent(i);
            double alpha_j = h_ch1_HSD_sys->GetBinContent(j) - h_ch1_HSD->GetBinContent(j);
            V_PHYS_SYS_CH1.at((i - 1) + (j - 1) * N_BINS_CH1) += alpha_i * alpha_j;
        }
    }
    */

    for(Int_t i = 1; i <= N_BINS_CH1; ++ i)
    {
        double alpha_i = h_ch1_HSD_sys->GetBinContent(i) - h_ch1_HSD->GetBinContent(i);
        alpha_coeff_ch1.at(i - 1) = alpha_i;
    }

    for(Int_t i = 1; i <= N_BINS_CH1; ++ i)
    {
        double delta = h_ch1_HSD_sys->GetBinContent(i) - h_ch1_HSD->GetBinContent(i);
        chi2_SYS_CH1.at(i - 1) += std::abs(delta);
    }
    

    ///////////////////////////////////////////////////////////////////////////
    // ch10
    ///////////////////////////////////////////////////////////////////////////

    /*
    for(Int_t l = 1; l <= N_BINS_CH10; ++ l)
    {
        for(Int_t k = 1; k <= N_BINS_CH10; ++ k)
        {
            for(Int_t j = 1; j <= N_BINS_CH10; ++ j)
            {
                for(Int_t i = 1; i <= N_BINS_CH10; ++ i)
                {
                    double alpha_i = h_ch10_HSD_sys->GetBinContent(i, j) - h_ch10_HSD->GetBinContent(i, j);
                    double alpha_k = h_ch10_HSD_sys->GetBinContent(k, l) - h_ch10_HSD->GetBinContent(k, l);
                    int output_index = (i - 1)
                                     + N_BINS_CH10 * (j - 1)
                                     + N_BINS_CH10 * N_BINS_CH10 * (k - 1)
                                     + N_BINS_CH10 * N_BINS_CH10 * N_BINS_CH10 * (l - 1);
                    V_PHYS_SYS_CH10.at(output_index) += alpha_i * alpha_k;
                }
            }
        }
    }
    */

    for(Int_t j = 1; j <= N_BINS_CH10; ++ j)
    {
        for(Int_t i = 1; i <= N_BINS_CH10; ++ i)
        {
            double alpha_i = h_ch10_HSD_sys->GetBinContent(i, j) - h_ch10_HSD->GetBinContent(i, j);
            alpha_coeff_ch10.at((i - 1) + N_BINS_CH10 * (j - 1)) = alpha_i;
        }
    }

    for(Int_t j = 1; j <= N_BINS_CH10; ++ j)
    {
        for(Int_t i = 1; i <= N_BINS_CH10; ++ i)
        {
            double delta = h_ch10_HSD_sys->GetBinContent(i, j) - h_ch10_HSD->GetBinContent(i, j);
            chi2_SYS_CH10.at((i - 1) + N_BINS_CH10 * (j - 1)) += std::abs(delta);
        }
    }

}


void reinit_all(const double N_events)
{


    flag_rebuild_ch0_HSD = true;
    flag_rebuild_ch0_SSD = true;
    flag_rebuild_ch0_HSD_sys = false;

    flag_rebuild_ch1_HSD = true;
    flag_rebuild_ch1_SSD = true;
    flag_rebuild_ch1_HSD_sys = false;

    flag_rebuild_ch10_HSD = true;
    flag_rebuild_ch10_SSD = true;
    flag_rebuild_ch10_HSD_sys = false;

    reset_systematic_arrays();
    rebuild_histograms(N_events);



    flag_rebuild_ch0_HSD = false;
    flag_rebuild_ch0_SSD = false;
    flag_rebuild_ch0_HSD_sys = true;

    flag_rebuild_ch1_HSD = false;
    flag_rebuild_ch1_SSD = false;
    flag_rebuild_ch1_HSD_sys = true;

    flag_rebuild_ch10_HSD = false;
    flag_rebuild_ch10_SSD = false;
    flag_rebuild_ch10_HSD_sys = true;

    gSYS_efficiency = 5.55e-2;
    rebuild_histograms(N_events);
    gSYS_efficiency = 0.0;
    systematic_arrays_add(alpha_coeff_ch0_sys1,
                          alpha_coeff_ch1_sys1,
                          alpha_coeff_ch10_sys1);
    

/*
    flag_rebuild_ch0_HSD = false;
    flag_rebuild_ch0_SSD = false;
    flag_rebuild_ch0_HSD_sys = true;

    flag_rebuild_ch1_HSD = false;
    flag_rebuild_ch1_SSD = false;
    flag_rebuild_ch1_HSD_sys = true;

    flag_rebuild_ch10_HSD = false;
    flag_rebuild_ch10_SSD = false;
    flag_rebuild_ch10_HSD_sys = true; 

    gSYS_coherent_energy_shift = 0.2e-2;
    rebuild_histograms(N_events);
    gSYS_coherent_energy_shift = 0.0;
    systematic_arrays_add(alpha_coeff_ch0_sys2,
                          alpha_coeff_ch1_sys2,
                          alpha_coeff_ch10_sys2); 
*/

    flag_rebuild_ch0_HSD = true;
    flag_rebuild_ch0_SSD = true;
    flag_rebuild_ch0_HSD_sys = false;

    flag_rebuild_ch1_HSD = true;
    flag_rebuild_ch1_SSD = true;
    flag_rebuild_ch1_HSD_sys = false;

    flag_rebuild_ch10_HSD = true;
    flag_rebuild_ch10_SSD = true;
    flag_rebuild_ch10_HSD_sys = false;
}


int main()
{

    gStyle->SetLabelFont(43);
    gStyle->SetLabelSize(15);
    gStyle->SetTitleFont(43);
    gStyle->SetTitleSize(20);

    load_psf_data();
    
    //rebuild_histograms();

    //TCanvas *c = new TCanvas("c", "c");
    //h_ch10_HSD->Draw("colz");
    //h_ch10_SSD->Draw("colz");

    

    ///////////////////////////////////////////////////////////////////////////
    // construct systematics
    ///////////////////////////////////////////////////////////////////////////
/*
    reset_systematic_arrays();

    flag_rebuild_ch0_HSD = false;
    flag_rebuild_ch0_SSD = false;
    flag_rebuild_ch0_HSD_sys = true;

    flag_rebuild_ch1_HSD = false;
    flag_rebuild_ch1_SSD = false;
    flag_rebuild_ch1_HSD_sys = true;

    flag_rebuild_ch10_HSD = false;
    flag_rebuild_ch10_SSD = false;
    flag_rebuild_ch10_HSD_sys = true;

    gSYS_efficiency = 5.55e-2;
    rebuild_histograms();
    systematic_arrays_add();
    
    gSYS_efficiency = 0.0;
    gSYS_coherent_energy_shift = 0.2e-2;
    rebuild_histograms();
    systematic_arrays_add(); 

    gSYS_coherent_energy_shift = 0.0;

    flag_rebuild_ch0_HSD = true;
    flag_rebuild_ch0_SSD = true;
    flag_rebuild_ch0_HSD_sys = false;

    flag_rebuild_ch1_HSD = true;
    flag_rebuild_ch1_SSD = true;
    flag_rebuild_ch1_HSD_sys = false;

    flag_rebuild_ch10_HSD = true;
    flag_rebuild_ch10_SSD = true;
    flag_rebuild_ch10_HSD_sys = false;
*/
//    rebuild_histograms();
    // already done

    



    // draw total E for 3, 4 bins

    int N_BINS_MIN = 1;
    int N_BINS_MAX = 100;
    TGraph *g_chi2_ch0 = new TGraph(N_BINS_MAX + 1 - N_BINS_MIN);
    TGraph *g_chi2_ch1 = new TGraph(N_BINS_MAX + 1 - N_BINS_MIN);
    TGraph *g_chi2_ch10 = new TGraph(N_BINS_MAX + 1 - N_BINS_MIN);
    int i = 0;
    for(int N_BINS = N_BINS_MIN; N_BINS <= N_BINS_MAX; ++ N_BINS)
    {
        N_BINS_CH0 = N_BINS;
        N_BINS_CH1 = N_BINS;
        N_BINS_CH10 = N_BINS;

        reinit_all(N2e);
        double chi2_ch0 = calc_chi2(h_ch0_SSD, h_ch0_HSD, alpha_coeff_ch0_sys1, alpha_coeff_ch0_sys2);
        double chi2_ch1 = calc_chi2(h_ch1_SSD, h_ch1_HSD, alpha_coeff_ch1_sys1, alpha_coeff_ch1_sys2);
        double chi2_ch10 = 0.0; // calc_chi2(h_ch10_SSD, h_ch10_HSD, alpha_coeff_ch10_sys1, alpha_coeff_ch10_sys2);

        TString fns_ch0;
        fns_ch0.Form("CH0_%d_bins", N_BINS_CH0);
        draw_hist_ratio_pull(h_ch0_SSD, h_ch0_HSD, &chi2_SYS_CH0, fns_ch0);

        TString fns_ch1;
        fns_ch1.Form("CH1_%d_bins", N_BINS_CH0);
        draw_hist_ratio_pull(h_ch1_SSD, h_ch1_HSD, &chi2_SYS_CH1, fns_ch1);

        g_chi2_ch0->SetPoint(i, (double)N_BINS, chi2_ch0);
        g_chi2_ch1->SetPoint(i, (double)N_BINS, chi2_ch1);
        g_chi2_ch10->SetPoint(i, (double)N_BINS, chi2_ch10);

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



    if(0)
    {
        N_BINS_CH0 = 4;
        reinit_all(N2e);
        std::cout << "chi2=" << calc_chi2(h_ch0_SSD, h_ch0_HSD, alpha_coeff_ch0_sys1, alpha_coeff_ch0_sys2) << std::endl;
        //rebuild_histograms();
//        h_ch0_HSD->Scale(N_events / h_ch0_HSD->Integral());
//        h_ch0_SSD->Scale(N_events / h_ch0_SSD->Integral());
//        h_ch1_HSD->Scale(2.0 * N_events / h_ch1_HSD->Integral());
//        h_ch1_SSD->Scale(2.0 * N_events / h_ch1_SSD->Integral());
//        h_ch10_HSD->Scale(N_events / h_ch10_HSD->Integral());
//        h_ch10_SSD->Scale(N_events / h_ch10_SSD->Integral());
        draw_hist_ratio_pull(h_ch0_SSD, h_ch0_HSD, &chi2_SYS_CH0, "CH0_4bins");
    }




    #if 0
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
            
            const int N_events = N2e;
            rebuild_histograms(N_events);
            //h_ch0_HSD->Scale(N_events / h_ch0_HSD->Integral());
            //h_ch0_SSD->Scale(N_events / h_ch0_SSD->Integral());
            //h_ch1_HSD->Scale(2.0 * N_events / h_ch1_HSD->Integral());
            //h_ch1_SSD->Scale(2.0 * N_events / h_ch1_SSD->Integral());
            //h_ch10_HSD->Scale(N_events / h_ch10_HSD->Integral());
            //h_ch10_SSD->Scale(N_events / h_ch10_SSD->Integral());

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
    #endif

    
    Int_t N_BINS_CH0 = 40;
    Int_t N_BINS_CH1 = 40;
    Int_t N_BINS_CH10 = 40;
    //rebuild_histograms(N2e);


    // draw chi2 as a function of number of events
    #if 0
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
    #endif


    return 0;

}
