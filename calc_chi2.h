#ifndef CALC_CHI2_H
#define CALC_CHI2_H



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
        if(content_model <= VERYSMALL) continue;

        delta.push_back(content_data - content_model);
        sigma.push_back(error_model);
        
//        std::cout << "i=" << i << " data=" << content_data << " model=" << content_model << std::endl;
    }
    TMatrixD *matrix = new TMatrixD(delta.size(), delta.size());

    {
        int j = 0;
        for(int jj = 1; jj <= h_model->GetNbinsX(); ++ jj)
        //for(int j = 0; j < matrix->GetNrows(); ++ j)
        {
            if(h_model->GetBinContent(jj) <= VERYSMALL) continue; 

            int i = 0;
            for(int ii = 1; ii <= h_model->GetNbinsX(); ++ ii)
            //for(int i = 0; i < matrix->GetNcols(); ++ i)
            {
                if(h_model->GetBinContent(ii) <= VERYSMALL) continue;

                double c = 0.0;
                if(i == j)
                {
                    double sigma_i = sigma.at(i);
                    double sigma_j = sigma.at(j);
                    c += sigma_i * sigma_j;
                }

                if(g_CHI2_SYSZERO == true)
                {
                    // ignore all systematic terms
                }
                else
                {
                    if((g_CHI2_OFFDIAGZERO == true) && (i != j))
                    {
                        // ignore off diagonal terms
                    }
                    else
                    {
                        double alpha_i = sys1.at(ii - 1);
                        double alpha_j = sys1.at(jj - 1);
                        //if(i == j)
                        {
                        c += alpha_i * alpha_j;
                        }
                        alpha_i = sys2.at(ii - 1);
                        alpha_j = sys2.at(jj - 1);
                        //if(i == j)
                        {
                        c += alpha_i * alpha_j;
                        }
                    }
                }

                //std::cout << "j=" << j << " i=" << i << " " << c << std::endl;
                matrix->operator[](j).operator[](i) = c;

                #if 0
                if(i == j)
                {
                //    std::cout << "diagonal: i=" << i << " " << c << std::endl;
                    std::cout << "diagonal, error composition: i=" << i
                              << " sigma=" << sigma.at(i) * sigma.at(j)
                              << " alpha=" << sys1.at(ii - 1) * sys1.at(jj - 1)
                              << " + " << sys2.at(ii - 1) * sys2.at(jj - 1)
                              << " | model=" << h_model->GetBinContent(ii)
                              << std::endl;
                }
                #endif

                ++ i;
            }

            ++ j;
        }

    }
    matrix->Invert();

    double chi2 = 0.0;
    for(int j = 0; j < matrix->GetNcols(); ++ j)
    {
        for(int i = 0; i < matrix->GetNrows(); ++ i)
        {
            double v_i = delta.at(i);
            double V_ij = matrix->operator[](j).operator[](i);
            double v_j = delta.at(j);
            //std::cout << "i=" << i << " j=" << j << " v_i=" << v_i << " V_ij=" << V_ij << " v_j=" << v_j << std::endl;
            chi2 += v_i * V_ij * v_j;
            /*
            if(i == j)
            {
                std::cout << "chi2=" << chi2
                          << " delta_i=" << delta.at(i)
                          << " delta_j=" << delta.at(j)
                          << std::endl;
            }
            */
        }
    }

    return chi2;
}



#endif // CALC_CHI2_H
