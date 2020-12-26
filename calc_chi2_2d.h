#ifndef CALC_CHI2_2D_H
#define CALC_CHI2_2D_H



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

double calc_chi2(
    TH2D* h_data,
    TH2D* h_model,
    const std::vector<double> &sys1,
    const std::vector<double> &sys2)
{
    std::cout << __func__ << " (2d)" << std::endl;

    std::vector<double> delta; // d - m
    std::vector<double> sigma; // sigma i

    int count = 0;
    for(Int_t j = 1; j <= h_model->GetNbinsY(); ++ j)
    {
        for(Int_t i = 1; i <= h_model->GetNbinsX(); ++ i)
        {
            double content_model = h_model->GetBinContent(i, j);
            double content_data = h_data->GetBinContent(i, j);
            double error_model = std::sqrt(content_model);

            //std::cout << i << " " << j << " > " << content_model << std::endl;
            
            if(content_model <= VERYSMALL) continue;
            ++ count;

            sigma.push_back(error_model);
            delta.push_back(content_data - content_model);
        }
    }
    //std::cout << "count=" << count << std::endl;
    TMatrixD *matrix = new TMatrixD(count, count);

    /*
    int count_sqrt = 0;
    for(Int_t i = 1; i <= h_model->GetNbinsX(); ++ i)
    {
        double content_model = h_model->GetBinContent(i, i);
        if(content_model <= VERYSMALL) continue;
        ++ count_sqrt;
    }
    std::cout << "count_sqrt=" << count_sqrt << std::endl;
    */

    {
        int index2 = 0;

        int l = 0;
        for(int ll = 1; ll <= h_model->GetNbinsY(); ++ ll)
        {
            int k = 0;
            for(int kk = 1; kk <= h_model->GetNbinsX(); ++ kk)
            {
                if(h_model->GetBinContent(kk, ll) <= VERYSMALL)
                {
                    //std::cout << "continue" << std::endl;
                    continue;
                }

                int index = 0;

                int j = 0;
                for(int jj = 1; jj <= h_model->GetNbinsY(); ++ jj)
                //for(int j = 0; j < matrix->GetNrows(); ++ j)
                {
                    int i = 0;
                    for(int ii = 1; ii <= h_model->GetNbinsX(); ++ ii)
                    //for(int i = 0; i < matrix->GetNcols(); ++ i)
                    {
                        if(h_model->GetBinContent(ii, jj) <= VERYSMALL)
                        {
                            //std::cout << "continue" << std::endl;
                            continue;
                        }

                        double c = 0.0;
                        if((i == k) && (j == l))
                        {
                            //double sigma_i = sigma.at((ii - 1) + h_model->GetNbinsX() * (jj - 1));
                            //double sigma_j = sigma.at((kk - 1) + h_model->GetNbinsX() * (ll - 1));
                            double sigma_i = sigma.at(index);
                            double sigma_j = sigma.at(index);
                            c += sigma_i * sigma_j;
                        }

                        if(g_CHI2_SYSZERO == true)
                        {
                            // ignore all systematic terms
                            std::cout << "IGNORE: SYSZERO" << std::endl;
                        }
                        else
                        {
                            if((g_CHI2_OFFDIAGZERO == true) && (index != index2))
                            {
                                // ignore off diagonal terms
                                std::cout << "IGNORE: OFFDIAG" << std::endl;
                            }
                            else
                            {
                                int index_i = (ii - 1) + h_model->GetNbinsX() * (jj - 1);
                                int index_j = (kk - 1) + h_model->GetNbinsX() * (ll - 1);
                                double alpha_i = sys1.at(index_i);
                                double alpha_j = sys1.at(index_j);
                                c += alpha_i * alpha_j;
                                alpha_i = sys2.at(index_i);
                                alpha_j = sys2.at(index_j);
                                c += alpha_i * alpha_j;
                                std::cout << "ACCEPT" << std::endl;
                            }
                        }

                        //std::cout << "ii=" << ii << " jj=" << jj << " kk=" << kk << " ll=" << ll << " i=" << i << " j=" << j << std::endl;
                        //std::cout << "j=" << j << " i=" << i << " " << c << std::endl;
                        //std::cout << l + count_sqrt * k << " "  << i + count_sqrt * j << std::endl;
                        //std::cout << index2 << " "  << index << std::endl;
                        matrix->operator[](index2).operator[](index) = c;

                        ++ i;

                        ++ index;
                    }

                    ++ j;
                }

                ++ k;

                ++ index2;
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


#endif // CALC_CHI2_2D_H
