
template <typename T , typename M>
void estimate_coef(std::vector<T> indep_var,std::vector<T> dep_var ,M &B_1 ,M &B_0 )
{
    M X_mean = alg_math::Math_Mean(indep_var);
    M Y_mean = alg_math::Math_Mean(dep_var);
    M SS_xy = 0;
    M SS_xx = 0;
    int n = indep_var.size();
    {
        std::vector<T> temp;
        temp  = alg_math::vect_Multiply(indep_var, dep_var);
        SS_xy = std::accumulate(temp.begin(), temp.end(), 0);
        SS_xy = SS_xy - n * X_mean * Y_mean;
    }
    {
        std::vector<T> temp;
        temp = alg_math::vect_Multiply(indep_var, indep_var);
        SS_xx = std::accumulate(temp.begin(), temp.end(), 0);
        SS_xx = SS_xx - n * X_mean * X_mean;
    }
    std::cout<< "SS_xy : " << SS_xy <<std::endl;
    std::cout<< "SS_xx : " << SS_xx <<std::endl;
    std::cout<< "X_mean : " << X_mean <<std::endl;
    std::cout<< "Y_mean : " << Y_mean <<std::endl;
    B_1 =  SS_xy / SS_xx;
    B_0 =  Y_mean - B_1 * X_mean ;
}


template <typename X, typename Z>
class Linear_Regression
{
    private:
    X b_0;
    X b_1;
    
    public:
    void fit(std::vector<std::vector<Z>> dataset)
    {
        estimate_coef(dataset[0],dataset[1] ,b_1 ,b_0 );
    }
    Z predict(const Z & test_data)
    {
        return b_0 + (b_1 * test_data); 
    }
};
