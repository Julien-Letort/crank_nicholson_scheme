#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

class BlackScholesSolver {
private:
    int N_;               // Nombre de points en espace
    int M_;               // Nombre de points en temps
    double S0_;           // Valeur maximale du sous-jacent
    double Smax_;
    double T_;            // Temps jusqu'à la maturité
    double date_;
    double K_;            // Prix d'exercice
    double sigma_;        // Volatilité
    double r_;            // Taux sans risque
    bool isCall_;         // Option call (true) ou put (false)
    bool isAmerican_;     // Option américaine (true) ou européenne (false)

    double deltaS_;       // Pas en espace
    double deltaT_;       // Pas en temps

    vector<vector<double>> output;    // Grille des prix de l'option
    vector<double> A, B, C;
    vector<double> tridiag_inf_1, tridiag_diag_1, tridiag_sup_1, tridiag_inf_2, tridiag_diag_2, tridiag_sup_2;

    vector<double> tridiag_matmult(vector<double> X, double cond1, double cond2){
        vector<double> Y(N_+1);
        for (int j = 1; j<N_; j++){
            Y[j] = X[j]*tridiag_diag_2[j-1] + X[j-1]*tridiag_inf_2[j-1] + X[j+1]*tridiag_sup_2[j-1];
        }
        Y[0] = cond1;
        Y[N_] = cond2;
        return Y;
    }

    vector<double> algo_thomas(vector<double> Y){
        double coeffs_1[N_-1];
        double coeffs_2[N_-1];
        vector<double> X(N_+1);
        X[0] = Y[0];
        X[N_] = Y[N_];
        vector<double> tmp_(N_-1);
        for (int i=0; i<N_-1; i++){
            tmp_[i] = Y[i+1];
        }
        tmp_[0] = tmp_[0] - tridiag_inf_1[0]*X[0];
        tmp_[N_-2] = tmp_[N_-2] - tridiag_sup_1[N_-2]*X[N_];

        coeffs_1[0] = tridiag_sup_1[0]/tridiag_diag_1[0];
        coeffs_2[0] = tmp_[0]/tridiag_diag_1[0];
        for (int i = 1; i < N_-1; i++){
            coeffs_1[i] = tridiag_sup_1[i]/(tridiag_diag_1[i] - tridiag_inf_1[i]*coeffs_1[i-1]);
            coeffs_2[i] = (tmp_[i] - tridiag_inf_1[i]*coeffs_2[i-1])/(tridiag_diag_1[i] - tridiag_inf_1[i]*coeffs_1[i-1]);
        }
        X[N_-1] = coeffs_2[N_-2];
        for (int k = N_-2; k>0; k--){
            X[k] = coeffs_2[k-1] - coeffs_1[k-1]*X[k+1];
        }
        return X;
    }

    void initialize_matrices() {
        for (int j=1; j<N_; j++){
            A[j-1] = (deltaT_/4)*( pow((sigma_*j),2.) - r_*j );
            tridiag_inf_1[j-1] = -A[j-1];
            tridiag_inf_2[j-1] = A[j-1];
            B[j-1] = (-deltaT_/2)*( pow((sigma_*j),2.) + r_ );
            tridiag_diag_1[j-1] = 1. - B[j-1];
            tridiag_diag_2[j-1] = 1. + B[j-1];
            C[j-1] = (deltaT_/4)*( pow((sigma_*j),2.) + r_*j );
            tridiag_sup_1[j-1] = -C[j-1];
            tridiag_sup_2[j-1] = C[j-1];
        }
    }

    double black_scholes_analytical(double S, double K, double T, double r, double sigma, bool isCall) {
        if (T == 0.0) {
            return isCall ? max(S - K, 0.0) : max(K - S, 0.0);
        }

        double d1 = (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
        double d2 = d1 - sigma * sqrt(T);
        double call_price = S * 0.5 * (1.0 + erf(d1 / sqrt(2.0))) - K * exp(-r * T) * 0.5 * (1.0 + erf(d2 / sqrt(2.0)));

        return isCall ? call_price : call_price + K * exp(-r * T) - S;
    }

public:
    BlackScholesSolver(int N, int M, double S0, double T, double K, double sigma, double r, bool isCall, bool isAmerican,double t=0)
        : N_(N), M_(M), S0_(S0), T_(T), K_(K), sigma_(sigma), r_(r), isCall_(isCall), isAmerican_(isAmerican),Smax_(2*max(S0,K)*exp(r*T)),
          deltaS_(Smax_ / N), deltaT_(T / M), output(M+1),
          A(N - 1), B(N - 1), C(N - 1), tridiag_inf_1(N-1), tridiag_diag_1(N-1), tridiag_sup_1(N-1), tridiag_inf_2(N-1), tridiag_diag_2(N-1), tridiag_sup_2(N-1),date_(t){}
    void solve_test(){
        if (isAmerican_){
            solve_amer();
        }
        else{
            solve_euro();
        }
    }
    void solve_euro(){
        initialize_matrices();
        vector<double> tmp(N_+1);
        for (int j = 0; j<N_+1; j++){
            tmp[j] = isCall_ ? max(j*deltaS_-K_, double(0.)) : max(K_ - j*deltaS_, double(0.));
        }
        output[M_] = tmp;
        for (int i = 0; i <M_+1 ;i++){
            output[i].resize(N_+1);
            output[i][0] = isCall_ ? 0 : K_*exp(r_*(T_ - i*deltaT_));
            output[i][N_] = isCall_ ? (Smax_-K_)*exp(r_*(T_ - i*deltaT_)) : 0;
        }
        for (int i =M_-1; i>=0; i--){
            tmp = tridiag_matmult(output[i+1], output[i][0], output[i][N_]);
            output[i] = algo_thomas(tmp);
        }
    }

    void solve_amer(){
        initialize_matrices();
        vector<double> tmp(N_+1);
        for (int j = 0; j<N_+1; j++){
            tmp[j] =isCall_ ? max(j*deltaS_-K_, double(0.)) : max(K_ - j*deltaS_, double(0.));
        }
        output[M_] = tmp;
        for (int i = 0; i <M_+1 ;i++){
            output[i].resize(N_+1);
            output[i][0] = isCall_ ? 0 : K_*exp(r_*(T_ - i*deltaT_));
            output[i][N_] = isCall_ ? (Smax_-K_)*exp(r_*(T_ - i*deltaT_)) : 0;
        }
        for (int i =M_-1; i>=0; i--){
            tmp = tridiag_matmult(output[i+1], output[i][0], output[i][N_]);
            output[i] = algo_thomas(tmp);
            for (int j = 0; j<N_+1; j++){
                output[i][j] = isCall_ ? max(output[i][j], double( (j*deltaS_-K_))) : max(output[i][j], double( (K_-j*deltaS_)));
            }
        }
    }

    void calculate_greeks() {
        cout << "Prix de l'option et ses Grecs (t=0) :\n";
        int idx_t = static_cast<int>(date_ / deltaT_);
        int idx_S = static_cast<int>(S0_ / deltaS_);
        vector<vector<double>> delta(M_+1);
        vector<vector<double>> gamma(M_+1);
        vector<vector<double>> theta(M_);
        for (int i =0; i<M_+1; i++){
            delta[i].resize(N_);
            gamma[i].resize(N_-1);
            for (int j = 0; j<N_; j++){
                delta[i][j] = (output[i][j+1] - output[i][j])/deltaS_;
            }
            for (int j = 0; j<N_-1; j++){
                gamma[i][j] = (output[i][j+2] - 2*output[i][j+1] + output[i][j])/(deltaS_*deltaS_);
            }
        }
        for (int i =0; i<M_; i++){
            theta[i].resize(N_+1);
            for (int j = 0; j<N_+1; j++){
                theta[i][j] = (output[i+1][j] - output[i][j])/deltaT_;
            }
        }
        double delta_r =0.0001;
        double delta_sigma=0.0001;
        BlackScholesSolver solver_rdiff(N_, M_, S0_, T_, K_, sigma_, r_ + delta_r, isCall_, isAmerican_, date_);
        BlackScholesSolver solver_sigmadiff(N_, M_, S0_, T_, K_, sigma_ + delta_sigma, r_, isCall_, isAmerican_, date_);
        solver_rdiff.solve_test();
        solver_sigmadiff.solve_test();
        double rho = (solver_rdiff.output[idx_t][idx_S]- output[idx_t][idx_S])/delta_r;
        double vega = (solver_sigmadiff.output[idx_t][idx_S]- output[idx_t][idx_S])/delta_sigma;
        cout << "S = " << setw(6) << S0_
             << " | P = " << setw(8) << output[idx_t][idx_S]
             << " | Delta = " << setw(8) << delta[idx_t][idx_S]
             << " | Gamma = " << setw(8) << gamma[idx_t][idx_S] 
             << " | Theta journalier= " << setw(8) << theta[idx_t][idx_S]/365
             << " | rho = " << setw(8) << rho
             << " | vega = " << setw(8) << vega             
             << endl;

        cout << "\nTableau des valeurs de l'option en fonction du sous-jacent :\n";
        cout << "S\tPrix de l'option\tDelta de l'option\n";
        for (int i = 0; i <= N_; ++i) {
            cout << i * deltaS_ << "\t" << output[idx_t][i] <<  "\t"<< delta[idx_t][i] << endl;
        }
        
        
    }

    void compare_with_analytical() {
        cout << "Comparaison avec la solution analytique :\n";
        for (int i = 0; i <= N_; ++i) {
            double analytical = black_scholes_analytical(i*deltaS_, K_, T_, r_, sigma_, isCall_);
            cout << "S = " << setw(6) << i*deltaS_
                 << " | P_num = " << setw(8) << output[i][0]
                 << " | P_ana = " << setw(8) << analytical << endl;
        }
    }
};

int main(int argc, char* argv[]) {
    if (argc != 11) {
        cerr << "Usage: test3.exe N M S0 T t K sigma r IsCall IsAmerican" << endl;
        return 1;
    }

    int N = atof(argv[1]);
    int M = atof(argv[2]);
    double S0 = atof(argv[3]);
    double T = atof(argv[4]);
    double t = atof(argv[5]);
    double K = atof(argv[6]);
    double sigma = atof(argv[7]);
    double r = atof(argv[8]);
    bool isCall = (atof(argv[9]) == 1);
    bool isAmerican = (atof(argv[10]) == 1);

    BlackScholesSolver solver(N, M, S0, T, K, sigma, r, isCall, isAmerican, t);
    solver.solve_test();
    solver.calculate_greeks();

    return 0;
}
