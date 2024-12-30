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
    double T_;            // Temps jusqu'à la maturité
    double K_;            // Prix d'exercice
    double sigma_;        // Volatilité
    double r_;            // Taux sans risque
    bool isCall_;         // Option call (true) ou put (false)
    bool isAmerican_;     // Option américaine (true) ou européenne (false)

    double deltaS_;       // Pas en espace
    double deltaT_;       // Pas en temps

    vector<double> S_;    // Grille des prix du sous-jacent
    vector<vector<double>> P_;    // Grille des prix de l'option
    vector<double> alpha_, beta_, gamma_;
    vector<double> lower_, diag_, upper_;

    void initialize_matrices() {
        for (int i = 0; i <= N_; ++i) {
            S_[i] = i * deltaS_;
            P_[i][M_] = isCall_ ? max(S_[i] - K_, 0.0) : max(K_ - S_[i], 0.0);
        }

        for (int j = 0; j <= M_; ++j) {
            double t = T_ * (M_ - j) / M_;
            P_[0][j] = isCall_ ? 0.0 : K_ * exp(-r_ * t);
            P_[N_][j] = isCall_ ? (S_[N_] - K_) * exp(-r_ * t) : 0.0;
        }

        for (int i = 1; i <= N_; ++i) {
            alpha_[i] = -0.25 * deltaT_ * ((sigma_ * sigma_) * (i * i) - r_ * i);
            beta_[i] = 0.5 * deltaT_ * ((sigma_ * sigma_) * (i * i) + r_);
            gamma_[i] = -0.25 * deltaT_ * ((sigma_ * sigma_) * (i * i) + r_ * i);
        }

        for (int i = 1; i < N_; ++i) {
            lower_[i - 1] = alpha_[i];
            diag_[i - 1] = 1.0 + beta_[i];
            upper_[i - 1] = gamma_[i];
        }
    }

    // Résolution d'un système tridiagonal (algorithme de Thomas)
    vector<double> solveTridiagonal(const vector<double>& lower, const vector<double>& diag, const vector<double>& upper, const vector<double>& rhs) {
        vector<double> c_prime(N_ - 1, 0.0), d_prime(N_ - 1, 0.0), x(N_ - 1, 0.0);

        // Vérification des valeurs initiales
        if (diag[0] == 0) {
            cerr << "Erreur: Division par zéro dans le solveur tridiagonal.\n";
            exit(1);
        }

        // Forward elimination
        c_prime[0] = upper[0] / diag[0];
        d_prime[0] = rhs[0] / diag[0];
        for (int i = 1; i < N_ - 1; ++i) {
            double denom = diag[i] - lower[i - 1] * c_prime[i - 1];
            if (denom == 0) {
                cerr << "Erreur: Division par zéro dans le solveur tridiagonal.\n";
                exit(1);
            }
            c_prime[i] = upper[i] / denom;
            d_prime[i] = (rhs[i] - lower[i - 1] * d_prime[i - 1]) / denom;
        }

        // Backward substitution
        x[N_ - 2] = d_prime[N_ - 2];
        for (int i = N_ - 3; i >= 0; --i) {
            x[i] = d_prime[i] - c_prime[i] * x[i + 1];
        }

        return x;
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

    void handle_american_exercise(int t) {
        if (isAmerican_) {
            for (int i = 1; i < N_; ++i) {
                double intrinsic_value = isCall_ ? max(S_[i] - K_, 0.0) : max(K_ - S_[i], 0.0);
                P_[i][t] = max(P_[i][t], intrinsic_value);
            }
        }
    }

public:
    BlackScholesSolver(int N, int M, double S0, double T, double K, double sigma, double r, bool isCall, bool isAmerican)
        : N_(N), M_(M), S0_(S0), T_(T), K_(K), sigma_(sigma), r_(r), isCall_(isCall), isAmerican_(isAmerican),
          deltaS_(2 * S0 / N), deltaT_(T / M), S_(N_ + 1), P_(N + 1, vector<double>(M + 1, 0.0)),
          alpha_(N + 1), beta_(N + 1), gamma_(N + 1), lower_(N_ - 1), diag_(N_ - 1), upper_(N_ - 1) {}

    void solve() {
        // Initialisation des grilles
        initialize_matrices();
        for (int t = M_ - 1; t >= 0; --t) {
            vector<double> rhs(N_ - 1);

            // Construction des vecteurs de la matrice tridiagonale et du vecteur RHS
            for (int i = 1; i < N_ - 1; ++i) {
                rhs[i - 1] = -alpha_[i + 1] * P_[i][t + 1] + (1 - beta_[i + 1]) * P_[i + 1][t + 1] - gamma_[i + 1] * P_[i + 2][t + 1];
            }

            // Ajout des conditions aux limites
            rhs[0] += lower_[0] * P_[0][t];
            rhs[N_ - 2] += upper_[N_ - 2] * P_[N_][t];

            // Résolution du système tridiagonal
            vector<double> solution = solveTridiagonal(lower_, diag_, upper_, rhs);

            // Mise à jour des valeurs
            for (int i = 1; i < N_; ++i) {
                P_[i][t] = solution[i - 1];
            }

            // Gestion de l'exercice anticipé pour les options américaines
            handle_american_exercise(t);
        }
    }

    void calculate_greeks() {
        cout << "Prix de l'option et ses Grecs (t=0) :\n";
        int idx = static_cast<int>(S0_ / deltaS_);
        double delta = (P_[idx + 1][0] - P_[idx - 1][0]) / (2 * deltaS_);
        double gamma = (P_[idx + 1][0] - 2 * P_[idx][0] + P_[idx - 1][0]) / (deltaS_ * deltaS_);
        cout << "S = " << setw(6) << S_[idx]
             << " | P = " << setw(8) << P_[idx][0]
             << " | Delta = " << setw(8) << delta
             << " | Gamma = " << setw(8) << gamma << endl;
    }

    void compare_with_analytical() {
        cout << "Comparaison avec la solution analytique :\n";
        for (int i = 0; i <= N_; ++i) {
            double analytical = black_scholes_analytical(S_[i], K_, T_, r_, sigma_, isCall_);
            cout << "S = " << setw(6) << S_[i]
                 << " | P_num = " << setw(8) << P_[i][0]
                 << " | P_ana = " << setw(8) << analytical << endl;
        }
    }
};

int main() {
    // Paramètres
    int N = 200; // Augmenter la résolution spatiale
    int M = 200; // Augmenter la résolution temporelle
    double S0 = 100.0;
    double T = 1.0;
    double K = 100.0;
    double sigma = 0.2;
    double r = 0.05;
    bool isCall = true;
    bool isAmerican = true; // Gestion des options américaines

    // Solveur
    BlackScholesSolver solver(N, M, S0, T, K, sigma, r, isCall, isAmerican);
    solver.solve();
    solver.calculate_greeks();
    solver.compare_with_analytical();

    return 0;
}
