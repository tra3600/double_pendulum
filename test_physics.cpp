// Tests de validation physique : g++ -O2 -std=c++17 test_physics.cpp -o test_physics && ./test_physics
#include <cstdio>
#include <functional>
#include "pendulum.hpp"
using namespace dp;

static int fails = 0;
static void check(const char* name, bool ok, double val, double ref) {
    std::printf("[%s] %-62s %.3e (réf %.3e)\n", ok ? " OK " : "FAIL", name, val, ref);
    if (!ok) ++fails;
}
static State run(const Params& P, State s, double T, double h, State (*step)(const Params&, const State&, double)) {
    const int n = int(T / h + 0.5);
    for (int i = 0; i < n; ++i) s = step(P, s, h);
    return s;
}

int main() {
    Params P;
    State s0{2.0, 1.0, 0.0, 0.0};

    // 1. Hamiltonien (via p) == T+V (via ω) : cohérence Legendre
    State s{1.1, -0.7, 0.9, 1.3};
    check("H(θ,p) == T+V (transformée de Legendre)", std::fabs(hamiltonian(P, s) - energy(P, s)) < 1e-12,
          std::fabs(hamiltonian(P, s) - energy(P, s)), 0);

    // 2. fromMomenta inverse de momenta
    auto p = momenta(P, s);
    State r = fromMomenta(P, s[0], s[1], p[0], p[1]);
    check("ω -> p -> ω aller-retour", std::fabs(r[2] - s[2]) + std::fabs(r[3] - s[3]) < 1e-12,
          std::fabs(r[2] - s[2]) + std::fabs(r[3] - s[3]), 0);

    // 3. Équation d'Euler-Lagrange vérifiée par différences finies sur L
    //    d/dt(∂L/∂ω) = ∂L/∂θ  (sans frottement)
    {
        auto Lag = [&](const State& q) { return kinetic(P, q) - potential(P, q); };
        State q{0.8, -0.3, 0.5, 1.2};
        const double e = 1e-6, dt = 1e-6;
        auto dLdth = [&](int i) { State a = q, b = q; a[i] += e; b[i] -= e; return (Lag(a) - Lag(b)) / (2 * e); };
        State qn = rk4(P, q, dt), qp = rk4(P, q, -dt);
        double worst = 0;
        for (int i = 0; i < 2; ++i) {
            double pn = momenta(P, qn)[i], pp = momenta(P, qp)[i];
            worst = std::fmax(worst, std::fabs((pn - pp) / (2 * dt) - dLdth(i)));
        }
        check("Euler-Lagrange : dp/dt = ∂L/∂θ (différences finies)", worst < 1e-5, worst, 0);
    }

    // 4. Petits angles : fréquences propres analytiques (m1=m2, l1=l2 -> ω²=(2±√2)g/l)
    {
        Modes m = smallAngleModes(P);
        double ex_slow = (2 - std::sqrt(2.0)) * P.g, ex_fast = (2 + std::sqrt(2.0)) * P.g;
        check("mode lent  ω² = (2-√2) g/l", std::fabs(m.w2_slow - ex_slow) < 1e-9, m.w2_slow, ex_slow);
        check("mode rapide ω² = (2+√2) g/l", std::fabs(m.w2_fast - ex_fast) < 1e-9, m.w2_fast, ex_fast);
        // Simulation en mode propre lent : θ1 = a, θ2 = √2 a (vecteur propre) -> θ1(t) = a cos(ω t)
        const double a = 1e-3, w = std::sqrt(m.w2_slow), T = 3.0;
        State q{a, std::sqrt(2.0) * a, 0, 0};
        q = run(P, q, T, 1e-3, rk4);
        double err = std::fabs(q[0] - a * std::cos(w * T)) / a;
        check("simulation vs solution analytique (mode lent)", err < 1e-3, err, 0);
    }

    // 5. Conservation de l'énergie : RK4, DP5 adaptatif, symplectique ; Euler diverge
    {
        const double E0 = energy(P, s0), T = 20.0;
        double drift;
        State a = run(P, s0, T, 1e-3, rk4);
        drift = std::fabs(energy(P, a) - E0) / std::fabs(E0);
        check("dérive d'énergie RK4  (h=1e-3, T=20 s)", drift < 1e-8, drift, 1e-8);
        double h = 0; State b = s0;
        for (int i = 0; i < 2000; ++i) b = advanceAdaptive(P, b, 0.01, h, 1e-12);
        drift = std::fabs(energy(P, b) - E0) / std::fabs(E0);
        check("dérive d'énergie DP5 adaptatif (tol=1e-12)", drift < 1e-9, drift, 1e-9);
        State c = run(P, s0, T, 1e-3, symplecticMidpoint);
        drift = std::fabs(energy(P, c) - E0) / std::fabs(E0);
        check("dérive d'énergie point-milieu symplectique", drift < 1e-4, drift, 1e-4);
        State d = run(P, s0, T, 1e-3, euler);
        double de = std::fabs(energy(P, d) - E0) / std::fabs(E0);
        check("Euler explicite : dérive énorme (attendu, contre-exemple)", de > 1e-2, de, 1e-2);
    }

    // 6. Ordre de convergence de RK4 (erreur ∝ h⁴) face à une référence DP5 serrée
    {
        double h = 0; State ref = s0;
        for (int i = 0; i < 100; ++i) ref = advanceAdaptive(P, ref, 0.01, h, 1e-13);  // T=1
        auto e = [&](double step) { State q = run(P, s0, 1.0, step, rk4); return std::fabs(q[0] - ref[0]) + std::fabs(q[1] - ref[1]); };
        double e1 = e(0.01), e2 = e(0.005), order = std::log2(e1 / e2);
        check("ordre de convergence de RK4 ≈ 4", order > 3.7 && order < 4.3, order, 4.0);
    }

    // 7. Tensions : Newton vérifié indépendamment par accélérations cartésiennes
    {
        State q{0.9, -0.4, 1.1, 0.7};
        const double dt = 1e-5;
        auto pos = [&](const State& z) { return cartesian(P, z); };
        State qp = rk4(P, q, dt), qm = rk4(P, q, -dt);
        Cart c0 = pos(q), c1 = pos(qp), c2 = pos(qm);
        double ax2 = (c1.x2 - 2 * c0.x2 + c2.x2) / (dt * dt), ay2 = (c1.y2 - 2 * c0.y2 + c2.y2) / (dt * dt);
        double ax1 = (c1.x1 - 2 * c0.x1 + c2.x1) / (dt * dt), ay1 = (c1.y1 - 2 * c0.y1 + c2.y1) / (dt * dt);
        Tension T = tensions(P, q);
        // m2 a2 = (0,-m2 g) + T2 (-sinθ2, cosθ2)
        double rx = P.m2 * ax2 - T.T2 * (-std::sin(q[1])), ry = P.m2 * ay2 + P.m2 * P.g - T.T2 * std::cos(q[1]);
        // m1 a1 = (0,-m1 g) + T1(-sinθ1, cosθ1) - T2(-sinθ2, cosθ2)
        double sx = P.m1 * ax1 - (T.T1 * -std::sin(q[0]) - T.T2 * -std::sin(q[1]));
        double sy = P.m1 * ay1 + P.m1 * P.g - (T.T1 * std::cos(q[0]) - T.T2 * std::cos(q[1]));
        double worst = std::fmax(std::fmax(std::fabs(rx), std::fabs(ry)), std::fmax(std::fabs(sx), std::fabs(sy)));
        check("2e loi de Newton vs tensions calculées", worst < 1e-3, worst, 0);
    }

    // 8. Frottement : dE/dt = -c(ω1²+ω2²) (théorème de l'énergie avec dissipation de Rayleigh)
    {
        Params Q = P; Q.c = 0.05;
        State q = s0;
        double E0 = energy(Q, q), diss = 0;
        const double h = 1e-3;
        for (int i = 0; i < 5000; ++i) {
            State n = rk4(Q, q, h);
            diss += Q.c * 0.5 * ((q[2] * q[2] + q[3] * q[3]) + (n[2] * n[2] + n[3] * n[3])) * h;
            q = n;
        }
        double err = std::fabs((E0 - energy(Q, q)) - diss);
        check("bilan d'énergie avec frottement : ΔE = -∫c|ω|²dt", err < 1e-5, err, 0);
    }

    // 9. Lyapounov : régulier (petit angle) ≈ 0 ; chaotique (grand angle) > 0
    {
        Lyapunov a; a.init(State{0.05, 0.05, 0, 0});
        Lyapunov b; b.init(State{2.5, 2.0, 0, 0});
        for (int i = 0; i < 40000; ++i) { a.step(P, 1e-3); b.step(P, 1e-3); }
        check("Lyapounov ≈ 0 pour petits angles (régulier)", a.value() < 0.3, a.value(), 0);
        check("Lyapounov > 1 pour grands angles (chaos)", b.value() > 1.0, b.value(), 1.0);
    }

    std::printf("\n%s (%d échec(s))\n", fails ? "ÉCHEC" : "TOUS LES TESTS PASSENT", fails);
    return fails ? 1 : 0;
}
