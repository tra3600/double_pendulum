// pendulum.hpp — Physique du pendule double (sans dépendance graphique)
//
// Mécanique analytique complète :
//   * Lagrangien       L = T - V
//   * Euler-Lagrange   d/dt(∂L/∂θ̇) - ∂L/∂θ = Q   (Q : forces non conservatives)
//   * Hamiltonien      H(θ, p), moments conjugués p = ∂L/∂θ̇
//   * Dissipation      fonction de Rayleigh  R = ½ c (θ̇1² + θ̇2²)
//   * Chaos            exposant de Lyapounov (algorithme de Benettin)
//   * Intégrateurs     Euler, RK4, RK45 Dormand-Prince adaptatif, point-milieu symplectique
//   * Théorie linéaire modes propres aux petits angles (validation analytique)
#pragma once
#include <array>
#include <cmath>

namespace dp {

constexpr double PI = 3.14159265358979323846;

// Paramètres physiques
struct Params {
    double m1 = 1.0, m2 = 1.0;   // masses (kg)
    double l1 = 1.0, l2 = 1.0;   // longueurs (m)
    double g  = 9.81;            // pesanteur (m/s²)
    double c  = 0.0;             // frottement visqueux articulaire (N·m·s)
};

// État : (θ1, θ2, ω1, ω2) — angles par rapport à la verticale descendante
using State = std::array<double, 4>;

// ---------------------------------------------------------------------------
// Énergies
//   T = ½(m1+m2) l1² ω1² + ½ m2 l2² ω2² + m2 l1 l2 ω1 ω2 cos(θ1-θ2)
//   V = -(m1+m2) g l1 cosθ1 - m2 g l2 cosθ2
// ---------------------------------------------------------------------------
inline double kinetic(const Params& P, const State& s) {
    const double d = s[0] - s[1];
    return 0.5 * (P.m1 + P.m2) * P.l1 * P.l1 * s[2] * s[2]
         + 0.5 * P.m2 * P.l2 * P.l2 * s[3] * s[3]
         + P.m2 * P.l1 * P.l2 * s[2] * s[3] * std::cos(d);
}
inline double potential(const Params& P, const State& s) {
    return -(P.m1 + P.m2) * P.g * P.l1 * std::cos(s[0])
           - P.m2 * P.g * P.l2 * std::cos(s[1]);
}
inline double energy(const Params& P, const State& s) { return kinetic(P, s) + potential(P, s); }

// Moments conjugués p_i = ∂L/∂θ̇_i
inline std::array<double, 2> momenta(const Params& P, const State& s) {
    const double d = s[0] - s[1], cd = std::cos(d);
    return {(P.m1 + P.m2) * P.l1 * P.l1 * s[2] + P.m2 * P.l1 * P.l2 * s[3] * cd,
            P.m2 * P.l2 * P.l2 * s[3] + P.m2 * P.l1 * P.l2 * s[2] * cd};
}
// Inverse : (θ, p) -> ω  (utile pour initialiser depuis l'espace des phases)
inline State fromMomenta(const Params& P, double t1, double t2, double p1, double p2) {
    const double d = t1 - t2, cd = std::cos(d);
    const double a = (P.m1 + P.m2) * P.l1 * P.l1, b = P.m2 * P.l1 * P.l2 * cd, e = P.m2 * P.l2 * P.l2;
    const double det = a * e - b * b;
    return {t1, t2, (e * p1 - b * p2) / det, (a * p2 - b * p1) / det};
}
// Hamiltonien H = T + V (identique à l'énergie, exprimé via p ; sert de contrôle croisé)
inline double hamiltonian(const Params& P, const State& s) {
    const auto p = momenta(P, s);
    const double d = s[0] - s[1], cd = std::cos(d);
    const double a = (P.m1 + P.m2) * P.l1 * P.l1, b = P.m2 * P.l1 * P.l2 * cd, e = P.m2 * P.l2 * P.l2;
    const double det = a * e - b * b;
    const double T = 0.5 * (e * p[0] * p[0] - 2 * b * p[0] * p[1] + a * p[1] * p[1]) / det;
    return T + potential(P, s);
}

// ---------------------------------------------------------------------------
// Équations d'Euler-Lagrange  M(θ) θ̈ = f(θ, θ̇)
//   M = [ (m1+m2) l1²        m2 l1 l2 cosδ ]
//       [ m2 l1 l2 cosδ      m2 l2²        ]
//   f1 = -m2 l1 l2 ω2² sinδ - (m1+m2) g l1 sinθ1 - c ω1
//   f2 = +m2 l1 l2 ω1² sinδ -  m2 g l2 sinθ2     - c ω2
//   (δ = θ1-θ2 ; les termes en ω² sont les forces centrifuge et de Coriolis)
// ---------------------------------------------------------------------------
inline State deriv(const Params& P, const State& s) {
    const double d = s[0] - s[1], sd = std::sin(d), cd = std::cos(d);
    const double a = (P.m1 + P.m2) * P.l1 * P.l1, b = P.m2 * P.l1 * P.l2 * cd, e = P.m2 * P.l2 * P.l2;
    const double f1 = -P.m2 * P.l1 * P.l2 * s[3] * s[3] * sd - (P.m1 + P.m2) * P.g * P.l1 * std::sin(s[0]) - P.c * s[2];
    const double f2 =  P.m2 * P.l1 * P.l2 * s[2] * s[2] * sd - P.m2 * P.g * P.l2 * std::sin(s[1])        - P.c * s[3];
    const double det = a * e - b * b;  // ≥ m1 m2 l1² l2² > 0 : jamais singulier
    return {s[2], s[3], (e * f1 - b * f2) / det, (a * f2 - b * f1) / det};
}

// ---------------------------------------------------------------------------
// Intégrateurs
// ---------------------------------------------------------------------------
inline State axpy(const State& y, double h, const State& k) {
    return {y[0] + h * k[0], y[1] + h * k[1], y[2] + h * k[2], y[3] + h * k[3]};
}

// Runge-Kutta d'ordre 4 (erreur globale O(h⁴))
inline State rk4(const Params& P, const State& y, double h) {
    const State k1 = deriv(P, y);
    const State k2 = deriv(P, axpy(y, h / 2, k1));
    const State k3 = deriv(P, axpy(y, h / 2, k2));
    const State k4 = deriv(P, axpy(y, h, k3));
    State r;
    for (int i = 0; i < 4; ++i) r[i] = y[i] + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    return r;
}

// Euler explicite (gardé pour comparaison : l'énergie diverge)
inline State euler(const Params& P, const State& y, double h) { return axpy(y, h, deriv(P, y)); }

// Point-milieu implicite en variables canoniques (θ, p) : schéma SYMPLECTIQUE.
// Équations de Hamilton :  θ̇ = ∂H/∂p = ω ,  ṗ = -∂H/∂θ = ∂L/∂θ
//   ṗ1 = -m2 l1 l2 ω1 ω2 sinδ - (m1+m2) g l1 sinθ1
//   ṗ2 = +m2 l1 l2 ω1 ω2 sinδ -  m2 g l2 sinθ2
// Le flux hamiltonien conserve le volume de l'espace des phases (Liouville) ;
// le schéma garde l'énergie bornée sur des temps très longs (pas de dérive).
inline std::array<double, 4> hamiltonRhs(const Params& P, double t1, double t2, double p1, double p2) {
    const State s = fromMomenta(P, t1, t2, p1, p2);
    const double d = t1 - t2, sd = std::sin(d);
    return {s[2], s[3],
            -P.m2 * P.l1 * P.l2 * s[2] * s[3] * sd - (P.m1 + P.m2) * P.g * P.l1 * std::sin(t1),
             P.m2 * P.l1 * P.l2 * s[2] * s[3] * sd - P.m2 * P.g * P.l2 * std::sin(t2)};
}
inline State symplecticMidpoint(const Params& P, const State& y, double h) {
    const auto p0 = momenta(P, y);
    std::array<double, 4> z0{y[0], y[1], p0[0], p0[1]}, z1 = z0;
    for (int it = 0; it < 30; ++it) {                    // point fixe sur z1 = z0 + h f((z0+z1)/2)
        std::array<double, 4> m, f;
        for (int i = 0; i < 4; ++i) m[i] = 0.5 * (z0[i] + z1[i]);
        f = hamiltonRhs(P, m[0], m[1], m[2], m[3]);
        double dz = 0;
        for (int i = 0; i < 4; ++i) { const double n = z0[i] + h * f[i]; dz = std::fmax(dz, std::fabs(n - z1[i])); z1[i] = n; }
        if (dz < 1e-14) break;
    }
    return fromMomenta(P, z1[0], z1[1], z1[2], z1[3]);
}

// Dormand-Prince 5(4) : un pas avec estimation d'erreur (pas adaptatif)
inline State dopri5(const Params& P, const State& y, double h, double& err) {
    static const double c2 = 1. / 5, c3 = 3. / 10, c4 = 4. / 5, c5 = 8. / 9;
    static const double a21 = 1. / 5;
    static const double a31 = 3. / 40, a32 = 9. / 40;
    static const double a41 = 44. / 45, a42 = -56. / 15, a43 = 32. / 9;
    static const double a51 = 19372. / 6561, a52 = -25360. / 2187, a53 = 64448. / 6561, a54 = -212. / 729;
    static const double a61 = 9017. / 3168, a62 = -355. / 33, a63 = 46732. / 5247, a64 = 49. / 176, a65 = -5103. / 18656;
    static const double b1 = 35. / 384, b3 = 500. / 1113, b4 = 125. / 192, b5 = -2187. / 6784, b6 = 11. / 84;
    static const double e1 = 71. / 57600, e3 = -71. / 16695, e4 = 71. / 1920, e5 = -17253. / 339200, e6 = 22. / 525, e7 = -1. / 40;
    (void)c2; (void)c3; (void)c4; (void)c5;
    auto comb = [&](double c1, const State& k1, double c2_, const State& k2, double c3_, const State& k3,
                    double c4_, const State& k4, double c5_, const State& k5) {
        State r;
        for (int i = 0; i < 4; ++i)
            r[i] = y[i] + h * (c1 * k1[i] + c2_ * k2[i] + c3_ * k3[i] + c4_ * k4[i] + c5_ * k5[i]);
        return r;
    };
    const State z{};
    const State k1 = deriv(P, y);
    const State k2 = deriv(P, comb(a21, k1, 0, z, 0, z, 0, z, 0, z));
    const State k3 = deriv(P, comb(a31, k1, a32, k2, 0, z, 0, z, 0, z));
    const State k4 = deriv(P, comb(a41, k1, a42, k2, a43, k3, 0, z, 0, z));
    const State k5 = deriv(P, comb(a51, k1, a52, k2, a53, k3, a54, k4, 0, z));
    const State k6 = deriv(P, comb(a61, k1, a62, k2, a63, k3, a64, k4, a65, k5));
    State y5;
    for (int i = 0; i < 4; ++i) y5[i] = y[i] + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
    const State k7 = deriv(P, y5);
    err = 0;
    for (int i = 0; i < 4; ++i)
        err = std::fmax(err, std::fabs(h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i])));
    return y5;
}

// Avance de dt total avec pas adaptatif (tolérance tol) ; h est conservé entre appels.
inline State advanceAdaptive(const Params& P, State y, double dt, double& h, double tol = 1e-10) {
    double t = 0;
    if (h <= 0) h = dt / 10;
    while (t < dt) {
        double step = std::fmin(h, dt - t), err;
        State yn = dopri5(P, y, step, err);
        if (err <= tol) { t += step; y = yn; }
        const double fac = std::fmin(5.0, std::fmax(0.2, 0.9 * std::pow(tol / (err + 1e-300), 0.2)));
        h = step * fac;
    }
    return y;
}

// Normalisation des angles dans ]-π, π] (n'affecte pas la physique : tout est 2π-périodique)
inline double wrap(double a) {
    a = std::fmod(a + PI, 2 * PI);
    if (a <= 0) a += 2 * PI;
    return a - PI;
}

// ---------------------------------------------------------------------------
// Positions cartésiennes (origine au pivot, y vers le HAUT)
// ---------------------------------------------------------------------------
struct Cart { double x1, y1, x2, y2; };
inline Cart cartesian(const Params& P, const State& s) {
    const double x1 = P.l1 * std::sin(s[0]), y1 = -P.l1 * std::cos(s[0]);
    return {x1, y1, x1 + P.l2 * std::sin(s[1]), y1 - P.l2 * std::cos(s[1])};
}
// Vitesses cartésiennes
inline Cart velocity(const Params& P, const State& s) {
    const double vx1 = P.l1 * s[2] * std::cos(s[0]), vy1 = P.l1 * s[2] * std::sin(s[0]);
    return {vx1, vy1, vx1 + P.l2 * s[3] * std::cos(s[1]), vy1 + P.l2 * s[3] * std::sin(s[1])};
}
// Tensions des fils (via Newton) — T1, T2 ≥ 0 tant que le fil reste tendu
struct Tension { double T1, T2; };
inline Tension tensions(const Params& P, const State& s) {
    const State a = deriv(P, s);
    // Accélération de m2 : a2 = a1 + l2(θ̈2 e_t2 - ω2² e_r2). T2 = m2 (g cosθ2 + l1·... ) -> on projette sur le fil 2
    const double d = s[0] - s[1];
    const double T2 = P.m2 * (P.g * std::cos(s[1]) + P.l2 * s[3] * s[3]
                              + P.l1 * (a[2] * std::sin(d) + s[2] * s[2] * std::cos(d)));
    // Équilibre de m1 le long du fil 1 : T1 = T2 cosδ + m1(g cosθ1 + l1 ω1²) + ...
    const double T1 = T2 * std::cos(d) + P.m1 * (P.g * std::cos(s[0]) + P.l1 * s[2] * s[2]);
    return {T1, T2};
}

// ---------------------------------------------------------------------------
// Théorie linéaire (petits angles) : pulsations propres, solution analytique.
//   (K - ω² M0) u = 0 avec M0 = M(0), K = diag((m1+m2) g l1, m2 g l2)
// ---------------------------------------------------------------------------
struct Modes { double w2_slow, w2_fast; };
inline Modes smallAngleModes(const Params& P) {
    const double a = (P.m1 + P.m2) * P.l1 * P.l1, b = P.m2 * P.l1 * P.l2, e = P.m2 * P.l2 * P.l2;
    const double k1 = (P.m1 + P.m2) * P.g * P.l1, k2 = P.m2 * P.g * P.l2;
    // det(K - λM0) = 0  ->  (a e - b²) λ² - (a k2 + e k1) λ + k1 k2 = 0
    const double A = a * e - b * b, B = -(a * k2 + e * k1), C = k1 * k2;
    const double D = std::sqrt(B * B - 4 * A * C);
    return {(-B - D) / (2 * A), (-B + D) / (2 * A)};
}

// ---------------------------------------------------------------------------
// Exposant de Lyapounov maximal (Benettin) : on suit une trajectoire voisine à
// distance d0, on mesure la croissance de l'écart sur un pas, puis on renormalise.
// ---------------------------------------------------------------------------
struct Lyapunov {
    State ref, twin;
    double sum = 0, T = 0, d0 = 1e-8;
    void init(const State& s) {
        ref = s; twin = s; twin[0] += d0; sum = 0; T = 0;
    }
    static double dist(const State& a, const State& b) {
        double q = 0; for (int i = 0; i < 4; ++i) q += (a[i] - b[i]) * (a[i] - b[i]);
        return std::sqrt(q);
    }
    void step(const Params& P, double h) {
        ref = rk4(P, ref, h); twin = rk4(P, twin, h);
        const double d = dist(ref, twin);
        sum += std::log(d / d0); T += h;
        const double r = d0 / d;                       // renormalisation
        for (int i = 0; i < 4; ++i) twin[i] = ref[i] + (twin[i] - ref[i]) * r;
    }
    double value() const { return T > 0 ? sum / T : 0; }
};

}  // namespace dp
