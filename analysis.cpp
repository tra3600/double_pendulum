// analysis.cpp — génère les figures de physique (SVG, aucune dépendance)
//   g++ -O2 -std=c++17 analysis.cpp -o analysis && ./analysis figures/
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <vector>
#include "pendulum.hpp"
using namespace dp;

// ---------- mini bibliothèque SVG ----------
struct Svg {
    std::ostringstream o; int W, H;
    Svg(int w, int h) : W(w), H(h) {
        o << "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 " << w << ' ' << h << "' width='" << w << "' height='" << h
          << "' font-family='Helvetica,Arial,sans-serif'>\n<rect width='100%' height='100%' fill='#0f1420'/>\n";
    }
    void text(double x, double y, const std::string& s, int sz = 13, const char* col = "#cfd6e6", const char* anchor = "start") {
        o << "<text x='" << x << "' y='" << y << "' font-size='" << sz << "' fill='" << col << "' text-anchor='" << anchor << "'>" << s << "</text>\n";
    }
    void line(double x1, double y1, double x2, double y2, const char* col, double w = 1, const char* extra = "") {
        o << "<line x1='" << x1 << "' y1='" << y1 << "' x2='" << x2 << "' y2='" << y2 << "' stroke='" << col << "' stroke-width='" << w << "' " << extra << "/>\n";
    }
    void circle(double x, double y, double r, const char* fill) { o << "<circle cx='" << x << "' cy='" << y << "' r='" << r << "' fill='" << fill << "'/>\n"; }
    void polyline(const std::vector<std::pair<double, double>>& p, const char* col, double w = 1.5, double op = 1) {
        o << "<polyline fill='none' stroke='" << col << "' stroke-width='" << w << "' stroke-opacity='" << op << "' stroke-linejoin='round' points='";
        for (auto& q : p) o << q.first << ',' << q.second << ' ';
        o << "'/>\n";
    }
    void save(const std::string& f) { o << "</svg>\n"; std::ofstream(f) << o.str(); std::printf("  écrit %s\n", f.c_str()); }
};
static std::string rgb(double r, double g, double b) {
    char buf[32]; std::snprintf(buf, 32, "#%02x%02x%02x", int(255 * std::clamp(r, 0., 1.)), int(255 * std::clamp(g, 0., 1.)), int(255 * std::clamp(b, 0., 1.))); return buf;
}
static std::string viridis(double t) {  // approx. de la palette viridis
    t = std::clamp(t, 0., 1.);
    static const double c[5][3] = {{0.267, 0.005, 0.329}, {0.231, 0.322, 0.546}, {0.128, 0.567, 0.551}, {0.369, 0.789, 0.383}, {0.993, 0.906, 0.144}};
    double x = t * 4; int i = std::min(3, int(x)); double f = x - i;
    return rgb(c[i][0] + f * (c[i + 1][0] - c[i][0]), c[i][1] + f * (c[i + 1][1] - c[i][1]), c[i][2] + f * (c[i + 1][2] - c[i][2]));
}
struct Axes {  // repère de tracé
    double x0, y0, w, h, xmin, xmax, ymin, ymax;
    double X(double x) const { return x0 + (x - xmin) / (xmax - xmin) * w; }
    double Y(double y) const { return y0 + h - (y - ymin) / (ymax - ymin) * h; }
};
static void frame(Svg& s, const Axes& a, const std::string& title, const std::string& xl, const std::string& yl, int nx = 5, int ny = 5, bool logy = false) {
    s.text(a.x0, a.y0 - 10, title, 15, "#ffffff");
    s.line(a.x0, a.y0 + a.h, a.x0 + a.w, a.y0 + a.h, "#6b7694");
    s.line(a.x0, a.y0, a.x0, a.y0 + a.h, "#6b7694");
    char b[32];
    for (int i = 0; i <= nx; ++i) {
        double v = a.xmin + (a.xmax - a.xmin) * i / nx; std::snprintf(b, 32, "%.3g", v);
        s.line(a.X(v), a.y0, a.X(v), a.y0 + a.h, "#252d42", 1); s.text(a.X(v), a.y0 + a.h + 16, b, 11, "#9aa4bf", "middle");
    }
    for (int i = 0; i <= ny; ++i) {
        double v = a.ymin + (a.ymax - a.ymin) * i / ny; std::snprintf(b, 32, logy ? "1e%.0f" : "%.3g", v);
        s.line(a.x0, a.Y(v), a.x0 + a.w, a.Y(v), "#252d42", 1); s.text(a.x0 - 6, a.Y(v) + 4, b, 11, "#9aa4bf", "end");
    }
    s.text(a.x0 + a.w / 2, a.y0 + a.h + 34, xl, 12, "#cfd6e6", "middle");
    s.o << "<text transform='translate(" << a.x0 - 44 << ',' << a.y0 + a.h / 2 << ") rotate(-90)' font-size='12' fill='#cfd6e6' text-anchor='middle'>" << yl << "</text>\n";
}

// ---------- 1. trajectoire de la masse 2 colorée par le temps + instantané ----------
static void figTrajectory(const std::string& dir, const Params& P) {
    State s{PI * 0.75, PI * 0.5, 0, 0};
    const int n = 6000; const double h = 0.002;   // 12 s
    std::vector<Cart> c; c.reserve(n);
    for (int i = 0; i < n; ++i) { c.push_back(cartesian(P, s)); s = rk4(P, s, h); }
    Svg g(640, 640); const double R = P.l1 + P.l2, sc = 280 / R, cx = 320, cy = 300;
    g.text(24, 30, "Trajectoire de la masse 2 (Lagrange, RK4)", 16, "#ffffff");
    g.text(24, 50, "θ1=135°, θ2=90°, t∈[0,12 s] — couleur = temps", 12, "#9aa4bf");
    for (int i = 1; i < n; ++i) {
        std::string col = viridis(double(i) / n);
        g.o << "<line x1='" << cx + sc * c[i - 1].x2 << "' y1='" << cy - sc * c[i - 1].y2 << "' x2='" << cx + sc * c[i].x2 << "' y2='" << cy - sc * c[i].y2
            << "' stroke='" << col << "' stroke-width='1.6'/>\n";
    }
    const Cart& k = c[n * 2 / 5];  // instantané
    g.line(cx, cy, cx + sc * k.x1, cy - sc * k.y1, "#e8ecf7", 3); g.line(cx + sc * k.x1, cy - sc * k.y1, cx + sc * k.x2, cy - sc * k.y2, "#e8ecf7", 3);
    g.circle(cx, cy, 5, "#8892ad"); g.circle(cx + sc * k.x1, cy - sc * k.y1, 9, "#5aa9ff"); g.circle(cx + sc * k.x2, cy - sc * k.y2, 9, "#ff6b6b");
    g.save(dir + "/trajectoire.svg");
}

// ---------- 2. dérive d'énergie : Euler / RK4 / symplectique / DP5 ----------
static void figEnergy(const std::string& dir, const Params& P) {
    State s0{2.0, 1.0, 0, 0};
    struct M { const char* nom; const char* col; std::function<State(State)> step; std::vector<std::pair<double, double>> pts; };
    const double h = 5e-3, T = 60; const int n = int(T / h);
    double ha = 0;
    std::vector<M> ms = {
        {"Euler explicite (ordre 1)", "#ff6b6b", [&](State s) { return euler(P, s, h); }, {}},
        {"Point-milieu symplectique (ordre 2)", "#f5c542", [&](State s) { return symplecticMidpoint(P, s, h); }, {}},
        {"RK4 (ordre 4)", "#5aa9ff", [&](State s) { return rk4(P, s, h); }, {}},
        {"Dormand-Prince 5(4) adaptatif", "#4fd18b", [&](State s) { return advanceAdaptive(P, s, h, ha, 1e-12); }, {}}};
    const double E0 = energy(P, s0);
    Axes a{80, 70, 640, 340, 0, T, -16, 1};
    Svg g(780, 500); frame(g, a, "Dérive relative de l'énergie |E(t)-E0|/|E0| (h = 5 ms)", "t (s)", "log10", 6, 6, true);
    for (auto& m : ms) {
        State s = s0;
        for (int i = 1; i <= n; ++i) {
            s = m.step(s);
            if (i % 20 == 0) { double d = std::fabs(energy(P, s) - E0) / std::fabs(E0) + 1e-17; m.pts.push_back({a.X(i * h), a.Y(std::clamp(std::log10(d), -16., 1.))}); }
        }
        g.polyline(m.pts, m.col, 1.8);
    }
    int yy = 440; for (auto& m : ms) { g.line(100, yy - 4, 125, yy - 4, m.col, 3); g.text(132, yy, m.nom, 12); yy += 18; if (yy > 490) break; }
    g.save(dir + "/derive_energie.svg");
}

// ---------- 3. effet papillon + exposant de Lyapounov ----------
static void figButterfly(const std::string& dir, const Params& P) {
    Svg g(780, 520);
    Axes a{80, 70, 640, 330, 0, 30, -9, 3};
    frame(g, a, "Effet papillon : écart entre deux pendules (δθ1 = 1e-8 rad)", "t (s)", "log10 ‖δ‖", 6, 6, true);
    struct Cas { const char* nom; State s; const char* col; };
    Cas cas[] = {{"petits angles (θ=0.05)", {0.05, 0.05, 0, 0}, "#4fd18b"}, {"chaotique (θ1=2.5, θ2=2.0)", {2.5, 2.0, 0, 0}, "#ff6b6b"}};
    int yy = 450; double lam = 0;
    for (auto& c : cas) {
        State A = c.s, B = c.s; B[0] += 1e-8;
        Lyapunov L; L.init(c.s);
        std::vector<std::pair<double, double>> pts; const double h = 1e-3;
        for (int i = 1; i <= 30000; ++i) {
            A = rk4(P, A, h); B = rk4(P, B, h); L.step(P, h);
            if (i % 100 == 0) { double d = Lyapunov::dist(A, B); pts.push_back({a.X(i * h), a.Y(std::clamp(std::log10(d + 1e-30), -9., 3.))}); }
        }
        g.polyline(pts, c.col, 1.8);
        char b[128]; std::snprintf(b, 128, "%s — λ ≈ %.2f s⁻¹  (horizon de prédiction ~ 1/λ ≈ %.2f s)", c.nom, L.value(), 1 / std::max(L.value(), 1e-9));
        g.line(100, yy - 4, 125, yy - 4, c.col, 3); g.text(132, yy, b, 12); yy += 20; lam = L.value();
    }
    // pente théorique λ pour le cas chaotique
    g.line(a.X(0), a.Y(-8), a.X(9 / std::max(lam, .01)), a.Y(-8 + 9 / std::max(lam, .01) * lam / std::log(10.0) * 1.0), "#ffffff", 1, "stroke-dasharray='5,4'");
    g.save(dir + "/effet_papillon.svg");
}

// ---------- 4. portraits de phase (θ1, ω1) régulier vs chaotique ----------
static void figPhase(const std::string& dir, const Params& P) {
    Svg g(900, 440);
    struct Cas { const char* nom; State s; const char* col; double th, om; };
    Cas cas[] = {{"E faible : mouvement quasi-périodique (tore KAM)", {0.5, 0.3, 0, 0}, "#4fd18b", 1.0, 4},
                 {"E élevée : mouvement chaotique", {2.5, 2.0, 0, 0}, "#ff6b6b", PI, 14}};
    for (int k = 0; k < 2; ++k) {
        Axes a{70.0 + k * 440, 60, 340, 300, -cas[k].th, cas[k].th, -cas[k].om, cas[k].om};
        frame(g, a, cas[k].nom, "θ1 (rad)", "ω1 (rad/s)", 4, 4);
        State s = cas[k].s; std::vector<std::pair<double, double>> pts;
        for (int i = 0; i < 60000; ++i) { s = rk4(P, s, 2e-3); if (i % 3 == 0) { double t = cas[k].th < 2 ? s[0] : wrap(s[0]); pts.push_back({a.X(t), a.Y(s[2])}); } }
        for (auto& p : pts) g.o << "<circle cx='" << p.first << "' cy='" << p.second << "' r='0.9' fill='" << cas[k].col << "' fill-opacity='0.55'/>\n";
    }
    g.save(dir + "/portrait_phase.svg");
}

// ---------- 5. section de Poincaré : θ1 = 0, ω1 > 0 -> (θ2, p2) à énergie fixée ----------
static void figPoincare(const std::string& dir, const Params& P) {
    const double E = -15.0;          // énergie fixée (J)
    Svg g(640, 560); Axes a{80, 60, 500, 440, -PI, PI, -7, 7};
    frame(g, a, "Section de Poincaré (θ1=0, ω1>0), E = -15 J (îlots KAM + mer chaotique)", "θ2 (rad)", "p2 (kg·m²/s)", 6, 6);
    int col = 0;
    for (double th2 : {0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, -0.4, -0.8, -1.2, -1.6, -2.0, -2.4}) {
        // conditions : θ1=0, ω2 donné par th2/p2 ; on résout ω1 > 0 pour atteindre E
        double w2 = 0.0; State s{0, th2, 0, w2};
        double V = potential(P, s);
        // énergie cinétique requise : E - V ; avec θ1=0 : T = ½(m1+m2)l1²ω1² + ½ m2 l2² ω2² + m2 l1 l2 ω1 ω2 cosδ ; on fixe ω2 = 0.3 th2 puis on résout ω1
        w2 = 0.3 * th2; s[3] = w2;
        double A = 0.5 * (P.m1 + P.m2) * P.l1 * P.l1, B = P.m2 * P.l1 * P.l2 * w2 * std::cos(-th2), C = 0.5 * P.m2 * P.l2 * P.l2 * w2 * w2 + V - E;
        double disc = B * B - 4 * A * C; if (disc < 0) continue;
        s[2] = (-B + std::sqrt(disc)) / (2 * A); if (s[2] <= 0) continue;
        std::string colr = viridis((col++ % 9) / 8.0);
        double prev = s[0];
        for (int i = 0; i <= 400000; ++i) {
            State n = rk4(P, s, 1e-3);
            if (prev < 0 && wrap(n[0]) >= 0 && n[2] > 0 && std::fabs(wrap(n[0])) < 0.5) {  // croisement θ1 : - -> +
                double f = -prev / (n[0] - prev); State q = s; for (int j = 0; j < 4; ++j) q[j] = s[j] + f * (n[j] - s[j]);
                auto p = momenta(P, q);
                g.circle(a.X(wrap(q[1])), a.Y(p[1]), 1.3, colr.c_str());
            }
            prev = wrap(n[0]); s = n; s[0] = wrap(s[0]);
        }
    }
    g.save(dir + "/poincare.svg");
}

// ---------- 6. carte du chaos : temps avant le premier retournement de la masse 2 ----------
static void figChaosMap(const std::string& dir, const Params& P) {
    const int N = 180; const double T = 10, h = 4e-3; const int steps = int(T / h);
    Svg g(640, 600); Axes a{80, 60, 480, 480, -PI, PI, -PI, PI};
    frame(g, a, "Carte du chaos : temps avant le 1er retournement de la masse 2", "θ1(0)", "θ2(0)", 4, 4);
    const double cs = a.w / N;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            State s{-PI + (i + 0.5) * 2 * PI / N, -PI + (j + 0.5) * 2 * PI / N, 0, 0};
            // Critère de retournement : |θ2 déroulé| > π
            double tflip = T; double th2 = s[1], prev = s[1];
            double e = energy(P, s);
            // Pour un état dont l'énergie < énergie minimale de retournement, jamais de retournement : échappe au calcul
            (void)e;
            for (int k = 0; k < steps; ++k) {
                s = rk4(P, s, h);
                th2 += s[1] - prev; prev = s[1];
                if (std::fabs(th2) > PI) { tflip = (k + 1) * h; break; }
            }
            g.o << "<rect x='" << a.x0 + i * cs << "' y='" << a.y0 + a.h - (j + 1) * cs << "' width='" << cs + 0.5 << "' height='" << cs + 0.5 << "' fill='" << viridis(std::sqrt(tflip / T)) << "'/>\n";
        }
    g.text(a.x0, a.y0 + a.h + 60, "jaune = jamais retourné en 10 s (régulier) ; violet = retournement immédiat ; mélange fin = chaos", 11, "#9aa4bf");
    g.save(dir + "/carte_chaos.svg");
}

// ---------- 7. modes propres aux petits angles ----------
static void figModes(const std::string& dir, const Params& P) {
    Svg g(780, 410);
    Modes m = smallAngleModes(P);
    Axes a{80, 60, 640, 250, 0, 6, -1.3, 1.3};
    frame(g, a, "Petits angles : superposition de 2 modes propres (analytique) vs simulation", "t (s)", "θ1 / amplitude", 6, 4);
    const double A = 0.01; State s{A, 0, 0, 0};
    const double w1 = std::sqrt(m.w2_slow), w2 = std::sqrt(m.w2_fast);
    // vecteurs propres pour m1=m2,l1=l2 : lent (1,√2), rapide (1,-√2) ; décomposition de (A,0) : A/2 sur chaque
    std::vector<std::pair<double, double>> ana, sim;
    const double h = 1e-3;
    for (int i = 0; i <= 6000; ++i) {
        double t = i * h;
        if (i % 10 == 0) {
            ana.push_back({a.X(t), a.Y(A / 2 * (std::cos(w1 * t) + std::cos(w2 * t)) / A)});
            sim.push_back({a.X(t), a.Y(s[0] / A)});
        }
        s = rk4(P, s, h);
    }
    g.polyline(ana, "#ffffff", 3.5, 0.45); g.polyline(sim, "#5aa9ff", 1.4);
    char b[160]; std::snprintf(b, 160, "ω_lent=%.3f rad/s, ω_rapide=%.3f rad/s  (formules : ω² = (2∓√2) g/l)", w1, w2);
    g.text(80, 374, b, 12); g.text(80, 392, "blanc épais = théorie linéaire ; bleu = intégration RK4 des équations non linéaires complètes", 12, "#9aa4bf");
    g.save(dir + "/modes_propres.svg");
}

int main(int argc, char** argv) {
    std::string dir = argc > 1 ? argv[1] : "figures";
    std::string cmd = "mkdir -p " + dir; if (std::system(cmd.c_str())) return 1;
    Params P;  // m1=m2=1 kg, l1=l2=1 m, g=9.81
    std::printf("Génération des figures dans %s/\n", dir.c_str());
    figTrajectory(dir, P);
    figModes(dir, P);
    figEnergy(dir, P);
    figButterfly(dir, P);
    figPhase(dir, P);
    figPoincare(dir, P);
    figChaosMap(dir, P);
    return 0;
}
