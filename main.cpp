// Pendule double — simulation interactive (mécanique de Lagrange + SDL2)
//
//   make && ./double_pendulum [l1 l2 theta1_deg theta2_deg [m1 m2 [frottement]]]
//
//   Entrées : longueurs (m), angles initiaux par rapport à la verticale (degrés), masses (kg), frottement c.
//   Touches : ESPACE pause · R reset · T jumeau perturbé (effet papillon) · D frottement on/off
//             ↑/↓ vitesse du temps · C efface les traînées · S capture (capture.bmp) · ÉCHAP quitter
//   Souris  : cliquer-glisser une masse pour fixer les conditions initiales (met en pause).
#include <SDL2/SDL.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <string>
#include "pendulum.hpp"
using namespace dp;

namespace {

constexpr int W = 1100, H = 700;
constexpr double STEP = 1.0 / 600.0;          // pas d'intégration RK4 (s)
constexpr double TWIN_EPS = 1e-6;             // perturbation du jumeau (rad)

// ---------------- primitives de dessin ----------------
void fillCircle(SDL_Renderer* r, int cx, int cy, int rad) {
    for (int dy = -rad; dy <= rad; ++dy) {
        const int dx = int(std::sqrt(double(rad * rad - dy * dy)));
        SDL_RenderDrawLine(r, cx - dx, cy + dy, cx + dx, cy + dy);
    }
}
void thickLine(SDL_Renderer* r, double x1, double y1, double x2, double y2, int w) {
    const double dx = x2 - x1, dy = y2 - y1, L = std::hypot(dx, dy);
    if (L < 1e-9) return;
    const double nx = -dy / L, ny = dx / L;
    for (int i = -w / 2; i <= w / 2; ++i)
        SDL_RenderDrawLineF(r, float(x1 + nx * i), float(y1 + ny * i), float(x2 + nx * i), float(y2 + ny * i));
}
void setCol(SDL_Renderer* r, int R, int G, int B, int A = 255) { SDL_SetRenderDrawColor(r, Uint8(R), Uint8(G), Uint8(B), Uint8(A)); }

// Police 3x5 pour les chiffres et quelques symboles (HUD numérique sans SDL_ttf)
const char* GLYPH[] = {
    "111101101101111", "010110010010111", "111001111100111", "111001111001111", "101101111001001",
    "111100111001111", "111100111101111", "111001001001001", "111101111101111", "111101111001111",  // 0-9
    "000000000000010", "000000111000000", "010111010111010",                                        // . - +
};
void drawGlyph(SDL_Renderer* r, int x, int y, int g, int s) {
    for (int i = 0; i < 15; ++i)
        if (GLYPH[g][i] == '1') { SDL_Rect q{x + (i % 3) * s, y + (i / 3) * s, s, s}; SDL_RenderFillRect(r, &q); }
}
void drawNumber(SDL_Renderer* r, int x, int y, double v, int dec, int s) {
    char b[32]; std::snprintf(b, 32, "%.*f", dec, v);
    for (const char* p = b; *p; ++p) {
        int g = *p == '.' ? 10 : *p == '-' ? 11 : *p == '+' ? 12 : (*p >= '0' && *p <= '9') ? *p - '0' : -1;
        if (g >= 0) drawGlyph(r, x, y, g, s);
        x += 4 * s;
    }
}

struct Trail {
    std::deque<std::pair<float, float>> pts;
    size_t cap = 2500;
    void add(float x, float y) { pts.emplace_back(x, y); if (pts.size() > cap) pts.pop_front(); }
    void draw(SDL_Renderer* r, int R, int G, int B) const {
        const size_t n = pts.size();
        for (size_t i = 1; i < n; ++i) {
            const double f = double(i) / n;                      // 0 = ancien, 1 = récent
            setCol(r, R, G, B, int(10 + 235 * f * f));
            SDL_RenderDrawLineF(r, pts[i - 1].first, pts[i - 1].second, pts[i].first, pts[i].second);
        }
    }
};

struct Sim {
    Params P;
    State s0{}, s{}, tw{};
    double t = 0, E0 = 0, speed = 1;
    bool paused = false, twin = true, damped = false;
    double cdamp = 0.15;
    Trail trail, trailTwin;
    std::deque<std::array<float, 3>> eh;                          // historique (T, V, E)
    std::deque<std::pair<float, float>> ph;                       // portrait de phase (θ1, ω1)
    Lyapunov ly;

    void reset() {
        s = s0; tw = s0; tw[0] += TWIN_EPS; t = 0; E0 = energy(P, s);
        trail.pts.clear(); trailTwin.pts.clear(); eh.clear(); ph.clear(); ly.init(s0);
    }
    void step() {
        Params Q = P; Q.c = damped ? cdamp : 0.0;
        s = rk4(Q, s, STEP); tw = rk4(Q, tw, STEP);
        if (!damped) ly.step(Q, STEP);
        t += STEP;
    }
};

}  // namespace

int main(int argc, char** argv) {
    Sim sim;
    double th1 = 120, th2 = -10;
    int maxFrames = -1; const char* shot = nullptr;
    std::string pos[8]; int np = 0;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--frames") && i + 1 < argc) maxFrames = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--screenshot") && i + 1 < argc) shot = argv[++i];
        else if (np < 8) pos[np++] = argv[i];
    }
    if (np >= 2) { sim.P.l1 = std::atof(pos[0].c_str()); sim.P.l2 = std::atof(pos[1].c_str()); }
    if (np >= 4) { th1 = std::atof(pos[2].c_str()); th2 = std::atof(pos[3].c_str()); }
    if (np >= 6) { sim.P.m1 = std::atof(pos[4].c_str()); sim.P.m2 = std::atof(pos[5].c_str()); }
    if (np >= 7) { sim.cdamp = std::atof(pos[6].c_str()); sim.damped = sim.cdamp > 0; }
    if (sim.P.l1 <= 0 || sim.P.l2 <= 0 || sim.P.m1 <= 0 || sim.P.m2 <= 0) { std::fprintf(stderr, "longueurs et masses doivent être > 0\n"); return 1; }
    sim.s0 = {th1 * PI / 180, th2 * PI / 180, 0, 0};
    sim.reset();

    if (SDL_Init(SDL_INIT_VIDEO) < 0) { std::fprintf(stderr, "SDL_Init: %s\n", SDL_GetError()); return 1; }
    SDL_Window* win = SDL_CreateWindow("Pendule double", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, W, H, SDL_WINDOW_SHOWN);
    if (!win) { std::fprintf(stderr, "SDL_CreateWindow: %s\n", SDL_GetError()); SDL_Quit(); return 1; }
    SDL_Renderer* ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!ren) ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_SOFTWARE);
    if (!ren) { std::fprintf(stderr, "SDL_CreateRenderer: %s\n", SDL_GetError()); SDL_DestroyWindow(win); SDL_Quit(); return 1; }
    SDL_SetRenderDrawBlendMode(ren, SDL_BLENDMODE_BLEND);

    const double ox = 380, oy = 270;                              // pivot (px)
    const double scale = 0.46 * std::min(W * 0.62, double(H)) / (sim.P.l1 + sim.P.l2);   // px / m
    auto toPx = [&](double x, double y) { return std::make_pair(ox + scale * x, oy - scale * y); };

    int drag = 0;                                                 // 1 = masse 1, 2 = masse 2
    bool quit = false; int frame = 0;
    Uint64 prev = SDL_GetPerformanceCounter(); double acc = 0, titleT = 0;

    while (!quit) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) quit = true;
            else if (e.type == SDL_KEYDOWN) {
                switch (e.key.keysym.sym) {
                    case SDLK_ESCAPE: quit = true; break;
                    case SDLK_SPACE: sim.paused = !sim.paused; break;
                    case SDLK_r: sim.reset(); break;
                    case SDLK_t: sim.twin = !sim.twin; break;
                    case SDLK_d: sim.damped = !sim.damped; sim.reset(); break;
                    case SDLK_c: sim.trail.pts.clear(); sim.trailTwin.pts.clear(); break;
                    case SDLK_UP: sim.speed = std::min(8.0, sim.speed * 1.5); break;
                    case SDLK_DOWN: sim.speed = std::max(0.05, sim.speed / 1.5); break;
                    case SDLK_s: {
                        SDL_Surface* sf = SDL_CreateRGBSurfaceWithFormat(0, W, H, 32, SDL_PIXELFORMAT_ARGB8888);
                        if (sf) { SDL_RenderReadPixels(ren, nullptr, sf->format->format, sf->pixels, sf->pitch); SDL_SaveBMP(sf, "capture.bmp"); SDL_FreeSurface(sf); std::puts("capture.bmp écrite"); }
                        break; }
                }
            } else if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT) {
                const Cart c = cartesian(sim.P, sim.s);
                auto p1 = toPx(c.x1, c.y1), p2 = toPx(c.x2, c.y2);
                const double d1 = std::hypot(e.button.x - p1.first, e.button.y - p1.second), d2 = std::hypot(e.button.x - p2.first, e.button.y - p2.second);
                if (d2 < 25 && d2 <= d1) drag = 2; else if (d1 < 25) drag = 1;
                if (drag) { sim.paused = true; }
            } else if (e.type == SDL_MOUSEBUTTONUP) {
                if (drag) { sim.s0 = sim.s; sim.s0[2] = sim.s0[3] = 0; sim.reset(); }
                drag = 0;
            } else if (e.type == SDL_MOUSEMOTION && drag) {
                const double mx = (e.motion.x - ox) / scale, my = -(e.motion.y - oy) / scale;
                if (drag == 1) sim.s[0] = std::atan2(mx, -my);
                else { const Cart c = cartesian(sim.P, sim.s); sim.s[1] = std::atan2(mx - c.x1, -(my - c.y1)); }
                sim.s[2] = sim.s[3] = 0; sim.tw = sim.s; sim.tw[0] += TWIN_EPS;
            }
        }

        // ---- temps ----
        const Uint64 now = SDL_GetPerformanceCounter();
        double dt = double(now - prev) / double(SDL_GetPerformanceFrequency()); prev = now;
        dt = std::min(dt, 0.05);
        if (maxFrames >= 0) dt = 1.0 / 60;                        // mode déterministe (tests)
        if (!sim.paused) {
            acc += dt * sim.speed;
            while (acc >= STEP) {
                sim.step(); acc -= STEP;
                const Cart c = cartesian(sim.P, sim.s), ct = cartesian(sim.P, sim.tw);
                static int k = 0;
                if (++k % 8 == 0) {
                    auto p = toPx(c.x2, c.y2), q = toPx(ct.x2, ct.y2);
                    sim.trail.add(float(p.first), float(p.second)); sim.trailTwin.add(float(q.first), float(q.second));
                    sim.eh.push_back({float(kinetic(sim.P, sim.s)), float(potential(sim.P, sim.s)), float(energy(sim.P, sim.s))});
                    if (sim.eh.size() > 700) sim.eh.pop_front();
                    sim.ph.push_back({float(wrap(sim.s[0])), float(sim.s[2])});
                    if (sim.ph.size() > 3000) sim.ph.pop_front();
                }
            }
        }

        // ---- rendu ----
        setCol(ren, 15, 20, 32); SDL_RenderClear(ren);
        // panneau gauche : pendule
        setCol(ren, 24, 31, 48); SDL_Rect left{10, 10, 740, 680}; SDL_RenderFillRect(ren, &left);
        sim.trail.draw(ren, 90, 169, 255);
        if (sim.twin) sim.trailTwin.draw(ren, 255, 107, 107);
        auto drawPend = [&](const State& st, int R, int G, int B, int rodA, bool main) {
            const Cart c = cartesian(sim.P, st);
            auto a = toPx(c.x1, c.y1), b = toPx(c.x2, c.y2);
            setCol(ren, R, G, B, rodA);
            thickLine(ren, ox, oy, a.first, a.second, main ? 3 : 2); thickLine(ren, a.first, a.second, b.first, b.second, main ? 3 : 2);
            const int r1 = int(6 + 4 * std::cbrt(sim.P.m1)), r2 = int(6 + 4 * std::cbrt(sim.P.m2));   // rayon ~ m^(1/3)
            fillCircle(ren, int(a.first), int(a.second), main ? r1 : 5); fillCircle(ren, int(b.first), int(b.second), main ? r2 : 5);
        };
        if (sim.twin) drawPend(sim.tw, 255, 107, 107, 150, false);
        drawPend(sim.s, 232, 236, 247, 255, true);
        setCol(ren, 136, 146, 173); fillCircle(ren, int(ox), int(oy), 5);

        // panneau droit haut : énergies T, V, E
        {
            SDL_Rect r{770, 10, 320, 210}; setCol(ren, 24, 31, 48); SDL_RenderFillRect(ren, &r);
            double lo = -1, hi = 1;
            for (auto& v : sim.eh) for (int k = 0; k < 3; ++k) { lo = std::min<double>(lo, v[k]); hi = std::max<double>(hi, v[k]); }
            auto Y = [&](double v) { return r.y + r.h - 12 - (v - lo) / (hi - lo) * (r.h - 24); };
            setCol(ren, 60, 70, 100); SDL_RenderDrawLineF(ren, r.x, float(Y(0)), r.x + r.w, float(Y(0)));
            const int cols[3][3] = {{245, 197, 66}, {79, 209, 139}, {255, 255, 255}};   // T jaune, V vert, E blanc
            for (int k = 0; k < 3; ++k) {
                setCol(ren, cols[k][0], cols[k][1], cols[k][2]);
                for (size_t i = 1; i < sim.eh.size(); ++i)
                    SDL_RenderDrawLineF(ren, float(r.x + 8 + (i - 1) * 0.44), float(Y(sim.eh[i - 1][k])), float(r.x + 8 + i * 0.44), float(Y(sim.eh[i][k])));
            }
            // dérive relative : barre
            const double drift = sim.eh.empty() ? 0 : std::fabs(sim.eh.back()[2] - sim.E0) / std::max(1e-12, std::fabs(sim.E0));
            setCol(ren, 245, 197, 66); SDL_Rect l1{r.x + 8, r.y + 8, 10, 4}; SDL_RenderFillRect(ren, &l1);
            setCol(ren, 79, 209, 139); SDL_Rect l2{r.x + 28, r.y + 8, 10, 4}; SDL_RenderFillRect(ren, &l2);
            setCol(ren, 255, 255, 255); SDL_Rect l3{r.x + 48, r.y + 8, 10, 4}; SDL_RenderFillRect(ren, &l3);
            setCol(ren, 207, 214, 230); drawNumber(ren, r.x + 150, r.y + 6, drift, 8, 2);
        }
        // panneau droit milieu : portrait de phase (θ1, ω1)
        {
            SDL_Rect r{770, 230, 320, 230}; setCol(ren, 24, 31, 48); SDL_RenderFillRect(ren, &r);
            setCol(ren, 60, 70, 100);
            SDL_RenderDrawLineF(ren, r.x + r.w / 2.f, float(r.y), r.x + r.w / 2.f, float(r.y + r.h));
            SDL_RenderDrawLineF(ren, float(r.x), r.y + r.h / 2.f, float(r.x + r.w), r.y + r.h / 2.f);
            double wm = 6; for (auto& p : sim.ph) wm = std::max<double>(wm, std::fabs(p.second));
            for (size_t i = 0; i < sim.ph.size(); ++i) {
                const double f = double(i) / sim.ph.size();
                setCol(ren, 90, 169, 255, int(40 + 215 * f));
                SDL_RenderDrawPointF(ren, float(r.x + r.w / 2.0 + sim.ph[i].first / PI * (r.w / 2 - 6)), float(r.y + r.h / 2.0 - sim.ph[i].second / wm * (r.h / 2 - 6)));
            }
        }
        // panneau droit bas : écart jumeau (log) + Lyapounov
        {
            SDL_Rect r{770, 470, 320, 220}; setCol(ren, 24, 31, 48); SDL_RenderFillRect(ren, &r);
            const double d = Lyapunov::dist(sim.s, sim.tw);
            const double lg = std::clamp(std::log10(d + 1e-30), -7.0, 1.0);                    // -7..1
            const double fr = (lg + 7) / 8;
            setCol(ren, 255, 107, 107); SDL_Rect bar{r.x + 20, r.y + 70, int((r.w - 40) * fr), 22}; SDL_RenderFillRect(ren, &bar);
            setCol(ren, 60, 70, 100); SDL_Rect fb{r.x + 20, r.y + 70, r.w - 40, 22}; SDL_RenderDrawRect(ren, &fb);
            setCol(ren, 207, 214, 230); drawNumber(ren, r.x + 20, r.y + 20, sim.t, 1, 3); drawNumber(ren, r.x + 170, r.y + 20, sim.ly.value(), 2, 3);
            drawNumber(ren, r.x + 20, r.y + 110, lg, 1, 3);
            drawNumber(ren, r.x + 20, r.y + 150, sim.speed, 2, 3);
        }
        SDL_RenderPresent(ren);

        titleT += dt;
        if (titleT > 0.25) {
            titleT = 0;
            const Tension tn = tensions(sim.P, sim.s);
            char title[256];
            std::snprintf(title, 256, "Pendule double | t=%.1fs x%.2f | E=%.3f J (T=%.2f V=%.2f) | lyapunov=%.2f/s | ecart jumeau=%.1e | T1=%.1fN T2=%.1fN%s%s",
                          sim.t, sim.speed, energy(sim.P, sim.s), kinetic(sim.P, sim.s), potential(sim.P, sim.s), sim.ly.value(),
                          Lyapunov::dist(sim.s, sim.tw), tn.T1, tn.T2, sim.damped ? " | frottement" : "", sim.paused ? " | PAUSE" : "");
            SDL_SetWindowTitle(win, title);
        }
        ++frame;
        if (maxFrames >= 0 && frame >= maxFrames) {
            if (shot) {
                SDL_Surface* sf = SDL_CreateRGBSurfaceWithFormat(0, W, H, 32, SDL_PIXELFORMAT_ARGB8888);
                if (sf) { SDL_RenderReadPixels(ren, nullptr, sf->format->format, sf->pixels, sf->pitch); SDL_SaveBMP(sf, shot); SDL_FreeSurface(sf); }
            }
            std::printf("t=%.3f  E=%.6f  dérive=%.2e  lyapounov=%.3f\n", sim.t, energy(sim.P, sim.s), std::fabs(energy(sim.P, sim.s) - sim.E0) / std::fabs(sim.E0), sim.ly.value());
            break;
        }
        if (maxFrames < 0 && !(SDL_GetWindowFlags(win) & SDL_WINDOW_SHOWN)) SDL_Delay(16);
    }
    SDL_DestroyRenderer(ren); SDL_DestroyWindow(win); SDL_Quit();
    return 0;
}
