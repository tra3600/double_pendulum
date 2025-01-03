#include <iostream>
#include <cmath>
//#include <SDL2/SDL.h>
//#include "SDL2/SDL2.h"
#include "../include/SDL.h"

// Constantes physiques
const double g = 9.81; // Accélération due à la gravité (m/s^2)

// Structure pour stocker les informations du pendule
struct Pendulum {
    double length;    // Longueur du pendule (m)
    double angle;     // Angle par rapport à la verticale (rad)
    double angularVelocity;  // Vitesse angulaire (rad/s)
    double angularAcceleration; // Accélération angulaire (rad/s^2)
};

// Fonction pour mettre à jour le mouvement du pendule double
void updatePendulum(Pendulum& p1, Pendulum& p2, double dt) {
    double num1 = -g * (2 * p1.angle + p2.angle) * sin(p1.angle);
    double num2 = -p2.angle * p2.angularVelocity * p2.angularVelocity * p2.length * sin(p1.angle - p2.angle);
    double num3 = -2 * sin(p1.angle - p2.angle) * p2.angularVelocity * p2.angularVelocity * p2.length * cos(p1.angle - p2.angle);
    double denom = p1.length * (2 - cos(2 * (p1.angle - p2.angle)));

    p1.angularAcceleration = (num1 + num2 + num3) / denom;

    num1 = 2 * sin(p1.angle - p2.angle) * (p1.angularVelocity * p1.angularVelocity * p1.length * cos(p1.angle - p2.angle) + g * cos(p1.angle));
    num2 = p2.angularVelocity * p2.angularVelocity * p2.length * cos(p1.angle - p2.angle);
    denom = p2.length * (2 - cos(2 * (p1.angle - p2.angle)));

    p2.angularAcceleration = (num1 + num2) / denom;

    p1.angularVelocity += p1.angularAcceleration * dt;
    p2.angularVelocity += p2.angularAcceleration * dt;
    p1.angle += p1.angularVelocity * dt;
    p2.angle += p2.angularVelocity * dt;
}

// Fonction pour dessiner le pendule double
void drawPendulum(SDL_Renderer* renderer, const Pendulum& p1, const Pendulum& p2, int originX, int originY) {
    int x1 = originX + p1.length * sin(p1.angle) * 100;
    int y1 = originY + p1.length * cos(p1.angle) * 100;
    int x2 = x1 + p2.length * sin(p2.angle) * 100;
    int y2 = y1 + p2.length * cos(p2.angle) * 100;

    SDL_SetRenderDrawColor(renderer, 255, 255, 255, SDL_ALPHA_OPAQUE);
    SDL_RenderDrawLine(renderer, originX, originY, x1, y1);
    SDL_RenderDrawLine(renderer, x1, y1, x2, y2);
    SDL_RenderFillCircle(renderer, x1, y1, 5);
    SDL_RenderFillCircle(renderer, x2, y2, 5);
}

int main(int argc, char* argv[]) {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        std::cerr << "Erreur lors de l'initialisation de SDL: " << SDL_GetError() << std::endl;
        return 1;
    }

    SDL_Window* window = SDL_CreateWindow("Simulation de Pendule Double", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, 800, 600, SDL_WINDOW_SHOWN);
    if (!window) {
        std::cerr << "Erreur lors de la création de la fenêtre: " << SDL_GetError() << std::endl;
        SDL_Quit();
        return 1;
    }

    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (!renderer) {
        std::cerr << "Erreur lors de la création du renderer: " << SDL_GetError() << std::endl;
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    // Initialisation des pendules
    Pendulum p1 = {1.0, M_PI / 4, 0.0, 0.0};
    Pendulum p2 = {1.0, M_PI / 4, 0.0, 0.0};

    bool quit = false;
    SDL_Event e;
    double dt = 0.01;

    while (!quit) {
        while (SDL_PollEvent(&e) != 0) {
            if (e.type == SDL_QUIT) {
                quit = true;
            }
        }

        // Mettre à jour le mouvement des pendules
        updatePendulum(p1, p2, dt);

        // Effacer l'écran
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, SDL_ALPHA_OPAQUE);
        SDL_RenderClear(renderer);

        // Dessiner les pendules
        drawPendulum(renderer, p1, p2, 400, 300);

        // Mettre à jour l'affichage
        SDL_RenderPresent(renderer);

        SDL_Delay(16); // Environ 60 FPS
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}