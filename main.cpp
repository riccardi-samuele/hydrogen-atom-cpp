#include "window.h"
#include "grid.h"
#include <chrono>

int main() {
    Atom H;
    std::vector<float> dots;
    for (int n : {32, 48, 64, 80, 96}) {  // the same box, finer and finer
        auto t0 = std::chrono::steady_clock::now();
        Grid g(n, 8);
        g.add(H);
        std::vector<double> psi = g.noise(1);
        double E = g.relax(psi);
        double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        printf("one point every %4.1f pm   energy: %.2f eV   (%.1f s)\n", g.h * 52.9177, E * EV, s);
        dots = cloud(g, psi, 100000);
    }

    Window win("hydrogen");
    while (win.open()) {
        win.begin(16);
        win.points(dots, 2, 0.25f);
        win.dot(H.x, H.y, H.z, 10, win.t().b);
        win.end();
    }
}
