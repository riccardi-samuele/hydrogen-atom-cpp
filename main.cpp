#include "window.h"
#include "grid.h"
#include <chrono>

int main() {
    Atom H;
    std::vector<double> h2, E;
    for (int n : {32, 48, 64, 80, 96}) {  // the same box, finer and finer
        auto t0 = std::chrono::steady_clock::now();
        Grid g(n, 8);
        g.add(H);
        std::vector<double> psi = g.noise(1);
        double e = g.relax(psi);
        double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        printf("one point every %4.1f pm   energy: %.2f eV   (%.1f s)\n", g.h * 52.9177, e * EV, s);
        h2.push_back(g.h * g.h); E.push_back(e * EV);
    }

    // The error shrinks with the square of the spacing, so energy against spacing^2 is a straight line.
    // A line through the three finest grids, followed down to spacing zero.
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (int i = 2; i < 5; i++) { sx += h2[i]; sy += E[i]; sxx += h2[i] * h2[i]; sxy += h2[i] * E[i]; }
    double slope = (3 * sxy - sx * sy) / (3 * sxx - sx * sx), E0 = (sy - slope * sx) / 3;
    printf("with no spacing at all: %.2f eV\n", E0);

    Window win("hydrogen");
    win.spin = false; win.yaw = win.pitch = 0;
    auto X = [](double v) { return float(-6 + v / 0.25 * 12); };     // spacing^2 from 0 to 0.25
    auto Y = [](double e) { return float((e + 13.25) / 0.4 * 3); };   // energy around -13.25 eV
    while (win.open()) {
        win.begin(16);
        for (int i = 0; i <= 300; i++) { double v = 0.26 * i / 300; win.dot(X(v), Y(E0 + slope * v), 0, 2, win.t().a, 0.5f); }
        for (int i = 0; i < 5; i++) win.dot(X(h2[i]), Y(E[i]), 0, 9, win.t().a);
        win.dot(X(0), Y(E0), 0, 12, win.t().b);
        win.end();
    }
}
