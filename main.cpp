#include "window.h"
#include "grid.h"

int main() {
    Atom H;
    Grid g(32, 8);  // 32 x 32 x 32 points, a box 850 pm wide
    g.add(H);
    printf("grid: %d^3 points, one every %.0f pm\n", g.n, g.h * 52.9177);

    std::vector<double> psi = g.noise(1);
    double E = g.relax(psi);
    printf("energy: %.2f eV\n", E * EV);

    std::vector<float> dots = cloud(g, psi, 100000);
    Window win("hydrogen");
    while (win.open()) {
        win.begin(16);
        win.points(dots, 2, 0.25f);
        win.dot(H.x, H.y, H.z, 10, win.t().b);
        win.end();
    }
}
