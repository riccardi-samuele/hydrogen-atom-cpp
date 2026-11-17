#include "window.h"
#include "shapes.h"
#include <GL/gl.h>
#include <chrono>

// all 14 states on one grid, relaxed and sorted
std::vector<Answers> solve(int n, double L) {
    auto t0 = std::chrono::steady_clock::now();
    Atom H;
    Grid g(n, L);
    g.add(H);
    std::vector<std::vector<double>> S;
    for (int s = 1; s <= 14; s++) S.push_back(g.noise(s));
    relax_all(g, S);
    std::vector<Answers> ans = sort_shapes(g, S);
    printf("one point every %.1f pm: %.0f s\n", g.h * 52.9177, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    return ans;
}

int main() {
    // The trick from the first atom: two grids, a straight line in spacing^2, followed down to zero.
    double ha = 2 * 24.0 / 96, hb = 2 * 24.0 / 120;
    std::vector<Answers> A = solve(96, 24), B = solve(120, 24), ans = B;
    for (size_t s = 0; s < ans.size(); s++) {
        size_t a = 0;  // the same shape on the coarse grid, found by its name (the order inside a shell can change)
        while (a < A.size() && A[a].name != B[s].name) a++;
        if (a == A.size()) { printf("the two grids disagree on %s\n", B[s].name.c_str()); return 1; }
        ans[s].E = (B[s].E * ha * ha - A[a].E * hb * hb) / (ha * ha - hb * hb);
        printf("%-7s %6.2f eV  %6.2f eV  -> no spacing: %6.2f eV\n", ans[s].name.c_str(), A[a].E * EV, B[s].E * EV, ans[s].E * EV);
    }

    // The red line: the electron drops from shell 3 to shell 2 and the energy leaves as light.
    // The brightest drop is 3d -> 2p. Wavelength = 1239.84 eV nm / energy.
    double e3 = 0, e2 = 0; int k3 = 0, k2 = 0;
    for (auto& a : ans) {
        if (a.name.rfind("3d", 0) == 0) { e3 += a.E; k3++; }
        if (a.name.rfind("2p", 0) == 0) { e2 += a.E; k2++; }
    }
    double drop = (e3 / k3 - e2 / k2) * EV, nm = 1239.84 / drop;
    printf("3d -> 2p: %.3f eV of light, %.1f nm\n", drop, nm);

    Window win("hydrogen");
    win.spin = false; win.yaw = win.pitch = 0;
    while (win.open()) {
        win.begin(16);
        float c[3];
        for (int i = 0; i < 1500; i++) {  // the rainbow, 380 to 750 nm, dim
            double w = 380 + 370.0 * i / 1500;
            Window::rainbow(w, c);
            for (int y = 0; y < 6; y++) win.dot(-9 + 18 * (w - 380) / 370, -1.5f + 0.5f * y, 0, 3, c, 0.15f);
        }
        Window::rainbow(nm, c);  // my line, bright
        for (int y = 0; y < 40; y++) win.dot(-9 + 18 * (nm - 380) / 370, -2.5f + 0.13f * y, 0, 4, c);
        win.end();
    }
}
