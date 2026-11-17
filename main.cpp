#include "window.h"
#include "atom.h"
#include <vector>

// Units: 1 length = 52.9 pm, 1 time = 0.0242 fs, electron mass and charge = 1.
const double PM = 52.9177, FS = 0.0241888;

int main() {
    Atom H;
    double x = 1, y = 0, vx = 0, vy = 1;  // I start the electron 53 pm away, at the speed of a circle
    printf("start: %.0f pm from the proton\n", x * PM);

    Window win("hydrogen");
    std::vector<float> trail;
    double t = 0, dt = 0.002;
    while (win.open()) {
        for (int i = 0; i < 25; i++) {  // the proton pulls the electron: a = -r / r^3
            double dx = x - H.x, dy = y - H.y, r = std::sqrt(dx * dx + dy * dy);
            vx -= dx / (r * r * r) * dt; vy -= dy / (r * r * r) * dt;
            x += vx * dt; y += vy * dt; t += dt;
        }
        trail.push_back(x); trail.push_back(y);
        if (trail.size() > 600) trail.erase(trail.begin(), trail.begin() + 2);
        if (win.frame % 300 == 0) printf("t = %.2f fs   distance: %.1f pm\n", t * FS, std::sqrt(x * x + y * y) * PM);

        win.begin(4);
        win.dot(H.x, H.y, H.z, 10, win.t().b);
        for (size_t i = 0; i < trail.size(); i += 2) win.dot(trail[i], trail[i + 1], 0, 2, win.t().a, float(i) / trail.size());
        win.dot(x, y, 0, 6, win.t().a);
        win.end();
    }
}
