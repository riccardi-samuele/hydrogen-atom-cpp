#include "window.h"
#include "atom.h"
#include <vector>

// Units: 1 length = 52.9 pm, 1 time = 0.0242 fs, electron mass and charge = 1.
const double PM = 52.9177, FS = 0.0241888;
const double C = 137.036;  // speed of light in these units

int main() {
    Atom H;
    double x = 1, y = 0, vx = 0, vy = 1;  // I start the electron 53 pm away, at the speed of a circle
    printf("start: %.0f pm from the proton\n", x * PM);

    // A charge that turns gives off light (Larmor, 1897): it loses 2/3 a^2 / c^3 energy every second.
    std::vector<float> path;  // where the electron is, every 0.0242 fs
    double t = 0, next = 0;
    long turns = 0;
    while (true) {
        double dx = x - H.x, dy = y - H.y, r = std::sqrt(dx * dx + dy * dy);
        if (r < 0.01) break;  // closer than half a picometer: it has fallen in
        double dt = 0.01 * r * std::sqrt(r);
        double ax = -dx / (r * r * r), ay = -dy / (r * r * r);
        double a2 = ax * ax + ay * ay, v2 = vx * vx + vy * vy;
        double P = 2.0 / 3.0 * a2 / (C * C * C);
        ax -= P * vx / v2; ay -= P * vy / v2;  // the lost energy slows it down
        vx += ax * dt; vy += ay * dt;
        double y0 = y;
        x += vx * dt; y += vy * dt; t += dt;
        if (y0 < 0 && y >= 0) turns++;
        if (t >= next) { path.push_back(x); path.push_back(y); next += 1; }
    }
    printf("fell into the proton after %.1f ps and %ld turns\n", t * FS / 1000, turns);

    Window win("hydrogen");
    size_t n = path.size() / 2, shown = 0;
    while (win.open()) {
        shown = std::min(n - 1, shown + n / 600);  // the whole fall in 10 seconds
        win.begin(4);
        win.dot(H.x, H.y, H.z, 10, win.t().b);
        for (size_t i = shown > 300 ? shown - 300 : 0; i < shown; i++) win.dot(path[2 * i], path[2 * i + 1], 0, 2, win.t().a, 0.5f);
        win.dot(path[2 * shown], path[2 * shown + 1], 0, 6, win.t().a);
        win.end(long(shown * FS));  // femtoseconds since the start
    }
}
