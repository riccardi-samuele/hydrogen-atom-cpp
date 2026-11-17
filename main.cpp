#include "window.h"
#include "shapes.h"
#include <GL/gl.h>
#include <chrono>

int main() {
    Atom H;
    Grid g(96, 24);  // 2500 pm wide
    g.add(H);
    auto t0 = std::chrono::steady_clock::now();
    std::vector<std::vector<double>> S;
    for (int s = 1; s <= 14; s++) S.push_back(g.noise(s));
    relax_all(g, S);
    std::vector<Answers> ans = sort_shapes(g, S);
    printf("%.0f s\n", std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    std::vector<std::vector<float>> shapes;
    for (size_t s = 0; s < S.size(); s++) {
        printf("%-7s energy: %6.2f eV   turn: %.2f   mirrors: %c%c%c\n", ans[s].name.c_str(), ans[s].E * EV, ans[s].L2,
               ans[s].px ? '-' : '+', ans[s].py ? '-' : '+', ans[s].pz ? '-' : '+');
        shapes.push_back(cloud(g, S[s], 40000));
    }

    Window win("hydrogen");
    win.spin = false; win.yaw = win.pitch = 0;
    while (win.open()) {
        win.begin(150);
        for (int s = 0; s < int(shapes.size()); s++) {  // rows of 1, 4 and 9
            int row = s < 1 ? 0 : s < 5 ? 1 : 2, col = s < 1 ? s : s < 5 ? s - 1 : s - 5, cols = row == 0 ? 1 : row == 1 ? 4 : 9;
            glPushMatrix();
            glTranslatef((col - (cols - 1) / 2.0f) * 38, 55 - row * 55, 0);
            glRotatef(20, 1, 0, 0);
            glRotatef(win.frame * 0.3f, 0, 1, 0);
            win.points(shapes[s], 1.5f, 0.12f);
            glPopMatrix();
        }
        win.end();
    }
}
