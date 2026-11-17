#include "window.h"
#include "grid.h"
#include <GL/gl.h>

int main() {
    Atom H;
    Grid g(80, 20);  // a bigger box, 2100 pm wide, for the next shapes
    g.add(H);
    std::vector<std::vector<float>> shapes;
    for (int s = 1; s <= 2; s++) {  // the second answer starts from different noise
        std::vector<double> psi = g.noise(s);
        double E = g.relax(psi);
        printf("state %d   energy: %.2f eV\n", s, E * EV);
        shapes.push_back(cloud(g, psi, 60000));
    }

    Window win("hydrogen");
    win.spin = false; win.yaw = win.pitch = 0;
    while (win.open()) {
        win.begin(40);
        for (int s = 0; s < int(shapes.size()); s++) {
            glPushMatrix();
            glTranslatef(-11 + 22 * s, 0, 0);
            glRotatef(20, 1, 0, 0);
            glRotatef(win.frame * 0.3f, 0, 1, 0);
            win.points(shapes[s], 2, 0.25f);
            glPopMatrix();
        }
        win.end();
    }
}
