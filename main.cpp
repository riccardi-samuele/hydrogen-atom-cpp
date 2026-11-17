#include "window.h"
#include "atom.h"

int main() {
    Window win("hydrogen");
    Atom H;
    while (win.open()) {
        win.begin(4);
        win.dot(H.x, H.y, H.z, 8, win.t().b);
        win.end();
    }
}
