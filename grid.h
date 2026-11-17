// The electron as a wave: one number for every point of a 3D grid.
// Schrodinger (1926): energy = how much the wave bends + where it sits in the pull of the proton.
// The wave relaxes: every step takes away a little of its energy, until only the lowest shape is left.
#pragma once
#include "atom.h"
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <cstdio>

const double EV = 27.2114;  // 1 energy unit in electronvolts

struct Grid {
    int n; double L, h;            // n points per side, from -L to +L
    std::vector<double> V;         // the pull of the nuclei at every point
    Grid(int n, double L) : n(n), L(L), h(2 * L / n), V(size_t(n) * n * n, 0) {}
    double x(int i) const { return i * h - L; }
    size_t id(int i, int j, int k) const { return (size_t(i) * n + j) * n + k; }

    void add(const Atom& a) {
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) for (int k = 0; k < n; k++) {
            double dx = x(i) - a.x, dy = x(j) - a.y, dz = x(k) - a.z;
            V[id(i, j, k)] -= a.Z / std::sqrt(dx * dx + dy * dy + dz * dz);
        }
    }

    // H psi: -1/2 (bending) + V psi, zero outside the box
    void H(const std::vector<double>& p, std::vector<double>& out) const {
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) for (int k = 0; k < n; k++) {
            size_t c = id(i, j, k);
            double s = -6 * p[c];
            if (i > 0) s += p[c - n * n]; if (i < n - 1) s += p[c + n * n];
            if (j > 0) s += p[c - n];     if (j < n - 1) s += p[c + n];
            if (k > 0) s += p[c - 1];     if (k < n - 1) s += p[c + 1];
            out[c] = -0.5 * s / (h * h) + V[c] * p[c];
        }
    }

    static double dot(const std::vector<double>& a, const std::vector<double>& b) {
        double s = 0;
        for (size_t c = 0; c < a.size(); c++) s += a[c] * b[c];
        return s;
    }
    static void normalize(std::vector<double>& p) { double s = std::sqrt(dot(p, p)); for (double& v : p) v /= s; }

    std::vector<double> noise(int seed) const {
        std::mt19937 g(seed); std::uniform_real_distribution<double> u(-1, 1);
        std::vector<double> p(V.size());
        for (double& v : p) v = u(g);
        normalize(p);
        return p;
    }

    // relax until the energy stops moving; returns the energy
    double relax(std::vector<double>& p) const {
        std::vector<double> Hp(p.size());
        double dt = h * h / 4, E = 1e9;
        for (int step = 1; step <= 100000; step++) {
            H(p, Hp);
            for (size_t c = 0; c < p.size(); c++) p[c] -= dt * Hp[c];
            normalize(p);
            if (step % 100 == 0) {
                H(p, Hp);
                double e = dot(p, Hp);
                if (step % 1000 == 0) printf("  step %5d   energy %.3f eV\n", step, e * EV);
                if (!(std::fabs(e - E) > 1e-7)) { printf("  stopped after %d steps\n", step); return e; }
                E = e;
            }
        }
        return E;
    }
};

// points where the electron is likely: each point of the grid gets dots in proportion to psi^2
std::vector<float> cloud(const Grid& g, const std::vector<double>& p, int count) {
    std::vector<double> sum(p.size());
    double total = 0;
    for (size_t c = 0; c < p.size(); c++) sum[c] = total += p[c] * p[c];
    std::mt19937 rng(7); std::uniform_real_distribution<double> u(0, 1);
    std::vector<float> out;
    if (!(total > 0)) return out;
    for (int m = 0; m < count; m++) {
        size_t c = std::lower_bound(sum.begin(), sum.end(), u(rng) * total) - sum.begin();
        if (c >= p.size()) c = p.size() - 1;
        int i = c / (g.n * g.n), j = c / g.n % g.n, k = c % g.n;
        out.insert(out.end(), {float(g.x(i) + g.h * (u(rng) - 0.5)), float(g.x(j) + g.h * (u(rng) - 0.5)),
                               float(g.x(k) + g.h * (u(rng) - 0.5)), float(p[c])});
    }
    return out;
}
