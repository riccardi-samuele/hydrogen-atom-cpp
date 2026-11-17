// States with the same energy can mix in any way, so the energy alone does not fix their shape.
// So I ask every state more questions, and keep the mixes that give clean answers:
//   does it look the same in a mirror?  (x -> -x, y -> -y, z -> -z)
//   how much does it turn?  (L^2 in total, Lz^2 around the z axis)
#pragma once
#include "grid.h"
#include <string>

// eigenvalues and eigenvectors of a small symmetric matrix (Jacobi, 1846)
void jacobi(std::vector<double>& A, std::vector<double>& Q, int m) {
    Q.assign(m * m, 0);
    for (int i = 0; i < m; i++) Q[i * m + i] = 1;
    for (int sweep = 0; sweep < 100; sweep++) {
        double off = 0;
        for (int p = 0; p < m; p++) for (int q = p + 1; q < m; q++) off += A[p * m + q] * A[p * m + q];
        if (off < 1e-24) return;
        for (int p = 0; p < m; p++) for (int q = p + 1; q < m; q++) {
            double apq = A[p * m + q];
            if (apq == 0) continue;
            double th = (A[q * m + q] - A[p * m + p]) / (2 * apq);
            double t = (th >= 0 ? 1 : -1) / (std::fabs(th) + std::sqrt(th * th + 1)), c = 1 / std::sqrt(t * t + 1), s = t * c;
            for (int k = 0; k < m; k++) { double a = A[k * m + p], b = A[k * m + q]; A[k * m + p] = c * a - s * b; A[k * m + q] = s * a + c * b; }
            for (int k = 0; k < m; k++) { double a = A[p * m + k], b = A[q * m + k]; A[p * m + k] = c * a - s * b; A[q * m + k] = s * a + c * b; }
            for (int k = 0; k < m; k++) { double a = Q[k * m + p], b = Q[k * m + q]; Q[k * m + p] = c * a - s * b; Q[k * m + q] = s * a + c * b; }
        }
    }
}

// new states = mixes of the old ones, with the weights in the columns of Q
void mix(std::vector<std::vector<double>>& S, const std::vector<double>& Q, const std::vector<int>& use) {
    int m = use.size();
    std::vector<std::vector<double>> T(m, std::vector<double>(S[0].size(), 0));
    for (int a = 0; a < m; a++) for (int b = 0; b < m; b++) {
        double c = Q[b * m + a];
        #pragma omp parallel for
        for (size_t i = 0; i < T[a].size(); i++) T[a][i] += c * S[use[b]][i];
    }
    for (int a = 0; a < m; a++) S[use[a]] = T[a];
}

// the matrix <a|op|b> for every pair of states in `use`
template <class Op> std::vector<double> ask(const std::vector<std::vector<double>>& S, const std::vector<int>& use, Op op) {
    int m = use.size();
    std::vector<std::vector<double>> O(m, std::vector<double>(S[0].size()));
    for (int a = 0; a < m; a++) op(S[use[a]], O[a]);
    std::vector<double> A(m * m);
    for (int a = 0; a < m; a++) for (int b = 0; b < m; b++) A[a * m + b] = Grid::dot(S[use[a]], O[b]);
    return A;
}

struct Answers { double E, L2, Lz2; int px, py, pz; std::string name; };

std::vector<Answers> sort_shapes(const Grid& g, std::vector<std::vector<double>>& S) {
    int n = g.n, m = S.size();
    auto mirror = [&](int axis) {
        return [&, axis](const std::vector<double>& p, std::vector<double>& o) {
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) for (int k = 0; k < n; k++)
                o[g.id(i, j, k)] = p[g.id(axis == 0 ? n - 1 - i : i, axis == 1 ? n - 1 - j : j, axis == 2 ? n - 1 - k : k)];
        };
    };
    // turning around an axis: (y d/dz - z d/dy) and the other two; L^2 = sum of their squares
    auto turn2 = [&](std::vector<int> axes) {
        return [&, axes](const std::vector<double>& p, std::vector<double>& o) {
            std::fill(o.begin(), o.end(), 0.0);
            std::vector<double> d(p.size()), dd(p.size());
            for (int ax : axes) {
                int a = (ax + 1) % 3, b = (ax + 2) % 3;  // turning around ax mixes the other two
                auto D = [&](const std::vector<double>& in, std::vector<double>& out) {
                    #pragma omp parallel for
                    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) for (int k = 0; k < n; k++) {
                        int c[3] = {i, j, k};
                        auto at = [&](int axis, int s) { int q[3] = {c[0], c[1], c[2]}; q[axis] += s;
                            return q[axis] < 0 || q[axis] >= n ? 0.0 : in[g.id(q[0], q[1], q[2])]; };
                        double da = (at(a, 1) - at(a, -1)) / (2 * g.h), db = (at(b, 1) - at(b, -1)) / (2 * g.h);
                        out[g.id(i, j, k)] = g.x(c[a]) * db - g.x(c[b]) * da;
                    }
                };
                D(p, d); D(d, dd);
                for (size_t c = 0; c < o.size(); c++) o[c] -= dd[c];
            }
        };
    };

    // 1. the energy first
    std::vector<int> all(m);
    for (int a = 0; a < m; a++) all[a] = a;
    std::vector<double> A = ask(S, all, [&](const std::vector<double>& p, std::vector<double>& o) { g.H(p, o); }), Q;
    jacobi(A, Q, m);
    std::vector<int> order(all);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return A[a * m + a] < A[b * m + b]; });
    std::vector<double> Qs(m * m);
    for (int a = 0; a < m; a++) for (int k = 0; k < m; k++) Qs[k * m + a] = Q[k * m + order[a]];
    mix(S, Qs, all);
    std::vector<double> E(m);
    for (int a = 0; a < m; a++) E[a] = A[order[a] * m + order[a]];

    // 2. inside every group of equal energy, the other questions pick the clean mixes
    std::vector<Answers> out(m);
    for (int start = 0, shell = 1; start < m; shell++) {
        std::vector<int> group;
        for (int a = start; a < m && (a == start || E[a] - E[a - 1] < 1 / EV); a++) group.push_back(a);  // a jump of 1 eV starts a new shell
        int k = group.size();
        auto L2 = ask(S, group, turn2({0, 1, 2})), Lz2 = ask(S, group, turn2({2}));
        auto Mx = ask(S, group, mirror(0)), My = ask(S, group, mirror(1)), Mz = ask(S, group, mirror(2));
        std::vector<double> M(k * k), W;
        for (int i = 0; i < k * k; i++) M[i] = 10 * L2[i] + Lz2[i] + 0.01 * Mx[i] + 0.02 * My[i] + 0.04 * Mz[i];
        jacobi(M, W, k);
        mix(S, W, group);
        L2 = ask(S, group, turn2({0, 1, 2})); Lz2 = ask(S, group, turn2({2}));
        Mx = ask(S, group, mirror(0)); My = ask(S, group, mirror(1)); Mz = ask(S, group, mirror(2));
        auto Hg = ask(S, group, [&](const std::vector<double>& p, std::vector<double>& o) { g.H(p, o); });
        for (int i = 0; i < k; i++) {
            Answers& r = out[group[i]];
            r = {Hg[i * k + i], L2[i * k + i], Lz2[i * k + i], Mx[i * k + i] < 0, My[i * k + i] < 0, Mz[i * k + i] < 0, ""};
            int l = std::lround((std::sqrt(1 + 4 * r.L2) - 1) / 2);
            std::string odd = std::string(r.px ? "x" : "") + (r.py ? "y" : "") + (r.pz ? "z" : "");
            r.name = std::to_string(shell) + "spdf"[l] + (l == 0 ? "" : l == 1 ? odd : !odd.empty() ? odd : r.Lz2 > 2 ? "x2-y2" : "z2");
        }
        start += k;
    }
    // same shell together, simplest shapes first
    std::vector<int> idx(all);
    std::stable_sort(idx.begin(), idx.end(), [&](int a, int b) {
        return out[a].name[0] != out[b].name[0] ? out[a].name[0] < out[b].name[0] : out[a].L2 < out[b].L2 - 0.5; });
    std::vector<std::vector<double>> S2; std::vector<Answers> o2;
    for (int a : idx) { S2.push_back(S[a]); o2.push_back(out[a]); }
    S = S2;
    return o2;
}
