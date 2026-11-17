// The window: zoom, drag to turn, three color themes, a frame counter. No physics in here.
#pragma once
#include <SDL2/SDL.h>
#include <GL/gl.h>
#include <cmath>
#include <cstdio>
#include <vector>

struct Theme { const char* name; float bg[3], a[3], b[3]; bool glow; };
const Theme THEMES[3] = {
    {"chalk", {0.10f, 0.11f, 0.12f}, {0.95f, 0.95f, 0.90f}, {0.95f, 0.55f, 0.25f}, true},
    {"night", {0.00f, 0.00f, 0.02f}, {0.35f, 0.65f, 1.00f}, {1.00f, 0.45f, 0.20f}, true},
    {"paper", {0.96f, 0.94f, 0.88f}, {0.10f, 0.12f, 0.20f}, {0.75f, 0.15f, 0.10f}, false},
};

struct Window {
    SDL_Window* sdl;
    int w = 1280, h = 720, frame = 0, theme = 0;
    float zoom = 1, yaw = 0.6f, pitch = 0.35f;
    bool spin = true, dragging = false;

    Window(const char* title) {
        SDL_Init(SDL_INIT_VIDEO);
        sdl = SDL_CreateWindow(title, 0, 0, w, h, SDL_WINDOW_OPENGL);
        SDL_GL_CreateContext(sdl);
        SDL_GL_SetSwapInterval(1);
        SDL_ShowCursor(SDL_DISABLE);
        printf("window %dx%d, theme: %s\n", w, h, THEMES[theme].name);
    }
    const Theme& t() const { return THEMES[theme]; }

    bool open() {  // handles keys and mouse, false when it is time to quit
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) return false;
            if (e.type == SDL_KEYDOWN) {
                SDL_Keycode k = e.key.keysym.sym;
                if (k == SDLK_ESCAPE) return false;
                if (k == SDLK_t) { theme = (theme + 1) % 3; printf("theme: %s\n", t().name); }
                if (k == SDLK_SPACE) spin = !spin;
            }
            if (e.type == SDL_MOUSEWHEEL) zoom *= e.wheel.y > 0 ? 1.1f : 1 / 1.1f;
            if (e.type == SDL_MOUSEBUTTONDOWN) dragging = true;
            if (e.type == SDL_MOUSEBUTTONUP) dragging = false;
            if (e.type == SDL_MOUSEMOTION && dragging) { yaw += e.motion.xrel * 0.01f; pitch += e.motion.yrel * 0.01f; }
        }
        return true;
    }

    // 3D view of a box `size` wide, centered on the proton
    void begin(float size) {
        if (spin) yaw += 0.004f;
        glViewport(0, 0, w, h);
        glClearColor(t().bg[0], t().bg[1], t().bg[2], 1);
        glClear(GL_COLOR_BUFFER_BIT);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, t().glow ? GL_ONE : GL_ONE_MINUS_SRC_ALPHA);
        glMatrixMode(GL_PROJECTION); glLoadIdentity();
        float a = float(w) / h, n = 0.1f * size;
        glFrustum(-a * n * 0.5f, a * n * 0.5f, -n * 0.5f, n * 0.5f, n, 10 * size);
        glMatrixMode(GL_MODELVIEW); glLoadIdentity();
        glTranslatef(0, 0, -1.3f * size / zoom);
        glRotatef(pitch * 57.3f, 1, 0, 0);
        glRotatef(yaw * 57.3f, 0, 1, 0);
    }

    void dot(float x, float y, float z, float px, const float* c, float alpha = 1) {
        glPointSize(px);
        glEnable(GL_POINT_SMOOTH);
        glBegin(GL_POINTS); glColor4f(c[0], c[1], c[2], alpha); glVertex3f(x, y, z); glEnd();
    }

    // a number drawn with seven lines per digit, like a microwave clock
    void number(long v, float x, float y, float s) {
        static const unsigned char seg[10] = {0x3f, 0x06, 0x5b, 0x4f, 0x66, 0x6d, 0x7d, 0x07, 0x7f, 0x6f};
        static const float L[7][4] = {{0,2,1,2},{1,2,1,1},{1,1,1,0},{0,0,1,0},{0,0,0,1},{0,1,0,2},{0,1,1,1}};
        char buf[32]; snprintf(buf, sizeof buf, "%ld", v);
        glLineWidth(2);
        glBegin(GL_LINES);
        glColor4f(t().a[0], t().a[1], t().a[2], 0.6f);
        for (int i = 0; buf[i]; i++) {
            if (buf[i] < '0' || buf[i] > '9') continue;
            for (int g = 0; g < 7; g++) if (seg[buf[i] - '0'] >> g & 1) {
                glVertex2f(x + i * 1.6f * s + L[g][0] * s, y + L[g][1] * s);
                glVertex2f(x + i * 1.6f * s + L[g][2] * s, y + L[g][3] * s);
            }
        }
        glEnd();
    }

    void end() {
        glMatrixMode(GL_PROJECTION); glLoadIdentity(); glOrtho(0, w, 0, h, -1, 1);
        glMatrixMode(GL_MODELVIEW); glLoadIdentity();
        number(frame, 20, 20, 10);
        SDL_GL_SwapWindow(sdl);
        static Uint32 last = 0;  // never faster than 60 frames a second
        while (SDL_GetTicks() - last < 16) SDL_Delay(1);
        last = SDL_GetTicks();
        frame++;
    }
};
