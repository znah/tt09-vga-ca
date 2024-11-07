#include <SDL.h>
#include "Vproject.h"
#include "verilated.h"

const int PAD = 48;
const int W = 640+PAD*2;
const int H = 480+PAD*2;
struct Bits {
    uint8_t r1:1, g1:1, b1:1, vsync:1, r0:1, g0:1, b0:1, hsync:1;
};

Uint32 pixels[W * H];
Vproject mod;


uint32_t bits2channel(uint8_t lo, uint8_t hi, int ch) {
    return (hi*170+lo*85) << (ch*8);
}

int main(int argc, char* argv[]) {
    SDL_Init(SDL_INIT_VIDEO);

    SDL_Window* window = SDL_CreateWindow("VGA CA demo", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, W, H, SDL_WINDOW_SHOWN);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    SDL_Texture* texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, W, H);

    bool quit = false;
    SDL_Event event;

    mod.rst_n = 0;
    mod.eval();
    mod.rst_n = 1;
    mod.eval();

    Bits prev = *(Bits*)(&mod.uo_out);
    int x=0, y=0;
    while (!quit) {
        while (SDL_PollEvent(&event) != 0) {
            if (event.type == SDL_QUIT) {
                quit = true;
            }
            switch( event.key.keysym.sym ) {
                case SDLK_ESCAPE: quit = true; break;
            }            
        }
        while(1) {
            mod.clk = !mod.clk;
            mod.eval();
            Bits v = *(Bits*)(&mod.uo_out);
            if (prev.hsync && !v.hsync) {
                x = 0;
                y += 1;
            }
            bool vsync = prev.vsync && !v.vsync;
            if (vsync) {
                x = 0;
                y = 0;
            }
            if (x<W && y<H) {
                pixels[y*W + x] = bits2channel(v.b0, v.b1, 0) | 
                    bits2channel(v.g0, v.g1, 1) | bits2channel(v.r0, v.r1, 2);
            }
            prev = v;

            x += mod.clk;
            if (vsync) break;
        }

        SDL_UpdateTexture(texture, nullptr, pixels, W * sizeof(Uint32));
        SDL_RenderClear(renderer);
        SDL_RenderCopy(renderer, texture, nullptr, nullptr);
        SDL_RenderPresent(renderer);

        SDL_Delay(5);  // Add a small delay to avoid using 100% CPU
    }

    SDL_DestroyTexture(texture);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}

