import PIL.Image
import numpy as np
import vcdvcd

vcd = vcdvcd.VCDVCD('tb.vcd')
output = vcd['tb.uo_out[7:0]']

pixel_x, pixel_y = 0, 0
PAD = 48*2
W, H = 640+PAD, 480+PAD
screen = np.zeros([H, W], np.uint8)
x, y = 0, 0

prev = 0
for time, clk in vcd['tb.clk'].tv:
    if clk != '0':
        continue
    val = int(output[time].replace('x', '1'), 2)
    if x<W and y<H:
        screen[y,x] = val
    x += 1
    if prev&0x80 and not val&0x80: #hsync
        x, y = 0, y+1
    if prev&0x08 and not val&0x08: #vsync
        x, y = 0, 0
    prev = val

bits = np.unpackbits(screen, bitorder='little').reshape(H, W, 8)
rgb = bits[...,:3]*170 + bits[...,4:7]*85
PIL.Image.fromarray(rgb).save('screen.png')
