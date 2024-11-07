# SPDX-FileCopyrightText: © 2024 Tiny Tapeout
# SPDX-License-Identifier: Apache-2.0

import cocotb
from cocotb.clock import Clock
from cocotb.triggers import ClockCycles, FallingEdge
# import numpy as np


@cocotb.test()
async def test_project(dut):
    dut._log.info("Start")

    # Set the clock period to 1/26 us (25 MHz)
    clock = Clock(dut.clk, 1/25, units="us")
    cocotb.start_soon(clock.start())

    # Reset
    dut._log.info("Reset")
    dut.ena.value = 1
    dut.ui_in.value = 0
    dut.uio_in.value = 0
    dut.rst_n.value = 0
    await ClockCycles(dut.clk, 10)
    dut.rst_n.value = 1

    pixel_x, pixel_y = 0, 0
    PAD = 48*2
    W, H = 640+PAD, 480+PAD
    ROW_CLOCKS = 800
    #screen = np.zeros([H, W], np.uint8)

    async def hsync():
        nonlocal pixel_x, pixel_y
        while True:
            await FallingEdge(dut.hsync)
            pixel_x = 0
            pixel_y += 1
            print('.', end='', flush=True)
    async def vsync():
        nonlocal pixel_x, pixel_y
        while True:
            await FallingEdge(dut.vsync)
            pixel_x = 0
            pixel_y = 0
    cocotb.start_soon(hsync())
    cocotb.start_soon(vsync())

    dut._log.info("Test project behavior")

    print(dir(dut.uo_out))
    for i in range(ROW_CLOCKS*4*8):
        await FallingEdge(dut.clk)
        # if pixel_x < W and pixel_y<H:
        #     screen[pixel_y,pixel_x] = dut.uo_out
        pixel_x += 1
    #bits = np.unpackbits(screen, bitorder='little').reshape(H, W, 8)
    #rgb = bits[...,:3]*170 + bits[...,4:7]*85
    #print(rgb.shape, bits.shape)
    
    # Set the input values you want to test
    #dut.ui_in.value = 20
    #dut.uio_in.value = 30

    # Wait for one clock cycle to see the output values
    await ClockCycles(dut.clk, 1)

    # The following assersion is just an example of how to check the output values.
    # Change it to match the actual expected output of your module:
    # assert dut.uo_out.value == 50

    # Keep testing the module by changing the input values, waiting for
    # one or more clock cycles, and asserting the expected output values.
