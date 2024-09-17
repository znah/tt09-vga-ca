/*
 * Copyright (c) 2024 Uri Shaked
 * SPDX-License-Identifier: Apache-2.0
 */

`default_nettype none

module tt_um_vga_example(
  input  wire [7:0] ui_in,    // Dedicated inputs
  output wire [7:0] uo_out,   // Dedicated outputs
  input  wire [7:0] uio_in,   // IOs: Input path
  output wire [7:0] uio_out,  // IOs: Output path
  output wire [7:0] uio_oe,   // IOs: Enable path (active high: 0=input, 1=output)
  input  wire       ena,      // always 1 when the design is powered, so you can ignore it
  input  wire       clk,      // clock
  input  wire       rst_n     // reset_n - low to reset
);

  // VGA signals
  wire hsync;
  wire vsync;
  wire [1:0] R;
  wire [1:0] G;
  wire [1:0] B;
  wire video_active;
  wire [9:0] pix_x;
  wire [9:0] pix_y;

  // TinyVGA PMOD
  assign uo_out = {hsync, B[0], G[0], R[0], vsync, B[1], G[1], R[1]};

  // Unused outputs assigned to 0.
  assign uio_out = 0;
  assign uio_oe  = 0;

  // Suppress unused signals warning
  wire _unused_ok = &{ena, ui_in, uio_in};

  hvsync_generator hvsync_gen(
    .clk(clk),
    .reset(~rst_n),
    .hsync(hsync),
    .vsync(vsync),
    .display_on(video_active),
    .hpos(pix_x),
    .vpos(pix_y)
  );
  
  parameter C = 64;
  parameter W = 10;
  reg [C-1:0] state [W];
  
  wire init_row = video_active & (pix_y==0);
  wire init_val = pix_x==320;
  wire serial = init_row ? init_val : (
     (state[0][0] ^ state[0][1] )
  );

  integer i;
  always @(posedge clk) begin
    if (video_active) begin
      for (i=0; i<W-1; ++i) begin
        state[i] <= {state[i+1][0], state[i][C-1:1]};
      end
      state[W-1] <= {serial, state[W-1][C-1:1]};
    end
  end

  
  wire [1:0] a = {serial,serial};

  assign R = video_active ? a : 2'b00;
  assign G = video_active ? a : 2'b00;
  assign B = video_active ? a : 2'b00;
  
endmodule
