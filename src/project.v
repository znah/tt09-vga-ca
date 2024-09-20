/*
 * Copyright (c) 2024 Alexander Mordvintsev
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
  
  parameter WIDTH = 640;
  parameter HEIGHT = 480;
  parameter GRID_W = 100;
  parameter logCELL_SIZE = 2;
  parameter CELL_SIZE = 1<<logCELL_SIZE;
  parameter PAD_LEFT = (WIDTH-GRID_W*CELL_SIZE)/2;
  
  wire [9:0] x = pix_x-PAD_LEFT;
  wire [7:0] cell_x = x[9:logCELL_SIZE];
  wire step = x[logCELL_SIZE-1];

  parameter L = GRID_W/2-1;
  `define REG(name) reg [L:0] name [2];
  `define SHIFT(data) data[1] <= {data[1][L-1:0], data[0][L]}; data[0][L:1] <= data[0][L-1:0];
  `define HEAD(data) data[0][0]
  `define TAIL(data,i) data[1][L-(i)]
  `REG(cells);
  `REG(next_cells);
  reg left;
  wire center = `TAIL(cells, 0);
  wire right = `TAIL(cells, 1);

  parameter [5:0] rule_colors [16] = '{6'b001011, 6'b101100, 6'b101001, 6'b001110,  6'b101010, 6'b101110, 6'b110001, 6'b100110, 
                                      6'b101010, 6'b101110, 6'b110001, 6'b100110,  6'b101010, 6'b101110, 6'b110001, 6'b100110};
  parameter [7:0] rules[16] = '{30, 18, 161, 150, 22, 110, 54, 122,     118, 57, 18, 165, 90, 180, 60, 146};
  reg [10:0] row_count;
  wire [3:0] rule_idx = row_count[10:7];
  wire [7:0] rule = rules[rule_idx];
  wire rule_cell = rule[{left,center,right}];
  
  wire copy_row = (pix_y&(CELL_SIZE-1)) != 0;  
  wire new_cell = copy_row ? center : rule_cell;

  wire in_grid = cell_x < GRID_W && video_active;
  reg init = 1;
  always @(negedge step) begin
    if (!rst_n) begin
      init <= 1;
    end
    if (in_grid) begin
      left <= `TAIL(cells, 0);
      `SHIFT(cells);
      if (pix_y == 0) begin
        if (init) begin
          `HEAD(cells) <= cell_x == GRID_W/2; // seed
         end else begin
          `HEAD(cells) <= `TAIL(next_cells, 0);
          `SHIFT(next_cells);
         end
      end else begin
        `HEAD(cells) <= new_cell;
      end
      if (pix_y == CELL_SIZE) begin
        `SHIFT(next_cells);
        `HEAD(next_cells) <= new_cell;
        init <= 0;
      end
    end
    if (cell_x == 0 && pix_y<=HEIGHT && pix_y%CELL_SIZE==0) begin // rule switching
      row_count <= row_count+(pix_y==HEIGHT ? 1-HEIGHT/CELL_SIZE : 1);
    end
  end

  wire c = `HEAD(cells)&in_grid;
  wire [5:0] color = c ? rule_colors[rule_idx] : 6'b000000;
  assign R = color[5:4];
  assign G = color[3:2];
  assign B = color[1:0];
endmodule
