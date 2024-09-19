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

  reg [GRID_W-1:0] cells;
  reg [GRID_W-1:0] next_cells;
  reg left;
  wire center = cells[GRID_W-1];
  wire right = cells[GRID_W-2];

  parameter RULE_BIT = 11;

  reg [RULE_BIT:0] row_count;
  wire in_grid = cell_x < GRID_W && video_active;
  wire copy_row = (pix_y&(CELL_SIZE-1)) != 0;
  wire rule30 = left ^ (center | right);  // 30
  wire rule110 = ((left|center) ^ (left&center&right)); // 110
  wire rule_sel = row_count[RULE_BIT];
  wire rule_cell = rule_sel ? rule110 : rule30;
  wire new_cell = copy_row ? center : rule_cell;

  reg init = 1;
  always @(negedge step) begin
    if (!rst_n) begin
      init <= 1;
    end
    if (in_grid) begin
      left <= cells[GRID_W-1];
      cells[GRID_W-1:1] <= cells[GRID_W-2:0];
      if (pix_y == 0) begin
        if (init) begin
          cells[0] <= cell_x == GRID_W/2; // seed
         end else begin
          cells[0] <= next_cells[GRID_W-1];
          next_cells <= {next_cells[GRID_W-2:0], next_cells[GRID_W-1]};
         end
      end else begin
        cells[0] <= new_cell;
      end
      if (pix_y == CELL_SIZE) begin
        next_cells <= {next_cells[GRID_W-2:0], new_cell};
        init <= 0;
      end
    end
    if (cell_x == 0 && video_active) begin // rule switching
      row_count <= row_count+1+(pix_y==HEIGHT-1 ? CELL_SIZE-HEIGHT : 0);
    end
  end

  wire c = cells[0]&in_grid;
  wire [5:0] color = c?(rule_sel?6'b001011:6'b101100):6'b000000;
  assign R = color[5:4];
  assign G = color[3:2];
  assign B = color[1:0];
endmodule
