/*
 * Copyright (c) 2024 Alexander Mordvintsev
 * SPDX-License-Identifier: Apache-2.0
 */

`default_nettype none

module tt_um_znah_vga_ca(
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
  
  parameter logCELL_SIZE = 2;
  parameter CELL_SIZE = 1<<logCELL_SIZE;
  parameter WIDTH = 640;
  parameter HEIGHT = 480;
  parameter GRID_W = 160;
  parameter PAD_LEFT = (WIDTH-GRID_W*CELL_SIZE)/2;
  
  wire [9:0] x = pix_x-PAD_LEFT;
  wire [7:0] cell_x = x[9:logCELL_SIZE];

  wire [logCELL_SIZE-1:0] fract_x = x[logCELL_SIZE-1:0];
  wire [logCELL_SIZE-1:0] fract_y = pix_y[logCELL_SIZE-1:0];
  
  // parameter L = GRID_W/4-1;
  // `define REG(name) (* mem2reg *) reg [L:0] name [4]
  // `define SHIFT(data) data[3] <= {data[3][L-1:0], data[2][L]}; \
  //                     data[2] <= {data[2][L-1:0], data[1][L]}; \
  //                     data[1] <= {data[1][L-1:0], data[0][L]}; \
  //                     data[0][L:1] <= data[0][L-1:0];
  // `define HEAD(data) data[0][0]
  // `define TAIL(data,i) data[3][L-(i)]
  parameter L = GRID_W-1;
  `ifdef SIM
    `define BUF(name) assign name``_buf = name
  `else
    `define BUF(name) sky130_fd_sc_hd__dlygate4sd3_1 name``buf_[L:0] ( .X(name``_buf), .A(name) )
  `endif
  `define REG(name) (* mem2reg *) reg[L:0] name; wire[L:0] name``_buf ; `BUF(name)

  `define SHIFT(data) data[L:1] <= data``_buf[L-1:0]
  `define HEAD(data) data[0]
  `define TAIL(data,i) data``_buf[L-(i)]

  `REG(cells);
  `REG(next_cells);
  reg left;
  wire center = `TAIL(cells, 0);
  wire right = `TAIL(cells, 1);

  reg [10:0] row_count;
  wire [2:0] i = row_count[10:8];
  reg [7:0] rules [0:7];
  initial begin
        rules[0] = 30;
        rules[1] = 110;
        rules[2] = 22;
        rules[3] = 73;
        rules[4] = 90;
        rules[5] = 146;
        rules[6] = 105;
        rules[7] = 102;
  end  

  wire [7:0] rule = rules[i];
  wire [5:0] rule_color = rule[6:1];
  
  wire seed_cell = cell_x == GRID_W/2;
  wire first_row_cell = row_count==0 ? seed_cell : `TAIL(next_cells, 0);
  wire rule_cell = fract_y==0 ? rule[{left,center,right}] : center;
  wire new_cell = pix_y==0 ? first_row_cell : rule_cell;

  wire in_grid = cell_x < GRID_W && video_active;
  wire row_end = pix_x == WIDTH;
  wire reset = ~rst_n;
  always @(posedge clk) begin
    if (reset) begin
      row_count <= 0;
      `HEAD(cells) <= 0;
    end else if (row_end && pix_y<HEIGHT && &fract_y) begin
      row_count <= row_count+1;
    end else if (row_end && pix_y==HEIGHT) begin
      row_count <= row_count-HEIGHT/CELL_SIZE+1;
    end

    if (~reset && in_grid && fract_x==0) begin
      left <= `TAIL(cells, 0);
      `SHIFT(cells);
      `HEAD(cells) <= new_cell;
      if (pix_y == 0) begin
        `SHIFT(next_cells);
      end else if (pix_y == CELL_SIZE) begin
        `SHIFT(next_cells);
        `HEAD(next_cells) <= new_cell;
      end
    end

  end

  wire c = `HEAD(cells) & in_grid;
  wire [5:0] color = c ? rule_color : 6'b000000;

  assign R = color[5:4];
  assign G = color[3:2];
  assign B = color[1:0];
endmodule