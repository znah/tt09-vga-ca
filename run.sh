# -Wall -O2
CFLAGS=`sdl2-config --cflags`
LDFLAGS=`sdl2-config --libs`

clear
SRC="main.cpp src/project.v src/hvsync_generator.v"
verilator --cc --exe --build -j 0 -O2 $SRC -o demo -CFLAGS $CFLAGS -DSIM && obj_dir/demo