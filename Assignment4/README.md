The command below executes faster than make.
icpc -diag-disable=10441 -O3 -Wall -Werror -o main flux.cc tendency.cc main.cc
tar -cvf hw4.tar flux.cc flux.h tendency.cc tendency.h main.cc plot_h.ncl