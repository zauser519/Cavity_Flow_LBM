// The solver is based on Hiroshi Watanabe's 2D Gray-Scott reaction diffusion
// code available at:
// https://github.com/kaityo256/sevendayshpc/tree/master/day5

#include "../../gray-scott/simulation/gray-scott.h"

#include <algorithm>
#include <mpi.h>
#include <random>
#include <stdexcept> // runtime_error
#include <vector>
