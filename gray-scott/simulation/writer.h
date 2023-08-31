#ifndef __WRITER_H__
#define __WRITER_H__

#include <adios2.h>
#include <mpi.h>

#include "../../gray-scott/simulation/gray-scott.h"
#include "../../gray-scott/simulation/settings.h"

#define LID 2
#define WALL 1
#define FLUID 0

#define nx 100
#define ny 100
#define nz 100
#define direc 19

class Writer
{
public:
    Writer(const Settings &settings, adios2::IO io);
    void open(const std::string &fname);
    void write(double u[nz][ny][nx], double v[nz][ny][nx], double w[nz][ny][nx], double rho[nz][ny][nx], int step);
    void close();

protected:
    Settings settings;

    adios2::IO io;
    adios2::Engine writer;
    adios2::Variable<double> var_u;
    adios2::Variable<double> var_v;
    adios2::Variable<double> var_w;
    adios2::Variable<double> var_rho;
    adios2::Variable<int> var_step;
};

#endif
