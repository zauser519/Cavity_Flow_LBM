#ifndef __WRITER_H__
#define __WRITER_H__

#include <mpi.h>

#include "../../gray-scott/simulation/gray-scott.h"
#include "../../gray-scott/simulation/settings.h"

class Writer
{
public:
    Writer(const Settings &settings, const GrayScott &sim);
    void open(const std::string &fname, bool append);
    void write(int step, const GrayScott &sim);
    void close();

protected:
    Settings settings;

    double var_u;
    double var_v;
    int var_step;
};

#endif
