#include "../../gray-scott/simulation/writer.h"

Writer::Writer(const Settings &settings, adios2::IO io)
: settings(settings), io(io)
{

    var_u = io.DefineVariable<double>("U", {nz, ny, nx}, {0, 0, 0}, {nz, ny, nx});
    var_v = io.DefineVariable<double>("V", {nz, ny, nx}, {0, 0, 0}, {nz, ny, nx});
    var_w = io.DefineVariable<double>("W", {nz, ny, nx}, {0, 0, 0}, {nz, ny, nx});
    var_rho = io.DefineVariable<double>("rho", {nz, ny, nx}, {0, 0, 0}, {nz, ny, nx});

    var_step = io.DefineVariable<int>("step");
}

void Writer::open(const std::string &fname)
{
    adios2::Mode mode = adios2::Mode::Write;
    writer = io.Open(fname, mode);
}

void Writer::write(double u[nz][ny][nx], double v[nz][ny][nx], double w[nz][ny][nx], double rho[nz][ny][nx], int step)
{

    writer.BeginStep();
    writer.Put<int>(var_step, &step);
    writer.Put<double>(var_u, (double*)u);
    writer.Put<double>(var_v, (double*)v);
    writer.Put<double>(var_w, (double*)w);
    writer.Put<double>(var_rho, (double*)rho);
    writer.EndStep();

}

void Writer::close() { writer.Close(); }
