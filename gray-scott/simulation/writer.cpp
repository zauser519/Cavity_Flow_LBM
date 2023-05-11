#include "../../gray-scott/simulation/writer.h"
 
#include <string>

#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>

Writer::Writer(const Settings &settings, const GrayScott &sim, adios2::IO io)
: settings(settings), io(io)
{
    var_u =
        io.DefineVariable<double>("U", {settings.L, settings.L, settings.L},
                                  {sim.offset_z, sim.offset_y, sim.offset_x},
                                  {sim.size_z, sim.size_y, sim.size_x});

    var_v =
        io.DefineVariable<double>("V", {settings.L, settings.L, settings.L},
                                  {sim.offset_z, sim.offset_y, sim.offset_x},
                                  {sim.size_z, sim.size_y, sim.size_x});
    var_step = io.DefineVariable<int>("step");
}

void Writer::Wopen(const std::string &fname, bool append)
{
    const char *path="/home/gp.sc.cc.tohoku.ac.jp/tseng/ADIOS2/Test/Tutorial/VE/share/adios2-examples/gray-scott/gs.bp/data.0";
    adios2::Mode mode = adios2::Mode::Write;
    if (append)
    {
        mode = adios2::Mode::Append;
    }
    writer = io.Open(fname, mode);
}

void Writer::Wwrite(int step, const GrayScott &sim, int fd)
{
    if (!sim.size_x || !sim.size_y || !sim.size_z)
    {
        return;
    }

    std::vector<double> u = sim.u_noghost();
    std::vector<double> v = sim.v_noghost();

    //if (!append)
    //{
        //@@@@@
        lseek(fd, step * (sizeof(u.size() * sizeof(double) + sizeof(v.size() * sizeof(double) + sizeof(int)))), SEEK_SET);
        write(fd, &step, sizeof(int));
        write(fd, u.data(), sizeof(u.size() * sizeof(double)));
        write(fd, v.data(), sizeof(v.size() * sizeof(double)));
    //}
    //else
    //{
        //@@@@@
        //write(open("gs.bp", O_APPEND | O_WRONLY | O_TRUNC, 0644), &step, sizeof(int));
        //write(open("gs.bp", O_APPEND | O_WRONLY | O_TRUNC, 0644), &var_u, sizeof(u.data()));
        //write(open("gs.bp", O_APPEND | O_WRONLY | O_TRUNC, 0644), &var_v, sizeof(v.data()));
    //}

}

void Writer::Wclose(int fd) 
{ 
    //@@@@@@@
    //writer.Close(); 
    //if (append)
    //{
        //@@@@@
        close(fd);
    //}
    //else
    //{
        //@@@@@
    //    close(open("gs.bp", O_APPEND | O_WRONLY, 0644));
    //}
    
}
