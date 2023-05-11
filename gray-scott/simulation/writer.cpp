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
     //Open thre files as original ADIOS2 did, include data.0 as simulated result, md.0 and md.idx as adios2 setting files
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

     //Write file into data.0
     //pointer
     lseek(fd, step * (sizeof(u.size() * sizeof(double) + sizeof(v.size() * sizeof(double) + sizeof(int)))), SEEK_SET);
     //Write simulation data
     write(fd, &step, sizeof(int));
     write(fd, u.data(), sizeof(u.size() * sizeof(double)));
     write(fd, v.data(), sizeof(v.size() * sizeof(double)));

}

void Writer::Wclose(int fd) 
{ 
    //Close file (ADIOS2's function, cause rewrite data.0)
    //writer.Close();
    //POSIX close file
    close(fd);
    
}
