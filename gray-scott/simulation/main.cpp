#include <fstream>
#include <iostream>
#include <sstream>
#include <fcntl.h>
#include <unistd.h>
#include <vector>

#include <adios2.h>
#include <mpi.h>
#include <cstdlib>

#include "../../gray-scott/common/timer.hpp"
#include "../../gray-scott/simulation/gray-scott.h"
#include "../../gray-scott/simulation/restart.h"
#include "../../gray-scott/simulation/writer.h"


void print_io_settings(const adios2::IO &io)
{
    std::cout << "Simulation writes data using engine type:              "
              << io.EngineType() << std::endl;
    auto ioparams = io.Parameters();
    std::cout << "IO parameters:  " << std::endl;
    for (const auto &p : ioparams)
    {
        std::cout << "    " << p.first << " = " << p.second << std::endl;
    }
}

void print_settings(const Settings &s, int restart_step)
{
    std::cout << "grid:             " << s.L << "x" << s.L << "x" << s.L
              << std::endl;
    if (restart_step > 0)
    {
        std::cout << "restart:          from step " << restart_step
                  << std::endl;
    }
    else
    {
        std::cout << "restart:          no" << std::endl;
    }
    std::cout << "steps:            " << s.steps << std::endl;
    std::cout << "plotgap:          " << s.plotgap << std::endl;
    std::cout << "F:                " << s.F << std::endl;
    std::cout << "k:                " << s.k << std::endl;
    std::cout << "dt:               " << s.dt << std::endl;
    std::cout << "Du:               " << s.Du << std::endl;
    std::cout << "Dv:               " << s.Dv << std::endl;
    std::cout << "noise:            " << s.noise << std::endl;
    std::cout << "output:           " << s.output << std::endl;
    std::cout << "adios_config:     " << s.adios_config << std::endl;
}

void print_simulator_settings(const GrayScott &s)
{
    std::cout << "process layout:   " << s.npx << "x" << s.npy << "x" << s.npz
              << std::endl;
    std::cout << "local grid size:  " << s.size_x << "x" << s.size_y << "x"
              << s.size_z << std::endl;
}


int main(int argc, char **argv)
{
    int provided, fd;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    int rank, procs, wrank;

    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    const unsigned int color = 1;
    MPI_Comm comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &comm);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &procs);

    ////offset is here
    //Calculate the chunk size
    int n=100;
    char *envval = NULL;

    if (argc < 2)
    {
        if (rank == 0)
        {
            std::cerr << "Too few arguments" << std::endl;
            std::cerr << "Usage: gray-scott settings.json" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    //Total data size
    if (NULL != (envval = getenv("DATASIZE")))
    {
        n = atoi(envval);
        if (n < 100)
            n= 100;
    };


    Settings settings = Settings::from_json(argv[1]);

    GrayScott sim(settings, comm);
    sim.init();

    adios2::ADIOS adios(settings.adios_config, comm);
    adios2::IO io_main = adios.DeclareIO("SimulationOutput");
    adios2::IO io_ckpt = adios.DeclareIO("SimulationCheckpoint");

    int restart_step = 0;

    //@@@@@@@@@@
    Writer writer_main(settings, sim, io_main);
    writer_main.Wopen(settings.output, (restart_step > 0));
    
    if (rank == 0)
    {
        //ftruncate(fd,(sizeof(sim.u_noghost().data())+sizeof(sim.v_noghost().data())+sizeof(int))*n);
    
#ifdef DEBUG
    printf("Total data size: %d\n", n);
#endif
        print_io_settings(io_main);
        std::cout << "========================================" << std::endl;
        print_settings(settings, restart_step);
        print_simulator_settings(sim);
        std::cout << "========================================" << std::endl;
    }

#ifdef ENABLE_TIMERS
    Timer timer_total;
    Timer timer_compute;
    Timer timer_write;

    std::ostringstream log_fname;
    log_fname << "gray_scott_pe_" << rank << ".log";

    std::ofstream log(log_fname.str());
    log << "step\ttotal_gs\tcompute_gs\twrite_gs" << std::endl;
#endif
    
    
    fd = open("/home/gp.sc.cc.tohoku.ac.jp/tseng/ADIOS2/Test/Tutorial/VE/share/adios2-examples/gray-scott/gs.bp/data.0", O_CREAT | O_WRONLY, 0644);
    for (int it = restart_step; it < settings.steps;)
    {
#ifdef ENABLE_TIMERS
        MPI_Barrier(comm);
        timer_total.start();
        timer_compute.start();
#endif

        sim.iterate();
        it++;

#ifdef ENABLE_TIMERS
        timer_compute.stop();
        MPI_Barrier(comm);
        timer_write.start();
#endif

        //Write every plotgap
        if (it % settings.plotgap == 0)
        {
            if (rank == 0)
            {
                std::cout << "Simulation at step " << it
                          << " writing output step     "
                          << it / settings.plotgap << std::endl;
            }

            //Write is here
            writer_main.Wwrite(it, sim, fd);
        }


#ifdef ENABLE_TIMERS
        double time_write = timer_write.stop();
        double time_step = timer_total.stop();
        MPI_Barrier(comm);

        log << it << "\t" << timer_total.elapsed() << "\t"
            << timer_compute.elapsed() << "\t" << timer_write.elapsed()
            << std::endl;
#endif
    }

    //Close the file
    writer_main.Wclose(fd);

#ifdef ENABLE_TIMERS
    log << "total\t" << timer_total.elapsed() << "\t" << timer_compute.elapsed()
        << "\t" << timer_write.elapsed() << std::endl;

    log.close();
#endif

    MPI_Finalize();
}
