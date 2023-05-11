#include "../../gray-scott/simulation/writer.h"
 
#include <string>
#include <iostream>
#include <fstream>

void define_bpvtk_attribute(const Settings &s)
{
    auto lf_VTKImage = [](const Settings &s) {
        const std::string extent = "0 " + std::to_string(s.L) + " " + "0 " +
                                   std::to_string(s.L) + " " + "0 " +
                                   std::to_string(s.L);

        const std::string imageData = R"(
        <?xml version="1.0"?>
        <VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">
          <ImageData WholeExtent=")" + extent +
                                      R"(" Origin="0 0 0" Spacing="1 1 1">
            <Piece Extent=")" + extent +
                                      R"(">
              <CellData Scalars="U">
                  <DataArray Name="U" />
                  <DataArray Name="V" />
                  <DataArray Name="TIME">
                    step
                  </DataArray>
              </CellData>
            </Piece>
          </ImageData>
        </VTKFile>)";

        //@@@@@@@ 
        //io.DefineAttribute<std::string>("vtk.xml", imageData);
        // NONE BEING USED
        //struct imageData
        //{
        //    std::string vtk.xml;
        //};
    };

    if (s.mesh_type == "image")
    {
        lf_VTKImage(s);
    }
    else if (s.mesh_type == "structured")
    {
        throw std::invalid_argument(
            "ERROR: mesh_type=structured not yet "
            "   supported in settings.json, use mesh_type=image instead\n");
    }
}

Writer::Writer(const Settings &settings, const GrayScott &sim)
: settings(settings)
{
    //@@@@@@@@@@@@@@@ 
    //io.DefineAttribute<double>("F", settings.F);
    //io.DefineAttribute<double>("k", settings.k);
    //io.DefineAttribute<double>("dt", settings.dt);
    //io.DefineAttribute<double>("Du", settings.Du);
    //io.DefineAttribute<double>("Dv", settings.Dv);
    //io.DefineAttribute<double>("noise", settings.noise);
    // NONE BEING USED
    //int n();
    //struct F
    //{
    //    double n() {return settings.F;}
    //};
    //struct k
    //{
    //    double n=settings.k;
    //};   
    //struct dt
    //{
    //    double n=settings.dt;
    //};
    //struct Du
    //{
    //    double n=settings.Du;
    //};
    //struct Dv
    //{
    //    double n=settings.Dv;
    //};
    //struct noise
    //{
    //    double n=settings.noise;
    //};

    // define VTK visualization schema as an attribute
    //if (!settings.mesh_type.empty())
    //{
    //    define_bpvtk_attribute(settings);
    //}

    //@@@@@@@@@
    //var_u =
    //    io.DefineVariable<double>("U", {settings.L, settings.L, settings.L},
    //                              {sim.offset_z, sim.offset_y, sim.offset_x},
    //                              {sim.size_z, sim.size_y, sim.size_x});

    //var_v =
    //    io.DefineVariable<double>("V", {settings.L, settings.L, settings.L},
    //                              {sim.offset_z, sim.offset_y, sim.offset_x},
    //                              {sim.size_z, sim.size_y, sim.size_x});

    //struct var_u
    //{
    //    std::string U, 
    //    int L={settings.L, settings.L, settings.L},
    //    double offset={sim.offset_z, sim.offset_y, sim.offset_x},
    //    double size={sim.size_z, sim.size_y, sim.size_x};
    //};

    //struct var_v{
    //    std::string V;
    //    int L={settings.L, settings.L, settings.L};
    //    double offset={sim.offset_z, sim.offset_y, sim.offset_x};
    //    double size={sim.size_z, sim.size_y, sim.size_x};
    //};

    //@@@@@@@@
    //    if (settings.adios_memory_selection)
    //{
    //    var_u.SetMemorySelection(
    //        {{1, 1, 1}, {sim.size_z + 2, sim.size_y + 2, sim.size_x + 2}});
    //    var_v.SetMemorySelection(
    //        {{1, 1, 1}, {sim.size_z + 2, sim.size_y + 2, sim.size_x + 2}});
    //}
    
    // ??
    //if (settings.adios_memory_selection)
    //{
    //    var_u.offset={1, 1, 1};
    //    var_u.size={sim.size_z + 2, sim.size_y + 2, sim.size_x + 2};
    //    var_v.offset={1, 1, 1};
    //    var_v.size={sim.size_z + 2, sim.size_y + 2, sim.size_x + 2};
    //}
    //@@@@@@@
    //var_step = io.DefineVariable<int>("step");
    //int var_step=step;
}

void Writer::open(const std::string &fname, bool append)
{
    //@@@@@
    //adios2::Mode mode = adios2::Mode::Write;
    //if (append)
    //{
    //    mode = adios2::Mode::Append;
    //}
    //writer = io.Open(fname, mode);

    std::fstream file;
    file.open(fname,std::ios::ate);
    if (append)
    {
        file.open(fname,std::ios::app);
    }
}

void Writer::write(int step, const GrayScott &sim)
{
    if (!sim.size_x || !sim.size_y || !sim.size_z)
    {
        //@@@@@
        //writer.BeginStep();
        //writer.EndStep();

        
        return;
    }

    if (settings.adios_memory_selection)
    {
        const std::vector<double> &u = sim.u_ghost();
        const std::vector<double> &v = sim.v_ghost();

        //@@@@@@@@@@@@
        //writer.Put<int>(var_step, &step);
        //writer.Put<double>(var_u, u.data());
        //writer.Put<double>(var_v, v.data());

    }
    else if (settings.adios_span)
    {
        //@@@@@@@
        //
        //writer.Put<int>(var_step, &step);
        //
        // provide memory directly from adios buffer
        //std::span u_span = writer.Put<double>(var_u);
        //std::span v_span = writer.Put<double>(var_v);

        // populate spans
        //sim.u_noghost(u_span.data());
        //sim.v_noghost(v_span.data());
        
    }
    else
    {
        std::vector<double> u = sim.u_noghost();
        std::vector<double> v = sim.v_noghost();

        //@@@@@@
        //writer.Put<int>(var_step, &step);
        //writer.Put<double>(var_u, u.data());
        //writer.Put<double>(var_v, v.data());
    }
}

void Writer::close() 
{ 
    //@@@@@@@
    //writer.Close(); 
    std::fstream file;
    file.close();
}
