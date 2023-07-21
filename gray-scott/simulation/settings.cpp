#include "../../gray-scott/simulation/settings.h"

#include <fstream>

#include "../../gray-scott/simulation/json.hpp"

void to_json(nlohmann::json &j, const Settings &s)
{
    j = nlohmann::json{{"nx", s.nx},
                       {"ny", s.ny},
                       {"nz", s.nz},
                       {"direc", s.direc},
                       {"steps", s.steps},
                       {"plotgap", s.plotgap},
                       {"LID", s.LID},
                       {"WALL", s.WALL},
                       {"FLUID", s.FLUID},
                       {"output", s.output},
                       {"adios_config", s.adios_config}};
}

void from_json(const nlohmann::json &j, Settings &s)
{
    j.at("nx").get_to(s.nx);
    j.at("ny").get_to(s.ny);
    j.at("nz").get_to(s.nz);
    j.at("direc").get_to(s.direc);
    j.at("steps").get_to(s.steps);
    j.at("plotgap").get_to(s.plotgap);
    j.at("LID").get_to(s.LID);
    j.at("WALL").get_to(s.WALL);
    j.at("FLUID").get_to(s.FLUID);
    j.at("output").get_to(s.output);
    j.at("adios_config").get_to(s.adios_config);
}

Settings::Settings()
{
    nx = 100;
    ny = 100;
    nz = 100;
    direc = 19;
    steps = 20000;
    plotgap = 200;
    LID = 2;
    WALL = 1;
    FLUID = 0;
    output = "foo.bp";
    adios_config = "adios2.xml";
}

Settings Settings::from_json(const std::string &fname)
{
    std::ifstream ifs(fname);
    nlohmann::json j;

    ifs >> j;

    return j.get<Settings>();
}
