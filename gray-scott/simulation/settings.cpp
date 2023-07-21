#ifndef __SETTINGS_H__
#define __SETTINGS_H__

#include <string>

struct Settings
{
    int nx;
    int ny;
    int nz;
    int direc;
    int steps;
    int plotgap;
    int LID;
    int WALL;
    int FLUID;
    std::string output;
    std::string adios_config;

    Settings();
    static Settings from_json(const std::string &fname);
};

#endif
