
//------------------------------------------------------------------------------
// cavity lbm3d for VE       by sxc    2022 11 30

// 2022 12 15  modify

// 2023 2 14  modify

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <omp.h>
#include <ftrace.h> // For VE??

#define LID 2
#define WALL 1
#define FLUID 0

#define nx 100
#define ny 100
#define nz 100
#define direc 19

#define NUM_THREADS 8 //

#ifndef max
#define max(x, y) (((x) <= (y)) ? (y) : (x))
#endif

#ifndef min
#define min(x, y) (((x) >= (y)) ? (y) : (x))
#endif

#define abs(x) (((x) < 0) ? (-x) : (x)) /* NOTICE: when using expression a+b, abs((a+b)) should be used instead of abs(a+b) */

double cpu_time();

void init_geo(int wall[nz][ny][nx]);

void init(int wall[nz][ny][nx], double u[nz][ny][nx], double v[nz][ny][nx], double w[nz][ny][nx], double rho[nz][ny][nx], double f[direc][nz][ny][nx], double ft[direc][nz][ny][nx],
          double rho_0, double u_0);

void stream(double f[direc][nz][ny][nx], double ft[direc][nz][ny][nx]);

void collision_bc(int wall[nz][ny][nx], double u[nz][ny][nx], double v[nz][ny][nx], double w[nz][ny][nx], double rho[nz][ny][nx], double f[direc][nz][ny][nx], double ft[direc][nz][ny][nx], double rho_0, double u_0,
                  double ow, double ow1);

void data_output(int wall[nz][ny][nx], double u[nz][ny][nx], double v[nz][ny][nx], double w[nz][ny][nx], double rho[nz][ny][nx], int atime);

extern void omp_set_num_threads(int num_threads); // Set number of threads
extern int omp_get_num_threads(void);             // Get number of threads
extern int omp_get_max_threads(void);             // Get upper bounds on number of threads
extern int omp_get_thread_num(void);              // Get thread number

// Beginning of user-specified region
extern int ftrace_region_begin(const char *id);
// End of user-specified region
extern int ftrace_region_end(const char *id);

int main(int argc, char *argv[])
{

  // define variables
  int atime, time_max, output;

  double Re, rho_0, u_0, L_ref, visc, tau, omega, ow, ow1;

  int wall[nz][ny][nx];

  double ft[direc][nz][ny][nx];

  double f[direc][nz][ny][nx], u[nz][ny][nx], v[nz][ny][nx], w[nz][ny][nx], rho[nz][ny][nx];

  //--------------------------------------------------------------------------------------------------

  atime = 0;
  time_max = 300000;
  output = 300000;
  Re = 100;
  rho_0 = 1.0;
  u_0 = 0.1;

  L_ref = 1.0 * (nx - 1);         /* 1.0 is the length of lattice */
  visc = u_0 * L_ref / Re;        /* viscosity derived from Re */
  tau = (6.0 * visc + 1.0) / 2.0; /* tau derived from viscosity */
  omega = 1.0 / tau;
  ow1 = (1.0f - omega);
  ow = omega;

  //-------------------------------------------------------------------------------------------------------
  // geo define

  init_geo(wall);
  //-------------------------------------------------------------------------------------------------------

  // initialization

  init(wall, u, v, w, rho, f, ft, rho_0, u_0);

  //-------------------------------------------------------------------------------------------------------
  // loop
  // stream  collision   BC

  for (atime = 0; atime <= time_max; atime++)
  {

    stream(f, ft);

    collision_bc(wall, u, v, w, rho, f, ft, rho_0, u_0, ow, ow1);

    if (atime % output == 0)
    {
      data_output(wall, u, v, w, rho, atime);
    }

  } //   --------  end loop ---------------------

  //-------------------------------------------------------------------------------------------------------

  return 0;
} /***** end main *****/

void init_geo(int wall[nz][ny][nx])
{

  for (int z = 0; z < nz; z++)
  {
    for (int y = 0; y < ny; y++)
    {
      for (int x = 0; x < nx; x++)
      {

        if (x == nx - 1 || x == 0 || y == 0 || y == ny - 1 || z == 0) // solid
        {
          wall[z][y][x] = 1;
        }
        else if (z == nz - 1) // top
        {
          wall[z][y][x] = 2;
        }
        else // fluid
        {
          wall[z][y][x] = 0;
        }
      }
    }
  }
}

void init(int wall[nz][ny][nx], double u[nz][ny][nx], double v[nz][ny][nx], double w[nz][ny][nx], double rho[nz][ny][nx], double f[direc][nz][ny][nx], double ft[direc][nz][ny][nx],
          double rho_0, double u_0)
{

  double uu, vv, ww, rhorho, square;

  for (int z = 0; z < nz; z++)
  {
    for (int y = 0; y < ny; y++)
    {
      for (int x = 0; x < nx; x++)
      {

        rho[z][y][x] = rho_0;
        u[z][y][x] = 0.0;
        v[z][y][x] = 0.0;
        w[z][y][x] = 0.0;

        if (wall[z][y][x] == 2)
        {
          u[z][y][x] = u_0; /* the top is moving lid with u_0 */
        }
      }
    }
  }

  for (int z = 0; z < nz; z++)
  {
    for (int y = 0; y < ny; y++)
    {
      for (int x = 0; x < nx; x++)
      {
        uu = u[z][y][x];
        vv = v[z][y][x];
        ww = w[z][y][x];

        rhorho = rho[z][y][x];
        square = uu * uu + vv * vv + ww * ww;

        f[0][z][y][x] = rhorho / 3. * (1. - 1.5 * square);
        f[1][z][y][x] = rhorho / 18. * (1. + 3. * uu + 4.5 * uu * uu - 1.5 * square);
        f[2][z][y][x] = rhorho / 18. * (1. - 3. * uu + 4.5 * uu * uu - 1.5 * square);
        f[3][z][y][x] = rhorho / 18. * (1. + 3. * vv + 4.5 * vv * vv - 1.5 * square);
        f[4][z][y][x] = rhorho / 18. * (1. - 3. * vv + 4.5 * vv * vv - 1.5 * square);
        f[5][z][y][x] = rhorho / 18. * (1. + 3. * ww + 4.5 * ww * ww - 1.5 * square);
        f[6][z][y][x] = rhorho / 18. * (1. - 3. * ww + 4.5 * ww * ww - 1.5 * square);
        f[7][z][y][x] = rhorho / 36. * (1. + 3. * (uu + vv) + 4.5 * (uu + vv) * (uu + vv) - 1.5 * square);
        f[8][z][y][x] = rhorho / 36. * (1. - 3. * (uu - vv) + 4.5 * (uu - vv) * (uu - vv) - 1.5 * square);
        f[9][z][y][x] = rhorho / 36. * (1. + 3. * (uu - vv) + 4.5 * (uu - vv) * (uu - vv) - 1.5 * square);
        f[10][z][y][x] = rhorho / 36. * (1. - 3. * (uu + vv) + 4.5 * (uu + vv) * (uu + vv) - 1.5 * square);
        f[11][z][y][x] = rhorho / 36. * (1. + 3. * (uu + ww) + 4.5 * (uu + ww) * (uu + ww) - 1.5 * square);
        f[12][z][y][x] = rhorho / 36. * (1. - 3. * (uu - ww) + 4.5 * (uu - ww) * (uu - ww) - 1.5 * square);
        f[13][z][y][x] = rhorho / 36. * (1. + 3. * (uu - ww) + 4.5 * (uu - ww) * (uu - ww) - 1.5 * square);
        f[14][z][y][x] = rhorho / 36. * (1. - 3. * (uu + ww) + 4.5 * (uu + ww) * (uu + ww) - 1.5 * square);
        f[15][z][y][x] = rhorho / 36. * (1. + 3. * (vv + ww) + 4.5 * (vv + ww) * (vv + ww) - 1.5 * square);
        f[16][z][y][x] = rhorho / 36. * (1. - 3. * (vv - ww) + 4.5 * (vv - ww) * (vv - ww) - 1.5 * square);
        f[17][z][y][x] = rhorho / 36. * (1. + 3. * (vv - ww) + 4.5 * (vv - ww) * (vv - ww) - 1.5 * square);
        f[18][z][y][x] = rhorho / 36. * (1. - 3. * (vv + ww) + 4.5 * (vv + ww) * (vv + ww) - 1.5 * square);
      }
    }
  }

  for (int i = 0; i < direc; i++)
  {
    for (int z = 0; z < nz; z++)
    {
      for (int y = 0; y < ny; y++)
      {
        for (int x = 0; x < nx; x++)
        {
          ft[i][z][y][x] = f[i][z][y][x];
        }
      }
    }
  }

} /***** end init *****/

//-------  stream ----------------------------------
void stream(double f[direc][nz][ny][nx], double ft[direc][nz][ny][nx])
{
  int x, y, z;
  int jxm, jxp, jym, jyp, jzm, jzp; // 这几个值是共享变量
                                    // int nxyz;
  // int temp;
  omp_set_num_threads(NUM_THREADS);

#ifdef FTRACE
  (void)ftrace_region_begin("stream_all");
#endif /* FTRACE */

#pragma omp parallel for private(jxp, jxm, jyp, jym, jzm, jzp)

  for (int z = 0; z < nz; z++) // unvectorized
  {
    for (int y = 0; y < ny; y++) // unvectorized
    {
      for (int x = 0; x < nx; x++)
      {

        jxp = (x + 1) % nx;
        jxm = (x - 1 + nx) % nx;
        jyp = (y + 1) % ny;
        jym = (y - 1 + ny) % ny;
        jzp = (z + 1) % nz;
        jzm = (z - 1 + nz) % nz;

        //  stream
        ft[0][z][y][x] = f[0][z][y][x];
        ft[1][z][y][x] = f[1][z][y][jxm];
        ft[2][z][y][x] = f[2][z][y][jxp];
        ft[3][z][y][x] = f[3][z][jym][x];
        ft[4][z][y][x] = f[4][z][jyp][x];
        ft[5][z][y][x] = f[5][jzm][y][x];
        ft[6][z][y][x] = f[6][jzp][y][x];
        ft[7][z][y][x] = f[7][z][jym][jxm];
        ft[8][z][y][x] = f[8][z][jym][jxp];
        ft[9][z][y][x] = f[9][z][jyp][jxm];
        ft[10][z][y][x] = f[10][z][jyp][jxp];
        ft[11][z][y][x] = f[11][jzm][y][jxm];
        ft[12][z][y][x] = f[12][jzm][y][jxp];
        ft[13][z][y][x] = f[13][jzp][y][jxm];
        ft[14][z][y][x] = f[14][jzp][y][jxp];
        ft[15][z][y][x] = f[15][jzm][jym][x];
        ft[16][z][y][x] = f[16][jzm][jyp][x];
        ft[17][z][y][x] = f[17][jzp][jym][x];
        ft[18][z][y][x] = f[18][jzp][jyp][x];
      }
    }
  }
#ifdef FTRACE
  (void)ftrace_region_end("stream_all");
#endif /* FTRACE */

  // 三维cavity的streaming步骤优化太复杂了，暂时先用之前简约的，后续确有需要，再简化。
}

//-------  collision and BC ----------------------------------
void collision_bc(int wall[nz][ny][nx], double u[nz][ny][nx], double v[nz][ny][nx], double w[nz][ny][nx], double rho[nz][ny][nx], double f[direc][nz][ny][nx], double ft[direc][nz][ny][nx],
                  double rho_0, double u_0, double ow, double ow1)
{

  double uu, vv, ww, rhorho, square, nuu, nvv, nww, nrhorho, nsquare;
  double us, vs, ws, rhos, rhow, US;

  int nyz;

  omp_set_num_threads(NUM_THREADS);

#ifdef FTRACE
  (void)ftrace_region_begin("collision_BC_top");
#endif /* FTRACE */
#pragma omp parallel for private(uu, vv, ww, rhorho, square, nuu, nvv, nww, nrhorho, nsquare)

  for (int y = 0; y < ny; y++)
  {
    for (int x = 0; x < nx; x++) //  vectorized loop.
    {

      uu = u[(nz - 1)][y][x];
      vv = v[(nz - 1)][y][x];
      ww = w[(nz - 1)][y][x];

      rhorho = rho[(nz - 1)][y][x];
      square = uu * uu + vv * vv + ww * ww;

      nuu = u[(nz - 2)][y][x];
      nvv = v[(nz - 2)][y][x];
      nww = w[(nz - 2)][y][x];

      nrhorho = rho[(nz - 2)][y][x];
      nsquare = nuu * nuu + nvv * nvv + nww * nww;

      f[0][nz - 1][y][x] = rhorho / 3. * (1. - 1.5 * square) + (f[0][nz - 2][y][x] - nrhorho / 3. * (1. - 1.5 * nsquare));

      f[1][nz - 1][y][x] = rhorho / 18. * (1. + 3. * uu + 4.5 * uu * uu - 1.5 * square) + (f[1][nz - 2][y][x] - nrhorho / 18. * (1. + 3. * nuu + 4.5 * nuu * nuu - 1.5 * nsquare));
      f[2][nz - 1][y][x] = rhorho / 18. * (1. - 3. * uu + 4.5 * uu * uu - 1.5 * square) + (f[2][nz - 2][y][x] - nrhorho / 18. * (1. - 3. * nuu + 4.5 * nuu * nuu - 1.5 * nsquare));
      f[3][nz - 1][y][x] = rhorho / 18. * (1. + 3. * vv + 4.5 * vv * vv - 1.5 * square) + (f[3][nz - 2][y][x] - nrhorho / 18. * (1. + 3. * nvv + 4.5 * nvv * nvv - 1.5 * nsquare));
      f[4][nz - 1][y][x] = rhorho / 18. * (1. - 3. * vv + 4.5 * vv * vv - 1.5 * square) + (f[4][nz - 2][y][x] - nrhorho / 18. * (1. - 3. * nvv + 4.5 * nvv * nvv - 1.5 * nsquare));
      f[5][nz - 1][y][x] = rhorho / 18. * (1. + 3. * ww + 4.5 * ww * ww - 1.5 * square) + (f[5][nz - 2][y][x] - nrhorho / 18. * (1. + 3. * nww + 4.5 * nww * nww - 1.5 * nsquare));
      f[6][nz - 1][y][x] = rhorho / 18. * (1. - 3. * ww + 4.5 * ww * ww - 1.5 * square) + (f[6][nz - 2][y][x] - nrhorho / 18. * (1. - 3. * nww + 4.5 * nww * nww - 1.5 * nsquare));

      f[7][nz - 1][y][x] = rhorho / 36. * (1. + 3. * (uu + vv) + 4.5 * (uu + vv) * (uu + vv) - 1.5 * square) + (f[7][nz - 2][y][x] - nrhorho / 36. * (1. + 3. * (nuu + nvv) + 4.5 * (nuu + nvv) * (nuu + nvv) - 1.5 * nsquare));
      f[8][nz - 1][y][x] = rhorho / 36. * (1. - 3. * (uu - vv) + 4.5 * (uu - vv) * (uu - vv) - 1.5 * square) + (f[8][nz - 2][y][x] - nrhorho / 36. * (1. - 3. * (nuu - nvv) + 4.5 * (nuu - nvv) * (nuu - nvv) - 1.5 * nsquare));
      f[9][nz - 1][y][x] = rhorho / 36. * (1. + 3. * (uu - vv) + 4.5 * (uu - vv) * (uu - vv) - 1.5 * square) + (f[9][nz - 2][y][x] - nrhorho / 36. * (1. + 3. * (nuu - nvv) + 4.5 * (nuu - nvv) * (nuu - nvv) - 1.5 * nsquare));
      f[10][nz - 1][y][x] = rhorho / 36. * (1. - 3. * (uu + vv) + 4.5 * (uu + vv) * (uu + vv) - 1.5 * square) + (f[10][nz - 2][y][x] - nrhorho / 36. * (1. - 3. * (nuu + nvv) + 4.5 * (nuu + nvv) * (nuu + nvv) - 1.5 * nsquare));

      f[11][nz - 1][y][x] = rhorho / 36. * (1. + 3. * (uu + ww) + 4.5 * (uu + ww) * (uu + ww) - 1.5 * square) + (f[11][nz - 2][y][x] - nrhorho / 36. * (1. + 3. * (nuu + nww) + 4.5 * (nuu + nww) * (nuu + nww) - 1.5 * nsquare));
      f[12][nz - 1][y][x] = rhorho / 36. * (1. - 3. * (uu - ww) + 4.5 * (uu - ww) * (uu - ww) - 1.5 * square) + (f[12][nz - 2][y][x] - nrhorho / 36. * (1. - 3. * (nuu - nww) + 4.5 * (nuu - nww) * (nuu - nww) - 1.5 * nsquare));
      f[13][nz - 1][y][x] = rhorho / 36. * (1. + 3. * (uu - ww) + 4.5 * (uu - ww) * (uu - ww) - 1.5 * square) + (f[13][nz - 2][y][x] - nrhorho / 36. * (1. + 3. * (nuu - nww) + 4.5 * (nuu - nww) * (nuu - nww) - 1.5 * nsquare));
      f[14][nz - 1][y][x] = rhorho / 36. * (1. - 3. * (uu + ww) + 4.5 * (uu + ww) * (uu + ww) - 1.5 * square) + (f[14][nz - 2][y][x] - nrhorho / 36. * (1. - 3. * (nuu + nww) + 4.5 * (nuu + nww) * (nuu + nww) - 1.5 * nsquare));

      f[15][nz - 1][y][x] = rhorho / 36. * (1. + 3. * (vv + ww) + 4.5 * (vv + ww) * (vv + ww) - 1.5 * square) + (f[15][nz - 2][y][x] - nrhorho / 36. * (1. + 3. * (nvv + nww) + 4.5 * (nvv + nww) * (nvv + nww) - 1.5 * nsquare));
      f[16][nz - 1][y][x] = rhorho / 36. * (1. - 3. * (vv - ww) + 4.5 * (vv - ww) * (vv - ww) - 1.5 * square) + (f[16][nz - 2][y][x] - nrhorho / 36. * (1. - 3. * (nvv - nww) + 4.5 * (nvv - nww) * (nvv - nww) - 1.5 * nsquare));
      f[17][nz - 1][y][x] = rhorho / 36. * (1. + 3. * (vv - ww) + 4.5 * (vv - ww) * (vv - ww) - 1.5 * square) + (f[17][nz - 2][y][x] - nrhorho / 36. * (1. + 3. * (nvv - nww) + 4.5 * (nvv - nww) * (nvv - nww) - 1.5 * nsquare));
      f[18][nz - 1][y][x] = rhorho / 36. * (1. - 3. * (vv + ww) + 4.5 * (vv + ww) * (vv + ww) - 1.5 * square) + (f[18][nz - 2][y][x] - nrhorho / 36. * (1. - 3. * (nvv + nww) + 4.5 * (nvv + nww) * (nvv + nww) - 1.5 * nsquare));

      rho[(nz - 1)][y][x] = rho[(nz - 2)][y][x];
      // u[(nz - 1) * ny * nx + y * nx + x] = u_0;
      // v[(nz - 1) * ny * nx + y * nx + x] = 0.0;
      // w[(nz - 1) * ny * nx + y * nx + x] = 0.0;
    }
  }

#ifdef FTRACE
  (void)ftrace_region_end("collision_BC_top");
#endif /* FTRACE */

  //====================================================================================================

#ifdef FTRACE
  (void)ftrace_region_begin("collision_BC_left");
#endif /* FTRACE */

#pragma omp parallel for
  for (int z = 0; z < nz - 1; z++) // vectorized loop.
  {
    for (int y = 0; y < ny; y++) // vectorized loop.
    {

      // rho[nyz] = rho_0;
      // u[nyz] = 0.0;
      // v[nyz] = 0.0;
      // w[nyz] = 0.0;

      f[0][z][y][0] = ft[0][z][y][0];
      f[1][z][y][0] = ft[2][z][y][0];
      f[2][z][y][0] = ft[1][z][y][0];
      f[3][z][y][0] = ft[4][z][y][0];
      f[4][z][y][0] = ft[3][z][y][0];
      f[5][z][y][0] = ft[6][z][y][0];
      f[6][z][y][0] = ft[5][z][y][0];
      f[7][z][y][0] = ft[10][z][y][0];
      f[8][z][y][0] = ft[9][z][y][0];
      f[9][z][y][0] = ft[8][z][y][0];
      f[10][z][y][0] = ft[7][z][y][0];
      f[11][z][y][0] = ft[14][z][y][0];
      f[12][z][y][0] = ft[13][z][y][0];
      f[13][z][y][0] = ft[12][z][y][0];
      f[14][z][y][0] = ft[11][z][y][0];
      f[15][z][y][0] = ft[18][z][y][0];
      f[16][z][y][0] = ft[17][z][y][0];
      f[17][z][y][0] = ft[16][z][y][0];
      f[18][z][y][0] = ft[15][z][y][0];
    }
  }
#ifdef FTRACE
  (void)ftrace_region_end("collision_BC_left");
#endif /* FTRACE */

  //====================================================================================================

#ifdef FTRACE
  (void)ftrace_region_begin("collision_BC_right");
#endif /* FTRACE */

#pragma omp parallel for //  right solid

  for (int z = 0; z < nz - 1; z++) // vectorized loop.
  {
    for (int y = 0; y < ny; y++) // vectorized loop.
    {
      //  nyz = nx * (ny * z + y) + nx - 1;

      // rho[nyz] = rho_0;
      // u[nyz] = 0.0;
      // v[nyz] = 0.0;
      // w[nyz] = 0.0;

      f[0][z][y][nx - 1] = ft[0][z][y][nx - 1];
      f[1][z][y][nx - 1] = ft[2][z][y][nx - 1];
      f[2][z][y][nx - 1] = ft[1][z][y][nx - 1];
      f[3][z][y][nx - 1] = ft[4][z][y][nx - 1];
      f[4][z][y][nx - 1] = ft[3][z][y][nx - 1];
      f[5][z][y][nx - 1] = ft[6][z][y][nx - 1];
      f[6][z][y][nx - 1] = ft[5][z][y][nx - 1];
      f[7][z][y][nx - 1] = ft[10][z][y][nx - 1];
      f[8][z][y][nx - 1] = ft[9][z][y][nx - 1];
      f[9][z][y][nx - 1] = ft[8][z][y][nx - 1];
      f[10][z][y][nx - 1] = ft[7][z][y][nx - 1];
      f[11][z][y][nx - 1] = ft[14][z][y][nx - 1];
      f[12][z][y][nx - 1] = ft[13][z][y][nx - 1];
      f[13][z][y][nx - 1] = ft[12][z][y][nx - 1];
      f[14][z][y][nx - 1] = ft[11][z][y][nx - 1];
      f[15][z][y][nx - 1] = ft[18][z][y][nx - 1];
      f[16][z][y][nx - 1] = ft[17][z][y][nx - 1];
      f[17][z][y][nx - 1] = ft[16][z][y][nx - 1];
      f[18][z][y][nx - 1] = ft[15][z][y][nx - 1];
    }
  }
#ifdef FTRACE
  (void)ftrace_region_end("collision_BC_right");
#endif /* FTRACE */

  //====================================================================================================

#ifdef FTRACE
  (void)ftrace_region_begin("collision_BC_bottom");
#endif /* FTRACE */

#pragma omp parallel for
  for (int y = 0; y < ny; y++) // vectorized loop.
  {
    for (int x = 0; x < nx; x++) // vectorized loop.
    {
      //  rho[nx * ny * 0 + nx * y + x] = rho_0;
      //  u[nx * ny * 0 + nx * y + x] = 0.0;
      //  v[nx * ny * 0 + nx * y + x] = 0.0;
      //  w[nx * ny * 0 + nx * y + x] = 0.0;

      f[0][0][y][x] = ft[0][0][y][x];
      f[1][0][y][x] = ft[2][0][y][x];
      f[2][0][y][x] = ft[1][0][y][x];
      f[3][0][y][x] = ft[4][0][y][x];
      f[4][0][y][x] = ft[3][0][y][x];
      f[5][0][y][x] = ft[6][0][y][x];
      f[6][0][y][x] = ft[5][0][y][x];
      f[7][0][y][x] = ft[10][0][y][x];
      f[8][0][y][x] = ft[9][0][y][x];
      f[9][0][y][x] = ft[8][0][y][x];
      f[10][0][y][x] = ft[7][0][y][x];
      f[11][0][y][x] = ft[14][0][y][x];
      f[12][0][y][x] = ft[13][0][y][x];
      f[13][0][y][x] = ft[12][0][y][x];
      f[14][0][y][x] = ft[11][0][y][x];
      f[15][0][y][x] = ft[18][0][y][x];
      f[16][0][y][x] = ft[17][0][y][x];
      f[17][0][y][x] = ft[16][0][y][x];
      f[18][0][y][x] = ft[15][0][y][x];
    }
  }

#ifdef FTRACE
  (void)ftrace_region_end("collision_BC_bottom");
#endif /* FTRACE */

  //====================================================================================================

#ifdef FTRACE
  (void)ftrace_region_begin("collision_BC_front");
#endif /* FTRACE */

#pragma omp parallel for
  for (int z = 0; z < nz - 1; z++) // vectorized loop.
  {
    for (int x = 0; x < nx; x++) // vectorized loop.
    {
      //  rho[nx * ny * z + nx * 0 + x] = rho_0;
      //  u[nx * ny * z + nx * 0 + x] = 0.0;
      //  v[nx * ny * z + nx * 0 + x] = 0.0;
      //   w[nx * ny * z + nx * 0 + x] = 0.0;

      f[0][z][0][x] = ft[0][z][0][x];
      f[1][z][0][x] = ft[2][z][0][x];
      f[2][z][0][x] = ft[1][z][0][x];
      f[3][z][0][x] = ft[4][z][0][x];
      f[4][z][0][x] = ft[3][z][0][x];
      f[5][z][0][x] = ft[6][z][0][x];
      f[6][z][0][x] = ft[5][z][0][x];
      f[7][z][0][x] = ft[10][z][0][x];
      f[8][z][0][x] = ft[9][z][0][x];
      f[9][z][0][x] = ft[8][z][0][x];
      f[10][z][0][x] = ft[7][z][0][x];
      f[11][z][0][x] = ft[14][z][0][x];
      f[12][z][0][x] = ft[13][z][0][x];
      f[13][z][0][x] = ft[12][z][0][x];
      f[14][z][0][x] = ft[11][z][0][x];
      f[15][z][0][x] = ft[18][z][0][x];
      f[16][z][0][x] = ft[17][z][0][x];
      f[17][z][0][x] = ft[16][z][0][x];
      f[18][z][0][x] = ft[15][z][0][x];
    }
  }

#ifdef FTRACE
  (void)ftrace_region_end("collision_BC_front");
#endif /* FTRACE */

  //====================================================================================================

#ifdef FTRACE
  (void)ftrace_region_begin("collision_BC_back");
#endif /* FTRACE */

#pragma omp parallel for
  for (int z = 0; z < nz - 1; z++) // vectorized loop.
  {
    for (int x = 0; x < nx; x++) // vectorized loop.
    {
      //   rho[nx * ny * z + nx * (ny - 1) + x] = rho_0;
      //    u[nx * ny * z + nx * (ny - 1) + x] = 0.0;
      //   v[nx * ny * z + nx * (ny - 1) + x] = 0.0;
      //     w[nx * ny * z + nx * (ny - 1) + x] = 0.0;

      f[0][z][ny - 1][x] = ft[0][z][ny - 1][x];
      f[1][z][ny - 1][x] = ft[2][z][ny - 1][x];
      f[2][z][ny - 1][x] = ft[1][z][ny - 1][x];
      f[3][z][ny - 1][x] = ft[4][z][ny - 1][x];
      f[4][z][ny - 1][x] = ft[3][z][ny - 1][x];
      f[5][z][ny - 1][x] = ft[6][z][ny - 1][x];
      f[6][z][ny - 1][x] = ft[5][z][ny - 1][x];
      f[7][z][ny - 1][x] = ft[10][z][ny - 1][x];
      f[8][z][ny - 1][x] = ft[9][z][ny - 1][x];
      f[9][z][ny - 1][x] = ft[8][z][ny - 1][x];
      f[10][z][ny - 1][x] = ft[7][z][ny - 1][x];
      f[11][z][ny - 1][x] = ft[14][z][ny - 1][x];
      f[12][z][ny - 1][x] = ft[13][z][ny - 1][x];
      f[13][z][ny - 1][x] = ft[12][z][ny - 1][x];
      f[14][z][ny - 1][x] = ft[11][z][ny - 1][x];
      f[15][z][ny - 1][x] = ft[18][z][ny - 1][x];
      f[16][z][ny - 1][x] = ft[17][z][ny - 1][x];
      f[17][z][ny - 1][x] = ft[16][z][ny - 1][x];
      f[18][z][ny - 1][x] = ft[15][z][ny - 1][x];
    }
  }

#ifdef FTRACE
  (void)ftrace_region_end("collision_BC_back");
#endif /* FTRACE */
//----------------------------------------------------------
#ifdef FTRACE
  (void)ftrace_region_begin("collision_BC_fluid");
#endif /* FTRACE */
#pragma omp parallel for private(us, vs, ws, rhos, rhow, US)

  for (int z = 1; z < nz - 1; z++)
  {
    for (int y = 1; y < ny - 1; y++) // vectorized loop.
    {
      for (int x = 1; x < nx - 1; x++) // vectorized loop.
      {

        rho[z][y][x] = ft[0][z][y][x] + ft[1][z][y][x] +
                       ft[2][z][y][x] + ft[3][z][y][x] + ft[4][z][y][x] + ft[5][z][y][x] + ft[6][z][y][x] +
                       ft[7][z][y][x] + ft[8][z][y][x] + ft[9][z][y][x] + ft[10][z][y][x] + ft[11][z][y][x] +
                       ft[12][z][y][x] + ft[13][z][y][x] + ft[14][z][y][x] + ft[15][z][y][x] + ft[16][z][y][x] +
                       ft[17][z][y][x] + ft[18][z][y][x];

        u[z][y][x] = (ft[1][z][y][x] + ft[7][z][y][x] + ft[9][z][y][x] + ft[11][z][y][x] + ft[13][z][y][x] - ft[2][z][y][x] - ft[8][z][y][x] - ft[10][z][y][x] - ft[12][z][y][x] - ft[14][z][y][x]) / rho[z][y][x];

        v[z][y][x] = (ft[3][z][y][x] + ft[7][z][y][x] + ft[8][z][y][x] + ft[15][z][y][x] + ft[17][z][y][x] - ft[4][z][y][x] - ft[9][z][y][x] - ft[10][z][y][x] - ft[16][z][y][x] - ft[18][z][y][x]) / rho[z][y][x];

        w[z][y][x] = (ft[5][z][y][x] + ft[11][z][y][x] + ft[12][z][y][x] + ft[15][z][y][x] + ft[16][z][y][x] - ft[6][z][y][x] - ft[13][z][y][x] - ft[14][z][y][x] - ft[17][z][y][x] - ft[18][z][y][x]) / rho[z][y][x];

        us = u[z][y][x];
        vs = v[z][y][x];
        ws = w[z][y][x];
        rhos = rho[z][y][x];

        rhow = ow * rhos;
        US = -1.5f * (us * us + vs * vs + ws * ws);

        f[0][z][y][x] = ow1 * ft[0][z][y][x] + rhow / 3.0f * (1.0f + US);
        f[1][z][y][x] = ow1 * ft[1][z][y][x] + rhow / 18.0f * (1.0f + 3.0f * us + 4.5f * us * us + US);
        f[2][z][y][x] = ow1 * ft[2][z][y][x] + rhow / 18.0f * (1.0f - 3.0f * us + 4.5f * us * us + US);
        f[3][z][y][x] = ow1 * ft[3][z][y][x] + rhow / 18.0f * (1.0f + 3.0f * vs + 4.5f * vs * vs + US);
        f[4][z][y][x] = ow1 * ft[4][z][y][x] + rhow / 18.0f * (1.0f - 3.0f * vs + 4.5f * vs * vs + US);
        f[5][z][y][x] = ow1 * ft[5][z][y][x] + rhow / 18.0f * (1.0f + 3.0f * ws + 4.5f * ws * ws + US);
        f[6][z][y][x] = ow1 * ft[6][z][y][x] + rhow / 18.0f * (1.0f - 3.0f * ws + 4.5f * ws * ws + US);

        f[7][z][y][x] = ow1 * ft[7][z][y][x] + rhow / 36.0f * (1.0f + 3.0f * (us + vs) + 4.5f * (us + vs) * (us + vs) + US);
        f[8][z][y][x] = ow1 * ft[8][z][y][x] + rhow / 36.0f * (1.0f - 3.0f * (us - vs) + 4.5f * (us - vs) * (us - vs) + US);
        f[9][z][y][x] = ow1 * ft[9][z][y][x] + rhow / 36.0f * (1.0f + 3.0f * (us - vs) + 4.5f * (us - vs) * (us - vs) + US);
        f[10][z][y][x] = ow1 * ft[10][z][y][x] + rhow / 36.0f * (1.0f - 3.0f * (us + vs) + 4.5f * (us + vs) * (us + vs) + US);

        f[11][z][y][x] = ow1 * ft[11][z][y][x] + rhow / 36.0f * (1.0f + 3.0f * (us + ws) + 4.5f * (us + ws) * (us + ws) + US);
        f[12][z][y][x] = ow1 * ft[12][z][y][x] + rhow / 36.0f * (1.0f - 3.0f * (us - ws) + 4.5f * (us - ws) * (us - ws) + US);
        f[13][z][y][x] = ow1 * ft[13][z][y][x] + rhow / 36.0f * (1.0f + 3.0f * (us - ws) + 4.5f * (us - ws) * (us - ws) + US);
        f[14][z][y][x] = ow1 * ft[14][z][y][x] + rhow / 36.0f * (1.0f - 3.0f * (us + ws) + 4.5f * (us + ws) * (us + ws) + US);

        f[15][z][y][x] = ow1 * ft[15][z][y][x] + rhow / 36.0f * (1.0f + 3.0f * (vs + ws) + 4.5f * (vs + ws) * (vs + ws) + US);
        f[16][z][y][x] = ow1 * ft[16][z][y][x] + rhow / 36.0f * (1.0f - 3.0f * (vs - ws) + 4.5f * (vs - ws) * (vs - ws) + US);
        f[17][z][y][x] = ow1 * ft[17][z][y][x] + rhow / 36.0f * (1.0f + 3.0f * (vs - ws) + 4.5f * (vs - ws) * (vs - ws) + US);
        f[18][z][y][x] = ow1 * ft[18][z][y][x] + rhow / 36.0f * (1.0f - 3.0f * (vs + ws) + 4.5f * (vs + ws) * (vs + ws) + US);
      }
    }
  }

#ifdef FTRACE
  (void)ftrace_region_end("collision_BC_fluid");
#endif /* FTRACE */

} /***** end collision *****/

void data_output(int wall[nz][ny][nx], double u[nz][ny][nx], double v[nz][ny][nx], double w[nz][ny][nx], double rho[nz][ny][nx], int atime)
{
  int kk1, kk2, kk3, kk4, kk5, kk6, kk7;
  int kk8, kk9, kk10, kk11, kk12, kk13, kk14, kk15, kk16;
  int x, y, z;
  char data_name[200];

  FILE *file_dat;

  kk1 = atime;
  kk2 = kk1 / 10000000;
  kk3 = kk1 - kk2 * 10000000;
  kk4 = kk3 / 1000000;
  kk5 = kk3 - kk4 * 1000000;
  kk6 = kk5 / 100000;
  kk7 = kk5 - kk6 * 100000;
  kk8 = kk7 / 10000;
  kk9 = kk7 - kk8 * 10000;
  kk10 = kk9 / 1000;
  kk11 = kk9 - kk10 * 1000;
  kk12 = kk11 / 100;
  kk13 = kk11 - kk12 * 100;
  kk14 = kk13 / 10;
  kk15 = kk13 - kk14 * 10;
  kk16 = kk15 / 1;

  sprintf(data_name, "%d%d%d%d%d%d%d%d.txt", kk2, kk4, kk6, kk8, kk10, kk12, kk14, kk16);

  file_dat = fopen(data_name, "w");

  fputs("VARIABLES=X,Y,Z,RHO,U,V,W\n", file_dat);
  fprintf(file_dat, "ZONE I=%6d, J=%6d, K=%6d, F=POINT\n", nx, ny, nz);

  for (int z = 0; z < nz; z++)
  {
    for (int y = 0; y < ny; y++)
    {
      for (int x = 0; x < nx; x++)
      {

        fprintf(file_dat, "%lf %lf  %lf  %20.6e %20.6e %20.6e %20.6e \n",
                (double)x, (double)y, (double)z, rho[z][y][x], u[z][y][x], v[z][y][x], w[z][y][x]);
      }
    }
  }

  fclose(file_dat); /**end of ordinary flie output **/
}

double cpu_time()
{

  struct timeval tm;
  double t;
  static int base_sec = 0, base_usec = 0;

  gettimeofday(&tm, NULL);

  if (base_sec == 0 && base_usec == 0)
  {
    base_sec = tm.tv_sec;
    base_usec = tm.tv_usec;
    t = 0.0;
  }
  else
  {
    t = (double)(tm.tv_sec - base_sec) +
        ((double)(tm.tv_usec - base_usec)) / 1.0e6;
  }

  return t;

} /***** end cpu_time *****/
