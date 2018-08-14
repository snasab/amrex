#ifndef _AMREX_PARAMS_H_
#define _AMREX_PARAMS_H_

#include <AMReX.H>
#include <AMReX_Particles.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

struct TestParams {
  int max_step;
  Real dt;
  int nx;
  int ny;
  int nz;
  int max_grid_size;
  int nppc;
  bool verbose;
};


#endif // _AMREX_PARAMS_H_
