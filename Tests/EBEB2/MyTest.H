#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <string>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

class MyTest
{
public:
    MyTest ();
    ~MyTest ();

    void test ();

private:
    void readParameters ();
    void initGrids ();
    void initializeEB ();
    void initializeEB2 ();
    void initData ();

    int n_cell = 64;
    int max_grid_size = 32;
    int max_coarsening_level = 5;

    amrex::Geometry geom;
    amrex::BoxArray grids;
    amrex::DistributionMapping dmap;
//    std::unique_ptr<amrex::EBFArrayBoxFactory> factory;
};

#endif