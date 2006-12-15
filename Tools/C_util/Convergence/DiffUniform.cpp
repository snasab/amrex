
//
// $Id: DiffUniform.cpp,v 1.9 2006-12-15 00:55:42 almgren Exp $
//

#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
using std::ios;
using std::set_new_handler;

#include <unistd.h>

#include "WritePlotFile.H"
#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"
#include "VisMF.H"

#ifndef NDEBUG
#include "TV_TempWrite.H"
#endif

#include "AVGDOWN_F.H"

#define GARBAGE 666.e+40

static
void
PrintUsage (const char* progName)
{
    std::cout << '\n';
    std::cout << "This utility computes the difference of two plotfiles " << std::endl
         << "which are at resolutions differing by a power of 2.   " << std::endl
         << "The finer solution is averaged down to the coarse grid" << std::endl
         << "resolution and then the difference is taken.          " << std::endl
         << std::endl;
    std::cout << "Usage:" << '\n';
    std::cout << progName << '\n';
    std::cout << "    infile  = inputFileName" << '\n';
    std::cout << "    exact   = exactFileName" << '\n';
    std::cout << "    outfile = outputFileName" << '\n';
    std::cout << "   [-help]" << '\n';
    std::cout << "   [-verbose]" << '\n';
    std::cout << '\n';
    exit(1);
}

int
finestLevelCoveringDomain(const AmrData& amrData);

IntVect
getRefRatio(const Box& crse,
	    const Box& fine);

bool
amrDatasHaveSameDerives(const AmrData& amrd1,
			const AmrData& amrd2);
int
main (int   argc,
      char* argv[])
{
    if (argc == 1)
        PrintUsage(argv[0]);

    BoxLib::Initialize(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
        PrintUsage(argv[0]);

    FArrayBox::setFormat(FABio::FAB_IEEE_32);
    //
    // Scan the arguments.
    //
    std::string iFileDir, iFile, eFile, oFile, oFileDir;

    bool verbose = false;
    if (pp.contains("verbose"))
    {
        verbose = true;
        AmrData::SetVerbose(true);
    }
    pp.query("infile", iFile);
    if (iFile.empty())
        BoxLib::Abort("You must specify `infile'");

    pp.query("exact", eFile);
    if (eFile.empty())
        BoxLib::Abort("You must specify `exact' file");

    pp.query("outfile", oFile);
    if (oFile.empty())
        BoxLib::Abort("You must specify `outfile'");

    DataServices::SetBatchMode();
    FileType fileType(NEWPLT);
    
    DataServices dataServicesC(iFile, fileType);
    DataServices dataServicesF(eFile, fileType);

    if (!dataServicesC.AmrDataOk() || !dataServicesF.AmrDataOk())
        //
        // This calls ParallelDescriptor::EndParallel() and exit()
        //
        DataServices::Dispatch(DataServices::ExitRequest, NULL);

    AmrData& amrDataC = dataServicesC.AmrDataRef();
    AmrData& amrDataF = dataServicesF.AmrDataRef();

    BL_ASSERT(amrDatasHaveSameDerives(amrDataC,amrDataF));
    int exact_level = finestLevelCoveringDomain(amrDataF);
    if (exact_level < 0)
    {
	std::cout << "Exact data does not contain a level covering the domain" << '\n';
	BoxLib::Abort();
    }
    if (verbose)
	std::cout << "Using level = " << exact_level << " in 'exact' file" << '\n';
	
    int finestLevel = amrDataC.FinestLevel();
    
    //
    // Compute the error
    //
    Array<MultiFab*> error(finestLevel+1);
    
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
    {
        const BoxArray& crseBA = amrDataC.boxArray(iLevel);
	int nComp              = amrDataC.NComp();
	const Box& domainC     = amrDataC.ProbDomain()[iLevel];
	const Box& domainF     = amrDataF.ProbDomain()[exact_level];
	IntVect refine_ratio   = getRefRatio(domainC, domainF);
	if (refine_ratio == IntVect())
	    BoxLib::Error("Cannot find refinement ratio from data to exact");
        std::cout << "Ratio for level " << iLevel << " is " << refine_ratio << std::endl;

	error[iLevel] = new MultiFab(crseBA, nComp, 0);
	error[iLevel]->setVal(GARBAGE);

	for (int iComp=0; iComp<nComp; ++iComp)
	{
//          MultiFab& exact = amrDataF.GetGrids(finestLevel,iComp);
            MultiFab& exact = amrDataF.GetGrids(exact_level,iComp);
            const BoxArray& exactBA = exact.boxArray();
            const BoxArray crseBA = ::BoxArray(exactBA).coarsen(refine_ratio);
            MultiFab aveExact(crseBA,1,0,Fab_allocate);
            int nc = exact.nComp();
            for (MFIter amfi(aveExact); amfi.isValid(); ++amfi)
            {
                const Box& crseBox = amfi.validbox();
                const int* a_lo = aveExact[amfi].loVect();
                const int* a_hi = aveExact[amfi].hiVect();
                const int* e_lo = exact[amfi].loVect();
                const int* e_hi = exact[amfi].hiVect();
                FORT_CV_AVGDOWN(aveExact[amfi].dataPtr(),
                                ARLIM(a_lo), ARLIM(a_hi), &nc,
                                exact[amfi].dataPtr(),
                                ARLIM(e_lo), ARLIM(e_hi),
                                crseBox.loVect(), crseBox.hiVect(),
                                refine_ratio.getVect());
            }

            // Copy result of coarsening into error as temporary storage
            error[iLevel]->copy(aveExact,0,iComp,nc);

            // Subtract coarse data from coarsened exact data
            MultiFab& data = amrDataC.GetGrids(iLevel,iComp);
            BL_ASSERT(data.boxArray() == error[iLevel]->boxArray());
            for (MFIter dmfi(data); dmfi.isValid(); ++dmfi)
            {
                FArrayBox err((*error[iLevel])[dmfi]);
                err.minus(data[dmfi],0,iComp,1);
            }
        }
    }

    WritePlotFile(error, amrDataC, oFile, verbose);
    
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
	delete error[iLevel];

    //
    // This calls ParallelDescriptor::EndParallel() and exit()
    //
    DataServices::Dispatch(DataServices::ExitRequest, NULL);
}


int
finestLevelCoveringDomain(const AmrData& amr_data)
{
    // Find the finest level covering the entire domain.  Return
    // -1 if there isn't one suitable
    int finest_level = amr_data.FinestLevel();
    const Array<Box>& domain_array = amr_data.ProbDomain();

    for (int iLevel=finest_level; iLevel>=0; --iLevel)
    {
	const BoxArray& ba = amr_data.boxArray(iLevel);
	BoxDomain bd;
	bd.add(BoxList(ba));
	BoxDomain complement = BoxLib::complementIn(domain_array[iLevel],bd);
	if (complement.isEmpty())
	    return iLevel;
    }
    return -1;
}

IntVect
getRefRatio(const Box& crse,
	    const Box& fine)
{
    // Compute refinement ratio between crse and fine boxes, return invalid
    // IntVect if there is none suitable
    ParmParse pp("");
    Array<int> rr_in(BL_SPACEDIM,-1);
    int Nrr = pp.countval("ref_ratio",Nrr);
    BL_ASSERT(Nrr==0 || Nrr==BL_SPACEDIM || Nrr==1);
    if (Nrr>0) 
    {
        pp.queryarr("ref_ratio",rr_in,0,Nrr);
        return IntVect(rr_in);
    }

    IntVect ref_ratio;
    for (int i=0; i<BL_SPACEDIM; ++i)
	ref_ratio[i] = fine.size()[i] / crse.size()[i];

    // Check results
    Box test1 = ::Box(fine).coarsen(ref_ratio);
    Box test2 = ::Box(test1).refine(ref_ratio);
    if (test1 != crse  ||  test2 != fine)
	ref_ratio = IntVect();
    return ref_ratio;
}

bool
amrDatasHaveSameDerives(const AmrData& amrd1,
			const AmrData& amrd2)
{
    const Array<std::string>& derives1 = amrd1.PlotVarNames();
    const Array<std::string>& derives2 = amrd2.PlotVarNames();
    int length = derives1.size();
    if (length != derives2.size())
	return false;
    for (int i=0; i<length; ++i)
	if (derives1[i] != derives2[i])
	    return false;
    return true;
}
    
