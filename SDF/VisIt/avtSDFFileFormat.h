/*****************************************************************************
 *
 * SDF (Self-Describing Format) VisIt reader
 * Copyright (c) 2010-2016, SDF Development Team
 *
 * Distributed under the terms of the BSD 3-clause License.
 * See the LICENSE file for details.
 *
 ****************************************************************************/

// ************************************************************************* //
//                            avtSDFFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_SDF_FILE_FORMAT_H
#define AVT_SDF_FILE_FORMAT_H

#include <avtSTMDFileFormat.h>

#include <vector>
#include "sdf.h"
#include "sdf_derived.h"
#include "stack_allocator.h"

class DBOptionsAttributes;
using std::vector;

// ****************************************************************************
//  Class: avtSDFFileFormat
//
//  Purpose:
//      Reads in SDF files as a plugin to VisIt.
//
//  Programmer: Keith Bennett
//  Creation:   Fri Oct 29 15:31:09 PST 2010
//
// ****************************************************************************

class avtSDFFileFormat : public avtSTMDFileFormat
{
  public:
                       avtSDFFileFormat(const char *, DBOptionsAttributes *);
    virtual           ~avtSDFFileFormat();

    //
    // This is used to return unconvention data -- ranging from material
    // information to information about block connectivity.
    //
    virtual void      *GetAuxiliaryData(const char *var, int domain,
                                        const char *type, void *args,
                                        DestructorFunction &);

    //
    // If you know the cycle number, overload this function.
    // Otherwise, VisIt will make up a reasonable one for you.
    //
    // virtual int         GetCycle(void);
    //

    virtual const char    *GetType(void)   { return "SDF"; };
    virtual void           FreeUpResources(void);

    virtual vtkDataSet    *GetMesh(int, const char *);
    virtual vtkDataArray  *GetVar(int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, const char *);

    virtual void           ActivateTimestep(void);
    virtual bool          HasInvariantMetaData(void) const { return false; };
    virtual bool          HasInvariantSIL(void) const      { return false; };
    virtual void           SetUpDomainConnectivity(void);
    int         GetCycle(void) { return step; }
    int         GetCycle(int) { return step; }
    double      GetTime(void) { return time; }
    double      GetTime(int) { return time; }
    int         GetCycleFromFilename(const char *f) const { return step; }
    double      GetTimeFromFilename(const char *f) const { return time; }
    bool        ReturnsValidCycle(void) const { return true; }
    bool        ReturnsValidTime(void) const { return true; }
/*
    virtual void           PopulateIOInformation(avtIOInformation& ioInfo);
    virtual void           SetDatabaseMetaData(avtDatabaseMetaData *md);

    void                   SetTimestep(int ts, int ns);


    void                   RegisterDatabaseMetaData(avtDatabaseMetaData *);

    void                  SetCache(avtVariableCache *);
    virtual bool          PerformsMaterialSelection(void) { return false; };
    virtual bool          HasVarsDefinedOnSubMeshes(void) { return false; };
    virtual bool          HasInvariantSIL(void) const      { return true; };
    virtual void          TurnMaterialSelectionOff(void);
    virtual void          TurnMaterialSelectionOn(const char *);

    virtual bool          CanCacheVariable(const char *) { return true; };

    bool                  CanDoStreaming(void)
                              { return canDoStreaming; };

    virtual void          RegisterVariableList(const char *,
                                          const std::vector<CharStrRef> &) {;};

    virtual void          RegisterDataSelections(
                              const std::vector<avtDataSelection_p>&,
                              std::vector<bool> *wasApplied) {;};

    void                  SetResultMustBeProducedOnlyOnThisProcessor(bool b)
                            { resultMustBeProducedOnlyOnThisProcessor = b; };


    virtual void          GetCycles(std::vector<int>&) { return; };
    virtual int           GetCycle(void) { return INVALID_CYCLE; };
    virtual int           GetCycle(int) { return INVALID_CYCLE; };
    virtual void          GetTimes(std::vector<double>&) { return; };
    virtual double        GetTime(void) { return INVALID_TIME; };
    virtual double        GetTime(int) { return INVALID_TIME; };

    //
    // These methods are designed so that we can distinguish between
    // the default implementation and one that is overridden in a plugin.
    // A plugin is expected to return INVALID_CYCLE or INVALID_TIME for any
    // situation in which it cannot return what it thinks is a valid
    // cycle/time from a filename. The default methods return something
    // slightly different. The reason we do this is that our guesses are
    // NOT to be trusted, but a plugin's guesses are. So, we need to know
    // the difference.
    //
    virtual int           GetCycleFromFilename(const char *f) const
                              { if (f[0] == '\0') return FORMAT_INVALID_CYCLE;
                                return GuessCycle(f); };
    virtual double        GetTimeFromFilename(const char *f) const
                              { if (f[0] == '\0') return FORMAT_INVALID_TIME;
                                return GuessTime(f); };

    void       AddMeshToMetaData(avtDatabaseMetaData *, std::string,
                                 avtMeshType, const double * = NULL, int = 1,
                                 int = 0, int = 3, int = 3);
    void       AddScalarVarToMetaData(avtDatabaseMetaData *, std::string,
                                      std::string, avtCentering,
                                      const double * = NULL,
                                      const bool = false);
    void       AddVectorVarToMetaData(avtDatabaseMetaData *, std::string,
                                      std::string, avtCentering, int = 3,
                                      const double * = NULL);
    void       AddTensorVarToMetaData(avtDatabaseMetaData *, std::string,
                                      std::string, avtCentering, int = 3);
    void       AddSymmetricTensorVarToMetaData(avtDatabaseMetaData *,
                              std::string, std::string, avtCentering, int = 3);
    void       AddMaterialToMetaData(avtDatabaseMetaData *, std::string,
                                     std::string,int,std::vector<std::string>);
    void       AddSpeciesToMetaData(avtDatabaseMetaData *, std::string,
                                    std::string, std::string, int,
                                    std::vector<int>,
                                    std::vector<std::vector<std::string> >);
    void       AddArrayVarToMetaData(avtDatabaseMetaData *,
                                     std::string, std::vector<std::string> &,
                                     std::string, avtCentering);
    void       AddArrayVarToMetaData(avtDatabaseMetaData *, std::string, int,
                                     std::string, avtCentering);

*/

  protected:
    // DATA MEMBERS
    int rank, ncpus, ndomains, step;
    int use_float, use_random, use_boundary, use_allboundary, use_ob_boundary;
    double time;
    comm_t comm;
    sdf_file_t *h;
    stack_handle_t *stack_handle;
    char *filename;
    bool gotMetadata;

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *);
    char *GetCompositeName(sdf_block_t *b);
    void *GetMaterial(const char *var, int domain);
    void *GetSpecies(const char *var, int domain);
    template <typename Real> void * GetMaterialType(sdf_block_t *, int);
    template <typename Real> void * GetSpeciesType(sdf_block_t *, int);
    vtkDataSet *GetCurve(int domain, sdf_block_t *b);
    sdf_block_t *GetArray(int, const char *);
    void OpenFile(int);
    void FillGhost(int domain, vtkDataSet *ds);
};


#endif
