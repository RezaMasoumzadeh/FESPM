
#ifndef __spm_h__
#define __spm_h__

#include <stdio.h>
#include <stdlib.h>

typedef int ErrorCode;

typedef enum {  DM_CONST, DM_NONLINEAR_CRITICALSLOPE, DM_NONLINEAR_POWERLAW,
                DM_UNDEFINED } FESLSPMDiffusionMode;
typedef enum { DTM_CONST, DTM_EXPLICIT } FESLSPMTimestepMode;

typedef enum { DIR_JMAX=0, DIR_IMAX, DIR_JMIN, DIR_IMIN } FESLSPMDirection;

typedef enum { FVT_ASCII, FVT_VTR, FVT_VTRB, FVT_UNDEFINED } FESLSPMFieldViewType;

typedef enum { SLINTERPOLATION_BILINEAR=0,
               SLINTERPOLATION_CSPLINE,
               SLINTERPOLATION_CSPLINE_MONOTONE,
               SLINTERPOLATION_UPWIND_CSPLINE,
               SLINTERPOLATION_UPWIND_CSPLINE_MONOTONE,
               SLINTERPOLATION_UNDEFINED } SLInterpolationType;

#define FESLSPM_IGNORE -1
#define FESLSPM_DECIDE -2

typedef struct _p_FESLSPM *FESLSPM;

typedef struct _p_FESPM *FESPM;

struct _p_FESLSPM {
    double xrange[2];
    double yrange[2];
    double time_end,dt,dt_const,time;
    int nx,ny,mx,my,nnodes,ncells;
    double dx,dy,dt_theta;
    int nt,nt_io,nt_max;
    FESLSPMDiffusionMode mode;
    FESLSPMTimestepMode dtmode;
    double *h_ic,*h,*h_c,*u_v,*u_h,*k0_n,*d_exp_n,*Sc_n; /* node */
    double *kappa; /* cell */
    double k0,d_exp,Sc;
    char bctype[4];
    double dt_min,dt_max;
    SLInterpolationType semilagrange_interpolation_type;

        /*
     * FEM diffusion module.
     */
    FESPM fespm;
};

#define FESLSPMCheckIsIgnore(a) ( (int)((a)) == FESLSPM_IGNORE ? 1 : 0 )

#include "spmgriditerator.h"

#define ERRQ(str) { printf("*** FESLSPM ERROR: -------------------------- Error Message --------------------------\n"\
                           "*** FESLSPM ERROR: %s\n"\
                           "*** FESLSPM ERROR: %s() line %d in %s\n"\
                           "*** FESLSPM ERROR: -------------------------- End of Error Message --------------------------\n",\
                           (str),__FUNCTION__,__LINE__,__FILE__); exit(EXIT_FAILURE); }

#define ERRQ2(str1,str2) { printf("*** FESLSPM ERROR: -------------------------- Error Message --------------------------\n"\
                                  "*** FESLSPM ERROR: %s %s\n"\
                                  "*** FESLSPM ERROR: %s() line %d in %s\n"\
                                  "*** FESLSPM ERROR: -------------------------- End of Error Message --------------------------\n",\
                                  (str1),(str2),__FUNCTION__,__LINE__,__FILE__); exit(EXIT_FAILURE); }


ErrorCode FESLSPMCreate(FESLSPM *s);
ErrorCode FESLSPMDestroy(FESLSPM *s);

ErrorCode FESLSPMSetGeometry(FESLSPM s,double x0,double x1,double y0,double y1);
ErrorCode FESLSPMSetResolution(FESLSPM s,int nx,int ny);

ErrorCode FESLSPMSetTimestepBounds(FESLSPM spm,double dt_min,double dt_max);
ErrorCode FESLSPMSetTimestep(FESLSPM s,double dt);
ErrorCode FESLSPMSetTimestepTheta(FESLSPM s,double theta);
ErrorCode FESLSPMSetMaxNumberSteps(FESLSPM s,int nt_max);
ErrorCode FESLSPMSetEndTime(FESLSPM s,double tend);
ErrorCode FESLSPMSetTimestepMode(FESLSPM s,FESLSPMTimestepMode mode);
ErrorCode FESLSPMSetTimestepIO(FESLSPM s,int nt);
ErrorCode FESLSPMSetConstantTimestep(FESLSPM s,double dt);
ErrorCode FESLSPMSetExplicitTimestep(FESLSPM s,double theta);

ErrorCode FESLSPMSetDiffusionMode(FESLSPM s,FESLSPMDiffusionMode mode);
ErrorCode FESLSPMSetConstantDiffusivity(FESLSPM s,double k0);
ErrorCode FESLSPMSetNonlinearDiffusivityCriticalSlope(FESLSPM s,double k0,double Sc,double d_exp);
ErrorCode FESLSPMSetNonlinearDiffusivityPowerLaw(FESLSPM s,double k0,double d_exp);

ErrorCode FESLSPMSetInitialHeight(FESLSPM s,double h0[]);
ErrorCode FESLSPMSetVerticalVelocity(FESLSPM s,double v[]);
ErrorCode FESLSPMSetHorizontalVelocity(FESLSPM s,double uw[]);

ErrorCode FESLSPMSetBoundaryConditions(FESLSPM s,const char bctype[]);
ErrorCode FESLSPMUpdate_EnforceBoundaryConditions(FESLSPM s,double field[]);

ErrorCode FESLSPMGetNodeCoordinate(FESLSPM s,int i,int j,double coor[]);
ErrorCode FESLSPMGetCentroidCoordinate(FESLSPM s,int i,int j,double coor[]);
ErrorCode FESLSPMGetVerticalVelocity(FESLSPM s,int *m,int *n,double **data);
ErrorCode FESLSPMGetHorizontalVelocity(FESLSPM s,int *m,int *n,double **data);
ErrorCode FESLSPMGetKappa(FESLSPM s,int *m,int *n,double **data);
ErrorCode FESLSPMGetHeight(FESLSPM s,int *m,int *n,double **data);

ErrorCode FESLSPMFieldView_ascii(FESLSPM s);
ErrorCode FESLSPMFieldView_VTRB(FESLSPM s,const char filename[]);
ErrorCode FESLSPMFieldView_VTR(FESLSPM s,const char filename[]);
ErrorCode FESLSPMView_Binary(FESLSPM s,const char name[]);
ErrorCode FESLSPMCheckpoint(FESLSPM s);
ErrorCode FESLSPMFieldView(FESLSPM s,FESLSPMFieldViewType type);
ErrorCode FESLSPMLoadFromCheckpoint(FESLSPM s,const char name[]);


ErrorCode FESLSPMComputeTimestep(FESLSPM s,double h[],double *dt);

ErrorCode FESLSPMView(FESLSPM s);
ErrorCode FESLSPMStepInitialize(FESLSPM s);

ErrorCode FESLSPMUpdate_SemiLagrangian(FESLSPM s,double dt,double hc[]);
ErrorCode FESLSPMStepAdvDiff(FESLSPM s);

ErrorCode FESLSPMNodeDataGetRange(FESLSPM s,double x[],int start,int stride,double range[]);
ErrorCode FESLSPMCellDataGetRange(FESLSPM s,double x[],int start,int stride,double range[]);
ErrorCode FESLSPMNodeDataSet(FESLSPM s,double x[],double val);
ErrorCode FESLSPMCellDataSet(FESLSPM s,double x[],double val);
ErrorCode FESLSPMNodeDataScale(FESLSPM s,double x[],double val);
ErrorCode FESLSPMCellDataScale(FESLSPM s,double x[],double val);
ErrorCode FESLSPMNodeDataShift(FESLSPM s,double x[],double val);
ErrorCode FESLSPMCellDataShift(FESLSPM s,double x[],double val);

ErrorCode FESLSPMCreateNodalField(FESLSPM s,int ncomponents,double *nodadata[]);
ErrorCode FESLSPMCreateCellField(FESLSPM s,int ncomponents,double *celldata[]);
ErrorCode FESLSPMDestroyNodalField(FESLSPM s,double **data);
ErrorCode FESLSPMDestroyCellField(FESLSPM s,double **data);
ErrorCode FESLSPMCellGetNodeValues(FESLSPM s,int ci,int cj,int ncomponents,double nodalfield[],double vertexvalues[]);
ErrorCode FESLSPMProjectNodalField(FESLSPM s,double nodedata[]);

ErrorCode FESLSPMSetSemiLagrangeInterpolationType(FESLSPM spm,SLInterpolationType type);

#endif
