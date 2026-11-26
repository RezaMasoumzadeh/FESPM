#ifndef __fespm_h__
#define __fespm_h__

#include "spm.h"

typedef struct _p_FESPM *FESPM;

/* Create the FEM diffusion module on the existing FESLSPM structured grid.
 * Set the desired diffusivity mode and parameters before calling this.
 * Supported modes: DM_CONST (sfcmod=4), DM_NONLINEAR_POWERLAW (5),
 * and DM_NONLINEAR_CRITICALSLOPE (6).
 */
ErrorCode FESPMCreate(FESLSPM spm, FESPM *fespm);

/* Destroy the FEM diffusion module. */
ErrorCode FESPMDestroy(FESPM *fespm);

/* Perform one FEM diffusion/uplift timestep after Semi-Lagrangian advection.
 * PDE: dh/dt = u_v + div(kappa grad h).
 * Uses explicit diffusion: M H^(n+1) = M H* + dt (M U_v - K(h*) H*).
 * For nonlinear modes, K is rebuilt from the supplied h_old at every step.
 */
ErrorCode FESPMDiffusionStep(FESPM fespm, FESLSPM spm, double h_old[]);

/* Estimate explicit-Euler stability limit from K q = lambda M q. */
ErrorCode FESPMGetStableTimestep(FESPM fespm, double *dt);

#endif
