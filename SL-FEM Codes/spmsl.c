#include <stdio.h>
#include <math.h>
#include "spm.h"

#define EPS 1.0e-8

/* Maple generated code */
ErrorCode matrixinverse_3x3(double A[3][3],double B[3][3])
{
	double t4,t6,t8,t10,t12,t14,t17;

	t4 = A[2][0] * A[0][1];
	t6 = A[2][0] * A[0][2];
	t8 = A[1][0] * A[0][1];
	t10 = A[1][0] * A[0][2];
	t12 = A[0][0] * A[1][1];
	t14 = A[0][0] * A[1][2];
	t17 = 0.1e1 / (t4 * A[1][2] - t6 * A[1][1] - t8 * A[2][2] + t10 * A[2][1] + t12 * A[2][2] - t14 * A[2][1]);

	B[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) * t17;
	B[0][1] = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]) * t17;
	B[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * t17;
	B[1][0] = -(-A[2][0] * A[1][2] + A[1][0] * A[2][2]) * t17;
	B[1][1] = (-t6 + A[0][0] * A[2][2]) * t17;
	B[1][2] = -(-t10 + t14) * t17;
	B[2][0] = (-A[2][0] * A[1][1] + A[1][0] * A[2][1]) * t17;
	B[2][1] = -(-t4 + A[0][0] * A[2][1]) * t17;
	B[2][2] = (-t8 + t12) * t17;

    return(0);
}

/* http://en.wikipedia.org/wiki/Spline_interpolation */
ErrorCode _cubicsplineinterp(double xi,double xn[],double yn[],double *yi)
{
    double A[3][3],B[3][3],k[3],b[3];
    double x1x02,x2x12;
    double qa[2],qb[2];
    int i;
    ErrorCode ierr;

    if (xi < xn[0]-EPS) { printf("cubicspline: out of range xn[0]\n"); }
    if (xi > xn[2]+EPS) { printf("cubicspline: out of range xn[2]\n"); }

    A[0][0] = 2.0 / (xn[1] - xn[0]);
    A[0][1] = 1.0 / (xn[1] - xn[0]);
    A[0][2] = 0.0;

    A[1][0] = A[0][1];
    A[1][1] = 2.0 * ( 1.0/(xn[1]-xn[0]) + 1.0/(xn[2]-xn[1]) );
    A[1][2] = 1.0/(xn[2]-xn[1]);

    A[2][0] = 0.0;
    A[2][1] = A[1][2];
    A[2][2] = 2.0 / (xn[2] - xn[1]);

    ierr = matrixinverse_3x3(A,B);

    x1x02 = (xn[1] - xn[0])*(xn[1] - xn[0]);
    x2x12 = (xn[2] - xn[1])*(xn[2] - xn[1]);

    b[0] = 3.0 * (yn[1] - yn[0])/x1x02;
    b[2] = 3.0 * (yn[2] - yn[1])/x2x12;
    b[1] = b[0] + b[2];

    k[0] = k[1] = k[2] = 0.0;
    for (i=0; i<3; i++) {
        k[0] += B[0][i] * b[i];
        k[1] += B[1][i] * b[i];
        k[2] += B[2][i] * b[i];
    }
    //printf("%1.6e %1.6e %1.6e\n",k[0],k[1],k[2]);

    qa[0] = k[0] * (xn[1] - xn[0]) - (yn[1] - yn[0]);
    qa[1] = k[1] * (xn[2] - xn[1]) - (yn[2] - yn[1]);

    qb[0] = -k[1] * (xn[1] - xn[0]) + (yn[1] - yn[0]);
    qb[1] = -k[2] * (xn[2] - xn[1]) + (yn[2] - yn[1]);

    //printf("a = (%1.6e , %1.6e)\n",qa[0],qa[1]);
    //printf("b = (%1.6e , %1.6e)\n",qb[0],qb[1]);

    if (xi < xn[1]) {
        double t = (xi - xn[0])/(xn[1] - xn[0]);

        *yi = (1.0-t) * yn[0] + t * yn[1] + t*(1.0-t)*(qa[0]*(1.0-t) + qb[0]*t);
    } else {
        double t = (xi - xn[1])/(xn[2] - xn[1]);

        *yi = (1.0-t) * yn[1] + t * yn[2] + t*(1.0-t)*(qa[1]*(1.0-t) + qb[1]*t);
    }

    return(0);
}

ErrorCode cubicspline_preparebasis(double xn[],double B[3][3])
{
    double A[3][3];
    double x1x02,x2x12;
    ErrorCode ierr;

    A[0][0] = 2.0 / (xn[1] - xn[0]);
    A[0][1] = 1.0 / (xn[1] - xn[0]);
    A[0][2] = 0.0;

    A[1][0] = A[0][1];
    A[1][1] = 2.0 * ( 1.0/(xn[1]-xn[0]) + 1.0/(xn[2]-xn[1]) );
    A[1][2] = 1.0/(xn[2]-xn[1]);

    A[2][0] = 0.0;
    A[2][1] = A[1][2];
    A[2][2] = 2.0 / (xn[2] - xn[1]);

    ierr = matrixinverse_3x3(A,B);
    return(0);
}

ErrorCode cubicspline_interpolate(double xi,double xn[],double B[3][3],double yn[],double *yi)
{
    double k[3],b[3];
    double x1x02,x2x12;
    double qa[2],qb[2];
    int i;

    if (xi < xn[0]-EPS) ERRQ("cubicspline: out of range xn[0]");
    if (xi > xn[2]+EPS) ERRQ("cubicspline: out of range xn[2]");

    x1x02 = (xn[1] - xn[0])*(xn[1] - xn[0]);
    x2x12 = (xn[2] - xn[1])*(xn[2] - xn[1]);

    b[0] = 3.0 * (yn[1] - yn[0])/x1x02;
    b[2] = 3.0 * (yn[2] - yn[1])/x2x12;
    b[1] = b[0] + b[2];

    k[0] = k[1] = k[2] = 0.0;
    for (i=0; i<3; i++) {
        k[0] += B[0][i] * b[i];
        k[1] += B[1][i] * b[i];
        k[2] += B[2][i] * b[i];
    }

    qa[0] = k[0] * (xn[1] - xn[0]) - (yn[1] - yn[0]);
    qa[1] = k[1] * (xn[2] - xn[1]) - (yn[2] - yn[1]);

    qb[0] = -k[1] * (xn[1] - xn[0]) + (yn[1] - yn[0]);
    qb[1] = -k[2] * (xn[2] - xn[1]) + (yn[2] - yn[1]);

    if (xi < xn[1]) {
        double t = (xi - xn[0])/(xn[1] - xn[0]);

        *yi = (1.0-t) * yn[0] + t * yn[1] + t*(1.0-t)*(qa[0]*(1.0-t) + qb[0]*t);
    } else {
        double t = (xi - xn[1])/(xn[2] - xn[1]);

        *yi = (1.0-t) * yn[1] + t * yn[2] + t*(1.0-t)*(qa[1]*(1.0-t) + qb[1]*t);
    }

    return(0);
}


ErrorCode SL_ComputeDeparturePoint_RK1(FESLSPM s,double dt,double xn[],double xd[])
{
    ErrorCode ierr;
    int ii,jj,nidx_ij;
    double vx[2];
    double xp[2],dtx,dty,dt_min,xboundary,yboundary;
    int outside_domain,bflag_x,bflag_y;


    /* fetch velocity  - no need to interpolate */
    ii = (int)( (xn[0] - s->xrange[0])/s->dx );
    jj = (int)( (xn[1] - s->yrange[0])/s->dy );

    nidx_ij = NodeIdx(s,ii,jj);

    vx[0] = s->u_h[2*nidx_ij+0];
    vx[1] = s->u_h[2*nidx_ij+1];

    xp[0] = xn[0];
    xp[1] = xn[1];

    /* note backwards advection */
    xp[0] = xp[0] - dt * vx[0];
    xp[1] = xp[1] - dt * vx[1];

    outside_domain = 0;
    bflag_x = bflag_y = 0;

    if (xp[0] < s->xrange[0]) {
        bflag_x = -1;
        outside_domain = 1;
    }
    if (xp[0] > s->xrange[1]) {
        bflag_x = 1;
        outside_domain = 1;
    }
    if (xp[1] < s->yrange[0]) {
        bflag_y = -1;
        outside_domain = 1;
    }
    if (xp[1] > s->yrange[1]) {
        bflag_y = 1;
        outside_domain = 1;
    }

    /*
     Solve for
     xboundary = xd - dtx vx
     yboundary = yd - dty vy
    */
    if (outside_domain == 1) {
        double xbound[2],ybound[2];

        switch (bflag_x) {
            case -1: /* left */
                xboundary = s->xrange[0];
                dt_min = (xn[0] - xboundary) / vx[0];
                xp[0] = xboundary;
                xp[1] = xn[1] - dt_min * vx[1];

                xbound[0] = xp[0];
                xbound[1] = xp[1];
                break;
            case  1: /* right */
                xboundary = s->xrange[1];
                dt_min = (xn[0] - xboundary) / vx[0];
                xp[0] = xboundary;
                xp[1] = xn[1] - dt_min * vx[1];

                xbound[0] = xp[0];
                xbound[1] = xp[1];
                break;
        }

        switch (bflag_y) {
            case -1: /* bottom */
                yboundary = s->yrange[0];
                dt_min = (xn[1] - yboundary) / vx[1];
                xp[0] = xn[0] - dt_min * vx[0];
                xp[1] = yboundary;

                ybound[0] = xp[0];
                ybound[1] = xp[1];
                break;
            case  1: /* top */
                yboundary = s->yrange[1];
                dt_min = (xn[1] - yboundary) / vx[1];
                xp[0] = xn[0] - dt_min * vx[0];
                xp[1] = yboundary;

                ybound[0] = xp[0];
                ybound[1] = xp[1];
                break;
        }

        if (bflag_x != 0) {
            if ( (xbound[0] >= s->xrange[0]) && (xbound[0] <= s->xrange[1]) ) {
                if ( (xbound[1] >= s->yrange[0]) && (xbound[1] <= s->yrange[1]) ) {
                    xp[0] = xbound[0];
                    xp[1] = xbound[1];
                }
            }
        }
        if (bflag_y != 0) {
            if ( (ybound[0] >= s->xrange[0]) && (ybound[0] <= s->xrange[1]) ) {
                if ( (ybound[1] >= s->yrange[0]) && (ybound[1] <= s->yrange[1]) ) {
                    xp[0] = ybound[0];
                    xp[1] = ybound[1];
                }
            }
        }

    }

    xd[0] = xp[0];
    xd[1] = xp[1];

    return(0);
}

ErrorCode SL_PointInsideDomain(FESLSPM s,double pos[],int *inside)
{
    *inside  = 1;
    if (pos[0] < s->xrange[0]-EPS) { *inside = 0; return(0); }
    if (pos[0] > s->xrange[1]+EPS) { *inside = 0; return(0); }
    if (pos[1] < s->yrange[0]-EPS) { *inside = 0; return(0); }
    if (pos[1] > s->yrange[1]+EPS) { *inside = 0; return(0); }

    return(0);
}

ErrorCode SL_ComputeVelocity(FESLSPM s,double pos[],double vel[])
{
    int ci,cj,cij;
    int nidx[4],k;
    double x0[2],xi[2],Ni[4];

    ci = (int)( (pos[0] - s->xrange[0])/s->dx );
    cj = (int)( (pos[1] - s->yrange[0])/s->dy );
    if (ci == s->mx) ci--;
    if (cj == s->my) cj--;

    cij = CellIdx(s,ci,cj);

    nidx[0] = NodeIdx(s,ci,  cj);
    nidx[1] = NodeIdx(s,ci+1,cj);
    nidx[2] = NodeIdx(s,ci,  cj+1);
    nidx[3] = NodeIdx(s,ci+1,cj+1);

    x0[0] = s->xrange[0] + ci * s->dx;
    x0[1] = s->yrange[0] + cj * s->dy;

    /*
     [x - x0]/dx = [xi - (-1)]/2
    */

    xi[0] = 2.0*(pos[0]-x0[0])/s->dx - 1.0;
    xi[1] = 2.0*(pos[1]-x0[1])/s->dy - 1.0;

    Ni[0] = 0.25 * (1.0 - xi[0]) * (1.0 - xi[1]);
    Ni[1] = 0.25 * (1.0 + xi[0]) * (1.0 - xi[1]);
    Ni[2] = 0.25 * (1.0 - xi[0]) * (1.0 + xi[1]);
    Ni[3] = 0.25 * (1.0 + xi[0]) * (1.0 + xi[1]);

    vel[0] = vel[1] = 0.0;
    for (k=0; k<4; k++) {
        vel[0] += Ni[k] * s->u_h[2*nidx[k]+0];
        vel[1] += Ni[k] * s->u_h[2*nidx[k]+1];
    }

    return(0);
}

ErrorCode SL_ComputeDeparturePoint_RK4(FESLSPM s,double dt,double xn[],double xd[],int *update_key)
{
    double k1[2],k2[2],k3[2],k4[2];
    int inside;
    double pos[2];
    ErrorCode ierr;

    *update_key = -1;

    pos[0] = xn[0];
    pos[1] = xn[1];
    ierr = SL_ComputeVelocity(s,pos,k1); /* t = t_k */
    k1[0] = -k1[0];
    k1[1] = -k1[1];

    pos[0] = xn[0] + 0.5 * dt * k1[0];
    pos[1] = xn[1] + 0.5 * dt * k1[1];
    ierr = SL_PointInsideDomain(s,pos,&inside);
    if (inside == 0) { *update_key = 0; return(0); }
    ierr = SL_ComputeVelocity(s,pos,k2); /* t = t_k - 0.5*dt */
    k2[0] = -k2[0];
    k2[1] = -k2[1];

    pos[0] = xn[0] + 0.5 * dt * k2[0];
    pos[1] = xn[1] + 0.5 * dt * k2[1];
    ierr = SL_PointInsideDomain(s,pos,&inside);
    if (inside == 0) { *update_key = 0; return(0); }
    ierr = SL_ComputeVelocity(s,pos,k3); /* t = t_k - 0.5*dt */
    k3[0] = -k3[0];
    k3[1] = -k3[1];

    pos[0] = xn[0] + dt * k3[0];
    pos[1] = xn[1] + dt * k3[1];
    ierr = SL_PointInsideDomain(s,pos,&inside);
    if (inside == 0) { *update_key = 0; return(0); }
    ierr = SL_ComputeVelocity(s,pos,k4); /* t = t_k - dt */
    k4[0] = -k4[0];
    k4[1] = -k4[1];

    xd[0] = xn[0] + (dt/6.0)*( k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0] );
    xd[1] = xn[1] + (dt/6.0)*( k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1] );
    ierr = SL_PointInsideDomain(s,xd,&inside);
    if (inside == 0) { *update_key = 0; return(0); }

    *update_key = 1;

    return(0);
}

ErrorCode SL_ComputeDeparturePoint(FESLSPM s,double dt,double xn[],double xd[])
{
    ErrorCode ierr;
    int update_key,point_inside;

    xd[0] = -100.0 * s->xrange[0];
    xd[1] = -100.0 * s->yrange[0];

    update_key = 0;
    ierr = SL_ComputeDeparturePoint_RK4(s,dt,xn,xd,&update_key);

    if (update_key == 0) { /* point left domain */
        ierr = SL_ComputeDeparturePoint_RK1(s,dt,xn,xd);
    }

    ierr = SL_PointInsideDomain(s,xd,&point_inside);
    if (point_inside == 0) {
        ERRQ("ComputeDeparturePoint: Failed... point outside domain");
    }


    return(0);
}

/*
[input]
  xd[] : single departure point
  h[] : array of all heights
[output]
  hd : interpolated height at departure point
*/
/* interpolate in y first, then x : arbitrary choice */
ErrorCode SL_InterpolateAPointCSpline(FESLSPM s,double xd[],double h[],double *hd)
{
    ErrorCode ierr;
    int ci,cj,cidx;
    int ni[3],nj[3],kk;
    double yj[3],xj[3];
    double invA[3][3],hinterp[3],hj[3];


    /* determine cell containing point */
    ci = (int)( (xd[0] - s->xrange[0])/s->dx );
    cj = (int)( (xd[1] - s->yrange[0])/s->dy );
    if (ci == s->mx) ci--;
    if (cj == s->my) cj--;

    cidx = CellIdx(s,ci,cj);

    /* determine 3 nj values to use for cspline interpolation */
    nj[0] = cj-1;
    nj[1] = cj;
    nj[2] = cj+1;

    if (nj[0] < 0) {
        nj[0] = 0;
        nj[1] = 1;
        nj[2] = 2;
    } else if (nj[2] >= s->ny) {
        nj[0] = s->ny-3;
        nj[1] = s->ny-2;
        nj[2] = s->ny-1;
    }

    /* determine 3 ni values to use for cspline interpolation */
    ni[0] = ci-1;
    ni[1] = ci;
    ni[2] = ci+1;

    if (ni[0] < 0) {
        ni[0] = 0;
        ni[1] = 1;
        ni[2] = 2;
    } else if (ni[2] >= s->nx) {
        ni[0] = s->nx-3;
        ni[1] = s->nx-2;
        ni[2] = s->nx-1;
    }

    for (kk=0; kk<3; kk++) {
        yj[kk] = s->yrange[0] + nj[kk] * s->dy;

        xj[kk] = s->xrange[0] + ni[kk] * s->dx;
    }

    if (xd[0] > xj[2]+EPS) { printf("xd>2 (%1.6e,%1.6e): %1.6e %1.6e %1.6e - %d %d %d\n",xd[0],xd[1],xj[0],xj[1],xj[2],ni[0],ni[1],ni[2]); exit(0); }
    if (xd[0] < xj[0]-EPS) { printf("xd<0 (%1.6e,%1.6e): %1.6e %1.6e %1.6e - %d %d %d\n",xd[0],xd[1],xj[0],xj[1],xj[2],ni[0],ni[1],ni[2]); exit(0); }

    if (xd[1] > yj[2]+EPS) { printf("yd>2 (%1.6e,%1.6e): %1.6e %1.6e %1.6e - %d %d %d\n",xd[0],xd[1],yj[0],yj[1],yj[2],nj[0],nj[1],nj[2]); exit(0); }
    if (xd[1] < yj[0]-EPS) { printf("yd<0 (%1.6e,%1.6e): %1.6e %1.6e %1.6e - %d %d %d\n",xd[0],xd[1],yj[0],yj[1],yj[2],nj[0],nj[1],nj[2]); exit(0); }

    /* interpolate three times to determin the points in the y direction */
    ierr = cubicspline_preparebasis(yj,invA);

    for (kk=0; kk<3; kk++) {
        int nidx = NodeIdx(s,ni[0],nj[kk]);

        hj[kk] = h[nidx];
    }
    ierr = cubicspline_interpolate(xd[1],yj,invA,hj,&hinterp[0]);

    for (kk=0; kk<3; kk++) {
        int nidx = NodeIdx(s,ni[1],nj[kk]);

        hj[kk] = h[nidx];
    }
    ierr = cubicspline_interpolate(xd[1],yj,invA,hj,&hinterp[1]);

    for (kk=0; kk<3; kk++) {
        int nidx = NodeIdx(s,ni[2],nj[kk]);

        hj[kk] = h[nidx];
    }
    ierr = cubicspline_interpolate(xd[1],yj,invA,hj,&hinterp[2]);

    /* interpolate once in the x direction */
    ierr = cubicspline_preparebasis(xj,invA);

    ierr = cubicspline_interpolate(xd[0],xj,invA,hinterp,hd);

    return(0);
}

/*
 Determines cubic-spline sample knots based on the direction of the flow field

 [input]
   xd[] : single departure point
   vd[] : velocity at departure point
   h[] : array of all heights
 [output]
   hd : interpolated height at departure point
 */
/* interpolate in y first, then x : arbitrary choice */
ErrorCode SL_InterpolateAPointUpwindCSpline(FESLSPM s,double xd[],double vd[],double h[],double *hd)
{
  ErrorCode ierr;
  int ci,cj,cidx;
  int ni[3],nj[3],kk;
  double yj[3],xj[3];
  double invA[3][3],hinterp[3],hj[3];


  /* determine cell containing point */
  ci = (int)( (xd[0] - s->xrange[0])/s->dx );
  cj = (int)( (xd[1] - s->yrange[0])/s->dy );
  if (ci == s->mx) ci--;
  if (cj == s->my) cj--;

  cidx = CellIdx(s,ci,cj);

  /*
     Determine 3 nj values to use for cspline interpolation
     biased towards direction where information is coming from
   */
  if (vd[1] > 0.0) {
    nj[0] = cj-1;
    nj[1] = cj;
    nj[2] = cj+1;
  } else {
    nj[0] = cj;
    nj[1] = cj+1;
    nj[2] = cj+2;
  }

  if (nj[0] < 0) {
    nj[0] = 0;
    nj[1] = 1;
    nj[2] = 2;
  } else if (nj[2] >= s->ny) {
    nj[0] = s->ny-3;
    nj[1] = s->ny-2;
    nj[2] = s->ny-1;
  }

  /* determine 3 ni values to use for cspline interpolation */
  if (vd[0] > 0.0) {
    ni[0] = ci-1;
    ni[1] = ci;
    ni[2] = ci+1;
  } else {
    ni[0] = ci;
    ni[1] = ci+1;
    ni[2] = ci+2;
  }

  if (ni[0] < 0) {
    ni[0] = 0;
    ni[1] = 1;
    ni[2] = 2;
  } else if (ni[2] >= s->nx) {
    ni[0] = s->nx-3;
    ni[1] = s->nx-2;
    ni[2] = s->nx-1;
  }

  for (kk=0; kk<3; kk++) {
    yj[kk] = s->yrange[0] + nj[kk] * s->dy;

    xj[kk] = s->xrange[0] + ni[kk] * s->dx;
  }

  if (xd[0] > xj[2]+EPS) { printf("xd>2 (%1.6e,%1.6e): %1.6e %1.6e %1.6e - %d %d %d\n",xd[0],xd[1],xj[0],xj[1],xj[2],ni[0],ni[1],ni[2]); exit(0); }
  if (xd[0] < xj[0]-EPS) { printf("xd<0 (%1.6e,%1.6e): %1.6e %1.6e %1.6e - %d %d %d\n",xd[0],xd[1],xj[0],xj[1],xj[2],ni[0],ni[1],ni[2]); exit(0); }

  if (xd[1] > yj[2]+EPS) { printf("yd>2 (%1.6e,%1.6e): %1.6e %1.6e %1.6e - %d %d %d\n",xd[0],xd[1],yj[0],yj[1],yj[2],nj[0],nj[1],nj[2]); exit(0); }
  if (xd[1] < yj[0]-EPS) { printf("yd<0 (%1.6e,%1.6e): %1.6e %1.6e %1.6e - %d %d %d\n",xd[0],xd[1],yj[0],yj[1],yj[2],nj[0],nj[1],nj[2]); exit(0); }

  /* interpolate three times to determin the points in the y direction */
  ierr = cubicspline_preparebasis(yj,invA);

  for (kk=0; kk<3; kk++) {
    int nidx = NodeIdx(s,ni[0],nj[kk]);

    hj[kk] = h[nidx];
  }
  ierr = cubicspline_interpolate(xd[1],yj,invA,hj,&hinterp[0]);

  for (kk=0; kk<3; kk++) {
    int nidx = NodeIdx(s,ni[1],nj[kk]);

    hj[kk] = h[nidx];
  }
  ierr = cubicspline_interpolate(xd[1],yj,invA,hj,&hinterp[1]);

  for (kk=0; kk<3; kk++) {
    int nidx = NodeIdx(s,ni[2],nj[kk]);

    hj[kk] = h[nidx];
  }
  ierr = cubicspline_interpolate(xd[1],yj,invA,hj,&hinterp[2]);

  /* interpolate once in the x direction */
  ierr = cubicspline_preparebasis(xj,invA);

  ierr = cubicspline_interpolate(xd[0],xj,invA,hinterp,hd);

  return(0);
}

/*
[input]
 xd[] : single departure point
 h[] : array of all heights
[output]
 hd : interpolated height at departure point
 */
/* interpolate in y first, then x : arbitrary choice */
ErrorCode SL_InterpolateAPointQ1(FESLSPM s,double xd[],double h[],double *hd)
{
    ErrorCode ierr;
    int ci,cj,cidx,k;
    int nidx[4];
    double hel[4],Ni[4],hint,x0[2],xi[2];

    /* determine cell containing point */
    ci = (int)( (xd[0] - s->xrange[0])/s->dx );
    cj = (int)( (xd[1] - s->yrange[0])/s->dy );
    if (ci == s->mx) ci--;
    if (cj == s->my) cj--;

    cidx = CellIdx(s,ci,cj);

    /* determine 4 nodes values to use for Q1 interpolation */
    nidx[0] = NodeIdx(s,ci  ,cj);
    nidx[1] = NodeIdx(s,ci+1,cj);
    nidx[2] = NodeIdx(s,ci  ,cj+1);
    nidx[3] = NodeIdx(s,ci+1,cj+1);

    x0[0] = s->xrange[0] + ci * s->dx;
    x0[1] = s->yrange[0] + cj * s->dy;

    xi[0] = 2.0*(xd[0]-x0[0])/s->dx - 1.0;
    xi[1] = 2.0*(xd[1]-x0[1])/s->dy - 1.0;

    Ni[0] = 0.25 * (1.0 - xi[0]) * (1.0 - xi[1]);
    Ni[1] = 0.25 * (1.0 + xi[0]) * (1.0 - xi[1]);
    Ni[2] = 0.25 * (1.0 - xi[0]) * (1.0 + xi[1]);
    Ni[3] = 0.25 * (1.0 + xi[0]) * (1.0 + xi[1]);

    for (k=0; k<4; k++) {
        hel[k] = h[nidx[k]];
    }
    hint = 0.0;
    for (k=0; k<4; k++) {
        hint += Ni[k] * hel[k];
    }
    *hd = hint;

    return(0);
}

/*
 [input]
 xd[] : single departure point
 h[] : array of all heights
 [output]
 hd : interpolated height at departure point
*/
/*
  Employs the method described in

  @inproceedings{paoli2006testing,
    title={Testing semi-Lagrangian schemes for two-phase flow applications},
    author={Paoli, R and Poinsot, T and Shariff, K},
    booktitle={Proceeding of the Summer program},
    year={2006},
    organization={Citeseer}
  }

  See equation (2.10)
  The essential idea is as follows
  [1] Compute the departure point x*
  [2] Compute a high-order approximation of \phi(x*) (using cubic splines), \phi_H(x*)
  [3] Compute the lower (\phi_{-}) / upper bounds (\phi_{+}) of \phi using the nodal values of \phi for the cell containint x*
  [4] Define \phi(x*) as
        \phi(x*) = \phi_{+} if \phi_H(x*) > \phi_{+}
        \phi(x*) = \phi_{-} if \phi_H(x*) < \phi_{-}
        \phi(x*) = \phi_H(x*) otherwise
*/
ErrorCode SL_InterpolateAPointQuasiMonotone(FESLSPM s,double xd[],double phi[],double *_phi_M)
{
  ErrorCode ierr;
  int ci,cj,cidx,k;
  int nidx[4];
  double phi_i,phi_H,phi_minus,phi_plus,phi_M;

  /* determine cell containing point */
  ci = (int)( (xd[0] - s->xrange[0])/s->dx );
  cj = (int)( (xd[1] - s->yrange[0])/s->dy );
  if (ci == s->mx) ci--;
  if (cj == s->my) cj--;

  cidx = CellIdx(s,ci,cj);

  /* determine 4 nodes values to use for Q1 interpolation */
  nidx[0] = NodeIdx(s,ci  ,cj);
  nidx[1] = NodeIdx(s,ci+1,cj);
  nidx[2] = NodeIdx(s,ci  ,cj+1);
  nidx[3] = NodeIdx(s,ci+1,cj+1);

  phi_minus = 1.0e20;
  phi_plus  = -1.0e20;
  for (k=0; k<4; k++) {
    phi_i = phi[ nidx[k] ];

    if (phi_i > phi_plus)  { phi_plus = phi_i; }
    if (phi_i < phi_minus) { phi_minus = phi_i; }
  }

  ierr = SL_InterpolateAPointCSpline(s,xd,phi,&phi_H);

  if (phi_H > phi_plus) {

    phi_M = phi_plus;

  } else if (phi_H < phi_minus) {

    phi_M = phi_minus;

  } else {

    phi_M = phi_H;

  }

  *_phi_M = phi_M;

  return(0);
}

ErrorCode SL_ApplyMonotoneLimiter(FESLSPM s,double xd[],double phi[],double phi_H,double *_phi_M)
{
  ErrorCode ierr;
  int ci,cj,cidx,k;
  int nidx[4];
  double phi_i,phi_minus,phi_plus,phi_M;

  /* determine cell containing point */
  ci = (int)( (xd[0] - s->xrange[0])/s->dx );
  cj = (int)( (xd[1] - s->yrange[0])/s->dy );
  if (ci == s->mx) ci--;
  if (cj == s->my) cj--;

  cidx = CellIdx(s,ci,cj);

  /* determine 4 nodes values to use for Q1 interpolation */
  nidx[0] = NodeIdx(s,ci  ,cj);
  nidx[1] = NodeIdx(s,ci+1,cj);
  nidx[2] = NodeIdx(s,ci  ,cj+1);
  nidx[3] = NodeIdx(s,ci+1,cj+1);

  phi_minus = 1.0e20;
  phi_plus  = -1.0e20;
  for (k=0; k<4; k++) {
    phi_i = phi[ nidx[k] ];
    if (phi_i > phi_plus)  { phi_plus = phi_i; }
    if (phi_i < phi_minus) { phi_minus = phi_i; }
  }

  if (phi_H > phi_plus) {
    phi_M = phi_plus;
  } else if (phi_H < phi_minus) {
    phi_M = phi_minus;
  } else {
    phi_M = phi_H;
  }
  *_phi_M = phi_M;

  return(0);
}

/*
[input]
  dt : timestep
  h[] : array of all heights
[output]
  h[] : array of all heights, corrected for horizontal advection
*/
ErrorCode FESLSPMUpdate_SemiLagrangian(FESLSPM s,double dt,double hc[])
{
    int nidx_ij,i;
    ErrorCode ierr;
    double xn[2],xd[2],hd;
    double *hnew;


    hnew = malloc(sizeof(double)*s->nnodes);

    NodalIteratorBegin(s)
    nidx_ij = NodeIdx(s,ni,nj);

    /* determine departure point */
    ierr = FESLSPMGetNodeCoordinate(s,ni,nj,xn);
    ierr = SL_ComputeDeparturePoint(s,dt,xn,xd);

    /* interpolate hc[] at the departure point */
    switch (s->semilagrange_interpolation_type) {
      case SLINTERPOLATION_BILINEAR:
        /* Note - if I had a discontinuity detector, I could switch between Q1 and C3 */
        ierr = SL_InterpolateAPointQ1(s,xd,hc,&hd);
        break;

      case SLINTERPOLATION_CSPLINE:
        ierr = SL_InterpolateAPointCSpline(s,xd,hc,&hd);
        break;

      case SLINTERPOLATION_CSPLINE_MONOTONE:
        ierr = SL_InterpolateAPointQuasiMonotone(s,xd,hc,&hd);
        break;

      case SLINTERPOLATION_UPWIND_CSPLINE:
      {
        double vd[2];

        /* compute velocity at departure point */
        ierr = SL_ComputeVelocity(s,xd,vd);

        /* do upwind selection of cspline points */
        ierr = SL_InterpolateAPointUpwindCSpline(s,xd,vd,hc,&hd);
      }
        break;

      case SLINTERPOLATION_UPWIND_CSPLINE_MONOTONE:
      {
        double vd[2],hd_H;

        ierr = SL_ComputeVelocity(s,xd,vd);
        ierr = SL_InterpolateAPointUpwindCSpline(s,xd,vd,hc,&hd_H);
        ierr = SL_ApplyMonotoneLimiter(s,xd,hc,hd_H,&hd);
      }
        break;

      default:
        ERRQ("A valid semi-Lagrangian interpolation type must be specified");
        break;
    }
    /* make update to h (delta h) */
    hnew[nidx_ij] = hd;

    NodalIteratorEnd(s)

    for (i=0; i<s->nnodes; i++) {
        hc[i] = hnew[i];
    }
    free(hnew);
    return(0);
}

/*
 dt = min(dx)/max(v)
 */
ErrorCode FESLSPMComputeTimestep_CFL(FESLSPM s,double uh[],double *dt)
{
  ErrorCode ierr;
  double sx,sy,ds,dt_min,vmax,vx,vy,vs;
  int nidx_ij;
  double CFL_max = 2.0;

  ds = s->dx;
  if (s->dy < ds) {
    ds = s->dy;
  }

  dt_min = 1.0e32;

  vmax = 0.0;
  NodalIteratorBegin(s)
  nidx_ij = NodeIdx(s,ni,nj);
  vx = uh[2*nidx_ij+0];
  vy = uh[2*nidx_ij+1];
  vs = vx*vx + vy*vy;
  if (vs > vmax) { vmax = vs; }
  NodalIteratorEnd(s)

  if (vmax < 1.0e-20) {
    *dt = 1.0e32;
  } else {
    *dt = CFL_max * (ds / sqrt(vmax));
  }

  return(0);
}
