#ifdef MAINFILE
#undef MAINFILE
#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>

#include "geo_head3mg.h"

/* Move markers by using Runge-Kutta method */
void movemark()
{
/* Counters */
int n1,n0,ynreset=0,ynres[1000],ynres1;
long int mm1,marknum1,m1,m2,m3,m4,p[8],wn1[500];
/* Val buffer */
double v[3][5],xold,yold,zold,mnu,mpb,mtk,eps1[500],vxyz1[10],markeii1,marksii1,markrii1,markdis1;
double vxwater,vywater,vzwater,dx,dy,dz,ival,dxgrid=0,dpdx,dpdy,dpdz,x,y,z,xkf,ykf,zkf,vxkoef,vykoef,vzkoef;
double topodens=2000.0;
/**/
/**/
/**/
/* Hydration front progress, melt extraction  */
/**/
/**/
/**/
/* Melt extraction */
if(timesum>0e+11 && nuend>3e+20)
	{
	meltextract();
	}
/*
*/
/**/
/* Topography surface calculation */
for (m3=0;m3<znumz;m3++)
for (m1=0;m1<xnumx;m1++)
	{
	for (m2=1;m2<ynumy;m2++)
		{
		m4=m3*xynumxy+m1*ynumy+m2;
		if(ro[m4]>topodens && ro[m4-1]<topodens)
			{
			bufv[2*xnumx*znumz+m3*xnumx+m1]=gy[m2-1]+(gy[m2]-gy[m2-1])*(topodens-ro[m4-1])/(ro[m4]-ro[m4-1]);
			break;
			}
		}
	}

/**/
/* Save number of markers */
marknum1=marknum;
/**/
/**/
/**/
/* Move markers */
#pragma omp parallel shared(ynres,n0,bufv,ynreset,markx,marky,markz,markk,markv,markt,markw,markd,marke,markex,markft,vx,vy,vz,pr,exx,eyy,ezz,exy,exz,eyz,sxx,syy,szz,sxy,sxz,syz) private(ynres1,n1,vxyz,xold,yold,zold,mm1,p,v,wn,eps,vxwater,vywater,vzwater,dx,dy,dz,dpdx,dpdy,dpdz,x,y,z,xkf,ykf,zkf,vxkoef,vykoef,vzkoef,m1,m2,m3,m4,mnu,mpb,mtk,ival,wn1,eps1,vxyz1,markeii1,marksii1,markrii1,markdis1) firstprivate(timestep,timesum,xnumx,ynumy,znumz,xynumxy,vyfluid,vymelt,nubeg,nuend,stp100,dxgrid,vdeep,zdeep,tdeep,xsize,ysize,zsize,outgrid,marknum,marknum1,gx,gy,gz,markll,marka0,marka1,markb0,markb1,marke0,marke1,marknu,markdh,markdv,markss,marks0,marks1,markmm,markn0,markn1)
{
/* Obtain cur thread number */
n1=omp_get_thread_num();
/* Obtain total number of threads */
if (n1==0) n0=omp_get_num_threads();
/* Reset marker reset identifier */
ynres1=0;
/*
*/
#pragma omp for schedule(static)
/*
*/
for (mm1=0;mm1<marknum;mm1++)
if((markx[mm1]>=0 && marky[mm1]>=0 && markz[mm1]>=0 && (double)(markx[mm1])<=xsize && (double)(marky[mm1])<=ysize && (double)(markz[mm1])<=zsize) || (outgrid!=1 && (markt[mm1]<50 || markt[mm1]>=100)))
	{
	/* Save initial coordinates */
	xold=markx[mm1];
	yold=marky[mm1];
	zold=markz[mm1];
	/**/
	/* Velocity, P, SIGii, EPSii invariants calc */
	vxyzcalc1(markx[mm1],marky[mm1],markz[mm1],vxyz1,wn1);
	epscalc1(markx[mm1],marky[mm1],markz[mm1],1,eps1,wn1);
	markeii1=pow(0.5*(eps1[4]*eps1[4]+eps1[5]*eps1[5]+eps1[6]*eps1[6])+eps1[7]*eps1[7]+eps1[8]*eps1[8]+eps1[9]*eps1[9],0.5);
	marksii1=pow(0.5*(eps1[54]*eps1[54]+eps1[55]*eps1[55]+eps1[56]*eps1[56])+eps1[57]*eps1[57]+eps1[58]*eps1[58]+eps1[59]*eps1[59],0.5);
	markdis1=eps1[60];
	mpb=eps1[10]*1e-5;
	mtk=(double)(markk[mm1]);
	/* Compute viscosity */
	if(markt[mm1]<50 && mtk>0)
		{
		/* Topography position */
		if(markt[mm1]>1 && markt[mm1]<20)
			{
			m4=2*xnumx*znumz+wn1[2]*xnumx+wn1[0];
			dx=(markx[mm1]-gx[wn1[0]])/(gx[wn1[0]+1]-gx[wn1[0]]);	
			dz=(markz[mm1]-gz[wn1[2]])/(gz[wn1[2]+1]-gz[wn1[2]]);	
			dy=(1.0-dx)*(1.0-dz)*bufv[m4]+dx*(1.0-dz)*bufv[m4+1]+(1.0-dx)*dz*bufv[m4+xnumx]+dx*dz*bufv[m4+xnumx+1];
			markto[mm1]=dy;
/*
{printf("A %ld %d %e %e %e %e    %e",mm1,markt[mm1],markx[mm1],marky[mm1],markz[mm1],mtk,markto[mm1]);getchar();}
*/
			}
/*
if(wn[2]==99 && wn[1]>20 && wn[0]>200) 
{
printf("%ld %d %e %e %e",mm1,markt[mm1],markeii,marksii,markv[mm1]);getchar();
}
*/
		/* Correct strain rate for numerical diffusion */
		markrii1=1.0;
		if(0==0 && markv[mm1]>0 && markx[mm1]>gx[3] && markx[mm1]<gx[xnumx-4] && marky[mm1]>gy[3] && marky[mm1]<gy[ynumy-4] && markz[mm1]>gz[3] && markz[mm1]<gz[znumz-4]) 
			{
/*
printf("%ld %d %e %e %e",mm1,markt[mm1],markeii,marksii,markv[mm1]);getchar();
*/
			markrii1*=pow(0.5*marksii1/markeii1/markv[mm1],0.87);
			}
		/* Viscosity save */
		mnu=viscalc(mtk,mpb,mm1,markt[mm1],markeii1,marksii1,markrii1,markdis1,wn1);
		markv[mm1]=mnu;
		/* Add plastic strain */
		if(marke[mm1]>0)
			{
			marke[mm1]+=timestep*markeii1;
/*
*/
			}
		/* Healing */
		if(marke[mm1]>0)
			{
			marke[mm1]-=timestep*healrate;
			if(marke[mm1]<0) marke[mm1]=0;
			}
		else
			{
			marke[mm1]+=timestep*healrate;
			if(marke[mm1]>0) marke[mm1]=0;
			}
        	/* Serpentinization of brittle crust and mantle faults at sub-surface */

        	if (0==0 && (markt[mm1]==9 || markt[mm1]==10) && marke[mm1]>dserp && marky[mm1]<yserp && markk[mm1]<tserp)
               		{
			markt[mm1]=13;
			}

		}
	/**/
	/**/
	/**/
	/* Water marker move */
	vxwater=vywater=vzwater=0;
	if(markt[mm1]>=50 && markt[mm1]<100) 
		{
		/* Water velocity */
		vywater=vyfluid; if(markd[mm1]>1100.0) vywater=vymelt;
		/* Fluid in rock */
		if(vyfluid>0 && (markk[mm1]==0 || markk[mm1]>298.0)) 
			{
			/* Horizontal,Vertical P-cell index */
			m1=wn1[0]; if(markx[mm1]>(gx[m1]+gx[m1+1])/2.0) m1+=1;
			if(m1<1) m1=1; if(m1>xnumx-2) m1=xnumx-2;
			m2=wn1[1]; if(marky[mm1]>(gy[m2]+gy[m2+1])/2.0) m2+=1;
			if(m2<1) m2=1; if(m2>ynumy-2) m2=ynumy-2;
			m3=wn1[2]; if(markz[mm1]>(gz[m3]+gz[m3+1])/2.0) m3+=1;
			if(m3<1) m3=1; if(m3>znumz-2) m3=znumz-2;
			/* Distances Calc */
			xkf=2.0/(gx[m1+1]-gx[m1-1]);
			ykf=2.0/(gy[m2+1]-gy[m2-1]);
			zkf=2.0/(gz[m3+1]-gz[m3-1]);
			/* Pressure gradients */
			x=xkf*(markx[mm1]-(gx[m1-1]+gx[m1])/2.0);
			y=xkf*(marky[mm1]-(gy[m2-1]+gy[m2])/2.0);
			z=xkf*(markz[mm1]-(gz[m3-1]+gz[m3])/2.0);
			/* P-SIGXX_EPSXX-Nodes num */
			/*    4----6 */
			/*   /|   /| */
			/*  / 5--/-7 */
			/* 0----2    */
			/* |    |    */
			/* 1----3    */
			p[0]=m3*xynumxy+m1*ynumy+m2;
			p[1]=p[0]+1;
			p[2]=p[0]+ynumy;
			p[3]=p[2]+1;
			p[4]=p[0]+xynumxy;
			p[5]=p[4]+1;
			p[6]=p[5]+ynumy;
			p[7]=p[6]+1;
			/* P derivatives */
			dpdx=xkf*((1.0-y)*(1.0-z)*(pr[p[2]]-pr[p[0]])+y*(1.0-z)*(pr[p[3]]-pr[p[1]])+(1.0-y)*z*(pr[p[6]]-pr[p[4]])+y*z*(pr[p[7]]-pr[p[5]]));
			dpdy=ykf*((1.0-x)*(1.0-z)*(pr[p[1]]-pr[p[0]])+x*(1.0-z)*(pr[p[3]]-pr[p[2]])+(1.0-x)*z*(pr[p[5]]-pr[p[4]])+x*z*(pr[p[7]]-pr[p[6]]));
			dpdz=zkf*((1.0-x)*(1.0-y)*(pr[p[4]]-pr[p[0]])+x*(1.0-y)*(pr[p[6]]-pr[p[2]])+(1.0-x)*y*(pr[p[5]]-pr[p[1]])+x*y*(pr[p[7]]-pr[p[3]]));
			/* Recalc velocity koefficients */
			vxkoef=(1000.0*GXKOEF-dpdx)/(2300.0*9.81);
			vykoef=(1000.0*GYKOEF-dpdy)/(2300.0*9.81);
			vzkoef=(1000.0*GZKOEF-dpdz)/(2300.0*9.81);
/*
printf("%ld %ld %ld    %ld %ld %ld   %e %e %e  %ld %d %e %e %e    %e %e %e   %e %e %e   %e %e %e",wn[0],wn[1],wn[2],m1,m2,m3,gx[m1],gy[m2],gz[m3],mm1,markt[mm1],markx[mm1],marky[mm1],markz[mm1],x,y,z,dpdx,dpdy,dpdz,vxkoef,vykoef,vzkoef);getchar();
*/
			if(vxkoef>2.0) vxkoef=2.0; if(vxkoef<-2.0) vxkoef=-2.0;
			if(vykoef>2.0) vykoef=2.0; if(vykoef<-2.0) vykoef=-2.0;
			if(vzkoef>2.0) vzkoef=2.0; if(vzkoef<-2.0) vzkoef=-2.0;
			/* Recalc velocity */
			vxwater=vywater*vxkoef;
			vzwater=vywater*vzkoef;
			vywater*=vykoef;
/*
printf("%ld %ld %ld %e %e %e  %ld %d %e %e %e    %e %e %e   %e %e %e   %e %e %e   %e %e %e ",m1,m2,m3,gx[m1],gy[m2],gz[m3],mm1,markt[mm1],markx[mm1],marky[mm1],markz[mm1],x,y,z,dpdx,dpdy,dpdz,vxkoef,vykoef,vzkoef,vxwater,vywater,vzwater);getchar();
*/
			}
		else
		/* Fluid in water */
			{
			vxwater=0;
			vywater=-ABSV(vywater);
			vzwater=0;
			}
		/**/
		}
/*
if(vywater>0) {printf("%e %e %e %e %e %e",tkpor,zmpor,vyfluid,vymelt,dmwamin,vywater);getchar();}
*/
	/**/
	/* Motion Calc ///////////////////////////////// */
	if (markmod==1)
		{
		/* Simple Velocity calc */
/*
		vxyzcalc1(markx[mm1],marky[mm1],markz[mm1],vxyz1,wn1);
*/
		v[0][0]=vxyz1[0]+vxwater; v[1][0]=vxyz1[1]+vywater; v[2][0]=vxyz1[2]+vzwater;
		}
	else
		{
		/* 4 Runge-Kutta koef calc */
/*
		vxyzcalc1(markx[mm1],marky[mm1],markz[mm1],vxyz1,wn1);
*/
		v[0][1]=vxyz1[0]+vxwater; v[1][1]=vxyz1[1]+vywater; v[2][1]=vxyz1[2]+vzwater;
		/**/
		vxyzcalc1(markx[mm1]+v[0][1]*timestep/2.0,marky[mm1]+v[1][1]*timestep/2.0,markz[mm1]+v[2][1]*timestep/2.0,vxyz1,wn1);
		v[0][2]=vxyz1[0]+vxwater; v[1][2]=vxyz1[1]+vywater; v[2][2]=vxyz1[2]+vzwater;
		/**/
		vxyzcalc1(markx[mm1]+v[0][2]*timestep/2.0,marky[mm1]+v[1][2]*timestep/2.0,markz[mm1]+v[2][2]*timestep/2.0,vxyz1,wn1);
		v[0][3]=vxyz1[0]+vxwater; v[1][3]=vxyz1[1]+vywater; v[2][3]=vxyz1[2]+vzwater;
		/**/
		vxyzcalc1(markx[mm1]+v[0][3]*timestep,marky[mm1]+v[1][3]*timestep,markz[mm1]+v[2][3]*timestep,vxyz1,wn1);
		v[0][4]=vxyz1[0]+vxwater; v[1][4]=vxyz1[1]+vywater; v[2][4]=vxyz1[2]+vzwater;
		/**/
		/* Vx,Vy,Vz calc after Runge-Kutta */
		v[0][0]=(v[0][1]+2.0*v[0][2]+2.0*v[0][3]+v[0][4])/6.0;
		v[1][0]=(v[1][1]+2.0*v[1][2]+2.0*v[1][3]+v[1][4])/6.0;
		v[2][0]=(v[2][1]+2.0*v[2][2]+2.0*v[2][3]+v[2][4])/6.0;
		}
	/**/
	/**/
	/**/
	/* Orthogonal motion only */
	if (outgrid==2)
		{
		if(markx[mm1]<0 || (double)(markx[mm1])>xsize) v[1][0]=v[2][0]=0;		
		if(marky[mm1]<0 || (double)(marky[mm1])>ysize) v[0][0]=v[2][0]=0;		
		if(markz[mm1]<0 || (double)(markz[mm1])>zsize) v[0][0]=v[1][0]=0;		
		}
	/**/
	/**/
	/**/
	/* Normal/Immobile markers */
	if(markt[mm1]<100)
		{
		/* Markers coming from the depth */
		if(marky[mm1]>zdeep && v[1][0]<0 && markk[mm1]<tdeep && markt[mm1]==10) markk[mm1]=tdeep;
		if(marky[mm1]>zdeep && v[1][0]<0 && markt[mm1]!=10) {markt[mm1]=10; markw[mm1]=0;}
		if(marky[mm1]>zdeep && markt[mm1]!=10 && markt[mm1]!=9) {markt[mm1]=9; markw[mm1]=0;}

/*
*/
		/* X,Y,Z calc after Runge-Kutta */
		markx[mm1]+=timestep*v[0][0]-dxgrid;
		marky[mm1]+=timestep*v[1][0];
		markz[mm1]+=timestep*v[2][0];
		/**/
		/* Markers coming from the depth */
		if((yold>vdeep && marky[mm1]<vdeep)|| (yold>zdeep && marky[mm1]<zdeep))
			{
			if(markk[mm1]<tdeep) markk[mm1]=tdeep;
			markt[mm1]=10;
			}
		/**/
		/* Out of grid marker reset */
		if(markx[mm1]<0 || marky[mm1]<0 || markz[mm1]<0 || (double)(markx[mm1])>xsize || (double)(marky[mm1])>ysize || (double)(markz[mm1])>zsize) 
			{
			markk[mm1]=0;
			markv[mm1]=0;
			}
		}
	else
		{
		/* Immobile markers */
		/* X,Y,Z calc after Runge-Kutta */
		markx[mm1]+=timestep*v[0][0]-dxgrid;
		marky[mm1]+=timestep*v[1][0];
		markz[mm1]+=timestep*v[2][0];
		/* Displacement calc */
		if(!ynres1)
			{
			dx=markx[mm1]-markk[mm1];
			dy=marky[mm1]-markd[mm1];
			dz=markz[mm1]-markw[mm1];
			ival=pow(dx*dx+dy*dy+dz*dz,0.5);
			/* Generate new markers, reset immobile marker positions y/n  */
			if(ival>stp100) ynres1=1;
			}
		}
	/**/
	/**/
	/**/
	/* Motion Calc ///////////////////////////////// */
	}
ynres[n1]=ynres1;
}
/* Check if there is a need for reset */
for(n1=0;n1<n0;n1++)
	{
	if(ynres[n1]) ynreset=1;
	}
/**/
/**/
/**/
/* Reset immobile markers */
if(ynreset)
	{
	if(printmod) printf("\n RESET immobile markers, Generate new markers OLD = %ld NEW = %ld \n",marknum,marknum1);
	for (mm1=0;mm1<marknum;mm1++)
	if(markt[mm1]>=100)
		{
		/* Create new normal marker */
		if(markx[mm1]>0 && marky[mm1]>0 && markz[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize && (double)(markz[mm1])<zsize)
			{
/*
printf("%e %e %e   %e %e %e   %e %e %e   %e %e %e",xold,yold,zold,markx[mm1],marky[mm1],markz[mm1],xnew,ynew,znew,v[0][0],v[1][0],v[2][0]);getchar();
*/
			markt[marknum1]=markt[mm1]-100;
			markx[marknum1]=markx[mm1];
			marky[marknum1]=marky[mm1];
			markz[marknum1]=markz[mm1];
			markk[marknum1]=0;
			markv[marknum1]=0;
			markft[marknum1]=0;
			markex[marknum1]=markex[mm1];
			marke[marknum1]=marke[mm1];
			markd[marknum1]=-1.0;
			markw[marknum1]=-1.0;
			markrr[marknum1]=0;
			markrn[marknum1]=0;
			/* Add aditional markers counter */
			marknum1++;
			}
		/* X,Y reset for immobile marker */
		markx[mm1]=markk[mm1];
		marky[mm1]=markd[mm1];
		markz[mm1]=markw[mm1];
		}
	}
/**/
/* New Mark num check */
if(printmod) printf("\n NEW Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1);
if(marknum1>MAXMRK) {printf("Space out in markx[]"); exit(0);}
/**/
/**/
/**/
/* Reset aditional markers */
mm1=0;
while(marknum1>marknum && mm1<marknum)
	{
	/* Reload marker */
	if(markt[mm1]<100 && (markx[mm1]<0 || marky[mm1]<0 || markz[mm1]<0 || (double)(markx[mm1])>xsize || (double)(marky[mm1])>ysize || (double)(markz[mm1])>zsize))
		{
		/* Decrease aditional markers counter */
		marknum1--;
		/* Type save */
		markt[mm1]=markt[marknum1];
		/* Temperature, Density, Water, Reset */
		markk[mm1]=0;
		markv[mm1]=0;
		markft[mm1]=0;
		markex[mm1]=markex[marknum1];
		marke[mm1]=marke[marknum1];
		markd[mm1]=-1.0;
		markw[mm1]=-1.0;
		markrr[mm1]=0;
		markrn[mm1]=0;
		/* X,Y reload */
		markx[mm1]=markx[marknum1];
		marky[mm1]=marky[marknum1];
		markz[mm1]=markz[marknum1];
		}
	/* Increase markers counter */
	mm1++;
	}
if(printmod) printf("\n Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1);
/* Set new marker number */
marknum=marknum1;
/**/
/* Clear melt extraction regions */
for (mm1=0;mm1<xnumx*znumz*5;mm1++)
	{
	bufv[mm1]=0;
	}
}
/* End Move markers by using Runge-Kutta method */




/* ro[],nu[] recalc after marker positions */
void ronurecalc()
{
/* Counters */
long int m1,m2,m3,m4;
long int m10,m20,m30;
int mm2,yn;
long int mm1,wn1[100];
double p_ga_in,p_pl_out,dx,dy,dz,ddx,ddy,ddz;
double mnu,mpb,mtk,mro,mcp,mht,mbb,mkt,mwa,dmwa,swt,swt1,ival;
double mro0,mcp0,mbb0;
double anu,aro,acp,akt,aht,abb,markwlev;
double wnu,wro,wcp,wkt,wht,wbb,dywa,markeii1,marksii1,markrii1,markdis1;
double eps1[100],wi1[100],ival2;
/**/
/**/
/* NEW_SFC*/
double local_sfc;
double dsx,dsz,dsd, dlwt;
double lambda[3]; // bathycentric coordinates
int  inl, inf;
/* NEW_SFC*/
/**/
/* Layering on sediments */
sedimnum++;
m1=(long int)(sedimnum/sedimcyc);
m2=((long int)(m1/2))*2;
if(m2==m1) yn=3; else yn=4;
/**/
/**/
/**/
/* ADD MARKERS TO THE v-CELLS ========================== */
/* Clear wt */
#pragma omp parallel for shared(sol0,sol1,ro,nu,nxx,nxy,nxz,nyz,cp,kt,ht,et,tk) private(m1) firstprivate(nodenum,nodenum2) schedule(static)
/*
*/
for (m1=0;m1<nodenum;m1++)
	{
	nu[nodenum+m1]=0;
	nxx[nodenum+m1]=0;
	nxy[nodenum+m1]=0;
	nxz[nodenum+m1]=0;
	nyz[nodenum+m1]=0;
	ro[nodenum+m1]=0;
	cp[nodenum+m1]=0;
	kt[nodenum+m1]=0;
	ht[nodenum+m1]=0;
	et[nodenum+m1]=0;
	tk[nodenum+m1]=0;
	rr[nodenum+m1]=0;
	sol0[m1]=0;
	sol1[m1]=0;
	sol0[nodenum+m1]=0;
	sol1[nodenum+m1]=0;
	sol0[nodenum2+m1]=0;
	sol1[nodenum2+m1]=0;
	}
/**/
/*
printf("%ld %e %e %e    %e %e %e",marknum,xsize,ysize,zsize,markx[0],marky[0],markz[0]); getchar();
*/
/* Add ro[] nu[] */
#pragma omp parallel for shared(sol0,sol1,ro,nu,nxx,rr,nxy,nxz,nyz,cp,pr,et,kt,tk,ht,markx,marky,markz,markt,markk,markd,markw,markto) private(local_sfc, dsx, dsz, dsd, dlwt,  lambda, inl, inf, p_ga_in,p_pl_out,marksii1,markrii1,markdis1,markeii1,ival2,wn1,wi1,eps1,m1,m2,m3,m4,m10,m20,m30,mm2,mm1,dx,dy,dz,ddx,ddy,ddz,mnu,mpb,mtk,mro,mcp,mht,mbb,mkt,mwa,dmwa,swt,swt1,ival,mro0,mcp0,mbb0,anu,aro,acp,akt,aht,abb,markwlev,wnu,wro,wcp,wkt,wht,wbb,dywa) firstprivate(yn,marknum,nodenum,xsize,ysize,zsize,eroslev,sedilev,waterlev,gx,gy,gz,markht,markkt,markkf,markcp,markkp,markn0,markn1,nubeg,nuend,zdeep,vdeep,nudeep,tdeep,dtdeep,drdeep,mnumy,mnumx,xynumxy,ynumy,nukoef,viscmod) schedule(static)
/*
*/
for (mm1=0;mm1<marknum;mm1++)
if(markx[mm1]>0 && marky[mm1]>0 && markz[mm1]>0 && markx[mm1]<xsize && marky[mm1]<ysize && markz[mm1]<zsize && markt[mm1]<50 && markk[mm1]>0)
	{
	/* Marker type */
	mm2=markt[mm1];
	/**/
	/* Check which surface model to use - only executed if no coupling with DAC external surface model NEW_DAC */
    if( (sfcmod==0) ){
        /* Erosion/sedimentation */
        if(mm2>1 && marky[mm1]<eroslev) mm2=markt[mm1]=0;
        if(mm2<2 && marky[mm1]>sedilev) mm2=markt[mm1]=yn;
        /* Water/air */
        if(mm2<2 && marky[mm1]>waterlev) mm2=markt[mm1]=1;
        if(mm2<2 && marky[mm1]<waterlev) mm2=markt[mm1]=0;
    }else{
        /*We use own surface*/

        /* check marker position relative to returned surface */
        /* since expensive, do only for shallower markers */
        if(    (nsfc>0)  && (marky[mm1] < 80000) ){


            //find sfc grid nodes that are left and in front
            //find sfc grid nodes that are left and in front
            inl = floor(markx[mm1]/xsstp);
            inf = floor(markz[mm1]/zsstp);

            // Add neighbour nodes elevation
            dlwt = 0;
            local_sfc = 0;
            //	Node 1
            dsx=markx[mm1]-(inl)*xsstp;
            dsz=markz[mm1]-(inf)*zsstp;
            dsd = (double)( sqrt( ( dsx*dsx + dsz*dsz ) ) );
            dlwt += 1/dsd;
            local_sfc += 1/dsd * sfc[inl*nsfcz+inf];
            //	Node 2
            dsx=markx[mm1]-(inl+1)*xsstp;
            dsz=markz[mm1]-(inf)*zsstp;
            dsd= (double)( sqrt( ( dsx*dsx + dsz*dsz ) ) );
            dlwt += 1/dsd;
            local_sfc += 1/dsd * sfc[(inl+1)*nsfcz+inf];
            //	Node 3
            dsx=markx[mm1]-(inl)*xsstp;
            dsz=markz[mm1]-(inf+1)*zsstp;
            dsd= (double)(sqrt ( ( dsx*dsx + dsz*dsz ) ) );
            dlwt += 1/dsd;
            local_sfc += 1/dsd * sfc[inl*nsfcz+inf+1];
            //	Node 4
            dsx=markx[mm1]-(inl+1)*xsstp;
            dsz=markz[mm1]-(inf+1)*zsstp;
            dsd= (double)(sqrt ( ( dsx*dsx + dsz*dsz ) ) );
            dlwt += 1/dsd;
            local_sfc += 1/dsd * sfc[(inl+1)*nsfcz+inf+1];

            local_sfc = waterlev - local_sfc/dlwt;

            if(isnan(local_sfc)){local_sfc=waterlev;}
            if(isinf(local_sfc)){local_sfc=waterlev;}

          /* do conversion */
            /* Erosion/sedimentation */
            if(mm2>1 && marky[mm1]< local_sfc) mm2=markt[mm1]=0;
            if(mm2<2 && marky[mm1]> local_sfc) mm2=markt[mm1]=yn;

           /* Water/air*/
            if(mm2<2 && marky[mm1]>waterlev) mm2=markt[mm1]=1;
            if(mm2<2 && marky[mm1]<waterlev) mm2=markt[mm1]=0;
        }


    }
    /* END OF SURFACE*/	/* P, T parameters calc */
	prcalc1((double)(markx[mm1]),(double)(marky[mm1]),(double)(markz[mm1]),eps1,wn1);
	mpb=eps1[10]*1e-5;
	mtk=(double)(markk[mm1]);
	/**/
	/**/
	/* Temperature reset for water/air */
/*
if(mm2<2) mtk=markk[mm1]=273.0; 
printf("k %ld %d %e %e ",mm1,mm2,mtk,mpb);getchar();
*/
	/**/
	/**/
	/**/
	/* Remove hot wet mantle */
/*
	if(mm2==11 && mtk>1273.15) mm2=markt[mm1]=10; 
*/
	/**/
	/* Rocks to rock+melt transformation */
	melting1(mtk,mpb,mm1,eps1);
	mm2=markt[mm1];
/*
*/
	/**/
	/* Marker Properties */
	/* Viscosity calc */
	mnu=markv[mm1];
	if(mnu<=0 || timestep<=0) 
		{
		epscalc1(markx[mm1],marky[mm1],markz[mm1],1,eps1,wn1);
		markeii1=pow(0.5*(eps1[4]*eps1[4]+eps1[5]*eps1[5]+eps1[6]*eps1[6])+eps1[7]*eps1[7]+eps1[8]*eps1[8]+eps1[9]*eps1[9],0.5);
		marksii1=pow(0.5*(eps1[54]*eps1[54]+eps1[55]*eps1[55]+eps1[56]*eps1[56])+eps1[57]*eps1[57]+eps1[58]*eps1[58]+eps1[59]*eps1[59],0.5);
		markdis1=eps1[60];
		markrii1=1.0;
		mnu=viscalc(mtk,mpb,mm1,mm2,markeii1,marksii1,markrii1,markdis1,wn1);
		markv[mm1]=mnu;
/*
*/
		}
	/**/
	mro0=mro=dencalc1(mtk,mpb,mm2,eps1);
	/* Mantle Depletion */
	if(mm2>=9 && mm2<=14 && markex[mm1]>0) mro0=mro*=1.0-0.04*markex[mm1];
	mbb0=mbb=eps1[20];
	mcp0=mcp=markcp[mm2];
	mkt=(markkt[mm2]+markkf[mm2]/(mtk+77.0))*exp(markkp[mm2]*mpb);
	/* Enhanced conductivity due to hydrothermal circulation (Gregg et al., 2009) */
/*
if(mm2>1 && mm2<20) {printf("%ld %d %e %e %e %e    %e",mm1,mm2,markx[mm1],marky[mm1],markz[mm1],mtk,markto[mm1]);getchar();}
*/
	if(timesum>0 && mm2>=7 && mm2<=13 && markto[mm1]>waterlev && (marky[mm1]-markto[mm1])<circdepth && mtk<circtemp)
		{
		mkt+=3.0*(nusseltno-1.0)*exp(0.75*(2.0-(mtk-273.15)/(circtemp-273.15)-(marky[mm1]-markto[mm1])/circdepth));
		}
	mht=markht[mm2];
/*
	if(mm2==7 && mtk>1000.0) mcp0=mcp*=2.0;
{printf("%ld %e %e %e",mm1,mtk,mro,mkt);getchar();}
if(mm2==9){printf("k %ld %d %e %e   %e    %e %e %e %e %e %e ",mm1,mm2,mpb,mtk,markk[mm1],mnu,mbb,mcp,mro,mkt,mht);getchar();}
if(m10==0){printf("k %ld %d    %e %e %e    %ld %ld %ld   %e %e   %e    %e %e %e %e %e %e ",mm1,mm2,markx[mm1],marky[mm1],markz[mm1],m10,m20,m30,mpb,mtk,markk[mm1],mnu,mbb,mcp,mro,mkt,mht);getchar();}
*/
	/* Melted rocks */
	if (mm2>20) 
		{
		if (meltpart0(mtk,mpb,mm1,mm2,eps1));
			{
			mro0=mro=eps1[23];
			mbb0=mbb=eps1[20];
/*
printf("%ld %d %e %e %e %e %e",mm1,mm2,mtk,mpb,mnu,eps1[21],val1[mm1]); getchar();
*/
			mcp0=mcp=eps1[25];
			mkt=eps1[26];
	                /* Basaltic melt viscosity */
        	        /* Effective NU calc check, Schott and Schmeling, 1995 */
			if((mm2==27 || mm2==28) && eps1[21]>0.001)
				{
                		ival=1e+15*exp(2.5+pow((1.0-eps1[21])/eps1[21],0.48)*(1.0-eps1[21]));
				if(ival<mnu) 
					{
					mnu=ival;
					if(mnu<markn0[mm2]) mnu=markn0[mm2]; 
					if(mnu>markn1[mm2]) mnu=markn1[mm2];
					if(mnu<nubeg) mnu=nubeg; 
					if(mnu>nuend) mnu=nuend;
					markv[mm1]=mnu;
					}
				
				}
			}
		}
/*
*/
	/**/
	/* Thermodynamic database use for ro, Cp */
	if (densimod>=2)
	if(mm2>1 && mm2!=5 && mm2!=6 && mm2!=25 && mm2!=26)
		{
		/* Compute TD variables */
		tdbasecalc1(mtk,mpb,mm2,mm1,eps1);
		mro=eps1[41];
		mwa=eps1[42];
		mcp=eps1[43];
		mbb=eps1[44];
		/**/
		/* Density && Water wt% save */
		if(1==0 || timesum<=1e+11 || mwa<=0 || markd[mm1]<=0) 
			{
			markw[mm1]=mwa;
			markd[mm1]=mro;
			}
		else
			{
			/* Recompute rock density on the basis of water density */
			wro=1050.0;
			dmwa=markw[mm1]-mwa;
			mro=mro/(1.0+dmwa*1e-2*(mro/wro-1.0));
			markd[mm1]=mro;
			}
		/**/
		/* Water only use mode */
		if(densimod==3)
			{
			mro=mro0;
			mbb=mbb0;
			mcp=mcp0;
			}
/*
{printf("TD1 %d %d %e %e   %e %e   %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1]);getchar();}
{printf("TD2 %d %d %e %e   %e %e   %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1]);getchar();}
*/
		}
	/**/
	/**/
	/**/
	/* Water/Air account ====================================*/
	if(waterlev>=0 && mm2<2)
		{
		/* Water/Air proportions definition */
		/* Water lewel definition relatively to the bottom of marker */
		dywa=(gy[wn[1]+1]-gy[wn[1]])/2.0/((double)(mnumy));
		markwlev=((double)(marky[mm1])-(waterlev-dywa))/(2.0*dywa);
/*
{printf("WAT %ld %d %e %e %e",mm1,mm2,marky[mm1],dywa,markwlev);getchar();}
*/
		if(markwlev<1.0 && markwlev>0)
			{
			/* Air Properties */
			anu=viscalc(mtk,mpb,mm1,0,markeii1,marksii1,markrii1,markdis1,wn1);
			/* Min,Max NU limitation */
			aro=dencalc1(mtk,mpb,0,eps1);
			abb=eps1[20];
			acp=markcp[0];
			akt=markkt[0];
			aht=markht[0];
			/**/
			/* Water Properties */
			wnu=viscalc(mtk,mpb,mm1,1,markeii1,marksii1,markrii1,markdis1,wn1);
			/* Min,Max NU limitation */
			wro=dencalc1(mtk,mpb,1,eps1);
			wbb=eps1[20];
			wcp=markcp[1];
			wkt=markkt[1];
			wht=markht[1];
			/**/
			/* Effective  properties recalculation */
			mnu=wnu*markwlev+anu*(1.0-markwlev);
			mro=wro*markwlev+aro*(1.0-markwlev);
			mbb=wbb*markwlev+abb*(1.0-markwlev);
			mcp=wcp*markwlev+acp*(1.0-markwlev);
			mkt=wkt*markwlev+akt*(1.0-markwlev);
			mht=wht*markwlev+aht*(1.0-markwlev);
			}
		}
	/* End Water/Air account ====================================*/
	/**/
	/* Limit viscosity of the grid */
/*
	if(mnu>nuend) mnu=nuend;
	if(mnu<nubeg) mnu=nubeg;
*/
	/* Limit viscosity in the middle */
/*
	if((nuend<9e+23 || timesum==0) && timesum<1e+11 && markx[mm1]>xsize/2.0-5e+4 && markx[mm1]<xsize/2.0+5e+4) 
		{
		if(mnu>nuend) mnu=nuend;
		mnu/=pow(100.0,1.0-ABSV(xsize/2.0-markx[mm1])/5e+4); 
		if(mnu<nubeg) mnu=nubeg;
		}
*/
	/* Limit viscosity */
	/* Lower boundary */
	/* Left,Right boundary */
	/* Increase density */
/*
	if(marky[mm1]>vdeep && mnu>nudeep) mnu=nudeep; 
	if(marky[mm1]>zdeep && markk[mm1]<tdeep-dtdeep) mro+=drdeep; 
	if(markx[mm1]<1e+4  && mnu>nudeep) mnu=nudeep; 
	if(markx[mm1]>xsize-1e+4  && mm2>1 && mnu<nuend && mtk<1000.0) mnu=nuend; 
	if(marky[mm1]>vdeep) mnu=nudeep; 
	if(marky[mm1]>vdeep && mnu>nudeep) mnu=nudeep; 
*/
	/**/
	/* Cell num calc for A-D Diagonal Corners of cur marker */
/*
	m10=m1serch(markx[mm1]);
	m20=m2serch(marky[mm1]);
	m30=m3serch(markz[mm1]);
*/
	m10=wn1[0];
	m20=wn1[1];
	m30=wn1[2];
	dx=gx[m10+1]-gx[m10];
	dy=gy[m20+1]-gy[m20];
	dz=gz[m30+1]-gz[m30];
	/**/
	/* Add nodes around current cell */
ival=0;
	for (m3=m30;m3<=m30+1;m3++)
		{
		/* Normalized Z-distance calc */
		ddz=(gz[m3]-markz[mm1])/dz;
		ddz=1.0-ABSV(ddz);
		for (m1=m10;m1<=m10+1;m1++)
			{
			/* Normalized X-distance calc */
			ddx=(gx[m1]-markx[mm1])/dx;
			ddx=1.0-ABSV(ddx);
			for (m2=m20;m2<=m20+1;m2++)
				{
				/* Normalized Y-distance calc */
				ddy=(gy[m2]-marky[mm1])/dy;
				ddy=1.0-ABSV(ddy);
				/* Wt calc */
				swt=ddx*ddy*ddz;
/*
ival+=swt;
if(m10==0 && m30==0) {printf("%ld  %ld %ld %ld  %e %e %e   %e %e %e   %e %e %e   %e %e",mm1,m1,m2,m3,markx[m1],marky[m2],markz[m3],dx,dy,dz,ddx,ddy,ddz,swt,ival);getchar();}
*/
				/**/
				/* Node num */
				m4=m3*xynumxy+m1*ynumy+m2;
				/**/
				/* Add viscosity for basic nodes */
				if(ddx>nukoef && ddy>nukoef && ddz>nukoef)
					{
					/* Add Viscosity */
					if(viscmod==0) ival2=mnu*swt;
					if(viscmod==1) ival2=log(mnu)*swt;
					if(viscmod==2) ival2=1.0/mnu*swt;
					nu[nodenum+m4]+=ival2;
					sol1[m4]+=swt;
					}
				/**/
if(0==0)
{
				/* Add viscosity for SIGxx, SIGyy, SIGzz (cell center) */
				if(m1>m10 && m2>m20 && m3>m30)
					{
					/* Add Viscosity */
					swt1=8.0*(0.5-ABSV(ddx-0.5))*(0.5-ABSV(ddy-0.5))*(0.5-ABSV(ddz-0.5));
					if(viscmod==0) ival2=mnu*swt1;
					if(viscmod==1) ival2=log(mnu)*swt1;
					if(viscmod==2) ival2=1.0/mnu*swt1;
					nxx[nodenum+m4]+=ival2;
					sol0[nodenum+m4]+=swt1;
					}
				/**/
				/* Add viscosity for SIGxy (z-edges) */
				if(m3==m30 && ddx>0.5 && ddy>0.5)
					{
					/* Add Viscosity */
					swt1=8.0*(ddx-0.5)*(ddy-0.5)*(0.5-ABSV(ddz-0.5));
					if(viscmod==0) ival2=mnu*swt1;
					if(viscmod==1) ival2=log(mnu)*swt1;
					if(viscmod==2) ival2=1.0/mnu*swt1;
					nxy[nodenum+m4]+=ival2;
					sol1[nodenum+m4]+=swt1;
					}
				/**/
				/* Add viscosity for SIGxz (y-edges) */
				if(m2==m20 && ddx>0.5 && ddz>0.5)
					{
					/* Add Viscosity */
					swt1=8.0*(ddx-0.5)*(0.5-ABSV(ddy-0.5))*(ddz-0.5);
					if(viscmod==0) ival2=mnu*swt1;
					if(viscmod==1) ival2=log(mnu)*swt1;
					if(viscmod==2) ival2=1.0/mnu*swt1;
					nxz[nodenum+m4]+=ival2;
					sol0[nodenum2+m4]+=swt1;
					}
				/**/
				/* Add viscosity for SIGyz (x-edges) */
				if(m1==m10 && ddy>0.5 && ddz>0.5)
					{
					/* Add Viscosity */
					swt1=8.0*(0.5-ABSV(ddx-0.5))*(ddy-0.5)*(ddz-0.5);
					if(viscmod==0) ival2=mnu*swt1;
					if(viscmod==1) ival2=log(mnu)*swt1;
					if(viscmod==2) ival2=1.0/mnu*swt1;
					nyz[nodenum+m4]+=ival2;
					sol1[nodenum2+m4]+=swt1;
					}
}
				/**/
				/* Add other properties */
				if(swt>0)
					{
					/* Add Properties */
					ival2=mro*swt;
					ro[nodenum+m4]+=ival2;
					ival2=mcp*mro*swt;
					cp[nodenum+m4]+=ival2;
					ival2=mkt*swt;
					kt[nodenum+m4]+=ival2;
					ival2=mht*swt;
					ht[nodenum+m4]+=ival2;
					ival2=mbb*swt;
					et[nodenum+m4]+=ival2;
					ival2=mtk*swt;
					tk[nodenum+m4]+=ival2;
					ival2=markrr[mm1]*swt;
if((mm2==7 || mm2==8 || mm2==27 || mm2==28) && ival2==0) ival2=5e-3*2.0/(2.0*asin(1.0));
					rr[nodenum+m4]+=ival2;
					sol0[m4]+=swt;
					}
/*
if(swt<0){printf("%ld  %e %e %e   %ld %ld %ld %ld  %e %e %e swt=%e ",mm1,dx,dy,dz,m1,m2,m3,m4,ddx,ddy,ddz,swt);getchar();}
if(swt<0){printf("%e %e  %e %e  %e %e ",gx[m1],gx[m1+1],gy[m2],gy[m2+1],gz[m3],gz[m3+1]);getchar();}
*/
				}
			}
		/* Use Zmax of marker corner */
		}
	}
/**/
/**/
/* Recalc ro[] nu[] */
#pragma omp parallel for shared(sol0,sol1,ro,nu,nxx,nxy,nxz,nyz,cp,kt,ht,et,tk) private(m3) firstprivate(nodenum,nodenum2,nubeg,nuend,viscmod) schedule(static)
/*
*/
for (m3=0;m3<nodenum;m3++)
	{
	/* Recompute viscosity for basic nodes */
	if(sol1[m3]>0)
		{
		nu[m3]=nu[nodenum+m3]/sol1[m3];
		if(viscmod==1) nu[m3]=exp(nu[m3]);
		if(viscmod==2) nu[m3]=1.0/nu[m3];
		}
	/**/
	/* Recompute viscosity SIGxx, SIGyy, SIGzz (cell center) */
	if(sol0[nodenum+m3]>0)
		{
		nxx[m3]=nxx[nodenum+m3]/sol0[nodenum+m3];
		if(viscmod==1) nxx[m3]=exp(nxx[m3]);
		if(viscmod==2) nxx[m3]=1.0/nxx[m3];
		}
	/**/
	/* Recompute viscosity for SIGxy (z-edges) */
	if(sol1[nodenum+m3]>0)
		{
		nxy[m3]=nxy[nodenum+m3]/sol1[nodenum+m3];
		if(viscmod==1) nxy[m3]=exp(nxy[m3]);
		if(viscmod==2) nxy[m3]=1.0/nxy[m3];
		}
	/**/
	/* Recompute viscosity for SIGxz (y-edges) */
	if(sol0[nodenum2+m3]>0)
		{
		nxz[m3]=nxz[nodenum+m3]/sol0[nodenum2+m3];
		if(viscmod==1) nxz[m3]=exp(nxz[m3]);
		if(viscmod==2) nxz[m3]=1.0/nxz[m3];
		}
	/**/
	/* Recompute viscosity for SIGyz (x-edges) */
	if(sol1[nodenum2+m3]>0)
		{
		nyz[m3]=nyz[nodenum+m3]/sol1[nodenum2+m3];
		if(viscmod==1) nyz[m3]=exp(nyz[m3]);
		if(viscmod==2) nyz[m3]=1.0/nyz[m3];
		}
	/**/
	/* Recompute other properties */
	if(sol0[m3]>0)
		{
		ro[m3]=ro[nodenum+m3]/sol0[m3];
		cp[m3]=cp[nodenum+m3]/sol0[m3]/ro[m3];
		kt[m3]=kt[nodenum+m3]/sol0[m3];
		ht[m3]=ht[nodenum+m3]/sol0[m3];
		et[m3]=et[nodenum+m3]/sol0[m3];
		tk[m3]=tk[nodenum+m3]/sol0[m3];
		rr[m3]=rr[nodenum+m3]/sol0[m3];
		}
	if(nu[m3]<nubeg) nu[m3]=nubeg;
	if(nu[m3]>nuend) nu[m3]=nuend;
	sol0[m3]=0;
	sol1[m3]=0;
	sol0[nodenum+m3]=0;
	sol1[nodenum+m3]=0;
	sol0[nodenum2+m3]=0;
	sol1[nodenum2+m3]=0;
	}
/**/
/**/
/* Set Boundary iconditions for T */
if (printmod) printf("\n AVERAGE TEMPERATURE CORRECTION FOR BOUNDARY CONDITIONS ...\n");
tkrecalc();
/*
*/
if (printmod) printf("AVERAGE TEMPERATURE OK!\n");
/**/
/**/
/**/
/* ADD MARKERS TO THE v-CELLS ========================== */
}
/* End ro[],nu[] recalc after marker positions */



/* Number of nearest left X-surface find */
long int m1serch(double x)
/* x - X coordinate */
{
/* Variables */
long int m1,m10=0,m11=xnumx-1;
/**/
/* Serch cycle */
do
	{
	m1=(m10+m11)/2;
	if (gx[m1]>x) m11=m1; else m10=m1;
	}
while((m11-m10)>1);
if(m10>xnumx-2) m10=xnumx-2;
/*
if(x<gx[m10] || x>gx[m10+1]) {printf("XXX %ld %ld %ld  %e %e  %e ",m10,m11,m1,gx[m10],gx[m11],x); getchar();}
*/
return m10;
}
/* Number of nearest left X-surface find */



/* Number of nearest upper Y-surface find */
long int m2serch(double y)
/* y - Y coordinate */
{
/* Variables */
long int m2,m20=0,m21=ynumy-1;
/**/
/* Serch cycle */
do
	{
	m2=(m20+m21)/2;
	if (gy[m2]>y) m21=m2; else m20=m2;
	}
while((m21-m20)>1);
if(m20>ynumy-2) m20=ynumy-2;
/*
if(y<gy[m20] || y>gy[m20+1]) {printf("YYY %ld %ld %ld  %e %e  %e ",m20,m21,m2,gy[m20],gy[m21],y); getchar();}
*/
return m20;
}
/* Number of nearest upper Y-surface find */



/* Number of nearest frontal Z-surface find */
long int m3serch(double z)
/* z - Z coordinate */
{
/* Variables */
long int m3,m30=0,m31=znumz-1;
/**/
/* Serch cycle */
do
	{
	m3=(m30+m31)/2;
	if (gz[m3]>z) m31=m3; else m30=m3;
	}
while((m31-m30)>1);
if(m30>znumz-2) m30=znumz-2;
/*
if(z<gz[m30] || z>gz[m30+1]) {printf("ZZZ %ld %ld %ld  %e %e  %e ",m30,m31,m3,gz[m30],gz[m31],z); getchar();}
*/
return m30;
}
/* Number of nearest frontal Z-surface find */




/* Nu calc after reological equation */
/* P-T-stress dependent rheology without/with brittle/ductile transition */
/* Reological equations */
/* Stress>SScr */
/* Power law dislocation creep: SSii={NU0*EEii*exp[(E+PV)/RT]}^(1/n) */
/* Effective viscosity: NU=1/2*{NU0*exp[(E+PV)/RT]}^(1/n)*EEii^[(1-n)/n] */
/* Stress<SScr */
/* Newtonian diffusion   creep: SSii=NU1*EEii*exp[(E+PV)/RT] */
/* Effective viscosity: NU=NU0/2*exp[(E+PV)/RT] */
/* NU1=NU0/SScr^(n-1) */
/* SScr - dislocation, diffusion transition stress */
/* SSii - second invariant of deviatoric stress tensor */
/* EEii - second invariant of strain rate tensor */
/* E - activation energy, J */
/* V - activation volume, J/bar */
/* R - gase constant 8.314 J/K */
/* Viscosity NU  calc after reological equations */
/* NU=SSii/(2*EEii) */
/* Brittle - Ductile transition */
/* sbrit=MINV(0.85e+5*pb,60e+6+0.6e+5*pb)*lambda;  (Schott & Scmeling, 1998) */
/* sbrit=MINV(0.667e+5*pb,51.2e+6+0.512e+5*pb)*lambda; (Brace & Kohlsstedt, 1980) */
double viscalc(double mtk, double mpb, long int mm1, int mm2, double markeii1, double marksii1, double markrii1, double markdis1, long int *wn1)
/* mtk - T, K */
/* mpb - P, bar */
/* x,y - XY location of point for Vx,Vy calc */
/* mm1 - Marker number */
/* mm2 - rock type */
{
/* Val buffer */
double y,e,n,rt=8.314*mtk,k1,e1,epsin,sduct,sbrit,nueff,smin,smax,nmin,nmax,strain,abrit,bbrit,nubrit,nunewt,nupowl;
/* Reological Eq par */
double lamb,mpb1,mnu2;
/* Counters */
long int m1;
int n1,n2;
/* Buffers */
double k1min,A,sig0,p,q,Pepsin,nupeierls,sigmin,sigmax,siginnew,sigin,mpf;
/* Bercovici flow law  */
int i,j;
double dissip,r,r0,dtsum,drlim,dsiginlim,dt,dt0;
double Adisl,Edisl,Vdisl,AA,Adiff,Ediff,Vdiff,m,BB;
double pival,phi1,phi2,eta,gammaI,Eggrow,Aggrow,Vggrow,GG,rscale,Gfac,GI,f0,fG,fI,Psi,DlnrDt;
double epsinpowl,siginpowl,ETA,sigin0,sigin1,sigin2,dsigin0,dsigin1,dsigin2;
/**/
/**/
/*
printf("%e %e",marky[mm1],mpb);getchar();
*/
/* Calc effective strain rate after second strain rate Tenzor invariant EEii=(1/2SUM(EPSik^2))^(1/2) */
epsin=markeii1;
sigin=marksii1;
/**/
/* Brittle "viscosity" calc: sbrit=A(eps)+LAM*B(eps)*P ------------------- */
/**/
/* Elimination of plastic strain outside of serpentine stability field */
/* Fluid pressure */
mpf=MAXV(mpb,(marky[mm1]-waterlev)*GYKOEF*1000.0*1e-5);
/* Depth, m */ 
y=mpf*3.0;
if(y>63000.0)
        {
        k1=1013.17699-0.060387633e-3*y-0.004289442e-6*y*y;
        }
else
        {
        k1=751.490422+6.00773668e-3*y-0.034690759e-6*y*y;
        }
/* No strain weakening outside of Serpentine stability */
/*
if(mtk>k1) marke[mm1]=0;
*/
/**/
/* Elimination of plastic strain at the high temperature */
/*
if(mtk>1000.0) marke[mm1]=0;
*/
/**/
/* Calc Second strain Tenzor invariant GAMMAii=(1/2SUM(EPSik^2))^(1/2) */
strain=ABSV(marke[mm1]);
/* Drop plastic strength on the melt pathway */
if(mm2>1 && mm2<20 && bufv[wn1[2]*xnumx+wn1[0]]>0 && marky[mm1]<bufv[wn1[2]*xnumx+wn1[0]] && marky[mm1]>bufv[wn1[2]*xnumx+wn1[0]]-dydike) strain+=lambmlt;
/*
*/
/**/
/* A,B coefficients calc depending on integral strain */
lamb=1.0; 
abrit=marka0[mm2]; 
bbrit=markb0[mm2]; 
if(strain>marke1[mm2])
	{
	abrit=marka1[mm2]; 
	bbrit=markb1[mm2];
	}
else
	{
	if(strain>marke0[mm2] && marke1[mm2]>marke0[mm2]) 
		{
		abrit=marka0[mm2]+(marka1[mm2]-marka0[mm2])*(strain-marke0[mm2])/(marke1[mm2]-marke0[mm2]);
		bbrit=markb0[mm2]+(markb1[mm2]-markb0[mm2])*(strain-marke0[mm2])/(marke1[mm2]-marke0[mm2]);
		}
	}
/* Fluid pressure subtract */
mpf=mpb*1e+5-(marky[mm1]-waterlev)*GYKOEF*1000.0;
/* Inverted value of Brittle "Mohr-Coulomb" viscosity calc */
nubrit=0;
if(epsin>0 && (bbrit>0 || abrit*lamb>0) && mm2>1)
	{
	/* Brittle strength calc */
	sbrit=abrit+bbrit*mpf*lamb;
	/* Negative pressure => mode 1 plasticity = open cracks */
	if(mpf<0) sbrit=abrit+mpf;
/*
*/
	/* Limit stress */
/*
	if (sbrit>marks1[mm2]) sbrit=marks1[mm2];
*/
	/* Brittle viscosity */
	if(sbrit>0)
		{
		nubrit=1.0/(0.5*sbrit/epsin/markrii1);
		}
	else
		{
		nubrit=1.0/markn0[mm2];
		sbrit=0;
		}
/*
if(mm2==9 || mm2==10) {printf("%ld %d  %e %e   %e %e %e   %e %e",mm1,mm2,mtk,mpb,mpf,epsin,markrii1,sbrit,nubrit);getchar();}
*/
	}
/* End Brittle "viscosity" calc: sbrit=A(eps)+LAM*B(eps)*P ------------------- */
/**/
/**/
/**/
/* Peierls plasticity-creep mechanism, data from Katayama and Karato, 2008 */
nupeierls=0;
if(1==0 && marks1[mm2]>0 && marks1[mm2]<1e+12 && mm2>1 && mm2<20 && sigin>1e+8 && mtk<1373.0 && epsin>0) 
	{
	n1=9;
	A=pow(10.0,7.8)*1e-12;/* Constant f(p,q), 1/s/MPa^2 */ 
	sig0=marks1[n1]; /* 9.1e+9 Dry Peierls stress at 0 K f(p,q), MPa Evans & Goetze, 1979 */
	if(0==1 && mm2!=9 && mm2!=10 && mm2!=14)
		{
		sig0=2.9e+9; /* Wet Peierls stress at 0 K f(p,q), MPa Katayama & Karato, 2008 */
		n1=11;
		}
	/* Using bisection */
	sigmin=1e+6;
	sigmax=sig0*0.9999;
	k1min=A*pow(sigmin,2.0)*exp(MAXV(-100.0,-(markdh[n1]+markdv[n1]*mpb)/rt*pow(1.0-sigmin/sig0,2.0)));
	for(n2=0;n2<15;n2++)
		{
		siginnew=(sigmax+sigmin)/2.0;
/*
		siginnew=pow(sigmax*sigmin,0.5);
*/
		k1=A*pow(siginnew,2.0)*exp(MAXV(-100,-(markdh[n1]+markdv[n1]*mpb)/rt*pow(1.0-siginnew/sig0,2.0)));
		if((k1<epsin && k1min<epsin) || (k1>epsin && k1min>epsin)) 
			{
			sigmin=siginnew;
			}
		else
			{
			sigmax=siginnew;
			}
/*
	p=1.0;
	q=2.0;
	k1=A*pow(siginnew,2.0)*exp(-(markdh[n1]+markdv[n1]*mpb)/rt*pow(1.0-pow(siginnew/sig0,p),q));
printf("%ld %d   %e %e   %e %e   %e %e",mm1,mm2,mtk,mpb,sigin,epsin,siginnew,1.0/nupeierls);getchar();
printf("%d   %ld %d   %e %e   %e %e   %e %e  %e %e     %e %e",n2,mm1,mm2,mtk,mpb,sigin,epsin,siginnew,k1,k1min,k1max,sigmin,sigmax);getchar();
*/
		}
		
	nupeierls=1.0/(0.5*siginnew/epsin);
/*
	if(nupeierls>1e-23 && (markx[mm1]<5e+4 || markx[mm1]>xsize-5e+4)) nupeierls=1e-23;
*/
	}
/*
printf("%ld %d  %e %e  %e %e  %e %e %e %e",mm1,mm2,x,y,mtk,mpb,nubrit,nunewt,nupowl,nueff);getchar();
*/
/**/
/**/
/**/
/* Ductile viscosity calc -------------------------------------------*/
/* Inverted value of newtonian NU set */
nunewt=0;
/**/
/* Inverted value of power-low NU set */
nupowl=0;
/**/
if(0==0) 
{
/* Check for the presence of ductile rheology */
if (marknu[mm2])
	{
	/* A)  Simple Newtonian rheology */
	/* Newtonian creep: SSii=NU0*2.0*EEii */
	/* Effective viscosity: NU=NU0 */
	/* Effective viscosity member in Stoks: NUs=NU */
	if(markdh[mm2]==0 && markdv[mm2]==0 && (markss[mm2]==0 || markmm[mm2]==1.0))
		{
		/* Inverted value of newtonian NU calc */
		nunewt=1.0/marknu[mm2];
		}
	/**/
	/**/
	/**/
	/* B)  P-T dependent, stress independent Newtonian rheology */
	/* Newtonian diffusion creep: SSii=NU0*EEii*exp[(E+PV)/RT] */
	/* Effective viscosity: NU=NU0*exp[(E+PV)/RT] */
	if((markdh[mm2]!=0 || markdv[mm2]!=0) && (markss[mm2]==0 || markmm[mm2]==1.0))
		{
		/* Inverted value of newtonian NU calc */
		e1=(markdh[mm2]+markdv[mm2]*mpb)/rt;
		if(e1>150.0) e1=150.0;
		/* Test creep Moresi & Solomatov (1995): SSii=NU0*exp[-a*(T-T0)] */
		if(markdh[mm2]<0 && markdv[mm2]>0) 
			{
			e1=markdh[mm2]*(mtk-markdv[mm2]);
			if(e1<-150.0) e1=-150.0;
			}
		/* Test creep Turkotte & Schubert(1982): SSii=NU0*exp[E/RTo(1-(T-T0)/T0)] */
		if(markdh[mm2]<0 && markdv[mm2]<0) 
			{
			e1=(-markdh[mm2])*(1.0-(mtk-(-markdv[mm2]))/(-markdv[mm2]))/8.314/(-markdv[mm2]);
			if(e1>150.0) e1=150.0;
			}
		nunewt=1.0/(marknu[mm2]*exp(e1));
		}
	/**/
	/**/
	/**/
	/* C)  P-T independent, stress dependent rheology without/with brittle/ductile transition */
	/* Stress>SScr */
	/* Power law creep: SSii={NU0*EEii}^(1/n) */
	/* Effective viscosity: NU=1/2*NU0^(1/n)*EEii^[(1-n)/n] */
	/* Effective viscosity member in Stoks: NUs=NU/n */
	/* Stress<SScr */
	/* Newtonian creep: SSii=NU1*EEii */
	/* Effective viscosity: NU=NU1/2 */
	/* Effective viscosity member in Stoks: NUs=NU */
	/* NU1=NU0/SScr^(n-1) */
	if((markdh[mm2]==0 && markdv[mm2]==0) && markss[mm2]!=0 && markmm[mm2]!=1.0)
		{
		/* Koef for stress independent creep NU1 calc */
		k1=marknu[mm2]/pow(markss[mm2],markmm[mm2]-1.0);
		/**/
		/* Inverted value of newtonian NU calc */
		nunewt=1.0/(0.5*k1);
		/**/
		/* Ductile power-low stress calc */
		sduct=pow(marknu[mm2]*epsin,1.0/markmm[mm2]);
		/**/
		/* Inverted value of power-low NU calc */
		if (epsin>0) nupowl=1.0/(0.5*sduct/epsin);
		}
	/**/
	/**/
	/**/
	/* D)  P-T-stress dependent rheology without/with brittle/ductile transition */
	/* Reological equations */
	/* Stress>SScr */
	/* Power law dislocation creep: SSii={NU0*EEii*exp[(E+PV)/RT]}^(1/n) */
	/* Effective viscosity: NU=1/2*{NU0*exp[(E+PV)/RT]}^(1/n)*EEii^[(1-n)/n] */
	/* Effective viscosity member in Stoks: NUs=NU/n */
	/* Stress<SScr */
	/* Newtonian diffusion   creep: SSii=NU1*EEii*exp[(E+PV)/RT] */
	/* Effective viscosity: NU=NU0/2*exp[(E+PV)/RT] */
	/* Effective viscosity member in Stoks: NUs=NU */
	/* NU1=NU0/SScr^(n-1) */
	if(marknu[mm2]>0 && markdh[mm2]>0 && (markdh[mm2]!=0 || markdv[mm2]!=0) && markss[mm2]!=0 && markmm[mm2]!=1.0)
		{
		/* T-P exponent for effective NU calc */
		e1=(markdh[mm2]+markdv[mm2]*mpb)/rt;
		if(e1>150.0) e1=150.0;
		e1=exp(e1);
		/**/
		/* Koef for stress independent creep NU1 calc */
		k1=marknu[mm2]/pow(markss[mm2],markmm[mm2]-1.0);
		/**/
		/* Inverted value of newtonian NU calc */
		nunewt=1.0/(0.5*k1*e1);
		/**/
		/* Ductile power-low stress calc */
		sduct=pow(marknu[mm2]*e1*epsin,1.0/markmm[mm2]);
		/**/
		/* Inverted value of power-low NU calc */
		if(epsin>0) nupowl=1.0/(0.5*sduct/epsin);
		}
	/* F)  P-T-stress dependent rheology of Mantle Karato & Wu (1993) */
	if(marknu[mm2]<0)
		{
		/* Diffusion viscosity calc for given grainsize, mm, Karato & Wu (1993) */
		mnu2=1.0; /* Grainsize, mm */ 
		nunewt=1.0/(0.5*pow(mnu2*1e-3/5e-10,2.5)*8e+10/marks1[mm2]*exp((markss[mm2]+mpb*markll[mm2])/rt));
		/* Dislocation viscosity calc, Karato & Wu (1993) */
		/* Effective viscosity1 calc */
		nupowl=0; 
		if(epsin>0) nupowl=1.0/(0.5*pow(epsin,1.0/markmm[mm2]-1.0)*8e+10/pow(-marknu[mm2],1.0/markmm[mm2])*exp((markdh[mm2]+markdv[mm2]*mpb)/rt/markmm[mm2]));
		}
	/* G)  P-T-stress dependent rheology of Mantle Bercovici (personal communication) */
	if(markdh[mm2]<0)
		{

/* Stress, strain rate and dissipation */
epsin=markeii1;
sigin=marksii1;
dissip=markdis1*markrn[mm1]; /* Ductile energy dissipation */
/* Grainsize parameters */
pival=2.0*asin(1.0);
r=r0=markrr[mm1];  /* Interface roughness or radius of curvature for pinned-state limit */
if(r0<=0) r=r0=markrr[mm1]=3e-3*2.0/pival;
dtsum=0; /* Cumulative time */
drlim=0.05; /* Limit for r change */
dsiginlim=1e-6; /* Limit for sigin change */
dt0=timestep; /* Timestep for grain growth */

/* Rheological parameters; using Hirth & Kohlstedt 2003 rheology data for olivine (dry only) */
/* Dislocation creep parameters */
n=markmm[mm2]; /* dislocation creep power-law exponent (could be n=3 also, which is what I often use for mathematical convenience) */
Adisl=marknu[mm2]; /* preexponential constant */ 
Edisl=-markdh[mm2]; /* 530kJ/mol activation energy for dislocation creep; uncertainty is ??4 */
Vdisl=markdv[mm2];  /* Hirth/Kohsltedt Table 2; for P=1 - 3GPa 14.e-6 most frequent, low as 6 high as 27 */
e1=(Edisl+mpb*Vdisl)/rt;
if(e1>150.0) e1=150.0;
e1=exp(-e1);
AA=Adisl*e1;  /* dislocation creep "compliance" */ 
/* Diffusion creep parameters */
m=3.0;  /* grain-size exponent */ 
Adiff=marks1[mm2]; /* preexponential constant */
Ediff=markss[mm2]; /* 375kJ/mol activation energy for diffusion creep uncertainty is ??50 */
Vdiff=markll[mm2]; /* activation volume 6.e-6 m^3/mol and uncertainty is ??4 */
e1=(Ediff+mpb*Vdiff)/rt;
if(e1>150.0) e1=150.0;
e1=exp(-e1);
BB=Adiff*e1; /* diffusion creep compliance */
BB=BB/pow(pival/2.0,m); /* the (pi/2)^m factor comes about b/c we're using r (pining body size */
	/* or interface roughness) instead of grainsize Rg and assuming a pinned state */
	/* where the Zener pinning factor  Z = 1- c(1-phi_i)(R_i/r)^2 ~ 0 */
        /* and thus R_i = r/sqrt(c(1-phi_i)); thus the avg is */
	/*	R = r/sqrt(c) * [phi1/sqrt(phi2) + phi2/sqrt(phi1)] */
        /* and using c=0.87 (from the Nature and EPSL paper etc), phi1=0.4 and phi2=0.6, */ 
	/* we get R=1.5707r, which  coincidentally is 1.5707 = pi/2 */    

/*
printf("A %ld %d %e %e %e   %e %e   %e %e   %e %e",mm1,mm2,epsin,sigin,dissip,mtk,mpb,dt0,r0,AA,BB); getchar();
*/

/* Update grainsize */
if(0==0)
{
j=0;
do
{
/* Grain growth parameters */
q=4.0; /* roughness coarsening exponent */ 
p=2.0; /* grain size coarsening exponent; we need this to scale the experimental growth rates */
phi1=0.4; phi2=1.0-phi1; /* volume fractions of Ol and Px; say Px is phase 1, Ol is phase 2; */
eta=3.0*phi1*phi2;  /* a simple function of how interface area density depends on phi1 and phi2 */
gammaI=1.0;  /* interface surface tension in Pa m */
Eggrow=300e+3; /* activation energy for grain growth; assuming standard for diffusion; 200 from Karato could be used also but not sure anyone including Shun believes it */
Aggrow=2e+4*pow(1e-6,p); /* pre-exponential factor; this is from Karato, repeated by Rozel; it should be redone one day;  not sure if there is something better */
Vggrow=Vdiff; /* activation volume; this is a total guess; assume it's similar to diff creep? */
e1=(Eggrow+mpb*Vggrow)/rt;
if(e1>150.0) e1=150.0;
e1=exp(-e1);
GG=Aggrow*e1;  /* grain-growth rate in units of s^(-1)*m^p */
/* interface coarsening is from comparison to Hiraga et al 2010 experiments and scaling arguments in Bercovici & Ricard 2013. */
rscale=1e-6; Gfac=100.0;   /* have used Gfac = 250 before; there is freedom/slop in adjusting this a bit */
GI=GG/Gfac*q/p*pow(r0,q-p);  /* interface coarsening GI is slower than GG, but also has different units this scaling is explained in Berco+Ricard 2013 */
/* damage partitioning fraction (how much deformational work goes into creating surface energy on interface) */
f0=1e-3; /* damage at 1000K; the Rozel paper suggest 2.e-1 which gives a HUGE damage number, */
          /* but it goes with the grain-growth activation  energy of 200.  for Eg=300 it is */
          /* usually much smaller, like 1e-3.  There is some freedome also in playing with this */
e1=2.0*pow(mtk/1000,2.9) - 2.0;
if(e1>150.0) e1=150.0;
e1=exp(-e1);
fG=f0*e1;  /* this T dependence is still from Rozel, assuming Eg=200; it needs to be adjusted for general Eg (will get this from Elvira) */
fI= fG;  /* assuming interface damage fraction is the same as the maximum for grain boundary damage. */
         /* In fact, if we include grain-boundary evolution eqns,  this fG is scaled again to vanish */
	 /* when we're in diffusion creep; but we don't include those eqns here. */
	 /* ALSO if we include the hysteresis effect (maybe one day) theh fI (interface damage) */
         /* will also be scaled to vanish but in dislocation creep; not worth doing yet */

/* Grain growth */
/* !!!! To look where you move relative to deformation regime boundary */
/*  Sometimes it's useful to get the field boundary grain-size, although this formula below */
/*   is really for the field-boundary interface roughness (b/c of how we redefined BB): */
/* Rf = (BB/AA/sigin^(n-1))^(1/m); */
/* Now for grain evolution laws; we're only going to use interface roughness in the pinned-state limit evolution equation (finally): */
/* !!!!!! Marker deformational work: */
/* Psi = 2*(AA*sigin^(n+1) + BB*sigin^2/r^m); */
/* !!!!!! Nodal deformational work: */ 
Psi = dissip;
/* DrDt = eta*GI/q/r0^(q-1) - fI*r0^2/gammaI/eta * Psi; */ 
DlnrDt = (eta*GI/q/pow(r0,q-1.0) - fI*r0*r0/gammaI/eta * Psi)/r0;
dt=dt0-dtsum; /* Remaining time for grain growth */
if(ABSV(DlnrDt*dt)>drlim) dt=dt*drlim/ABSV(DlnrDt*dt);
r=r0*exp(DlnrDt*dt);
if(r>1e-2) r=1e-2;
if(r<1e-9) r=1e-9;
/* Update r */
j=j+1;
dtsum=dtsum+dt;
r0=r;
}
while(dtsum<dt0 && j<50);
/* Update grainsize */
markrr[mm1]=r;
}

/*
printf("B %d %e %e %e",j,r,nunewt,nupowl); getchar();
*/

/* !!!!!! Viscosity calculation */
/* Diffusion creep viscosity */
nunewt=pow(r,m)/BB/2.0;
nupowl=0;
if(epsin>0)
{
/* Dislocation creep strain rate and viscosity */
/* Assuming that epsindisl=epsin */
epsinpowl=epsin;
siginpowl=pow(epsinpowl/AA,1.0/n); 
nupowl=siginpowl/2.0/epsinpowl;
if(0==0)
{
/* Resulting effective viscosity */
ETA=1.0/(1.0/nunewt+1.0/nupowl);
/* Viscosity check */
/*
if(ETA<markn0[mm2]) ETA=markn0[mm2]; 
if(ETA>markn1[mm2]) ETA=markn1[mm2];
if(ETA<nubeg) ETA=nubeg; 
if(ETA>nuend) ETA=nuend;
*/
/*
printf("C %e %e %e %e",r,nunewt,nupowl,ETA); getchar();
*/

/* Bisection limits */
/* First initial stress computing from viscosity */
sigin0=2.0*epsin*ETA;
sigin=sigin0;
/* Dislocation creep strain rate and viscosity */
epsinpowl=AA*pow(sigin,n); 
nupowl=sigin/2.0/epsinpowl;
/* Resulting effective viscosity */
ETA=1.0/(1.0/nunewt+1.0/nupowl);
/* Viscosity check */
/*
if(ETA<markn0[mm2]) ETA=markn0[mm2]; 
if(ETA>markn1[mm2]) ETA=markn1[mm2];
if(ETA<nubeg) ETA=nubeg; 
if(ETA>nuend) ETA=nuend;
*/
/* New stress */
sigin1=2.0*epsin*ETA;
/* First stress change */
dsigin0=sigin1-sigin0;
/* Second initial stress taken as sigin1 */
sigin=sigin1;
epsinpowl=AA*pow(sigin,n); 
nupowl=sigin/2.0/epsinpowl;
/* Resulting effective viscosity */
ETA=1.0/(1.0/nunewt+1/nupowl);
/* Viscosity check */
/*
if(ETA<markn0[mm2]) ETA=markn0[mm2]; 
if(ETA>markn1[mm2]) ETA=markn1[mm2];
if(ETA<nubeg) ETA=nubeg; 
if(ETA>nuend) ETA=nuend;
*/
/* New stress */
sigin2=2.0*epsin*ETA;
/* Second stress change */
dsigin1=sigin2-sigin1;

/*
printf("D %e %e %e   %e %e",sigin0,sigin1,sigin2,dsigin0,dsigin1); getchar();
*/

/* Bisection on the stress */
if(ABSV(dsigin1/sigin)>dsiginlim && ((dsigin0<0 && dsigin1>0) || (dsigin0>0 && dsigin1<0)))
for(i=1;i<50;i++)
    {
    /* Current stress value */
    sigin=(sigin0+sigin1)/2.0;
    /* Dislocation creep strain rate and viscosity */
    epsinpowl=AA*pow(sigin,n); 
    nupowl=sigin/2.0/epsinpowl;
    /* Resulting effective viscosity */
    ETA=1.0/(1.0/nunewt+1.0/nupowl);
    /* Viscosity check */
/*
    if(ETA<markn0[mm2]) ETA=markn0[mm2]; 
    if(ETA>markn1[mm2]) ETA=markn1[mm2];
    if(ETA<nubeg) ETA=nubeg; 
    if(ETA>nuend) ETA=nuend;
*/
    /* New stress */
    sigin2=2.0*epsin*ETA;
    /* Current stress change */
    dsigin2=sigin2-sigin;
    /* Exit iteration */
    if(ABSV(dsigin2/sigin)<dsiginlim) break;
    /* Changing bisection limits */
    if((dsigin2>0 && dsigin0>0) || (dsigin2<0 && dsigin0<0))
        {
        sigin0=sigin;
        }
    else
        {
        sigin1=sigin;
        }
    }
/*
printf("E %d %e %e %e %e  %e %e  %e",i,sigin,sigin0,sigin1,sigin2,dsigin0,dsigin2,nupowl); getchar();
*/
}
/* Inverted dislocation creep viscosity */
if(sigin>sbrit)
    {
    /* Dislocation creep viscosity corrected for yield stress */
    nupowl=2.0*AA*pow(sbrit,n-1);
    }
    else
    {
    nupowl=1.0/nupowl;
    }
/*
*/
}
nunewt=1.0/nunewt;
/*
*/

		}
	}
/* End Ductile viscosity calc -------------------------------------------*/
}
/**/
/**/
/**/
/* Inverted value of effective viscosity calc, check */
/*
nueff=1.0/nubeg;
*/
nueff=nunewt+nupowl;
/*
if(mm2==9 || mm2==10) {printf("F %ld %d %e %e %e %e %e",mm1,mm2,nunewt,nupowl,nueff,nupeierls,nubrit); getchar();}
*/
if(nupeierls>nueff) nueff=nupeierls;
if(mm2==7 || mm2==8 || mm2==9 || mm2==10)
	{
	if(mtk>1300.0)
		{
		/* Reset formation time */
		markft[mm1]=0;
		}
	else
		{
		/* Save formation time */
		if(markft[mm1]<=0) markft[mm1]=timesum;
		}
	}
if(nubrit>nueff)
	{
	nueff=nubrit;
	/* Mark plastic yeilding */
	if(marke[mm1]==0) marke[mm1]=1e-20;
	marke[mm1]=ABSV(marke[mm1]);
	}
else
	{
	/* Stop plastic yeilding */
	marke[mm1]=-ABSV(marke[mm1]);
	}
/* Update ductile viscosity ratio */
if(markdh[mm2]<0 && timesum>1e+11) markrn[mm1]=(nueff-nunewt)/nueff;
/*
printf("%ld %d  %e %e  %e %e  %e %e %e %e",mm1,mm2,x,y,mtk,mpb,nubrit,nunewt,nupowl,nueff);getchar();
*/
/**/
/* Inverted Viscosity check */
if(nueff<=0) nueff=1.0/markn1[mm2];
/**/
/* Viscosity calc */
nueff=1.0/nueff;
/*
printf("A %ld %d  %e %e %e   %e %e  %e %e",mm1,mm2,markx[mm1],marky[mm1],markz[mm1],mtk,mpb,epsin,nueff);getchar();
printf("A %ld %d %e %e %e",mm1,mm2,mtk,mpb,nueff);getchar();
*/
/**/
/* Viscosity check */
if(nueff<markn0[mm2]) nueff=markn0[mm2]; 
if(nueff>markn1[mm2]) nueff=markn1[mm2];
if(nueff<nubeg) nueff=nubeg; 
if(nueff>nuend) nueff=nuend;
/*
*/
/**/
/*
if(strain>1.5 && nubrit>0 && x>1000000.0 && y<400000.0) {printf("%ld %d  %e %e  %e %e  %e %e %e %e  %e %e %e",mm1,mm2,x,y,mtk,mpb,strain,abrit,bbrit,lamb,1.0/nubrit,epsin,nueff); getchar();}
*/
/* Return calculated viscosity */
/*
printf("B   %ld %d %e %e %e",mm1,mm2,mtk,mpb,nueff);getchar();
*/
return nueff;
}
/* Nu calc after reological equation */




/* Thermodynamic database use for ro, Cp */
void tdbasecalc1(double mtk, double mpb, int mm2, long int mm1, double *eps1)
{
/* TD Database variables,  dTK,dPB - TK, PB step for tabulation in TD database */
double H0,H1,H2,H3,R0,R1,R2,R3,W0,W1,W2,W3,n,e;
/* Val Buffers */
int n1,n2,mm3;
double mhh0,mhh1,mdhh,maa,mwa,dmwa,wro,mro,mcp,mbb;
long int m1=wn[0];
double sy1,e1;
/**/
/* Reset TD variables */
eps1[40]=eps1[41]=eps1[42]=eps1[43]=eps1[44]=eps1[45]=0;
/**/
/* TD base type */
switch (mm2)
	{
	/* Sediments */
	case 2:
	case 3:
	case 4: 
	case 15: 
	mm3=5; break;
	/* Molten Sediments */
	case 22:
	case 23:
	case 24: 
	case 35: 
	mm3=6; break;
	/* Basalts */
	case 7: 
	case 16:
	case 17:
	case 18:
	mm3=7; break;
	/* Molten Basalt */
	case 27: 
	case 36:
	case 37:
	case 38:
	mm3=8; break;
	/* Gabbro */
	case 8: 
	mm3=3; break;
	/* Molten Gabbro */
	case 28: 
	mm3=4; break;
	/* Dry peridotite */
	case 9:
	case 12:
	case 14:
	case 10: 
	mm3=0; break;
	/* Wet peridotite */
	case 13:
	case 11: 
	mm3=1; break;
	/* Molten peridotite */
	case 34: 
	mm3=2; break;
	/* Unknown type */
	default: {printf("Unknown rock type for TD database %d",mm2); exit(0);}
	}
/* ABCD-4Cell Number */
e=(mtk-tkmin)/tkstp;
if(e<0) e=0;
if(e>(double)(tknum-1)) e=(double)(tknum-1);
n=(mpb-pbmin)/pbstp;
if(n<0) n=0;
if(n>(double)(pbnum-1)) n=(double)(pbnum-1);
n1=(int)(e);
if(n1>tknum-2) n1=tknum-2;
n2=(int)(n);
if(n2>pbnum-2) n2=pbnum-2;
/* e,n Calc */
e=(e-(double)(n1));
n=(n-(double)(n2));
/* Ro H values */
/* 0 2 */
/* 1 3 */
R0=td[n1  ][n2  ][mm3][0]*1000.0;
R1=td[n1  ][n2+1][mm3][0]*1000.0;
R2=td[n1+1][n2  ][mm3][0]*1000.0;
R3=td[n1+1][n2+1][mm3][0]*1000.0;
H0=td[n1  ][n2  ][mm3][1]*1000.0*4.1837;
H1=td[n1  ][n2+1][mm3][1]*1000.0*4.1837;
H2=td[n1+1][n2  ][mm3][1]*1000.0*4.1837;
H3=td[n1+1][n2+1][mm3][1]*1000.0*4.1837;
W0=td[n1  ][n2  ][mm3][2];
W1=td[n1  ][n2+1][mm3][2];
W2=td[n1+1][n2  ][mm3][2];
W3=td[n1+1][n2+1][mm3][2];
/* Ro calc by interpolation */
mro=((R0*(1.0-n)+R1*n)*(1.0-e)+(R2*(1.0-n)+R3*n)*e);
/* Water wt% calc by interpolation */
mwa=((W0*(1.0-n)+W1*n)*(1.0-e)+(W2*(1.0-n)+W3*n)*e);
/* Add porocity fluid */
/* Erosion surface */
if(marks0[mm2]>0 && marky[mm1]<zmpor && mtk<tkpor) 
	{
	dmwa=marks0[mm2]*(tkpor-mtk)/(tkpor-273.15)*(zmpor-MAXV(marky[mm1],sedilev))/(zmpor-sedilev);
	mwa+=dmwa;
	wro=1050.0;
	mro=mro/(1.0+dmwa*1e-2*(mro/wro-1.0));
/*
if(sy1>10000.0){printf("TD1 %d %d %e %e   %e %e   %e %e",mm2,mm3,mtk-273.15,sy1,mwa,mro,dmwa,zmpor);getchar();}
if(sy1>10000.0){printf("TD2 %d %d %e %e   %e %e   %e %e",mm2,mm3,mtk-273.15,sy1,mwa,mro,dmwa,zmpor);getchar();}
*/
	}
/* Cp calc by interpolation */
mcp=((H2-H0)*(1.0-n)+(H3-H1)*n)/tkstp;
if(mcp<1e+2) mcp=1e+2; else if(mcp>5e+4) mcp=5e+4;
/* Effective adiabatic betta=1/V*dV/dT=ro/T*[-dH/dP+V] calc by interpolation */
mbb=(2.0/(R1+R0)-(H1-H0)/pbstp/1e+5)*(1.0-e)+(2.0/(R3+R2)-(H3-H2)/pbstp/1e+5)*e;
mbb*=mro/mtk;
if(mbb<-1e-2) mbb=-1e-2; else if(mbb>1e-2) mbb=1e-2;
/* Effective compressibility term alpha=1/ro*d(ro)/dP calc by interpolation */
maa=(2.0/(R1+R0)*(R1-R0)*(1.0-e)+2.0/(R3+R2)*(R3-R2)*e)/pbstp/1e+5;
if(maa<0) maa=0;
/* Activation enthalpy recalc using enthalpy changes */
/* Current Enthalpy */
mhh1=((H0*(1.0-n)+H1*n)*(1.0-e)+(H2*(1.0-n)+H3*n)*e);
/* Pmin Enthalpy */
mhh0=(td[n1][0 ][mm3][1]*(1.0-e) + td[n1+1][0 ][mm3][1]*e)*1000.0*4.1837;
/* Enthalpy Difference calc */
mdhh=(mhh1-mhh0);
/*
{printf("TD1 %d %d %e %e   %e %e %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,mgg,mro,mwa,mcp,mbb);getchar();}
{printf("TD1 %d %d %e %e   %e %e %e %e %e %e   %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,mgg,mro,mwa,mcp,mbb,maa,mdhh,mhh1,mhh0);getchar();}
eps1[47]=mhh1;
eps1[48]=mhh0;
*/
/* Save TD variables */
eps1[41]=mro;
eps1[42]=mwa;
eps1[43]=mcp;
eps1[44]=mbb;
eps1[45]=maa;
eps1[46]=mdhh;
}
/* Thermodynamic database use for ro, Cp */





/* Melt fraction, latent heat calculation */
void meltpart11(double mtk, double mpb, int mm2, double *eps1)
/* mtk - T, K */
/* mpb - P, bar */
/* x,y - XY location of point for Vx,Vy calc */
/* mm1 - mark number */
/* mm2 - mark type */
/* yn  - type of calculation: 0 - Ro, 1 - Nu, 2 - Cp, 3 - kt */
{
/* Val buffer */
double xmelt=0,hlatent=0,ival,mpg,Fcpxout,Tcpxout;
long int m1;
double ykm=mpb*3e-3,ts=0,tl=0,tll=0;
/**/
/**/
/**/
/* Calculate melt fraction using marker type */
if (ykm>0 && ykm*1000.0<zdeep)
switch(mm2)
	{
	/* Sediments: latent heat 300 kJ/kg (Bittner & Schmeling, 1995) */
	case 15:
	case 35:
	case 3:
	case 4:
	case 5:
	case 23:
	case 24:
	case 25:
	/* Wet Solidus Temperature, Johannes, 1985, Poli & Schmidt, 2002 */
	if (ykm<36.0) 
		{
		ts=889.0+536.6/(ykm+1.609)+18.21/(ykm+1.609)/(ykm+1.609);
		}
	else
		{
		ts=831.3+2.0*ykm;
		}
	/* Dry Granite Liquidus, Johannes, 1985 */
	tl=1262.0+3.0*ykm;
	hlatent=300000.0;
	break;
	/**/
	/* Basalt, Gabbro: latent heat 380 kJ/kg (Bittner & Schmeling, 1995) */
	case 6:
	case 26:
	case 36:
	case 37:
	case 38:
	case 16:
	case 17:
	case 18:
	case 7:
	case 8:
	case 27:
	case 28: 
	/* Wet solidus, Schmidt & Poli, 1998  */
/*
	if (ykm<48.0) 
		{
		ts=972.6-2111.0/(ykm+10.63)+70033.0/(ykm+10.63)/(ykm+10.63);
		}
	else
		{
		ts=935.4+0.1162*ykm+0.006937*ykm*ykm;
		}
*/
	/* Dry Toleitic Basalt solidus, Hess, 1989 */
	ts=1327.15+3.02*ykm;
	/* Dry Toleitic Basalt Liquidus, Hess, 1989 */
	tl=1423.15+3.5*ykm;
	hlatent=380000.0;
	break;
	/**/
	/* Dry Peridotite: latent heat 400 kJ/kg Turcotte & Schubert, 1982, p.171 */
	case 11:
	case 14:
	case 34:
	/* Wet solidus, Schmidt & Poli, 1998  */
	if (ykm<72.0) 
		{
		ts=1239.8+1493.0/(ykm+9.701);
		}
	else
		{
		ts=1266.3-0.3948*ykm+0.003893*ykm*ykm;
		}
	/* Dry Peridotite Liquidus, Hess, 1989 */
	tl=2073.15+3.8*ykm;
	hlatent=400000.0;
	break;
	/**/
	/* Wet Peridotite: latent heat 400 kJ/kg Turcotte & Schubert, 1982, p.171 */
	case  9:
	case 10:
	case 12:
	case 29:
	case 30:
	case 32:
	/* Dry mantle melting Katz et al., 2003 */
	mpg=mpb*1e-4;
	ts=273.15+1085.7+132.9*mpg-5.1*mpg*mpg;
	tll=273.15+1475.0+80.0*mpg-3.2*mpg*mpg;
	tl=273.15+1780.0+45.0*mpg-2.0*mpg*mpg;
	xmelt=0;
	if(mtk>=ts)
		{
		if(mtk>=tl)
			{
			xmelt=1.0;
			}
		else
			{
			/* 0.15 = 15 wt% of Cpx in the mantle */
			Fcpxout=0.15/(0.5+0.08*mpg);
			Tcpxout=pow(Fcpxout,1.0/1.5)*(tll-ts)+ts;
			if(mtk>Tcpxout)
				{
				/* Opx melting */
				xmelt=Fcpxout+(1.0-Fcpxout)*pow((mtk-Tcpxout)/(tl-Tcpxout),1.5);
				}
			else
				{
				/* Cpx melting */
				xmelt=pow((mtk-ts)/(tll-ts),1.5);
				}
			/* Compute artificial Tliquidus */
			tl=(mtk-ts)/xmelt+ts;
			}
		}
	hlatent=400000.0;
	break;
	/**/
	/* Other rocks - No melting */
	default:
	break;
	}
/**/
/* Melt fraction, latent heat calculation */
eps1[21]=eps1[22]=0;
if(tl)
	{
	/* Melt fraction calc, check */
/*
	if(ts>=tl-90.0) ts=tl-90.0;
*/
	xmelt=(mtk-ts)/(tl-ts);
	if(xmelt<0) xmelt=0;
	if(xmelt>1.0) xmelt=1.0;
	eps1[21]=xmelt;
	/* Latent heat calc */
	hlatent*=xmelt;
	eps1[22]=hlatent;
/*
if(xmelt>0 && (mm2==10 || mm2==30 || mm2==9 || mm2==29)) {printf("Meltpart1: %d %e %e  %e %e %e",mm2,mtk,mpb,ts,tl,xmelt);getchar();}
if(mm2<20 && xmelt) {printf("Meltpart1: %d %e %e  %e %e",mm2,mtk,mpb,xmelt,hlatent);getchar();}
*/
	}
/**/
}
/* Melt fraction, latent heat calculation */


/* Rock to rock+melt transformation */
void melting1(double mtk, double mpb, long int mm1, double *eps1)
/* mtk - T, K */
/* mpb - P, bar */
/* mm1 - mark number */
{
/* Marker type */
int mm2=(int)markt[mm1];
/**/
/**/
/* Melting related cahnge of the marker type */
/* Check marker type */
if (mm2>1)
if (mpb>0)
switch(mm2)
	{
	/* Sediments, upper crust */
	case 15:
	case 35:
	case 3:
	case 4:
	case 5:
	case 23:
	case 24:
	case 25:
	/* Basalt, Gabbro */
	case 6:
	case 7:
	case 8:
	case 16:
	case 17:
	case 18:
	case 36:
	case 37:
	case 38:
	case 27:
	case 28:
	case 26:
	meltpart11(mtk,mpb,mm2,eps1);
	if(eps1[21]>0 && mm2<20) {markt[mm1]+=20; marke[mm1]=0;}
	if(eps1[21]<=0 && mm2>20){markt[mm1]-=20; marke[mm1]=0;}
/*
if (mm2!=markt[mm1]){printf("Granite/Basalt %ld %d %d  %e %e  %e",mm1,mm2,markt[mm1],mtk,mpb,eps1[21]);getchar();}
*/
 	return;
	/**/
	/* Hydrated Peridotite */
	case 11:
	case 13:
	case 14:
	case 34:
	meltpart11(mtk,mpb,mm2,eps1);
	if(eps1[21]>0 && mm2<20) {markt[mm1]=34; marke[mm1]=0;}
	if(eps1[21]<=0 && mm2>20){markt[mm1]=14; marke[mm1]=0;}
/*
if (mm2!=markt[mm1]){printf("Peridotite %ld %d %d  %e %e  %e",mm1,mm2,markt[mm1],mtk,mpb,eps1[21]);getchar();}
*/
 	return;
	/* Dry Peridotite */
	case  9:
	case 10:
	case 12:
	case 29:
	case 30:
	case 32:
	meltpart11(mtk,mpb,mm2,eps1);
	if(eps1[21]>0 && mm2<20) {markt[mm1]+=20; marke[mm1]=0;}
	if(eps1[21]<=0 && mm2>20){markt[mm1]-=20; marke[mm1]=0;}
 	return;
	/* Others */
	default: return;
	}
}
/* Rock to rock+melt transformation */




/* Melt fraction, density, viscosity, heat capacity calculation */
double meltpart0(double mtk, double mpb, long int mm1, int mm2, double *eps1)
/* mtk - T, K */
/* mpb - P, bar */
/* mm1 - mark number */
/* mm2 - mark type */
{
/* Val buffer */
double xmelt=0,ival,p_ga_in,p_pl_out,dmpb,dmtk,epsin,sduct,nueff,smin,smax,nmin,nmax,cpadd=0,xmelt0;
long int m1;
int mm3;
/**/
/* Check marker type */
if (mm2>21)
	{
	/* Calculate melt fraction */
	meltpart11(mtk,mpb,mm2,eps1);
	xmelt0=eps1[21];if(xmelt0<0) xmelt=0;
	xmelt=eps1[21]-markex[mm1];if(xmelt<0) xmelt=0;
	/**/
	/* Solid rock type */
	mm3=mm2-20; if (mm2==34) mm3=11;
	/**/
	/* Standard adiabatic term: al=bro/(1+bro*(Tk-298.15)) */
	eps1[20]=(markbb[mm2]*xmelt+markbb[mm3]*(1.0-xmelt))/(1.0-(markbb[mm2]*xmelt+markbb[mm3]*(1.0-xmelt))*(mtk-298.15));
	/**/
	/* Density */
	/* Ro=ro0 */
	if (densimod==0) 
		{
		eps1[23]=markro[mm2]*xmelt+markro[mm3]*(1.0-xmelt);
		}
	/* Ro=ro0*(1-bro*(TK-298.15))*(1+aro*(Pkbar-0.001)) */
	else
		{
		if(mm2>=29 && mm2<=34 && markex[mm1]>0)
			{
			/* Mantle Depletion */
			eps1[23]=xmelt*markro[mm2]*(1.0-markbb[mm2]*(mtk-298.15))*(1.0+markaa[mm2]*(mpb-1.0)*1e-3)+(1.0-0.04*markex[mm1])*(1.0-xmelt)*markro[mm3]*(1.0-markbb[mm3]*(mtk-298.15))*(1.0+markaa[mm3]*(mpb-1.0)*1e-3);
			}
		else
			{
			eps1[23]=xmelt*markro[mm2]*(1.0-markbb[mm2]*(mtk-298.15))*(1.0+markaa[mm2]*(mpb-1.0)*1e-3)+(1.0-xmelt)*markro[mm3]*(1.0-markbb[mm3]*(mtk-298.15))*(1.0+markaa[mm3]*(mpb-1.0)*1e-3);
			}
		}
	/**/
	/**/
	/**/
	/* Heat capacity */
	eps1[25]=markcp[mm2]*xmelt+markcp[mm3]*(1.0-xmelt);
	/**/
	/* heat conductivity */
	eps1[26]=((markkt[mm2]+markkf[mm2]/(mtk+77.0))*exp(markkp[mm2]*mpb))*xmelt+((markkt[mm3]+markkf[mm3]/(mtk+77.0))*exp(markkp[mm3]*mpb))*(1.0-xmelt);
	/**/
	/* Additional melting adiabatic term, heat capacity */
	if(xmelt>0 && xmelt<1.0)
		{
		/* Melting adiabatic term: alm=-ro*(dHlat/dP)/T */
		/* Numerical differentiation */
		dmpb=mpb*0.001;
		meltpart11(mtk,mpb-dmpb,mm2,eps1);
		ival=eps1[22];
		meltpart11(mtk,mpb+dmpb,mm2,eps1);
		ival-=eps1[22];
		ival*=eps1[23]/(mtk*2.0*dmpb*1e+5);
		eps1[20]+=ival;
		/**/
		/* Melting heat capacity term: cpm=dHlat/dT */
		/* Numerical differentiation */
		dmtk=1.0;
		meltpart11(mtk+dmtk,mpb,mm2,eps1);
		ival=eps1[22];
		meltpart11(mtk-dmtk,mpb,mm2,eps1);
		ival-=eps1[22];
		ival/=2.0*dmtk;
		eps1[25]+=ival;
/*
printf("Meltpart: %ld %d %e %e %e  %e %e %e %e %e %e",mm1,mm2,mtk,mpb,xmelt,eps1[20],eps1[22],eps1[23],eps1[24],eps1[25],eps1[26]);getchar();
*/
		}
	return 1.0;
	}
eps1[20]=eps1[21]=eps1[22]=eps1[23]=eps1[24]=eps1[25]=eps1[26]=0;
return 0;
}
/* Rock to rock+melt transformation */




/* Calculation of Vx,Vy,Vz for current location by Linear Interpolation */
void vxyzcalc1(double x, double y, double z, double *vxyz1, long int *wn1)
/* x,y,z - XYZ location of point for Vx,Vy,Vz calc */
{
/* Vx,Vy-nodes num */
long int m1,m2,m3,m10,m20,m30,m40,m10min,m20min,m30min;
/* Relativ coord */
double dx,dy,dz,ddx,ddy,ddz,swt,ival;
/**/
/**/
/**/
/* Clear Vx,Vy,Vz */
vxyz1[0]=vxyz1[1]=vxyz1[2]=0;
/* Out of grid */
if(x<0) x=0; if(x>xsize) x=xsize;
if(y<0) y=0; if(y>ysize) y=ysize;
if(z<0) z=0; if(z>zsize) z=zsize;
/* Fing Cell indexes */
wn1[0]=m1=m1serch(x);
wn1[1]=m2=m2serch(y);
wn1[2]=m3=m3serch(z);
/**/
/* Vx Cell indexes */
m10min=m1;
if(m10min<0) m10min=0;
if(m10min>xnumx-2) m10min=xnumx-2;
m20min=m2;
if(y<(gy[m20min]+gy[m20min+1])/2.0) m20min-=1;
if(m20min<0) m20min=0;
if(m20min>ynumy-3) m20min=ynumy-3;
m30min=m3;
if(z<(gz[m30min]+gz[m30min+1])/2.0) m30min-=1;
if(m30min<0) m30min=0;
if(m30min>znumz-3) m30min=znumz-3;
/* Vx Cell dimensions */
dx= gx[m10min+1]-gx[m10min];
dy=(gy[m20min+2]-gy[m20min])/2.0;
dz=(gz[m30min+2]-gz[m30min])/2.0;
/* Interpolate from 8 nodes */
ival=0;
for (m30=m30min;m30<=m30min+1;m30++)
	{
	/* Normalized Z-distance calc */
	ddz=((gz[m30]+gz[m30+1])/2.0-z)/dz;
	if (m30==m30min) ddz=-ddz;
	ddz=1.0-ddz;
	for (m10=m10min;m10<=m10min+1;m10++)
		{
		/* Normalized X-distance calc */
		ddx=(gx[m10]-x)/dx;
		if (m10==m10min) ddx=-ddx;
		ddx=1.0-ddx;
		for (m20=m20min;m20<=m20min+1;m20++)
			{
			/* Normalized Y-distance calc */
			ddy=((gy[m20]+gy[m20+1])/2.0-y)/dy;
			if (m20==m20min) ddy=-ddy;
			ddy=1.0-ddy;
			/* Pos in vx[] */
			m40=m30*xynumxy+m10*ynumy+m20;
			/* Calc wt */
			swt=ddx*ddy*ddz;
/*
printf("Vx %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,vx[m40]);getchar();
*/
			/* Add Vx for current node */
			vxyz1[0]+=vx[m40]*swt;
			/* Add Wt for current node */
			ival+=swt;
			}
		}
	}
if(ABSV(ival-1.0)>1e-6){printf("Vx %e %e %e   %ld %ld %ld  %e %e",x,y,z,m10min,m20min,m30min,ival,vxyz[0]);getchar();}
/*
*/
/**/
/* Vy Cell indexes */
m10min=m1;
if(x<(gx[m10min]+gx[m10min+1])/2.0) m10min-=1;
if(m10min<0) m10min=0;
if(m10min>xnumx-3) m10min=xnumx-3;
m20min=m2;
if(m20min<0) m20min=0;
if(m20min>ynumy-2) m20min=ynumy-2;
m30min=m3;
if(z<(gz[m30min]+gz[m30min+1])/2.0) m30min-=1;
if(m30min<0) m30min=0;
if(m30min>znumz-3) m30min=znumz-3;
/* Vy Cell dimensions */
dx=(gx[m10min+2]-gx[m10min])/2.0;
dy= gy[m20min+1]-gy[m20min];
dz=(gz[m30min+2]-gz[m30min])/2.0;
/* Interpolate from 8 nodes */
ival=0;
for (m30=m30min;m30<=m30min+1;m30++)
	{
	/* Normalized Z-distance calc */
	ddz=((gz[m30]+gz[m30+1])/2.0-z)/dz;
	if (m30==m30min) ddz=-ddz;
	ddz=1.0-ddz;
	for (m10=m10min;m10<=m10min+1;m10++)
		{
		/* Normalized X-distance calc */
		ddx=((gx[m10]+gx[m10+1])/2.0-x)/dx;
		if (m10==m10min) ddx=-ddx;
		ddx=1.0-ddx;
		for (m20=m20min;m20<=m20min+1;m20++)
			{
			/* Normalized Y-distance calc */
			ddy=(gy[m20]-y)/dy;
			if (m20==m20min) ddy=-ddy;
			ddy=1.0-ddy;
			/* Pos in vy[] */
			m40=m30*xynumxy+m10*ynumy+m20;
			/* Calc wt */
			swt=ddx*ddy*ddz;
/*
printf("Vx %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,vx[m40]);getchar();
*/
			/* Add Vy for current node */
			vxyz1[1]+=vy[m40]*swt;
			/* Add Wt for current node */
			ival+=swt;
			}
		}
	}
if(ABSV(ival-1.0)>1e-6){printf("Vy %e %e %e   %ld %ld %ld  %e %e",x,y,z,m10min,m20min,m30min,ival,vxyz[1]);getchar();}
/*
*/
/**/
/* Vz Cell indexes */
m10min=m1;
if(x<(gx[m10min]+gx[m10min+1])/2.0) m10min-=1;
if(m10min<0) m10min=0;
if(m10min>xnumx-3) m10min=xnumx-3;
m20min=m2;
if(y<(gy[m20min]+gy[m20min+1])/2.0) m20min-=1;
if(m20min<0) m20min=0;
if(m20min>ynumy-3) m20min=ynumy-3;
m30min=m3;
if(m30min<0) m30min=0;
if(m30min>znumz-2) m30min=znumz-2;
/* Vz Cell dimensions */
dx=(gx[m10min+2]-gx[m10min])/2.0;
dy=(gy[m20min+2]-gy[m20min])/2.0;
dz= gz[m30min+1]-gz[m30min];
/* Interpolate from 8 nodes */
ival=0;
for (m30=m30min;m30<=m30min+1;m30++)
	{
	/* Normalized Z-distance calc */
	ddz=(gz[m30]-z)/dz;
	if (m30==m30min) ddz=-ddz;
	ddz=1.0-ddz;
	for (m10=m10min;m10<=m10min+1;m10++)
		{
		/* Normalized X-distance calc */
		ddx=((gx[m10]+gx[m10+1])/2.0-x)/dx;
		if (m10==m10min) ddx=-ddx;
		ddx=1.0-ddx;
		for (m20=m20min;m20<=m20min+1;m20++)
			{
			/* Normalized Y-distance calc */
			ddy=((gy[m20]+gy[m20+1])/2.0-y)/dy;
			if (m20==m20min) ddy=-ddy;
			ddy=1.0-ddy;
			/* Pos in vz[] */
			m40=m30*xynumxy+m10*ynumy+m20;
			/* Calc wt */
			swt=ddx*ddy*ddz;
/*
printf("Vz %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,vx[m40]);getchar();
*/
			/* Add Vx for current node */
			vxyz1[2]+=vz[m40]*swt;
			/* Add Wt for current node */
			ival+=swt;
			}
		}
	}
if(ABSV(ival-1.0)>1e-6){printf("Vz %e %e %e   %ld %ld %ld  %e %e",x,y,z,m10min,m20min,m30min,ival,vxyz[2]);getchar();}
/*
*/
/**/
}
/* Calculation of Vx,Vy,Vz for current location by Linear Interpolation */






/* Calculation of P, SIG and  EPSxx, EPSyy, EPSzz, EPSxy, EPSxz, EPSyz for current location by Linear Interpolation */
void epscalc1(double x, double y, double z, int yn, double *eps1, long int *wn1)
/* x,y,z - XYZ location of point for Vx,Vy,Vz calc */
/* yn - cell address available in wn[]  Y(1)/N(0) */
{
/* Vx,Vy-nodes num */
long int m1,m2,m3,m10,m20,m30,m40,m10min,m20min,m30min;
/* Relativ coord */
double dx,dy,dz,ddx,ddy,ddz,swt,ival;
/**/
/**/
/**/
/* Clear All EPS,SIG */
eps1[4]=eps1[5]=eps1[6]=eps1[7]=eps1[8]=eps1[9]=eps1[10]=0;
eps1[54]=eps1[55]=eps1[56]=eps1[57]=eps1[58]=eps1[59]=eps1[60]=0;
/* Out of grid */
if(x<0) x=0; if(x>xsize) x=xsize;
if(y<0) y=0; if(y>ysize) y=ysize;
if(z<0) z=0; if(z>zsize) z=zsize;
/* Fing Cell indexes */
if(yn==0)
	{
	wn1[0]=m1=m1serch(x);
	wn1[1]=m2=m2serch(y);
	wn1[2]=m3=m3serch(z);
	}
else
	{
	m1=wn1[0];
	m2=wn1[1];
	m3=wn1[2];
	}
/**/
/**/
/**/
/* EPSxx, EPSyy, EPSzz Cell indexes */
m10min=m1;
if(x>(gx[m10min]+gx[m10min+1])/2.0) m10min+=1;
if(m10min<1) m10min=1;
if(m10min>xnumx-2) m10min=xnumx-2;
m20min=m2;
if(y>(gy[m20min]+gy[m20min+1])/2.0) m20min+=1;
if(m20min<1) m20min=1;
if(m20min>ynumy-2) m20min=ynumy-2;
m30min=m3;
if(z>(gz[m30min]+gz[m30min+1])/2.0) m30min+=1;
if(m30min<1) m30min=1;
if(m30min>znumz-2) m30min=znumz-2;
/* EPS Cell dimensions */
dx=(gx[m10min+1]-gx[m10min-1])/2.0;
dy=(gy[m20min+1]-gy[m20min-1])/2.0;
dz=(gz[m30min+1]-gz[m30min-1])/2.0;
/* Interpolate from 8 nodes */
ival=0;
for (m30=m30min;m30<=m30min+1;m30++)
	{
	/* Normalized Z-distance calc */
	ddz=((gz[m30-1]+gz[m30])/2.0-z)/dz;
	if (m30==m30min) ddz=-ddz;
	ddz=1.0-ddz;
	if(ddz<0) ddz=0;
	if(ddz>1.0) ddz=1.0;
	for (m10=m10min;m10<=m10min+1;m10++)
		{
		/* Normalized X-distance calc */
		ddx=((gx[m10-1]+gx[m10])/2.0-x)/dx;
		if (m10==m10min) ddx=-ddx;
		ddx=1.0-ddx;
		if(ddx<0) ddx=0;
		if(ddx>1.0) ddx=1.0;
		for (m20=m20min;m20<=m20min+1;m20++)
			{
			/* Normalized Y-distance calc */
			ddy=((gy[m20-1]+gy[m20])/2.0-y)/dy;
			if (m20==m20min) ddy=-ddy;
			ddy=1.0-ddy;
			if(ddy<0) ddy=0;
			if(ddy>1.0) ddy=1.0;
			/* Pos in vx[] */
			m40=m30*xynumxy+m10*ynumy+m20;
			/* Calc wt */
			swt=ddx*ddy*ddz;
/*
printf("Vx %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,vx[m40]);getchar();
*/
			/* Add EPS for current node */
			eps1[4]+=exx[m40]*swt;
			eps1[5]+=eyy[m40]*swt;
			eps1[6]+=ezz[m40]*swt;
			/* Add SIG for current node */
			eps1[54]+=sxx[m40]*swt;
			eps1[55]+=syy[m40]*swt;
			eps1[56]+=szz[m40]*swt;
			/* Add DISSIPATION for current node */
			eps1[60]+=(sxx[m40]*exx[m40]+syy[m40]*eyy[m40]+szz[m40]*ezz[m40])*swt;
/*
*/
			/* Add P   for current node */
			eps1[10]+=pr[m40]*swt;
			/* Add Wt for current node */
			ival+=swt;
			}
		}
	}
/*
if(ABSV(ival-1.0)>1e-6){printf("EPS  %e %e %e   %ld %ld %ld  %e %e",x,y,z,m10min,m20min,m30min,ival,eps1[4]);getchar();}
*/
/**/
/* EPSxy Cell indexes */
m10min=m1;
if(m10min<1) m10min=1;
if(m10min>xnumx-3) m10min=xnumx-3;
m20min=m2;
if(m20min<1) m20min=1;
if(m20min>ynumy-3) m20min=ynumy-3;
m30min=m3;
if(z<(gz[m30min]+gz[m30min+1])/2.0) m30min-=1;
if(m30min<0) m30min=0;
if(m30min>znumz-3) m30min=znumz-3;
/* EPSxy Cell dimensions */
dx= gx[m10min+1]-gx[m10min];
dy= gy[m20min+1]-gy[m20min];
dz=(gz[m30min+2]-gz[m30min])/2.0;
/* Interpolate from 8 nodes */
ival=0;
for (m30=m30min;m30<=m30min+1;m30++)
	{
	/* Normalized Z-distance calc */
	ddz=((gz[m30]+gz[m30+1])/2.0-z)/dz;
	if (m30==m30min) ddz=-ddz;
	ddz=1.0-ddz;
	if(ddz<0) ddz=0;
	if(ddz>1.0) ddz=1.0;
	for (m10=m10min;m10<=m10min+1;m10++)
		{
		/* Normalized X-distance calc */
		ddx=(gx[m10]-x)/dx;
		if (m10==m10min) ddx=-ddx;
		ddx=1.0-ddx;
		if(ddx<0) ddx=0;
		if(ddx>1.0) ddx=1.0;
		for (m20=m20min;m20<=m20min+1;m20++)
			{
			/* Normalized Y-distance calc */
			ddy=(gy[m20]-y)/dy;
			if (m20==m20min) ddy=-ddy;
			ddy=1.0-ddy;
			if(ddy<0) ddy=0;
			if(ddy>1.0) ddy=1.0;
			/* Pos in vx[] */
			m40=m30*xynumxy+m10*ynumy+m20;
			/* Calc wt */
			swt=ddx*ddy*ddz;
/*
printf("Vx %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,vx[m40]);getchar();
*/
			/* Add EPSxy for current node */
			eps1[7]+=exy[m40]*swt;
			/* Add SIGxy for current node */
			eps1[57]+=sxy[m40]*swt;
			/* Add DISSIPATION for current node */
			eps1[60]+=2.0*sxy[m40]*exy[m40]*swt;
/*
*/
			/* Add Wt for current node */
			ival+=swt;
			}
		}
	}
/*
if(ABSV(ival-1.0)>1e-6){printf("EPSxy %e %e %e   %ld %ld %ld  %e %e",x,y,z,m10min,m20min,m30min,ival,eps1[7]);getchar();}
*/
/**/
/* EPSxz Cell indexes */
m10min=m1;
if(m10min<1) m10min=1;
if(m10min>xnumx-3) m10min=xnumx-3;
m20min=m2;
if(y<(gy[m20min]+gy[m20min+1])/2.0) m20min-=1;
if(m20min<0) m20min=0;
if(m20min>ynumy-3) m20min=ynumy-3;
m30min=m3;
if(m30min<1) m30min=1;
if(m30min>znumz-3) m30min=znumz-3;
/* EPSxz Cell dimensions */
dx= gx[m10min+1]-gx[m10min];
dy=(gy[m20min+2]-gy[m20min])/2.0;
dz= gz[m30min+1]-gz[m30min];
/* Interpolate from 8 nodes */
ival=0;
for (m30=m30min;m30<=m30min+1;m30++)
	{
	/* Normalized Z-distance calc */
	ddz=(gz[m30]-z)/dz;
	if (m30==m30min) ddz=-ddz;
	ddz=1.0-ddz;
	if(ddz<0) ddz=0;
	if(ddz>1.0) ddz=1.0;
	for (m10=m10min;m10<=m10min+1;m10++)
		{
		/* Normalized X-distance calc */
		ddx=(gx[m10]-x)/dx;
		if (m10==m10min) ddx=-ddx;
		ddx=1.0-ddx;
		if(ddx<0) ddx=0;
		if(ddx>1.0) ddx=1.0;
		for (m20=m20min;m20<=m20min+1;m20++)
			{
			/* Normalized Y-distance calc */
			ddy=((gy[m20]+gy[m20+1])/2.0-y)/dy;
			if (m20==m20min) ddy=-ddy;
			ddy=1.0-ddy;
			if(ddy<0) ddy=0;
			if(ddy>1.0) ddy=1.0;
			/* Pos in vz[] */
			m40=m30*xynumxy+m10*ynumy+m20;
			/* Calc wt */
			swt=ddx*ddy*ddz;
/*
printf("Vz %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,vx[m40]);getchar();
*/
			/* Add EPSxz for current node */
			eps1[8]+=exz[m40]*swt;
			/* Add SIGxz for current node */
			eps1[58]+=sxz[m40]*swt;
			/* Add DISSIPATION for current node */
			eps1[60]+=2.0*sxz[m40]*exz[m40]*swt;
/*
*/
			/* Add Wt for current node */
			ival+=swt;
			}
		}
	}
/*
if(ABSV(ival-1.0)>1e-6){printf("EPSxz %e %e %e   %ld %ld %ld  %e %e",x,y,z,m10min,m20min,m30min,ival,eps1[8]);getchar();}
*/
/**/
/* EPSyz Cell indexes */
m10min=m1;
if(x<(gx[m10min]+gx[m10min+1])/2.0) m10min-=1;
if(m10min<0) m10min=0;
if(m10min>xnumx-3) m10min=xnumx-3;
m20min=m2;
if(m20min<1) m20min=1;
if(m20min>ynumy-3) m20min=ynumy-3;
m30min=m3;
if(m30min<1) m30min=1;
if(m30min>znumz-3) m30min=znumz-3;
/* Vy Cell dimensions */
dx=(gx[m10min+2]-gx[m10min])/2.0;
dy= gy[m20min+1]-gy[m20min];
dz= gz[m30min+1]-gz[m30min];
/* Interpolate from 8 nodes */
ival=0;
for (m30=m30min;m30<=m30min+1;m30++)
	{
	/* Normalized Z-distance calc */
	ddz=(gz[m30]-z)/dz;
	if (m30==m30min) ddz=-ddz;
	ddz=1.0-ddz;
	if(ddz<0) ddz=0;
	if(ddz>1.0) ddz=1.0;
	for (m10=m10min;m10<=m10min+1;m10++)
		{
		/* Normalized X-distance calc */
		ddx=((gx[m10]+gx[m10+1])/2.0-x)/dx;
		if (m10==m10min) ddx=-ddx;
		ddx=1.0-ddx;
		if(ddx<0) ddx=0;
		if(ddx>1.0) ddx=1.0;
		for (m20=m20min;m20<=m20min+1;m20++)
			{
			/* Normalized Y-distance calc */
			ddy=(gy[m20]-y)/dy;
			if (m20==m20min) ddy=-ddy;
			ddy=1.0-ddy;
			if(ddy<0) ddy=0;
			if(ddy>1.0) ddy=1.0;
			/* Pos in vy[] */
			m40=m30*xynumxy+m10*ynumy+m20;
			/* Calc wt */
			swt=ddx*ddy*ddz;
/*
printf("Vx %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,vx[m40]);getchar();
*/
			/* Add EPSyz for current node */
			eps1[9]+=eyz[m40]*swt;
			/* Add SIGyz for current node */
			eps1[59]+=syz[m40]*swt;
			/* Add DISSIPATION for current node */
			eps1[60]+=2.0*syz[m40]*eyz[m40]*swt;
/*
*/
			/* Add Wt for current node */
			ival+=swt;
			}
		}
	}
if(ABSV(ival-1.0)>1e-6){printf("EPSyz %e %e %e   %ld %ld %ld  %e %e",x,y,z,m10min,m20min,m30min,ival,eps1[9]);getchar();}
if(eps1[60]<0){printf("DISSIP %e %e %e   %ld %ld %ld  %e %e",x,y,z,m10min,m20min,m30min,ival,eps1[60]);getchar();}
/*
*/
/**/
}
/* Calculation of EPSxx, EPSyy, EPSzz, EPSxy, EPSxz, EPSyz for current location by Linear Interpolation */






/* Calculation of TK current location by Linear Interpolation */
void tkcalc1(double x, double y, double z, double *eps1, long int *wn1, double *wi1)
/* x,y,z - XYZ location of point for Vx,Vy,Vz calc */
{
/* Vx,Vy-nodes num */
long int m1,m2,m3,m4,m10,m20,m30,m40;
int n1;
/* Relativ coord */
double dx,dy,dz,ddx,ddy,ddz,swt,ival;
/**/
/**/
/**/
/* Cell num calc for A-D Diagonal Corners of cur marker */
wn1[0]=m10=m1serch(x);
wn1[1]=m20=m2serch(y);
wn1[2]=m30=m3serch(z);
wi1[0]=dx=gx[m10+1]-gx[m10];
wi1[1]=dy=gy[m20+1]-gy[m20];
wi1[2]=dz=gz[m30+1]-gz[m30];
/**/
/* Clear TK */
eps1[0]=eps1[1]=eps1[2]=eps1[3]=0;
/* Add nodes around current cell */
n1=0;
ival=0;
for (m3=m30;m3<=m30+1;m3++)
	{
	/* Normalized Z-distance calc */
	ddz=(gz[m3]-z)/dz;
	ddz=1.0-ABSV(ddz);
	for (m1=m10;m1<=m10+1;m1++)
		{
		/* Normalized X-distance calc */
		ddx=(gx[m1]-x)/dx;
		ddx=1.0-ABSV(ddx);
		for (m2=m20;m2<=m20+1;m2++)
			{
			/* Normalized Y-distance calc */
			ddy=(gy[m2]-y)/dy;
			ddy=1.0-ABSV(ddy);
			/* Wt calc */
			wi1[3+n1]=swt=ddx*ddy*ddz;
			/**/
			/* Node num */
			wn1[3+n1]=m4=m3*xynumxy+m1*ynumy+m2;
/*
printf("%ld  %ld %ld %ld  %e %e %e   %e %e %e   %e %e %e   %e %d %e",m4,m1,m2,m3,x,y,z,dx,dy,dz,ddx,ddy,ddz,swt,n1,tk[m4]);getchar();
*/
			/**/
			/* Add TK for current node */
			eps1[0]+=swt*tk0[m4];
			eps1[1]+=swt*tk1[m4];
			eps1[2]+=swt*tk[m4];
			eps1[3]+=swt*tk2[m4];
			/* Add Wt */
			ival+=swt;
			/* Add counter */
			n1++;
			}
		}
	}
/*
if(eps1[2]<900.0){printf("T  %e %e %e   %ld %ld %ld  %e %e",x,y,z,m10,m20,m30,ival,eps1[2]);getchar();}
{printf("T  %e %e %e   %ld %ld %ld  %e %e",x,y,z,m10,m20,m30,ival,eps1[2]);getchar();}
if(ABSV(ival-1.0)>1e-6){printf("T  %e %e %e   %ld %ld %ld  %e %e",x,y,z,m10,m20,m30,ival,eps1[2]);getchar();}
*/
/**/
}
/* Calculation of TK for current location by Linear Interpolation */




/* Calculation of P for current location by Linear Interpolation */
void  prcalc1(double x, double y, double z, double *eps1, long int *wn1)
/* x,y,z - XYZ location of point for Vx,Vy,Vz calc */
{
/* Vx,Vy-nodes num */
long int m1,m2,m3,m10,m20,m30,m40,m10min,m20min,m30min;
/* Relativ coord */
double dx,dy,dz,ddx,ddy,ddz,swt,ival,mpb1;
/**/
/**/
/**/
/* Clear P */
eps1[10]=0;
/* Out of grid */
if(x<0) x=0; if(x>xsize) x=xsize;
if(y<0) y=0; if(y>ysize) y=ysize;
if(z<0) z=0; if(z>zsize) z=zsize;
/* Fing Cell indexes */
wn1[0]=m1=m1serch(x);
wn1[1]=m2=m2serch(y);
wn1[2]=m3=m3serch(z);
/**/
/* P Cell indexes */
m10min=m1;
if(x>(gx[m10min]+gx[m10min+1])/2.0) m10min+=1;
if(m10min<1) m10min=1;
if(m10min>xnumx-2) m10min=xnumx-2;
m20min=m2;
if(y>(gy[m20min]+gy[m20min+1])/2.0) m20min+=1;
if(m20min<1) m20min=1;
if(m20min>ynumy-2) m20min=ynumy-2;
m30min=m3;
if(z>(gz[m30min]+gz[m30min+1])/2.0) m30min+=1;
if(m30min<1) m30min=1;
if(m30min>znumz-2) m30min=znumz-2;
/* P Cell dimensions */
dx=(gx[m10min+1]-gx[m10min-1])/2.0;
dy=(gy[m20min+1]-gy[m20min-1])/2.0;
dz=(gz[m30min+1]-gz[m30min-1])/2.0;
/* Interpolate from 8 nodes */
ival=0;
for (m30=m30min;m30<=m30min+1;m30++)
	{
	/* Normalized Z-distance calc */
	ddz=((gz[m30-1]+gz[m30])/2.0-z)/dz;
	if (m30==m30min) ddz=-ddz;
	ddz=1.0-ddz;
	for (m10=m10min;m10<=m10min+1;m10++)
		{
		/* Normalized X-distance calc */
		ddx=((gx[m10-1]+gx[m10])/2.0-x)/dx;
		if (m10==m10min) ddx=-ddx;
		ddx=1.0-ddx;
		for (m20=m20min;m20<=m20min+1;m20++)
			{
			/* Normalized Y-distance calc */
			ddy=((gy[m20-1]+gy[m20])/2.0-y)/dy;
			if (m20==m20min) ddy=-ddy;
			ddy=1.0-ddy;
			/* Pos in vx[] */
			m40=m30*xynumxy+m10*ynumy+m20;
			/* Calc wt */
			swt=ddx*ddy*ddz;
/*
printf("Vx %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,vx[m40]);getchar();
*/
			/* Add P for current node */
			eps1[10]+=pr[m40]*swt;
			/* Add Wt for current node */
			ival+=swt;
			}
		}
	}
/* Lithostatic pressure */
/*
mpb1=(y-1e+4)*1e+5/3.0;
if(mpb1<0) mpb1=0;
if(eps1[10]<0) eps1[10]=0;
if(eps1[10]>2.0*mpb1) eps1[10]=2.0*mpb1;
*/
/**/
/*
if(ABSV(ival-1.0)>1e-6){printf("P  %e %e %e   %ld %ld %ld  %e %e",x,y,z,m10min,m20min,m30min,ival,vxyz[0]);getchar();}
*/
/**/
}
/* Calculation of P for current location by Linear Interpolation */


/* Calc ro for given P,T */
double dencalc1(double mtk, double mpb, int mm2, double *eps1)
/* mtk - T, K */
/* mpb - P, bar */
/* x,y - XY location of point for Vx,Vy calc */
/* mm2 - Rock number */
{
/* Val buffer */
double ival;
eps1[20]=0;
/**/
/* Constant density */
if (densimod==0) return markro[mm2];
/**/
/* Ro=ro0*(1-bro*(TK-298.15))*(1+aro*(Pkbar-0.001)) */
/* Adiabatic term: al=bro/(1-bro*(Tk-298.15)) */
ival=markro[mm2]*(1.0-markbb[mm2]*(mtk-298.15))*(1.0+markaa[mm2]*(mpb-1.0)*1e-3);
eps1[20]=markbb[mm2]/(1.0-markbb[mm2]*(mtk-298.15));
return  ival;
}
/* End Calc ro for given P,T */




/*1[21]-markex[mm1]-meltmin) Melt extraction, crust deposition */
void meltextract()
{
/* Variables */
long int m1,m2,m3,m4=xnumx*znumz,m5,m6,m7,m5min,m6min,m7min,m1cur,m2cur,m3cur;
long int mm1,wn1[500];
double mtk,mpb,meltsum=0,meltsum0,meltsum1[500],meltnum[10001],meltlim=2.0,crustsum=0,crustsum0,crustsum1[500],extractsum=0,extractsum0,extractsum1[500];
int yn,yn1,n1,n2,nn0,nn1;
double eps1[100],mycur,mymin,mymax,mycel,myavr;
/**/
printf("Crust growth BEG \n");
/**/
/* Clear melt extraction regions */
for (m1=0;m1<xnumx*znumz;m1++)
	{
	bufv[m1]=0;
	bufv[m1+m4]=ydike;
	bufv[m1+m4*2]=0;
	bufv[m1+m4*3]=ysize-stp100;
	}
/**/
/* Clear molten mantle counter */
for (m1=0;m1<=10000;m1++)
	{
	meltnum[m1]=0;
	}
/**/
/**/
/* Check molten markers */
#pragma omp parallel shared(nn0,bufv,meltsum1,markx,marky,markz,markk,markt,markex) private(nn1,meltsum0,mm1,m1,m2,m3,mpb,mtk,wn1,eps1) firstprivate(stp100,xsize,ysize,zsize,m4,meltmin,timesum,xnumx,marknum)
{
/* Obtain cur thread number */
nn1=omp_get_thread_num();
/* Obtain total number of threads */
if (nn1==0) nn0=omp_get_num_threads();
/* Reset melt counter */
meltsum0=0;
/*
*/
#pragma omp for schedule(static)
/*
*/
for (mm1=0;mm1<marknum;mm1++)
if(markx[mm1]>=0 && marky[mm1]>=0 && markz[mm1]>=0 && (double)(markx[mm1])<=xsize && (double)(marky[mm1])<=ysize && (double)(markz[mm1])<=zsize)
if(marky[mm1]<ysize-stp100)
{
/* Find magma chamber depth */
if(markt[mm1]==27 || markt[mm1]==28)
	{
	m1=m1serch(markx[mm1]);
	m3=m3serch(markz[mm1]);
	m2=m3*xnumx+m1+m4;
	if(bufv[m2]>marky[mm1]) 
		{
		bufv[m2]=marky[mm1];
		}
	}
/* Compute melt extraction */
if(markt[mm1]==9 || markt[mm1]==10 || markt[mm1]==29 || markt[mm1]==30) 
	{
	/* P, T parameters calc */
	prcalc1((double)(markx[mm1]),(double)(marky[mm1]),(double)(markz[mm1]),eps1,wn1);
	mpb=eps1[10]*1e-5;
	mtk=(double)(markk[mm1]);
	/* Melting degree calc */
	meltpart11(mtk,mpb,markt[mm1],eps1);
	/* Extract melts before crustal growth starts */
	if(eps1[21]-markex[mm1]>meltmin) 
		{
		/* Add global melt counter */
		meltsum0+=(eps1[21]-markex[mm1]-meltmin);
		/* Add local melt counter */
/*
		#pragma omp atomic
*/
		bufv[m4*2+wn1[2]*xnumx+wn1[0]]+=(eps1[21]-markex[mm1]-meltmin);
		/* Extract melt */
		if(timesum<=1e+11) 
			{
			markex[mm1]+=(eps1[21]-markex[mm1]-meltmin); 
			}
		}
	/* Find minimal molten mantle depth */
	if(eps1[21]>0 && marky[mm1]<bufv[m4*3+wn1[2]*xnumx+wn1[0]])
		{
		bufv[m4*3+wn1[2]*xnumx+wn1[0]]=marky[mm1];
		}
	}
}
/* Save melt counter */
meltsum1[nn1]=meltsum0;
}
/* Add melt counter */
for (n1=0;n1<nn0;n1++)
	{
	meltsum+=meltsum1[n1];
	}
/**/
/* Mark melt extraction regions */
for (m3=0;m3<znumz;m3++)
for (m1=1;m1<xnumx-1;m1++)
	{
	m2=m3*xnumx+m1+m4;
	if(bufv[m2]<ydike && bufv[m2]<bufv[m2-1] && bufv[m2]<bufv[m2+1]) bufv[m2-m4]=bufv[m2];
	}
/**/
/**/
/* Crustal growth */
if(timesum>1e+11 && meltsum>0)
{
/* Melt movement cycle counter */
n1=0;
/* Melt exchange cycle */
do
	{
	/* Melt exchange mark */
	yn=0;
	/* Move melt in the mantle */
	for (m3=0;m3<znumz;m3++)
	for (m1=1;m1<xnumx-1;m1++)
	if(bufv[m3*xnumx+m1+m4*2]>0)
		{
		/* Melt level in the current cell */
		m2=m3*xnumx+m1+m4*2;
		mycel=bufv[m2+m4]+bufv[m2]/((double)(mnumx*mnumy*mnumz))*ysize/((double)(ynumy-1));
		mymin=mycel;
		m7min=m2;
		m3cur=m3;
		m1cur=m1;
		m2cur=m2;
		/* Melt level exchange with the surrounding cells */
		do
			{
			yn1=0;
			for (m5=m3cur-1;m5<=m3cur+1;m5++)
			for (m6=m1cur-1;m6<=m1cur+1;m6++)
			if(m5>=0 && m5<znumz && m6>=0 && m6<xnumx)
				{
				/* Melt level in the current cell */
				m7=m5*xnumx+m6+m4*2;
				mycur=bufv[m7+m4]+bufv[m7]/((double)(mnumx*mnumy*mnumz))*ysize/((double)(ynumy-1));
				/* Find shallowest melt level */
				if(mycur<mymin)
					{
					mymin=mycur;
					m5min=m5;
					m6min=m6;
					m7min=m7;
					}
				}
			/* Equilibrate melt levels */
			if(m7min!=m2cur && mycel-mymin>1.0)
				{
/*
printf("A %ld %ld %e %e %e %e ",m2cur-m4*2,m7min-m4*2,bufv[m2cur],bufv[m7min],bufv[m2cur+m4],bufv[m7min+m4]);getchar();
*/
				/* Average melt level */
				myavr=0.5*(mycel+mymin);
				if(myavr<bufv[m2cur+m4])
					{
					/* All melt goes to one cell */
					bufv[m7min]+=bufv[m2cur];
					bufv[m2cur]=0;
					}
				else
					{
					/* Melt partitionning */
					bufv[m7min]+=0.5*(mycel-mymin)/(ysize/((double)(ynumy-1)))*((double)(mnumx*mnumy*mnumz));
					bufv[m2cur]-=0.5*(mycel-mymin)/(ysize/((double)(ynumy-1)))*((double)(mnumx*mnumy*mnumz));
					}
/*
printf("B %ld %ld %e %e %e %e ",m2cur-m4*2,m7min-m4*2,bufv[m2cur],bufv[m7min],bufv[m2cur+m4],bufv[m7min+m4]);getchar();
*/
				/* Melt exchange mark */
				yn=1;
				yn1=1;
				m1cur=m6min;
				m3cur=m5min;
				m2cur=m7min;
				mycel=bufv[m2cur+m4]+bufv[m2cur]/((double)(mnumx*mnumy*mnumz))*ysize/((double)(ynumy-1));
				mymin=mycel;
				}
			}
		while(yn1==1);
		}
	n1++;
	}
while(yn==1);
printf("MELT MOVEMENT CYCLES=%d \n",n1);
/**/
/**/
/* Define Mantle->Crust transmutation treshold */
#pragma omp parallel for shared(bufv,meltnum,markx,marky,markz,markk,markt) private(mm1,m2,m5,mpb,mtk,wn1,eps1,mycel) firstprivate(mnumx,mnumy,mnumz,stp100,xsize,ysize,zsize,m4,xnumx,marknum) schedule(static)
/*
*/
for (mm1=0;mm1<marknum;mm1++)
if(markt[mm1]==9 || markt[mm1]==10 || markt[mm1]==29 || markt[mm1]==30) 
if(markx[mm1]>=0 && marky[mm1]>=0 && markz[mm1]>=0 && (double)(markx[mm1])<=xsize && (double)(marky[mm1])<=ysize && (double)(markz[mm1])<=zsize)
if(marky[mm1]<ysize-stp100)
	{
	/* P, T parameters calc */
	prcalc1((double)(markx[mm1]),(double)(marky[mm1]),(double)(markz[mm1]),eps1,wn1);
	mpb=eps1[10]*1e-5;
	mtk=(double)(markk[mm1]);
	/* Melting degree calc */
	meltpart11(mtk,mpb,markt[mm1],eps1);
	/* Melt level in the current cell */
	m2=wn1[2]*xnumx+wn1[0]+m4*2;
	mycel=bufv[m2+m4]+bufv[m2]/((double)(mnumx*mnumy*mnumz))*ysize/((double)(ynumy-1));
	/* Mantle->Crust transmutation */
	if(eps1[21]>0 && bufv[m2]>0 && marky[mm1]>=bufv[m2+m4] && marky[mm1]<=mycel)
		{
		m5=(int)(bufv[m2]*10000.0); if(m5>10000) m5=10000;
/*
		#pragma omp atomic
*/
		meltnum[m5]+=1.0;
		}
	}
/* Recalc molten mantle counter */
mycur=0.375;
for (m1=9999;m1>=0;m1--)
	{
	meltnum[m1]+=meltnum[m1+1];
	if(meltnum[m1+1]<meltsum && meltnum[m1]>=meltsum)
		{
		/* Compute mantle transmutation melt amount */
		mycur=((double)(m1)+(meltnum[m1]-meltsum)/(meltnum[m1]-meltnum[m1+1]))/10000.0;	
		break;
		}
	}
printf("AVAILABLE MELT=%e  MELT TRESHOLD=%e \n",meltsum,mycur);
/**/
/**/
/* Mantle->crust transmutation */
#pragma omp parallel shared(bufv,crustsum1,markx,marky,markz,markk,markt) private(crustsum0,nn1,mm1,m2,mpb,mtk,wn1,eps1,mycel) firstprivate(mnumx,mnumy,mnumz,stp100,xsize,ysize,zsize,m4,xnumx,marknum)
{
/*
*/
/* Obtain cur thread number */
nn1=omp_get_thread_num();
/* Reset crust counter */
crustsum0=0;
/*
*/
if(meltsum>0)
{
#pragma omp for schedule(static)
/*
*/
for (mm1=0;mm1<marknum;mm1++)
if(markt[mm1]==9 || markt[mm1]==10 || markt[mm1]==29 || markt[mm1]==30) 
if(markx[mm1]>=0 && marky[mm1]>=0 && markz[mm1]>=0 && (double)(markx[mm1])<=xsize && (double)(marky[mm1])<=ysize && (double)(markz[mm1])<=zsize)
if(marky[mm1]<ysize-stp100)
	{
	/* P, T parameters calc */
	prcalc1((double)(markx[mm1]),(double)(marky[mm1]),(double)(markz[mm1]),eps1,wn1);
	mpb=eps1[10]*1e-5;
	mtk=(double)(markk[mm1]);
	/* Melting degree calc */
	meltpart11(mtk,mpb,markt[mm1],eps1);
	/* Melt level in the current cell */
	m2=wn1[2]*xnumx+wn1[0]+m4*2;
	mycel=bufv[m2+m4]+bufv[m2]/((double)(mnumx*mnumy*mnumz))*ysize/((double)(ynumy-1));
	/* Mantle->Crust transmutation */
	if(eps1[21]>0 && bufv[m2]>mycur && marky[mm1]>=bufv[m2+m4] && marky[mm1]<=mycel)
		{
		markt[mm1]=27;
		markex[mm1]=0;
		markrr[mm1]=0;
		markrn[mm1]=0;
		marke[mm1]=0;
		crustsum0+=1.0;
		}
	}
}
/* Save crust counter */
crustsum1[nn1]=crustsum0;
}
/* Add crust counter */
for (n1=0;n1<nn0;n1++)
	{
	crustsum+=crustsum1[n1];
	}
/**/
/* Extract melt from molten markers */
/* Mantle->crust transmutation */
#pragma omp parallel shared(bufv,extractsum1,markx,marky,markz,markk,markex,markt) private(extractsum0,nn1,mm1,mpb,mtk,wn1,eps1,mycel) firstprivate(mnumx,mnumy,mnumz,stp100,xsize,ysize,zsize,xnumx,meltmin,crustsum,meltsum,marknum)
{
/*
*/
/* Obtain cur thread number */
nn1=omp_get_thread_num();
/* Reset extract counter */
extractsum0=0;
if(meltsum>0 && crustsum>0)
{
#pragma omp for schedule(static)
/*
*/
for (mm1=0;mm1<marknum;mm1++)
if(markt[mm1]==9 || markt[mm1]==10 || markt[mm1]==29 || markt[mm1]==30) 
if(markx[mm1]>=0 && marky[mm1]>=0 && markz[mm1]>=0 && (double)(markx[mm1])<=xsize && (double)(marky[mm1])<=ysize && (double)(markz[mm1])<=zsize)
if(marky[mm1]<ysize-stp100)
	{
	/* P, T parameters calc */
	prcalc1((double)(markx[mm1]),(double)(marky[mm1]),(double)(markz[mm1]),eps1,wn1);
	mpb=eps1[10]*1e-5;
	mtk=(double)(markk[mm1]);
	/* Melting degree calc */
	meltpart11(mtk,mpb,markt[mm1],eps1);
	/* Extract melts */
	if(eps1[21]-markex[mm1]>meltmin) 
		{
		/* Add extracted melt counters */
		extractsum0+=(eps1[21]-markex[mm1]-meltmin)*crustsum/meltsum;
		/* Extract melt */
		markex[mm1]+=(eps1[21]-markex[mm1]-meltmin)*crustsum/meltsum; 
		}
	}
}
/* Save extract counte */
extractsum1[nn1]=extractsum0;
}
/* Add extract counters */
for (n1=0;n1<nn0;n1++)
	{
	extractsum+=extractsum1[n1];
	}
}
/**/
printf("AVAILABLE MELT=%e  EXTRACTED MELT=%e  NEW CRUST =%e \n",meltsum,extractsum,crustsum);
}
/* Melt extraction, crust deposition */



