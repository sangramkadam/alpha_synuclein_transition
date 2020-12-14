#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
//#include<mpi.h>


#define np 100000        //Number of particles
#define np_c 500          //Number of colloidal particles
#define Temp 1.0        //Temperature of the system
#define dt_mpc 0.1      //MPCD Time step
#define dt_md 0.001     //MD Time step
#define stepLimit 20000000//Maximum number of steps
#define rho 10.0        //Density of liquid
#define epsA 3.0        //Parameter of LJ potential
#define epsB 5.0        //Parameter of LJ potential
#define epsC 15.0        //Parameter of LJ potential
#define sigma_sc 0.5    //Parameter of LJ potential
#define sigma_cc 1.0    //Parameter of LJ potential
#define alpha 2.27      //Rotation angle
#define NINT(a) ((a) >= 0.0 ? (int)((a)+0.5) : (int)((a)-0.5))	//Nearest integer




#define Sqr(x)     ((x) * (x))
#define Cube(x)    ((x) * (x) * (x))
#define AllocMem(a, n, t)  a = (t *) malloc ((n) * sizeof (t))


//Random number generator
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
        int j;
        long k;
        static long idum2=123456789;
        static long iy=0;
        static long iv[NTAB];
        float temp;

        if (*idum <= 0) {
                if (-(*idum) < 1) *idum=1;
                else *idum = -(*idum);
                idum2=(*idum);
                for (j=NTAB+7;j>=0;j--) {
                        k=(*idum)/IQ1;
                        *idum=IA1*(*idum-k*IQ1)-k*IR1;
                        if (*idum < 0) *idum += IM1;
                        if (j < NTAB) iv[j] = *idum;
                }
                iy=iv[0];
        }
        k=(*idum)/IQ1;
        *idum=IA1*(*idum-k*IQ1)-k*IR1;
        if (*idum < 0) *idum += IM1;
        k=idum2/IQ2;
        idum2=IA2*(idum2-k*IQ2)-k*IR2;
        if (idum2 < 0) idum2 += IM2;
        j=iy/NDIV;
        iy=iv[j]-idum2;
        iv[j] = *idum;
        if (iy < 1) iy += IMM1;
        if ((temp=AM*iy) > RNMX) return RNMX;
        else return temp;
}



void Initialize(void);
void Position (void);
void Velocity (void);
void Force(void);
double Rangauss (void);
void Sphere(double *r1,double *r2,double *r3);
void Grid(double *x1,double *x2,double *x3,double a);
void PeriodicBoundary(double* x,double*y,double*z);
void ShiftGrid(double x1,double x2,double x3);
void RestoreGrid(double x1,double x2,double x3);
void CreateNighborsList(int nc,int nc_dim,int neigh[][27]);
void CreateLinkedList(int hoc[],int ll[],double a,int n,int nc,int nc_dim);
void SingleMPCDstep(void);
void SingleMDstep(int nc,int nc_dim,int neigh[][27],double rcell);
void VerletUpdate1(void);
void VerletUpdate2(void);
void Sample(double *ke,double* sumv);
void RadialDistribution(void);
void Diffusion(void);

//Global Variables

double x[np],y[np],z[np];                   //current positions of solvent particles
double xc[np_c],yc[np_c],zc[np_c];          //current positions of colloidal particles
double gxc[np_c],gyc[np_c],gzc[np_c];       //global positions of colloidal particles

double vx[np],vy[np],vz[np];                //velocities of solvent particles
double vxc[np_c],vyc[np_c],vzc[np_c];       //velocities of colloidal particles

double fx[np],fy[np],fz[np];                //force on solvent particles
double fxc[np_c],fyc[np_c],fzc[np_c];       //force on colloidal paticles

double box;                                 //box size
double mc;                                  //mass of colloidal particle
double rc_sc,rc_cc;                         //cutoff radius
double rc2_sc,rc2_cc;                       //cutoff radius squared
double rc6_sc,rc6_cc,sigma6_sc,sigma6_cc;
double ecut_sc,ecut_cc;                     //cutoff potential
double r_eq,ks,kappa;
double msd[10000];                          //mean square displacement
double RN[10000];                           //end to end polymer distance
double xo,yo,zo;                            //position at start of observation
double rate;			    //Rate of conversion per bead
int step=0,nstep=0;
int solventFlag[np],monomerFlag[np_c],na,nb,nb_total;
int stepCount=0,collision_frq;
long idum;
double t,tt;
double T;

FILE *fp;

////////////////////////////////////////////////////////////////

int main()
{
//fp=fopen("time","w");
//fclose(fp);
double cx,cy,cz;
//MPI_Init(NULL,NULL);
double start_time,end_time,current_time;    //Start and end time for program
//start_time=MPI_Wtime();	                    //Record start time
int world_rank;								//rank of the process in MPI_COMM_WORLD
int world_size;								//size of MPI_COMM_WORLD
char file[64];
double rnum;
//MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
//MPI_Comm_size(MPI_COMM_WORLD,&world_size);

idum=-98397324;
rate=0.0000001;

//srand(time(NULL)+world_rank);
int nc,nc_dim,i;
double rcell;

Initialize();
nc_dim=(int)(box/rc_sc);
nc=Cube(nc_dim);
rcell=box/(double)nc_dim;
int hoc[nc],neigh[nc][27];                 
CreateNighborsList(nc,nc_dim,neigh);

collision_frq=(int)(dt_mpc/dt_md);
t=0.0;
    while(stepCount<stepLimit)
    {
        if((stepCount)%5000==0)
        {
        sprintf(file,"nb%d_rate_e-7_eps3_5_15.xyz",nb_total);
        fp=fopen(file,"a");
        fprintf(fp,"%d\n%d\n",np_c,stepCount-1);
        for(i=0;i<np_c;i++)
        {
            fprintf(fp,"%d %lf %lf %lf\n",monomerFlag[i],xc[i],yc[i],zc[i]);
        }
        fclose(fp);
        }
    SingleMDstep(nc,nc_dim,neigh,rcell);
    t = stepCount * dt_md;
    stepCount++;
        if(stepCount%collision_frq==0)
        {
        SingleMPCDstep();
        }

        if(stepCount%100000==0)
        {    
        fp=fopen("time","a");
        fprintf(fp,"%d %d\t",world_rank,stepCount);
        fclose(fp);
        }

    }
    
//end_time=MPI_Wtime();
//fp=fopen("time","a");
//fprintf(fp,"run time for process %d is %lf\n",world_rank,end_time-start_time);
//fclose(fp);
//MPI_Finalize();
return 0;
}

///////////////////////////////////////////////////////////////




void Initialize(void)
{
int i;
double eps;
nb_total = 10;
na = np-nb;
eps = epsA;
/*
    for(i=0;i<np_c;i++)
    {                       
	if(i<nb)
	monomerFlag[i]=1;
	else		    //flag 0 A
    	monomerFlag[i]=0;   //flag 1 B
    }                       //flag 2 A which is nighbor of B
    */
    for(i=0;i<10000;i++)
    {
    msd[i]=0.0;
    }
box=pow((np/rho+((double)(np_c))*4.188790205*Cube(sigma_sc)),1.0/3.0);
mc=rho*4.188790205*Cube(sigma_sc);
rc_sc=pow(2.0,1.0/6.0)*sigma_sc;
rc_cc=2.5*sigma_cc;
rc2_sc=Sqr(rc_sc);
rc2_cc=Sqr(rc_cc);
sigma6_sc=Sqr(Cube(sigma_sc));
sigma6_cc=Sqr(Cube(sigma_cc));
rc6_sc=(sigma6_sc/(Cube(rc2_sc)));
rc6_cc=(sigma6_cc/(Cube(rc2_cc)));
ecut_sc=4.0*eps*rc6_sc*(rc6_sc-1.0);
ecut_cc=4.0*eps*rc6_cc*(rc6_cc-1.0);
r_eq = rc_cc+0.1;
ks = 20.0;


Position();
Velocity ();
Force();
}




void Position(void)
{
int i,j;
int xi,yi,zi;
double d,dd,xr,yr,zr;


//Positions of colloidal particles


//Cluster positions

j=8;
nb=0;
    for(i=0;i<np_c;i++)
    {
    zi=(int)(i/(j*j));
    yi=(int)((i-j*j*zi)/j);
    xi=(int)(i-j*j*zi-j*yi);
    xc[i]=(double)xi;
    yc[i]=(double)yi;
    zc[i]=(double)zi;
        if(xi<2 && yi<2 && zi<3 && nb<nb_total)
        {
        nb++;
        monomerFlag[i]=1;
        }
        else monomerFlag[i]=0;
    gxc[i]=xc[i];
    gyc[i]=yc[i];
    gzc[i]=zc[i];
    }

//Positions of solvent particles

    for(i=0;i<np;i++)
    {
    d=0.0;
    dd=box;
    while(d<Sqr(sigma_sc))
    {
    x[i]=rand()*box/RAND_MAX;
    y[i]=rand()*box/RAND_MAX;
    z[i]=rand()*box/RAND_MAX;
        for(j=0;j<np_c;j++)
        {
        if(j==0)
        {
        xr=x[i]-xc[j];
        yr=y[i]-yc[j];
        zr=z[i]-zc[j];
        xr=xr-box*NINT(xr/box);
        yr=yr-box*NINT(yr/box);
        zr=zr-box*NINT(zr/box);	
        d=xr*xr+yr*yr+zr*zr;
        }
        else
        {
        xr=x[i]-xc[j];
        yr=y[i]-yc[j];
        zr=z[i]-zc[j];
        xr=xr-box*NINT(xr/box);
        yr=yr-box*NINT(yr/box);
        zr=zr-box*NINT(zr/box);	
        dd=xr*xr+yr*yr+zr*zr;
        }
        if(dd<d)d=dd;
        }
    }
    
    }

}

//Function for Assigning random velocities with guassian distribution


void Velocity (void)
{
double sumvx,sumvy,sumvz,sumv2;     //sum of velocities
double fs;                          //scaling factor
double T;                           //Temperature
int i,j,k;

sumvx=0.0;
sumvy=0.0;
sumvz=0.0;
sumv2=0.0;
    for (i=0;i<np_c;i++)
    {
    vxc[i]=Rangauss();
    vyc[i]=Rangauss();
    vzc[i]=Rangauss();
    sumvx=sumvx+vxc[i];
    sumvy=sumvy+vyc[i];
    sumvz=sumvz+vzc[i];
    sumv2=sumv2 + (vxc[i]*vxc[i]+vyc[i]*vyc[i]+vzc[i]*vzc[i]);
    }
    for (i=0;i<np;i++)
    {
    vx[i]=Rangauss();
    vy[i]=Rangauss();
    vz[i]=Rangauss();
    sumvx=sumvx+vx[i];
    sumvy=sumvy+vy[i];
    sumvz=sumvz+vz[i];
    sumv2=sumv2 + (vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
    }
sumvx=sumvx/(np+np_c);                     //Average velocity
sumvy=sumvy/(np+np_c);
sumvz=sumvz/(np+np_c);
sumv2=sumv2/(np+np_c);
fs=sqrt(3.0*Temp/sumv2);            //scaling factor
sumv2=0.0;
  for (i=0;i<np_c;i++)                  //scaling of velocities
    {
    vxc[i]=(vxc[i]-sumvx)*fs;
    vyc[i]=(vyc[i]-sumvy)*fs;
    vzc[i]=(vzc[i]-sumvz)*fs;
    sumv2=sumv2+(vxc[i]*vxc[i]+vyc[i]*vyc[i]+vzc[i]*vzc[i]);
    }
    for (i=0;i<np;i++)                  //scaling of velocities
    {
    vx[i]=(vx[i]-sumvx)*fs;
    vy[i]=(vy[i]-sumvy)*fs;
    vz[i]=(vz[i]-sumvz)*fs;
    sumv2=sumv2+(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
    }	
T=sumv2/(3.0*(np));
}


//Function for Guassian distribution

double Rangauss (void)
{
double u,v,R2,r;
R2=0.0;
u=0.0;
v=0.0;
    while(R2==0||R2>=1.0)
    {
    u=rand();
    v=rand();
    u=u/RAND_MAX;
    v=v/RAND_MAX;
    u=2.0*u-1.0;
    v=2.0*v-1.0;
    R2=u*u+v*v;
    }
r=u*sqrt(-2.0*log(R2)/R2);
return r;
}

void Force(void)
{
int i;

    for(i=0;i<np;i++)                 //Initialize forces to zero
    {
    fx[i]=0.0;
    fy[i]=0.0;
    fz[i]=0.0;
    }
    for(i=0;i<np_c;i++)
    {
    fxc[i]=0.0;
    fyc[i]=0.0;
    fzc[i]=0.0;
    }
}

void Sphere(double *r1,double *r2,double *r3)
{
double x1,x2,z,phi;
x1=rand();
x2=rand();
x1=x1/RAND_MAX;
x2=x2/RAND_MAX;
z=2.0*x1-1.0;
phi=2.0*3.14159*x2;
*r1=sqrt(1-z*z)*cos(phi);
*r2=sqrt(1-z*z)*sin(phi);
*r3=z; 
}

void Grid(double *x1,double *x2,double *x3,double a)
{
*x1=rand();
*x2=rand();
*x3=rand();
*x1=*x1/RAND_MAX;
*x2=*x2/RAND_MAX;
*x3=*x3/RAND_MAX;
*x1=(*x1-0.5)*a;
*x2=(*x2-0.5)*a;
*x3=(*x3-0.5)*a;
}



void PeriodicBoundary(double* xx,double* yy,double* zz)
{
if(*xx<0.0)*xx=*xx+box;
if(*xx>box)*xx=*xx-box;
if(*yy<0.0)*yy=*yy+box;
if(*yy>box)*yy=*yy-box;
if(*zz<0.0)*zz=*zz+box;
if(*zz>box)*zz=*zz-box; 
}
void ShiftGrid(double x1,double x2,double x3)
{
int i;
    for(i=0;i<np;i++)
    {
    x[i]=x[i]+x1;                   
    y[i]=y[i]+x2;
    z[i]=z[i]+x3;
    PeriodicBoundary(&x[i],&y[i],&z[i]);
    }
}

void RestoreGrid(double x1,double x2,double x3)
{
int i;
    for(i=0;i<np;i++)
    {
    x[i]=x[i]-x1;                   
    y[i]=y[i]-x2;
    z[i]=z[i]-x3;
    PeriodicBoundary(&x[i],&y[i],&z[i]);
    }
}

void CreateNighborsList(int nc,int nc_dim,int neigh[][27])
{
int icell,jcell,ix,iy,iz,h,i,j,k,ii,jj,kk;
for(icell=0;icell<nc;icell++)
{
iz=(int) (icell/(nc_dim*nc_dim));
iy=(int)((icell-nc_dim*nc_dim*iz)/nc_dim);
ix=icell-nc_dim*nc_dim*iz-nc_dim*iy;
h=0;
for(k=iz-1;k<=iz+1;k++)
	{
	for(j=iy-1;j<=iy+1;j++)
		{
		for(i=ix-1;i<=ix+1;i++)
			{
			ii=i;
			jj=j;
			kk=k;
	
			if(ii<0)
			ii=ii+nc_dim;
			if(ii>nc_dim-1)
			ii=ii-nc_dim;
			if(jj<0)
			jj=jj+nc_dim;
			if(jj>nc_dim-1)
			jj=jj-nc_dim;
			if(kk<0)
			kk=kk+nc_dim;
			if(kk>nc_dim-1)
			kk=kk-nc_dim;
			jcell=ii+nc_dim*jj+nc_dim*nc_dim*kk;
			neigh[icell][h]=jcell;
			h++;
			}
		}
	}
}



}


void CreateLinkedList(int hoc[],int ll[],double a,int n,int nc,int nc_dim)
{
int cellx,celly,cellz,cell,i;
    for(i=0;i<nc;i++)               
    {	
    hoc[i]=-1;
    }

    for(i=0;i<n;i++)               
    {                               
    cellx=(int)(x[i]/a);
    celly=(int)(y[i]/a);
    cellz=(int)(z[i]/a);
    cell=cellx+nc_dim*celly+nc_dim*nc_dim*cellz;
    ll[i]=hoc[cell];
    hoc[cell]=i;	
    }

}






void SingleMPCDstep(void)
{
double vcmx,vcmy,vcmz;
double dvx,dvy,dvz;
double dv_dot_r;
double dvx_para,dvx_ortho;          //parallel and perpendicular- 
double dvy_para,dvy_ortho;          //-components of dv wrt r
double dvz_para,dvz_ortho;          // 
double x1,x2,x3;
double r1,r2,r3;
double a;
int nc,nc_dim;
int id,jcell,m;
a=1.0;
nc_dim=(int)(box/a);
nc=nc_dim*nc_dim*nc_dim;
a=box/(double)nc_dim;
int hoc[nc],ll[np];

Grid(&x1,&x2,&x3,a);
ShiftGrid(x1,x2,x3);
CreateLinkedList(hoc,ll,a,np,nc,nc_dim);

    for(jcell=0;jcell<nc;jcell++)
    {
    vcmx=0.0;
    vcmy=0.0;
    vcmz=0.0;
    m=0;
    id=hoc[jcell];
        while(id!=-1)
        {                           //calculate centre of mass velocity
        vcmx=vcmx+vx[id];
        vcmy=vcmy+vy[id];
        vcmz=vcmz+vz[id];
        m++;
        id=ll[id];
        }
    vcmx=vcmx/(double)m;
    vcmy=vcmy/(double)m;
    vcmz=vcmz/(double)m;
    Sphere(&r1,&r2,&r3);            //choose random vector on a sphere
    
    id=hoc[jcell];
        while(id!=-1)
        {
        dvx=vx[id]-vcmx;
        dvy=vy[id]-vcmy;
        dvz=vz[id]-vcmz;
        dv_dot_r=dvx*r1+dvy*r2+dvz*r3;
        dvx_para=dv_dot_r*r1;
        dvy_para=dv_dot_r*r2;
        dvz_para=dv_dot_r*r3;
        dvx_ortho=dvx-dvx_para;
        dvy_ortho=dvy-dvy_para;
        dvz_ortho=dvz-dvz_para;
        vx[id]=vcmx+dvx_ortho*cos(alpha)+dvx_para+sin(alpha)*(dvy_ortho*r3-dvz_ortho*r2);
        vy[id]=vcmy+dvy_ortho*cos(alpha)+dvy_para+sin(alpha)*(dvz_ortho*r1-dvx_ortho*r3);
        vz[id]=vcmz+dvz_ortho*cos(alpha)+dvz_para+sin(alpha)*(dvx_ortho*r2-dvy_ortho*r1);
        id=ll[id];
        }
    }

RestoreGrid(x1,x2,x3);

}








void SingleMDstep(int nc,int nc_dim,int neigh[][27],double rcell)
{
double xi,yi,zi;                            //position of i'th particle
double xr,yr,zr,r,r2,r2i,r6i;               //distance between particles
double en,ke,etot,sumv;
double eps;
int cellx,celly,cellz;                      //cell index of particles
int cell,jcell;
int ll[np];
int i,j,n,m,id,ff;
int hoc[nc];
CreateLinkedList(hoc,ll,rcell,np,nc,nc_dim);
en=0.0;
etot=0.0;
ke=0.0;
VerletUpdate1();

//Force calculation
    for(i=0;i<np;i++)
    {
    fx[i]=0.0;
    fy[i]=0.0;
    fz[i]=0.0;
    }
    for(i=0;i<np_c;i++)
    {
    fxc[i]=0.0;
    fyc[i]=0.0;
    fzc[i]=0.0;
    }

    for(n=0;n<np_c;n++)
    {
    cellx=(int)(xc[n]/rcell);
    celly=(int)(yc[n]/rcell);
    cellz=(int)(zc[n]/rcell);
    cell=cellx+nc_dim*celly+nc_dim*nc_dim*cellz;
    
        for(m=0;m<27;m++)
    	{
    	jcell=neigh[cell][m];
    	id=hoc[jcell];
    		while(id!=-1)
    		{
    		xr=xc[n]-x[id];
    		yr=yc[n]-y[id];
           	zr=zc[n]-z[id];
    		xr=xr-box*NINT(xr/box);
    		yr=yr-box*NINT(yr/box);
    		zr=zr-box*NINT(zr/box);	
    		r2=xr*xr+yr*yr+zr*zr;
    	   		if(r2<rc2_sc)
    			{

                    eps=epsA;
    	    	r2i=(Sqr(sigma_sc))/r2;
                r6i=r2i*r2i*r2i;
                ff=eps*48.0*r2i*r6i*(r6i-0.5)/Sqr(sigma_sc);
                fxc[n]=fxc[n]+ff*xr;
                fyc[n]=fyc[n]+ff*yr;
                fzc[n]=fzc[n]+ff*zr;
                fx[id]=fx[id]-ff*xr;
                fy[id]=fy[id]-ff*yr;
                fz[id]=fz[id]-ff*zr;
                en=en+4.0*eps*r6i*(r6i-1.0)-ecut_sc;	
    			}
    		id=ll[id];
    		}
    	}
    }

    for(i=0;i<np_c-1;i++)
    {    
    xi=xc[i];
    yi=yc[i];
    zi=zc[i];
        for(j=i+1;j<np_c;j++)
        {
            if(monomerFlag[i]%2==0 && monomerFlag[j]%2==0)
            eps = epsA;
	    else if(monomerFlag[i]==1 && monomerFlag[j]==1)
	    eps = epsC;
	    else
	    eps = epsB;
        xr=xi-xc[j];
        xr=xr-box*NINT(xr/box);
        yr=yi-yc[j];
        yr=yr-box*NINT(yr/box);
        zr=zi-zc[j];
        zr=zr-box*NINT(zr/box);
        r2=xr*xr+yr*yr+zr*zr;
	    if(r2<1.5*1.5)
	    {
                if(monomerFlag[i]==1 && monomerFlag[j]==0)
		{
		    if(ran2(&idum)<rate)
			    monomerFlag[j]=1;
		}
                if(monomerFlag[i]==0 && monomerFlag[j]==1)
		{
		    if(ran2(&idum)<rate)
			    monomerFlag[i]=1;
		}
	    
	    }
            if(r2<rc2_cc)          //???????? rc2 or rc2_cc
            {
            r2i=(sigma_cc*sigma_cc)/r2;
            r6i=r2i*r2i*r2i;
            ff=48.0*eps*r2i*r6i*(r6i-0.5)/Sqr(sigma_cc);
            fxc[i]=fxc[i]+ff*xr;
            fxc[j]=fxc[j]-ff*xr;
            fyc[i]=fyc[i]+ff*yr;
            fyc[j]=fyc[j]-ff*yr;
            fzc[i]=fzc[i]+ff*zr;
            fzc[j]=fzc[j]-ff*zr;
            en=en+4.0*eps*r6i*(r6i-1.0)-ecut_cc;
            }

        }
    }




VerletUpdate2();
Sample(&ke,&sumv);
en=en/((double)(np+np_c));
etot=en+ke;
/*
    if(stepCount%100==0)
    {
    fp=fopen("energy","a");
    fprintf(fp," %d %lf %lf %lf %lf %lf %d %d\n",stepCount,ke,en,etot,sumv,T,na,nb);
    fclose(fp);
    }
*/
}



void VerletUpdate1(void)
{
int i;
double xold,yold,zold;
    for(i=0;i<np;i++)
    {   
    x[i]=x[i]+dt_md*vx[i]+0.5*dt_md*dt_md*fx[i];
    y[i]=y[i]+dt_md*vy[i]+0.5*dt_md*dt_md*fy[i];
    z[i]=z[i]+dt_md*vz[i]+0.5*dt_md*dt_md*fz[i];
    vx[i]=vx[i]+0.5*dt_md*fx[i];
    vy[i]=vy[i]+0.5*dt_md*fy[i];
    vz[i]=vz[i]+0.5*dt_md*fz[i];
    PeriodicBoundary(&x[i],&y[i],&z[i]);
    }

    for(i=0;i<np_c;i++)
    {
    xold=xc[i];
    yold=yc[i];
    zold=zc[i];
    xc[i]=xc[i]+dt_md*vxc[i]+0.5*dt_md*dt_md*fxc[i]/mc;
    yc[i]=yc[i]+dt_md*vyc[i]+0.5*dt_md*dt_md*fyc[i]/mc;
    zc[i]=zc[i]+dt_md*vzc[i]+0.5*dt_md*dt_md*fzc[i]/mc;
    gxc[i]=gxc[i]+(xc[i]-xold);
    gyc[i]=gyc[i]+(yc[i]-yold);
    gzc[i]=gzc[i]+(zc[i]-zold);
    vxc[i]=vxc[i]+0.5*dt_md*fxc[i]/mc;
    vyc[i]=vyc[i]+0.5*dt_md*fyc[i]/mc;
    vzc[i]=vzc[i]+0.5*dt_md*fzc[i]/mc;
    PeriodicBoundary(&xc[i],&yc[i],&zc[i]);
    }

}

void VerletUpdate2(void)
{
int i;
    for(i=0;i<np;i++)
    {
    vx[i]=vx[i]+0.5*dt_md*fx[i];
    vy[i]=vy[i]+0.5*dt_md*fy[i];
    vz[i]=vz[i]+0.5*dt_md*fz[i];
    }
    for(i=0;i<np_c;i++)
    {
    vxc[i]=vxc[i]+0.5*dt_md*fxc[i]/mc;
    vyc[i]=vyc[i]+0.5*dt_md*fyc[i]/mc;
    vzc[i]=vzc[i]+0.5*dt_md*fzc[i]/mc;
    }
}



void Sample(double *ke,double* sumv)
{
double sumvx,sumvy,sumvz,sumv2;
int i;
sumvx=0.0;
sumvy=0.0;
sumvz=0.0;
sumv2=0.0;
    for(i=0;i<np;i++)
    {
    sumvx=sumvx+vx[i];
    sumvy=sumvy+vy[i];
    sumvz=sumvz+vz[i];
    sumv2=sumv2+vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
    }
    for(i=0;i<np_c;i++)
    {
    sumvx=sumvx+vxc[i]*mc;
    sumvy=sumvy+vyc[i]*mc;
    sumvz=sumvz+vzc[i]*mc;
    sumv2=sumv2+(vxc[i]*vxc[i]+vyc[i]*vyc[i]+vzc[i]*vzc[i]);
    }
T=sumv2/(3.0*((double)(np+np_c)));
*ke=0.5*sumv2/((double)(np+np_c));
	//printf("%lf %lf %lf %lf\n",sumvx,sumvy,sumvz,*ke);
*sumv=sqrt(sumvx*sumvx+sumvy*sumvy+sumvz*sumvz);

}


void RadialDistribution(void)
{
int i,j,k,ig;
int nbin = 100;         //total number of bins
double delg;            //bin size
double g[nbin];         //radial distribution function
double r;               //distance between two particles
double vb;              //volume between bin i+1 and i
double nid;             //number of ideal gas part. in vb
double xr,yr,zr,r2;

delg=box/(2*nbin);
    for(i=0;i<=nbin;i++)
    {
    g[i]=0.0;
    }
    
    for (i=0; i<np_c;i++)
    {
        for (j=0;j<np; j++)
        {
            if(solventFlag[j]==1)
            {    
            xr=xc[i]-x[j];
            xr=xr-box*NINT(xr/box);
            yr=yc[i]-y[j];
            yr=yr-box*NINT(yr/box);
            zr=zc[i]-z[j];
            zr=zr-box*NINT(zr/box);
            r=sqrt(xr*xr+yr*yr+zr*zr);
                if(r<(box/2.0))
                {
                ig=(int)(r/delg);
                g[ig]=g[ig]+2;
                }
            }
        }
    }
    
/*fp=fopen("num","r");
fscanf(fp,"%d",&k);
fclose(fp);

char file[64];
sprintf(file,"radial");
fp=fopen(file,"w");
    for(i=0;i<nbin;i++)
    {
    r=delg*(i+0.5);
    vb=((i+1)*(i+1)*(i+1)-i*i*i)*delg*delg*delg;
    nid=(4/3)*3.14159*vb;
    g[i]=g[i]/(np*nid);
    fprintf(fp,"%lf %lf  \n",r,g[i]);
    }
fclose(fp);
*/
}



void Diffusion(void)
{
int i;
double sum;
double xcm,ycm,zcm;
xcm = 0.0;
ycm = 0.0;
zcm = 0.0;
    for(i=0;i<np_c;i++)
    {   
    xcm = xcm+gxc[i];
    ycm = ycm+gyc[i];
    zcm = zcm+gzc[i];
    }
xcm = xcm/((double)np_c);
ycm = ycm/((double)np_c);
zcm = zcm/((double)np_c);
    if(stepCount%100000==0)
    {
    xo = xcm;
    yo = ycm;
    zo = zcm;    
    step=0;
    nstep++;
    }
    if(stepCount>=100000&&stepCount%10==0)
    {
    sum = (xcm-xo)*(xcm-xo) + (ycm-yo)*(ycm-yo) + (zcm-zo)*(zcm-zo);
    msd[step] = msd[step]+sum;
    step++;
    }
}



