#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define NINT(a) ((a) >= 0.0 ? (int)((a)+0.5) : (int)((a)-0.5))	//Nearest integer

int main()
{
int i,j,k;
int n;
n=500;
int s[n];
int t,nn,clust_count[n];
int count[100];
int count1[100];
int dd;
int na,nb;
double x[n],y[n],z[n];
double xr,yr,zr,rr,dr[n];
double xcm,ycm,zcm;
double rcm;
double rmax,vmax;
double rho,rho1;
double sigma,rc,rc2;
double avg;
double width;
double box;
FILE *ip,*op;
sigma=1.0;
rc=1.2*sigma;
rc2=rc*rc;
box=21.73074;
xcm=0.0;
ycm=0.0;
zcm=0.0;
na=0;
nb=0;
rmax=0.0;
ip=fopen("cluster.xyz","r");
  fscanf(ip,"%d\n%d\n",&n,&t);
    for(i=0;i<n;i++)
    {
    fscanf(ip,"%d %lf %lf %lf\n",&s[i],&x[i],&y[i],&z[i]);
    clust_count[i]=i;
    xcm=xcm+x[i];
    ycm=ycm+y[i];
    zcm=zcm+z[i];
        if(s[i]==0)na++;
        else nb++;
    }
xcm=xcm/(double)n;
ycm=ycm/(double)n;
zcm=zcm/(double)n;
    for(i=0;i<n;i++)
    {
    xr=xcm-x[i];
    yr=ycm-y[i];
    zr=zcm-z[i];
    rr=xr*xr+yr*yr+zr*zr;
    dr[i]=sqrt(rr);
        if(rmax<dr[i])
            rmax=dr[i];
    }
fclose(ip);
vmax=4.0*3.14*rmax*rmax*rmax/3.0;
width=0.25;
    for(i=0;i<100;i++)count[i]=0;
    for(i=0;i<100;i++)count1[i]=0;
    for(i=0;i<n;i++)
    {
        if(s[i]==0)
        {
        dd=(int)(dr[i]/width);
        count[dd]++;
        }
        else
        {
        dd=(int)(dr[i]/width);
        count1[dd]++;
        }
    }
    if(na+nb!=n)printf("ERROR\n");
    for(i=0;i<40;i++)
    {
    rr=width*(double)i+0.5*width;
    rho=((double)count[i])/(4.0*3.14*rr*rr*width);
    rho1=((double)count1[i])/(4.0*3.14*rr*rr*width);
    printf("%lf %lf %lf\n",(double)i*width+0.5*width,rho,rho1);
    }

return 0;
}
