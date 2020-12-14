#include<stdio.h>
#include<stdlib.h>
#define NINT(a) ((a) >= 0.0 ? (int)((a)+0.5) : (int)((a)-0.5))	//Nearest integer

int main()
{
int i,j,k;
int n;
n=500;
int s[n];
int t,nn,clust_count[n];
int count;
double x[n],y[n],z[n];
double xr,yr,zr,rr;
double sigma,rc,rc2;
double box;
FILE *ip,*op;
sigma=1.0;
rc=1.5*sigma;
rc2=rc*rc;
box=21.73074;
ip=fopen("snapshot_to_be_analyzed.xyz","r");
  fscanf(ip,"%d\n%d\n",&nn,&t);
    for(i=0;i<n;i++)
    {
    fscanf(ip,"%d %lf %lf %lf\n",&s[i],&x[i],&y[i],&z[i]);
    clust_count[i]=i;
    }
fclose(ip);

    for(i=0;i<n;i++)
    {
    	for(j=0;j<n;j++)
	{
	xr=x[i]-x[j];
	yr=y[i]-y[j];
	zr=z[i]-z[j];
	xr=xr-box*NINT(xr/box);
	yr=yr-box*NINT(yr/box);
	zr=zr-box*NINT(zr/box);
	rr=xr*xr+yr*yr+zr*zr;
	    if(rr<rc2)
	    {
	        if(clust_count[i]<clust_count[j])
		clust_count[j]=clust_count[i];
		else
		clust_count[i]=clust_count[j];
	    }
	}
    	for(j=0;j<n;j++)
	{
	xr=x[i]-x[j];
	yr=y[i]-y[j];
	zr=z[i]-z[j];
	xr=xr-box*NINT(xr/box);
	yr=yr-box*NINT(yr/box);
	zr=zr-box*NINT(zr/box);
	rr=xr*xr+yr*yr+zr*zr;
	    if(rr<rc2)
	    {
	        if(clust_count[i]<clust_count[j])
		clust_count[j]=clust_count[i];
		else
		clust_count[i]=clust_count[j];
	    }
	}

    }
count=0;
    for(i=0;i<n;i++)
    {
        if(clust_count[i]==0)count++;
        printf("%d %d\n",i,clust_count[i]);
    }
    op=fopen("cluster.xyz","w");
    fprintf(op,"%d\n%d\n",count,1);
    for(i=0;i<n;i++)
        if(clust_count[i]==0)
	{
	//if(s[i]==0)
	fprintf(op,"%d %lf %lf %lf\n",s[i],x[i],y[i],z[i]);
	//else
	//fprintf(op,"Np %lf %lf %lf\n",x[i],y[i],z[i]);
	}
    fclose(op);


return 0;
}
