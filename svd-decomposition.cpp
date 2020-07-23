#include<iostream>
#include<cmath>
#include <fstream>
#include <stdio.h>
using namespace std;


#define NR_END 1
#define FREE_ARG char*
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))

double **dmatrix(int nrl, int nrh, int ncl, int nch);
double *dvector(int nl, int nh);
void free_dvector(double *v, int nl, int nh);
double pythag(double a, double b);
void svdcmp(double **a, int m, int n, double w[], double **v);
double myround(double a);
double **dmatrix(int, int, int, int);
double *dvector(int, int);
void svdcmp(double **, int, int, double *, double **);
int main() {
    
    int i,j,k;
    double t,temp;
    
    int m, n,time,M,N,MN=1;
    fstream fin("input.txt");
    fstream fout("output.txt", ios::app);
    if ( ! fin) cout << "failed" <<endl;
    
    fin>>time;
    fout<<time<<endl;
    for(int p=0; p<time; p++)
    {
        double **a;
        double *w;
        double **u,**v;
        fin>>M>>N;
        if(N>M)
        {
            MN = N;
        }
        else MN = M;
        /* 矩阵均以M,N中的最大值来申请空间，避免越界 */
        double *t1 = new double[MN];
        double *t2 = new double[MN];
        a = dmatrix(1, MN, 1, MN);
        u = dmatrix(1, MN, 1, MN);
        w = dvector(1, MN);
        v = dmatrix(1, MN, 1, MN);
        
        for (int i = 1; i < M+1; i++)
            for (int j = 1; j < N+1; j++)
            fin >> a[i][j];
        
      //  printf("=== a ===\n");
      //  print_r(a, M, N);
    

        for (i = 1; i <= M; i++) {
            for (j = 1; j <= N; j++)
                u[i][j] = a[i][j];
        }
      
        svdcmp(u, MN, MN, w, v);
   
        for (i = 1; i <= N; i++) {
            for (j = i+1; j <= N; j++) {
                if (w[i] < w[j]) {
                    t = w[i];
                    w[i] = w[j];
                    w[j] = t;
                    for (k = 1; k <= M; k++) {
                        t1[k] = u[k][i];
                    }
                    for (k = 1; k <= M; k++) {
                        u[k][i] = u[k][j];
                    }
                    for (k = 1; k <= M; k++) {
                        u[k][j] = t1[k];
                    }
     
                    
                    for (k = 1; k <= N; k++) {
                        t2[k] = v[k][i];
                    }
                    for (k = 1; k <= N; k++) {
                        v[k][i] = v[k][j];
                    }
                    for (k = 1; k <= N; k++) {
                        v[k][j] = t2[k];
                    }
                }
            }
        }
//        printf("%d %d\n",M,M);
//        print_r(u, M, M);
        fout<<M<<' '<<M<<endl;
        for (int f = 1; f <= M; f++)
        {
            for (int g = 1; g <= M; g++)
            {
                if(u[f][g]==0)  fout  <<fixed<<setprecision(2)<<"0.00"<<" ";
                else
                fout  <<fixed<<setprecision(2)<<myround(u[f][g])<<" ";
            }
            fout<<endl;
        }
     
        double **W;
        W = dmatrix(1, MN, 1, MN);
        for (i = 1; i<= M; i++) {
            for (j =1; j <= N; j++) {
                if (i==j) {
                    W[i][j] = w[i];
                } else {
                    W[i][j] = 0.0;
                }
            }
        }
        fout<<M<<' '<<N<<endl;
//        printf("%d %d\n",M,N);
//        print_r(W, M, N);
        for (int f = 1; f <= M; f++) {
            for (int g = 1; g <= N; g++) {
                if(W[f][g]==0)  fout  <<fixed<<setprecision(2)<<"0.00"<<" ";
                else
                fout  <<fixed<<setprecision(2)<<myround(W[f][g])<<" ";
            }
            fout<<endl;
        }
        
        fout<<N<<' '<<N<<endl;
//        printf("%d %d\n",N,N);
//        print_r(v, N, N);
        for (int f = 1; f <= N; f++) {
            for (int g = 1; g <= N; g++) {
                if(v[f][g]==0)  fout  <<fixed<<setprecision(2)<<"0.00"<<" ";
                else
                fout  <<fixed<<setprecision(2)<<myround(v[f][g])<<" ";
            }
            fout<<endl;
        }
   
        delete  [] t1;
        delete  [] t2;
    }
    fin.close();
    fout.close();
    return 0;
}

double myround(double a)
{
    double diff,ans;
    if(a>0)
    {
        diff = a*100-round(a*100);
        if(diff > 0.5 || diff == 0.5 )
        {
            
            return(a+0.01);
        }
        //ans += (double)diff%100.0
        //ans += diff/100.0;
        else
            return(a);
    }
    else
    {
        diff = -a*100-round((-a)*100);
        if(diff > 0.5 || diff == 0.5 )
        {
            return(a-0.01);
        }
        //ans += (double)diff%100.0
        //ans += diff/100.0;
        else
            return(a);
    }
}


void svdcmp(double **a, int m, int n, double w[], double **v)
{
    int flag,i,its,j,jj,k,l,nm;
    double anorm,c,f,g,h,s,scale,x,y,z,*rv1;
 
    rv1=dvector(1,n);
    g=scale=anorm=0.0;
    for (i=1;i<=n;i++) {
        l=i+1;
        rv1[i]=scale*g;
        g=s=scale=0.0;
        if (i <= m) {
            for (k=i;k<=m;k++) scale += fabs(a[k][i]);
            if (scale) {
                for (k=i;k<=m;k++) {
                    a[k][i] /= scale;
                    s += a[k][i]*a[k][i];
                }
                f=a[i][i];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                a[i][i]=f-g;
                for (j=l;j<=n;j++) {
                    for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
                    f=s/h;
                    for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
                }
                for (k=i;k<=m;k++) a[k][i] *= scale;
            }
        }
        w[i]=scale *g;
        g=s=scale=0.0;
        if (i <= m && i != n) {
            for (k=l;k<=n;k++) scale += fabs(a[i][k]);
            if (scale) {
                for (k=l;k<=n;k++) {
                    a[i][k] /= scale;
                    s += a[i][k]*a[i][k];
                }
                f=a[i][l];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                a[i][l]=f-g;
                for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
                for (j=l;j<=m;j++) {
                    for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
                    for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
                }
                for (k=l;k<=n;k++) a[i][k] *= scale;
            }
        }
        anorm = DMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
    }
    for (i=n;i>=1;i--) {
        if (i < n) {
            if (g) {
                for (j=l;j<=n;j++)
                    v[j][i]=(a[i][j]/a[i][l])/g;
                for (j=l;j<=n;j++) {
                    for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
                    for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
                }
            }
            for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
        }
        v[i][i]=1.0;
        g=rv1[i];
        l=i;
    }
    for (i=IMIN(m,n);i>=1;i--) {
        l=i+1;
        g=w[i];
        for (j=l;j<=n;j++) a[i][j]=0.0;
        if (g) {
            g=1.0/g;
            for (j=l;j<=n;j++) {
                for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
                f=(s/a[i][i])*g;
                for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
            }
            for (j=i;j<=m;j++) a[j][i] *= g;
        } else for (j=i;j<=m;j++) a[j][i]=0.0;
        ++a[i][i];
    }
 
    for (k=n;k>=1;k--) {
        for (its=1;its<=30;its++) {
            flag=1;
            for (l=k;l>=1;l--) {
                nm=l-1;
                if ((double)(fabs(rv1[l])+anorm) == anorm) {
                    flag=0;
                    break;
                }
                if ((double)(fabs(w[nm])+anorm) == anorm) break;
            }
            if (flag) {
                c=0.0;
                s=1.0;
                for (i=l;i<=k;i++) {
                    f=s*rv1[i];
                    rv1[i]=c*rv1[i];
                    if ((double)(fabs(f)+anorm) == anorm) break;
                    g=w[i];
                    h=pythag(f,g);
                    w[i]=h;
                    h=1.0/h;
                    c=g*h;
                    s = -f*h;
                    for (j=1;j<=m;j++) {
                        y=a[j][nm];
                        z=a[j][i];
                        a[j][nm]=y*c+z*s;
                        a[j][i]=z*c-y*s;
                    }
                }
            }
            z=w[k];
            if (l == k) {
                if (z < 0.0) {
                    w[k] = -z;
                    for (j=1;j<=n;j++) v[j][k] = -v[j][k];
                }
                break;
            }
            if (its == 30) printf("no convergence in 30 svdcmp iterations\n");
            x=w[l];
            nm=k-1;
            y=w[nm];
            g=rv1[nm];
            h=rv1[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=pythag(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
            c=s=1.0;
            for (j=l;j<=nm;j++) {
                i=j+1;
                g=rv1[i];
                y=w[i];
                h=s*g;
                g=c*g;
                z=pythag(f,h);
                rv1[j]=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g = g*c-x*s;
                h=y*s;
                y *= c;
                for (jj=1;jj<=n;jj++) {
                    x=v[jj][j];
                    z=v[jj][i];
                    v[jj][j]=x*c+z*s;
                    v[jj][i]=z*c-x*s;
                }
                z=pythag(f,h);
                w[j]=z;
                if (z) {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                }
                f=c*g+s*y;
                x=c*y-s*g;
                for (jj=1;jj<=m;jj++) {
                    y=a[jj][j];
                    z=a[jj][i];
                    a[jj][j]=y*c+z*s;
                    a[jj][i]=z*c-y*s;
                }
            }
            rv1[l]=0.0;
            rv1[k]=f;
            w[k]=x;
        }
    }
    free_dvector(rv1,1,n);
}
double **dmatrix(int nrl, int nrh, int ncl, int nch)

{
    int i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
    double **m;
    m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
    m += NR_END;
    m -= nrl;
    m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
    m[nrl] += NR_END;
    m[nrl] -= ncl;
    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
    return m;
}
 
double *dvector(int nl, int nh)
{
    double *v;
    v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
    return v-nl+NR_END;
}
 
void free_dvector(double *v, int nl, int nh)
{
    free((FREE_ARG) (v+nl-NR_END));
}
 
double pythag(double a, double b)
{
    double absa,absb;
    absa=fabs(a);
    absb=fabs(b);
    if (absa > absb) return absa*sqrt(1.0+(absb/absa)*(absb/absa));
    else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
}

