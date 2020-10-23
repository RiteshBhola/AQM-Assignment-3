#include<stdlib.h>
#include<stdio.h>
#include<assert.h>
#include<math.h>
#include</home/ritesh/display.h>
double correlation(double *p,int N_t,int dt,float m)
{
  double sum=0.0;
  int j;
  for(int i=0;i<N_t;++i)
  {
    j=i+dt;
    if(j>N_t) j=(i+dt)%N_t;
    sum+=p[i]*p[j];
    
  }
  return(sum/N_t);
}
void init_genrand64(unsigned long long);
double genrand64_real2(void);
void main()
{ 
  init_genrand64(5785);
  int N_t=1200,N_sweep=100000;
  int N_sep=500;
  float h,m,w;
  m=w=.1;
  h=5;
  double temp;
  char fname1[200],fname2[200];
  //sprintf(fname1,"scalar_N_t%d_m%0.6f.dat",N_t,m);
  sprintf(fname2,"correlation_N_t%d_m%0.6f.dat",N_t,m);
  //FILE* f1=fopen(fname1,"w");
  FILE* f2=fopen(fname2,"w");
  
  double x2=0.0,x3=0.0,x4=0.0;
  
  int t_plu=0,t_min=0,i,t;
  float x_new=0.0, s_old=0.0, s_new=0.0;
  float* randm=(float*)malloc(2*N_t*sizeof(float));
  double* path=(double*)malloc(N_t*sizeof(double));
  int* index =(int*)malloc(N_t*sizeof(int));
  float accrate=0.0;
  for (i = 0; i < N_t; i += 1)
  {
    path[i]=0.0;
  }
  
  for (int k = 0; k < N_sweep; k += 1)
  {
    accrate=0.0;
/*    fdisp2("h",h,"accrate",accrate);*/
  for (i = 0; i < N_t; i += 1)
  {
    index[i]=(long long)((genrand64_real2())*N_t);
/*    disp1("index",index[i]);*/
  }
  for (i = 0; i < 2*N_t; i += 1)
  {
    randm[i]=genrand64_real2();
/*    fdisp1("randm",randm[i]);*/
  }
  
  for (i = 0; i < N_t; i += 1)
  {
    t=index[i];
    t_min=(t+ N_t -1)%N_t;
    t_plu=(t+1)%N_t;
    x_new=path[t] + h*(randm[i]-0.5);
/*    fdisp2("h",h,"rand",randm[t]-0.5);*/
    s_old= (m/2.0)*(path[t_plu]-path[t])*(path[t_plu]-path[t]) + (m/2.0)*(path[t]-path[t_min])*(path[t]-path[t_min]) + (m/2.0)*w*w*path[t]*path[t];
    s_new= (m/2.0)*(path[t_plu]-x_new)*(path[t_plu]-x_new) + (m/2.0)*(x_new-path[t_min])*(x_new-path[t_min]) + (m/2.0)*w*w*x_new*x_new;
    
    if(randm[N_t + i] < exp(-s_new + s_old))
    {
      path[t]=x_new;
      accrate=accrate+(1.0/N_t);
/*      fdisp2("h",h,"accrate",accrate);*/
    }
  }
  
  h=(h*accrate)/0.8;
/*  fdisp2("h",h,"accrate",accrate);*/
  x2=0.0;
  x3=0.0;
  /*if(k%N_sep==0&&k>N_sweep/10)
  {
    for (i = 0; i < N_t; i += 1)
  { 
    temp=path[i]*path[i];
    x2+=temp;
    x3+=path[i]*temp;
    x4+=temp*temp;
  }
  
    x2=(x2*m*m)/N_t;
    x3=(x3*m*m*m)/N_t;
    x4=(x4*m*m*m*m)/N_t;
    fprintf(f1,"%d %0.16f %0.16f %0.16f\n",k,x2,x3,x4);
  }*/
    for (int dt =1;dt<=80;++dt)
    {
      fprintf(f2,"%0.16f ",correlation(path,N_t,dt,m));
    }
      fprintf(f2,"\n");
  }
  //fclose(f1);
  fclose(f2);
}
