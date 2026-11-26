// Independent review implementation. Equations: Joshi & Rabin 2023 Table 2,
// q = 1000 C_Fe SAR(T) u(t). Geometry/data supplied separately from pinned Vitrify.
// Each undirected face contributes equal/opposite W/m^3. Explicit midpoint RK2.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <limits>
using namespace std;
double capacity(double t,int lab){
 if(lab==0) return (931-.49*t-.00091*t*t)*(1031+3.57*t-.014*t*t-.000061710*t*t*t);
 double cp = t< -120 ? 2402.3+7.6*t : t<=-119 ? 66747+545.24*t : 2678.8+11.73*t+.039*t*t;
 return (1078-.43*t)*cp;
}
double conductivity(double t,int lab){
 if(lab==0)return .39+.00052*t;
 if(t<=-116)return .440+.000816*t;
 if(t< -96)return .234-.001*t;
 if(t<=-77)return .374+.000459*t;
 if(t< -34)return .374+.0005*t;
 return .374+.000459*t;
}
// Primitive of polynomial rho(T)*cp(T), piecewise differences from -150 C.
double polyint(double t,const vector<double>&a){double v=0,p=t;for(size_t j=0;j<a.size();j++){v+=a[j]*p/(j+1);p*=t;}return v;}
double enthalpy(double t,int lab){
 vector<double>rho=lab==0?vector<double>{931,-.49,-.00091}:vector<double>{1078,-.43};
 auto part=[&](double lo,double hi,vector<double>cp){vector<double>a(rho.size()+cp.size()-1,0);for(size_t i=0;i<rho.size();i++)for(size_t j=0;j<cp.size();j++)a[i+j]+=rho[i]*cp[j];return polyint(hi,a)-polyint(lo,a);};
 if(lab==0)return part(-150,t,{1031,3.57,-.014,-.000061710});
 double v=part(-150,min(t,-120.),{2402.3,7.6});
 if(t> -120)v+=part(-120,min(t,-119.),{66747,545.24});
 if(t> -119)v+=part(-119,t,{2678.8,11.73,.039});return v;
}
int main(int argc,char**argv){
 if(argc!=2)return 2; ifstream in(argv[1]);int n,nf,clamp,scheme;double dx,dtmax,duration,h,env,u0,u1,u2;
 in>>n>>nf>>dx>>dtmax>>duration>>clamp>>scheme>>h>>env>>u0>>u1>>u2;
 vector<int>lab(n),exposed(n,6),a(nf),b(nf);vector<double>fe(n),T(n,-150),k(n),r(n),r2(n),mid(n),entry(n,-1),leave(n,-1);
 for(int i=0;i<n;i++)in>>lab[i]>>fe[i];for(int f=0;f<nf;f++){in>>a[f]>>b[f];exposed[a[f]]--;exposed[b[f]]--;}
 if(!in)return 3;
 int steps=ceil(duration/dtmax);double dt=duration/steps,vol=dx*dx*dx;
 auto duty=[&](double time){return time<=35?u0+(u1-u0)*time/35:u1+(u2-u1)*(time-35)/35;};
 auto rhs=[&](const vector<double>&temp,double time,vector<double>&out){
  double source=0,boundary=0;
  for(int i=0;i<n;i++){double eff=min(0.,temp[i]);if(clamp)eff=max(-80.,eff);double q=1000*fe[i]*(691-4.6*(eff+80))*duty(time);double loss=exposed[i]*h*(env-temp[i])/dx;source+=q*vol;boundary+=loss*vol;out[i]=q+loss;k[i]=conductivity(temp[i],lab[i]);}
  for(int f=0;f<nf;f++){int i=a[f],j=b[f];double flux=2/(1/k[i]+1/k[j])*(temp[j]-temp[i])/(dx*dx);out[i]+=flux;out[j]-=flux;}
  for(int i=0;i<n;i++)out[i]/=capacity(temp[i],lab[i]);return pair<double,double>(source,boundary);
 };
 double minrate=1e99,window=1e99,mintime=0,mintemp=0,windowtime=0,windowtemp=0,energy=0,bcenergy=0,minsumflux=0;int minlab=0;
 auto aggregate=[&](const vector<double>&temp,const vector<double>&rate,double time){for(int i=0;i<n;i++)if(lab[i]>0){double v=rate[i]*60;if(v<minrate){minrate=v;mintime=time;mintemp=temp[i];minlab=lab[i];}if(temp[i]>=-100&&temp[i]<=-40&&v<window){window=v;windowtime=time;windowtemp=temp[i];}}};
 for(int step=0;step<steps;step++){
  double time=step*dt;auto p=rhs(T,time,r);aggregate(T,r,time);
  if(scheme==2){for(int i=0;i<n;i++)mid[i]=T[i]+.5*dt*r[i];p=rhs(mid,time+.5*dt,r2);aggregate(mid,r2,time+.5*dt);}else r2=r;
  energy+=dt*p.first;bcenergy+=dt*p.second;
  for(int i=0;i<n;i++){double next=T[i]+dt*r2[i];if(lab[i]>0){if(T[i]< -100&&next>=-100&&entry[i]<0)entry[i]=time+dt*(-100-T[i])/(next-T[i]);if(T[i]< -40&&next>=-40&&leave[i]<0)leave[i]=time+dt*(-40-T[i])/(next-T[i]);}T[i]=next;}
 }
 rhs(T,duration,r);aggregate(T,r,duration);
 double cold=1e99,hot=-1e99,transit=1e99,H=0,iron=0;int missing=0;
 for(int i=0;i<n;i++){H+=enthalpy(T[i],lab[i])*vol;iron+=fe[i]*vol*1000;if(lab[i]>0){cold=min(cold,T[i]);hot=max(hot,T[i]);if(entry[i]<0||leave[i]<0)missing++;else transit=min(transit,3600/(leave[i]-entry[i]));}}
 cout<<setprecision(14)<<"{\"n\":"<<n<<",\"steps\":"<<steps<<",\"dt\":"<<dt<<",\"min_rate\":"<<minrate<<",\"min_time\":"<<mintime<<",\"min_temp\":"<<mintemp<<",\"min_label\":"<<minlab<<",\"window_rate\":"<<window<<",\"window_time\":"<<windowtime<<",\"window_temp\":"<<windowtemp<<",\"cold\":"<<cold<<",\"hot\":"<<hot<<",\"transit\":"<<transit<<",\"missing\":"<<missing<<",\"source_J\":"<<energy<<",\"boundary_J\":"<<bcenergy<<",\"enthalpy_J\":"<<H<<",\"energy_rel_error\":"<<(H-energy-bcenergy)/max(1.,abs(energy+bcenergy))<<",\"iron_g\":"<<iron<<",\"duty_s\":"<<35*(u0+2*u1+u2)/2<<",\"pass\":"<<((minrate>=50&&cold>=-35&&hot<=0)?"true":"false")<<"}\n";
}
