// C01 adaptation of Longevidence MIT independent solver; see LICENSE.
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

struct Minimum { double rate=1e99,time=0,temp=0; int cell=-1,label=-1; };
int main(int argc,char**argv){
 if(argc!=2)return 2; ifstream in(argv[1]); int n,nf,clamp,scheme; double dx,dtmax,duration,h,env,u0,u1,u2;
 in>>n>>nf>>dx>>dtmax>>duration>>clamp>>scheme>>h>>env>>u0>>u1>>u2;
 vector<int> lab(n),exposed(n,6),a(nf),b(nf); vector<double>fe(n),T(n,-150),k(n),r(n),r2(n),mid(n),prevT(n),prevR(n);
 for(int i=0;i<n;i++)in>>lab[i]>>fe[i]; for(int f=0;f<nf;f++){in>>a[f]>>b[f];exposed[a[f]]--;exposed[b[f]]--;}
 if(!in)return 3;
 int steps=ceil(duration/dtmax);double dt=duration/steps,vol=dx*dx*dx;
 auto duty=[&](double time){return time<=35?u0+(u1-u0)*time/35:u1+(u2-u1)*(time-35)/35;};
 auto rhs=[&](const vector<double>&temp,double time,vector<double>&out){double source=0,boundary=0;
  for(int i=0;i<n;i++){double eff=min(0.,temp[i]);if(clamp)eff=max(-80.,eff);double q=1000*fe[i]*(691-4.6*(eff+80))*duty(time);double loss=exposed[i]*h*(env-temp[i])/dx;source+=q*vol;boundary+=loss*vol;out[i]=q+loss;k[i]=conductivity(temp[i],lab[i]);}
  for(int f=0;f<nf;f++){int i=a[f],j=b[f];double flux=2/(1/k[i]+1/k[j])*(temp[j]-temp[i])/(dx*dx);out[i]+=flux;out[j]-=flux;}
  for(int i=0;i<n;i++)out[i]/=capacity(temp[i],lab[i]); return pair<double,double>(source,boundary);
 };
 double lows[4]={-1e9,-123,-123,-100}, highs[4]={1e9,0,-35,-40}, thresholds[5]={-123,-100,-40,-35,0};
 string names[4]={"all","glass","glass35","window"}; Minimum mins[2][4]; vector<vector<double>>cross(5,vector<double>(n,-1)); double prevtime=-1,energy=0,bcenergy=0;
 auto update=[&](Minimum& m,double rate,double time,double temp,int i){if(rate<m.rate){m={rate,time,temp,i,lab[i]};}};
 auto aggregate=[&](const vector<double>&temp,const vector<double>&rate,double time){
  for(int i=0;i<n;i++)if(lab[i]>0){
   for(int q=0;q<5;q++)if(prevtime>=0&&cross[q][i]<0&&prevT[i]<thresholds[q]&&temp[i]>=thresholds[q])cross[q][i]=prevtime+(time-prevtime)*(thresholds[q]-prevT[i])/(temp[i]-prevT[i]);
   for(int d=0;d<2;d++)if(d==0||lab[i]==2)for(int w=0;w<4;w++){
    if(temp[i]>=lows[w]&&temp[i]<=highs[w])update(mins[d][w],rate[i]*60,time,temp[i],i);
    if(prevtime>=0&&w>0&&temp[i]!=prevT[i])for(double edge:{lows[w],highs[w]})if((prevT[i]<edge&&temp[i]>=edge)||(prevT[i]>edge&&temp[i]<=edge)){
      double z=(edge-prevT[i])/(temp[i]-prevT[i]); update(mins[d][w],60*(prevR[i]+z*(rate[i]-prevR[i])),prevtime+z*(time-prevtime),edge,i);
    }
   }
  }
  prevT=temp;prevR=rate;prevtime=time;
 };
 for(int step=0;step<steps;step++){
  double time=step*dt;auto p=rhs(T,time,r);aggregate(T,r,time);
  for(int i=0;i<n;i++)mid[i]=T[i]+.5*dt*r[i];p=rhs(mid,time+.5*dt,r2);aggregate(mid,r2,time+.5*dt);
  energy+=dt*p.first;bcenergy+=dt*p.second;for(int i=0;i<n;i++)T[i]+=dt*r2[i];
 }
 rhs(T,duration,r);aggregate(T,r,duration);double H=0,iron=0;
 for(int i=0;i<n;i++){H+=enthalpy(T[i],lab[i])*vol;iron+=fe[i]*vol*1000;}
 cout<<setprecision(15)<<"{\"n\":"<<n<<",\"steps\":"<<steps<<",\"dt\":"<<dt<<",\"source_J\":"<<energy<<",\"boundary_J\":"<<bcenergy<<",\"enthalpy_J\":"<<H<<",\"energy_rel_error\":"<<(H-energy-bcenergy)/max(1.,abs(energy+bcenergy))<<",\"iron_g\":"<<iron<<",\"duty_s\":"<<35*(u0+2*u1+u2)/2;
 for(int d=0;d<2;d++){
  string pre=d==0?"specimen":"myocardium";double cold=1e99,hot=-1e99,transit=1e99;int missing=0;int misses[5]={0,0,0,0,0};
  for(int i=0;i<n;i++)if(lab[i]>0&&(d==0||lab[i]==2)){
    cold=min(cold,T[i]);hot=max(hot,T[i]);for(int q=0;q<5;q++)if(cross[q][i]<0)misses[q]++;
    if(cross[1][i]<0||cross[2][i]<0)missing++;else transit=min(transit,3600/(cross[2][i]-cross[1][i]));
  }
  cout<<",\""<<pre<<"\":{\"cold\":"<<cold<<",\"hot\":"<<hot<<",\"terminal_pass\":"<<((cold>=-35&&hot<=0)?"true":"false")<<",\"transit_100_40\":"<<transit<<",\"missing_100_40\":"<<missing<<",\"missing_crossings\":[";
  for(int q=0;q<5;q++)cout<<(q?",":"")<<misses[q];cout<<"]";
  for(int w=0;w<4;w++){auto m=mins[d][w];cout<<",\""<<names[w]<<"\":{\"rate\":"<<m.rate<<",\"time\":"<<m.time<<",\"temp\":"<<m.temp<<",\"cell\":"<<m.cell<<",\"label\":"<<m.label<<",\"crossing_123_time\":"<<(m.cell>=0?cross[0][m.cell]:-1)<<"}";}
  cout<<"}";
 }
 cout<<"}\n";
}
