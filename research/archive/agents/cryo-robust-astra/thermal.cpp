// Thermal correlations and discretization adapted from Marc Davici's Vitrify
// (Apache-2.0, upstream/LICENSE). Independent face-list implementation.
#include <algorithm>
#include <cmath>
#include <vector>
#include <limits>
extern "C" void solve(int n,int nf,const int *faces,const int *label,const double *iron,
 double dx,double dt,int steps,const double *duty,int clamp,double h,int chunk,double *out,double *end){
 std::vector<double> T(n,-150), k(n),cap(n),flux(n),entry(n,-1),exit(n,-1);
 std::vector<int> exposed(n,6);
 for(int f=0;f<2*nf;f++) exposed[faces[f]]--;
 double minrate=1e9,inrate=1e9, energy=0, glass=0, minT=0, maxT=0, sumT=0;
 int count=0;
 for(int i=0;i<n;i++) count+=label[i]>0;
 for(int step=0;step<steps;step++){
  double power=duty[step];
  for(int i=0;i<n;i++){
   double t=T[i];
   if(label[i]){
    double cp=t< -120?7.60*t+2402.3:t<=-119?545.24*t+66747:0.039*t*t+11.73*t+2678.8;
    cap[i]=(1078-.43*t)*cp;
    k[i]=t<=-116?.000816*t+.440:t< -96?-.001*t+.234:t<=-77?.000459*t+.374:t< -34?.0005*t+.374:.000459*t+.374;
   }else{
    cap[i]=(931-.49*t-.00091*t*t)*(1031+3.57*t-.014*t*t-.000061710*t*t*t);
    k[i]=.39+.00052*t;
   }
   double effective=std::min(0.,t);
   if(clamp) effective=std::max(effective,-80.);
   flux[i]=iron[i]*1000*(691-4.6*(effective+80))*power + exposed[i]*h*(22-t)/dx;
   energy+=flux[i]*dt*dx*dx*dx;
  }
  for(int f=0;f<nf;f++){
   int a=faces[2*f], b=faces[2*f+1];
   double q=(2*k[a]*k[b]/(k[a]+k[b]))*(T[b]-T[a])/(dx*dx);
   flux[a]+=q; flux[b]-=q;
  }
  minT=1e9; maxT=-1e9; sumT=0;
  for(int i=0;i<n;i++){
   double t=T[i], v=flux[i]/cap[i], next=t+dt*v;
   if(label[i]){
    if(60*v<minrate){minrate=60*v; out[1]=step*dt; out[2]=t; out[3]=i;}
    if(t>=-100 && t<=-40 && 60*v<inrate){inrate=60*v;out[5]=step*dt;out[6]=t;out[7]=i;}
    if(t< -100 && next>=-100 && entry[i]<0)entry[i]=step*dt+(-100-t)/v;
    if(t< -40 && next>=-40 && exit[i]<0)exit[i]=step*dt+(-40-t)/v;
    minT=std::min(minT,next);maxT=std::max(maxT,next);sumT+=next;
   }
   T[i]=next;
  }
  if(!glass && (step+1)%chunk==0 && sumT/count>=-123) glass=(sumT/count+150)*60/((step+1)*dt);
 }
 double transit=1e9;int missing=0;
 for(int i=0;i<n;i++){
  end[i]=T[i];
  if(label[i] && exit[i]<0)missing++;
  if(label[i] && exit[i]>=0)transit=std::min(transit,3600/(exit[i]-entry[i]));
 }
 out[0]=minrate;out[4]=inrate;out[8]=minT;out[9]=maxT;out[10]=sumT/count;
 out[11]=transit;out[12]=missing;out[13]=energy;out[14]=glass;
}
