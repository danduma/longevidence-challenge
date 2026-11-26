from pathlib import Path
s=Path('inputs/independent_original.cpp').read_text()
s=s.replace('if(argc!=2)return 2;', 'if(argc!=3)return 2;')
a=s.index(' double minrate=')
s=s[:a]+r'''
 vector<double> initial(n);
 vector<double> amin(n,1e99),at(n),aT(n),wmin(n,1e99),wt(n),wT(n),gmin(n,1e99),gt(n),gT(n);
 vector<vector<double>>cross(4,vector<double>(n,-1));
 vector<double>bounds={-123,-100,-40,-35};
 double energy=0,bcenergy=0; double minObserved=1e99;
 auto aggregate=[&](const vector<double>&temp,const vector<double>&rate,double time){
  for(int i=0;i<n;i++){
   double v=60*rate[i];minObserved=min(minObserved,v);
   if(v<amin[i]){amin[i]=v;at[i]=time;aT[i]=temp[i];}
   if(temp[i]>=-100&&temp[i]<=-40&&v<wmin[i]){wmin[i]=v;wt[i]=time;wT[i]=temp[i];}
   if(temp[i]>=-123&&temp[i]<=-35&&v<gmin[i]){gmin[i]=v;gt[i]=time;gT[i]=temp[i];}
  }
 };
 for(int step=0;step<steps;step++){
  double time=step*dt;auto p=rhs(T,time,r);aggregate(T,r,time);if(step==0)initial=r;
  if(scheme==2){for(int i=0;i<n;i++)mid[i]=T[i]+.5*dt*r[i];p=rhs(mid,time+.5*dt,r2);aggregate(mid,r2,time+.5*dt);}else r2=r;
  energy+=dt*p.first;bcenergy+=dt*p.second;
  for(int i=0;i<n;i++){
   double next=T[i]+dt*r2[i];
   for(int j=0;j<4;j++)if(T[i]<bounds[j]&&next>=bounds[j]&&cross[j][i]<0)cross[j][i]=time+dt*(bounds[j]-T[i])/(next-T[i]);
   T[i]=next;
  }
 }
 rhs(T,duration,r);aggregate(T,r,duration);
 double H=0,iron=0;for(int i=0;i<n;i++){H+=enthalpy(T[i],lab[i])*vol;iron+=fe[i]*vol*1000;}
 ofstream out(argv[2]);out<<setprecision(17);
 out<<"index,label,fe,terminal,initial_rate,all_rate,all_time,all_temp,window_rate,window_time,window_temp,glass_rate,glass_time,glass_temp,cross123,cross100,cross40,cross35\n";
 for(int i=0;i<n;i++)out<<i<<","<<lab[i]<<","<<fe[i]<<","<<T[i]<<","<<initial[i]*60<<","<<amin[i]<<","<<at[i]<<","<<aT[i]<<","<<wmin[i]<<","<<wt[i]<<","<<wT[i]<<","<<gmin[i]<<","<<gt[i]<<","<<gT[i]<<","<<cross[0][i]<<","<<cross[1][i]<<","<<cross[2][i]<<","<<cross[3][i]<<"\n";
 cout<<setprecision(17)<<"{\"n\":"<<n<<",\"steps\":"<<steps<<",\"dt\":"<<dt<<",\"source_J\":"<<energy<<",\"boundary_J\":"<<bcenergy<<",\"enthalpy_J\":"<<H<<",\"energy_rel_error\":"<<(H-energy-bcenergy)/max(1.,abs(energy+bcenergy))<<",\"iron_g\":"<<iron<<",\"duty_s\":"<<35*(u0+2*u1+u2)/2<<",\"minimum_observed_rate\":"<<minObserved<<"}\n";
}
'''
Path('domains.cpp').write_text(s)
