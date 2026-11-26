"""Check rates exactly at interpolated local -100/-40 C crossings.
Generate an instrumented copy of our solver; no change to reviewed artifacts.
"""
from pathlib import Path
import subprocess,json
from run_review import run
base=Path('independent.cpp').read_text()
base=base.replace('if(!in)return 3;', '''if(!in)return 3;
 vector<vector<int>> neighbors(n);for(int f=0;f<nf;f++){neighbors[a[f]].push_back(b[f]);neighbors[b[f]].push_back(a[f]);}
 double crossing_min=1e99;
''')
base=base.replace('for(int i=0;i<n;i++){double next=T[i]+dt*r2[i];', '''
  // The full next field must be available while evaluating an individual crossing.
  for(int i=0;i<n;i++)if(lab[i]>0){
   double next=T[i]+dt*r2[i];
   for(double threshold:{-100.,-40.})if(T[i]<threshold && next>=threshold){
    double fraction=(threshold-T[i])/(next-T[i]);double eff=min(0.,threshold);if(clamp)eff=max(-80.,eff);
    double heat=1000*fe[i]*(691-4.6*(eff+80))*duty(time+fraction*dt);
    double ki=conductivity(threshold,lab[i]);
    for(int j:neighbors[i]){double tj=T[j]+fraction*dt*r2[j];double kj=conductivity(tj,lab[j]);heat+=2/(1/ki+1/kj)*(tj-threshold)/(dx*dx);}
    heat+=exposed[i]*h*(env-threshold)/dx;
    crossing_min=min(crossing_min,60*heat/capacity(threshold,lab[i]));
   }
  }
  for(int i=0;i<n;i++){double next=T[i]+dt*r2[i];''')
base=base.replace('cout<<setprecision(14)<<', 'cout<<setprecision(14)<<"{\\\"crossing_rate\\\":"<<crossing_min<<","<<')
# The original JSON opens with {; remove that second brace.
base=base.replace('<<"{\\\"n\\\":', '<<"\\\"n\\\":')
Path('independent_crossing.cpp').write_text(base)
subprocess.run(['c++','-O3','-std=c++17','independent_crossing.cpp','-o','independent_crossing'],check=True)
rows=[]
for mode,mean in [('extrapolated',1.47),('clamped',1.83)]:
 for design in ['ramp','flat']:
  rows.append(run(.5,design,mean=mean,mode=mode,solver='independent_crossing'))
# A stronger heat-leak perturbation avoids interpreting the marginal h=5 endpoint.
rows.append(run(.5,'flat',mean=1.83,mode='clamped',h=10,env=-150,solver='independent_crossing'))
Path('results/crossing.json').write_text(json.dumps(rows,indent=2)+'\n')
