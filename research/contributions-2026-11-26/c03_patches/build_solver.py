from pathlib import Path
s=Path('snapshots/independent.cpp').read_text()
s=s.replace('if(argc!=2)return 2;', 'if(argc!=3)return 2;')
s=s.replace('int steps=ceil', 'vector<double> glass(n,-1), endglass(n,-1), postrate(n,1e99), winrate(n,1e99);\n int steps=ceil')
s=s.replace('double v=rate[i]*60;', 'double v=rate[i]*60; if(temp[i]>=-123&&temp[i]<=-35)postrate[i]=min(postrate[i],v); if(temp[i]>=-100&&temp[i]<=-40)winrate[i]=min(winrate[i],v);')
s=s.replace('if(lab[i]>0){if(T[i]< -100', 'if(lab[i]>0){if(T[i]<-123&&next>=-123&&glass[i]<0)glass[i]=time+dt*(-123-T[i])/(next-T[i]); if(T[i]<-35&&next>=-35&&endglass[i]<0)endglass[i]=time+dt*(-35-T[i])/(next-T[i]); if(T[i]< -100')
s=s.replace('cout<<setprecision(14)', 'ofstream vox(argv[2]);vox<<setprecision(17)<<"label,fe,terminal,glass,endglass,entry,leave,postrate,winrate\\n"; for(int i=0;i<n;i++)vox<<lab[i]<<","<<fe[i]<<","<<T[i]<<","<<glass[i]<<","<<endglass[i]<<","<<entry[i]<<","<<leave[i]<<","<<postrate[i]<<","<<winrate[i]<<"\\n";\n cout<<setprecision(14)')
Path('solver.cpp').write_text('// C03 adaptation: per-voxel terminal, crossing and restricted-rate diagnostics. Original MIT notice in snapshots/LONGEVIDENCE-LICENSE.\n'+s)
