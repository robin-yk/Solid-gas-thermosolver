/* Note 2b: fixed inventory, ideal vacancy and Ti3+ mixing, y = 4x. */
(function(root){
  'use strict';
  function solve(p) {
    var N=p.inventory, R=p.radius, T=p.temperature, A=p.A, xi=p.xi;
    var mean=N*79.866e-6/2, kt=8.617333262145e-5*T;
    if (![N,R,T,A,xi].every(Number.isFinite) || !(mean>0 && mean<.25 && R>0 && T>0 && A>=0 && xi>0))
      throw Error('Use an inventory between 0 and 6260 µmol O g⁻¹, positive radius, temperature and decay length, and non-negative stabilization.');
    var edge=Math.min(R,Math.max(20*xi,2)), edges=[Math.min(2,R)];
    for(var i=0;i<=2000;i++) edges.push(edge*i/2000);
    for(i=0;i<=200;i++) edges.push(edge+(R-edge)*i/200);
    edges=Array.from(new Set(edges)).sort((a,b)=>a-b);
    var depth=[], weights=[], energy=[];
    for(i=1;i<edges.length;i++) {
      var z=(edges[i-1]+edges[i])/2;
      depth.push(z); energy.push(-A*Math.exp(-z/xi));
      weights.push((Math.pow(R-edges[i-1],3)-Math.pow(R-edges[i],3))/Math.pow(R,3));
    }
    function occ(mu,E) {
      var lo=0,hi=.25;
      for(var j=0;j<56;j++) {
        var x=(lo+hi)/2, h=Math.log(x/(1-x))+2*Math.log(4*x/(1-4*x));
        if(E+kt*h<mu)lo=x; else hi=x;
      }
      return (lo+hi)/2;
    }
    var lo=-A-100*kt,hi=100*kt,mu;
    for(var n=0;n<55;n++) {
      mu=(lo+hi)/2;
      var sum=0; for(i=0;i<depth.length;i++)sum+=weights[i]*occ(mu,energy[i]);
      if(sum<mean)lo=mu;else hi=mu;
    }
    mu=(lo+hi)/2;
    var shell=0,total=0;
    for(i=0;i<depth.length;i++) {var q=weights[i]*occ(mu,energy[i]);total+=q;if(depth[i]<Math.min(2,R))shell+=q;}
    var span=Math.min(R,Math.max(10,10*xi)), profile=[];
    for(i=0;i<=400;i++){var d=span*i/400;profile.push([d,occ(mu,-A*Math.exp(-d/xi))]);}
    return {surface:occ(mu,-A),interior:occ(mu,-A*Math.exp(-R/xi)),mean:mean,shell:shell/mean,closure:total-mean,mu:mu,profile:profile,span:span};
  }
  root.VacancyDistribution={solve:solve};
  if(typeof module!=='undefined')module.exports={solve:solve};
})(typeof window!=='undefined'?window:globalThis);
