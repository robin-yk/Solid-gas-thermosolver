(function(){
  'use strict';
  var $=id=>document.getElementById(id), latest;
  function draw(){
    try {
      var p={inventory:+$('vdInventory').value,radius:+$('vdRadius').value,temperature:+$('vdTemperature').value+273.15,A:+$('vdA').value,xi:+$('vdXi').value};
      var r=VacancyDistribution.solve(p); latest=r;
      var K=FigKit,f=K.square(),a=f.pane;
      var X=K.lin(0,r.span,a.x0,a.x1),Y=K.lin(0,25,a.y1,a.y0);
      K.frame(f,[a.x0,a.x1],[a.y0,a.y1]);
      K.axisX(f,X,a.y1,[0,.25,.5,.75,1].map(v=>v*r.span),'Depth from surface (nm)',v=>String(+v.toFixed(2)));
      K.axisY(f,Y,a.x0,[0,5,10,15,20,25],'Vacant oxygen sites (%)',String);
      f.line(X(0),Y(100*r.mean),X(r.span),Y(100*r.mean),'#777777',1.5,'5,4');
      f.path(r.profile.map((q,i)=>(i?'L':'M')+X(q[0])+','+Y(100*q[1])).join(' '),'#0072B2',2.5);
      K.legend(f,a.x0+100,a.y0+20,[{col:'#0072B2',text:'Equilibrium profile'},{col:'#777777',text:'Uniform inventory',dash:'5,4'}]);
      $('vdFigure').innerHTML=f.done();
      $('vdResults').textContent='Surface: '+(100*r.surface).toFixed(2)+'% vacant oxygen sites. Interior: '+(100*r.interior).toFixed(3)+'%. Outermost '+Math.min(2,p.radius)+' nm: '+(100*r.shell).toFixed(1)+'% of all vacancies.';
      $('vdCheck').textContent='Inventory closure: '+Math.abs(r.closure/r.mean).toExponential(1)+' relative error (acceptance < 10⁻⁸).';
      $('vdError').textContent='';
    }catch(e){latest=null;$('vdError').textContent=e.message;$('vdFigure').innerHTML='';$('vdResults').textContent='';$('vdCheck').textContent='';}
  }
  function init(){
    ['vdInventory','vdRadius','vdTemperature','vdA','vdXi'].forEach(id=>$(id).addEventListener('change',draw));
    $('vdSample').addEventListener('change',function(){if(this.value){$('vdInventory').value=this.value;draw();}});
    $('vdInventory').addEventListener('change',()=>$('vdSample').value='');
    $('vdSvg').addEventListener('click',function(){if(!latest)return;var url=URL.createObjectURL(new Blob([$('vdFigure').innerHTML],{type:'image/svg+xml'})),a=document.createElement('a');a.href=url;a.download='vacancy-distribution.svg';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);});
    draw();
  }
  if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',init);else init();
})();
