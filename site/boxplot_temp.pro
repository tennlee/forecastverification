pro boxplot_temp

observed=[-1, 8,12,13,18,10,16,19,23,24]
forecast1=[ 5,10, 9,15,22,13,17,17,19,23]
forecast2=[ 4, 6,12,14,16,11,15,17,20,20]

n=n_elements(observed)
nbars=3

;Sort forecasts
ob=observed(sort(observed))
f1=forecast1(sort(forecast1))
f2=forecast2(sort(forecast2))

print,ob
print,f1
print,f2

s25=round(n/4.-.01)
s75=round(3.*n/4.-.01)
print,'s25,s75=',s25,s75

x=[-1,1,1,-1,-1]*0.2
y=[0,0,1,1,0]

minv=[min(ob),min(f1),min(f2)]
maxv=[max(ob),max(f1),max(f2)]
medv=[median(ob),median(f1),median(f2)]

q25=[ob[s25],f1[s25],f2[s25]]
q75=[ob[s75],f1[s75],f2[s75]]

;Plot boxes
plot,[0,nbars+1],[min(minv),max(maxv)],/nodata,xticks=4,/xminor,xticklen=.001, $
  xtickname=[' ','Observed','Forecast 1','Forecast 2',' '], $
  ytitle='Temperature (C)',/yminor

for i=0,nbars-1 do begin
  oplot,1+i+x,[q25[i],q25[i],q75[i],q75[i],q25[i]]
  oplot,1+i+x[0:1],[medv[i],medv[i]]
  oplot,1+[i,i],[minv[i],q25[i]]
  oplot,1+[i,i],[maxv[i],q75[i]]
  oplot,1+i+x[0:1],[minv[i],minv[i]]
  oplot,1+i+x[0:1],[maxv[i],maxv[i]]
  endfor

end