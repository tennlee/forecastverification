pro value,probcat,n,o,pclim,show=show

;Plot relative value as a function of cost/loss ratio
;For a probability forecast plots the envelope of curve maxima

;Input:   probcat      probability categories
;         n            # forecasts in each fcst probability category
;         o            # occurrences in each fcst probability category
;         pclim        climatological probability (base rate)
;Keywds:  show         show plot (default=don't)

if n_elements(show) eq 0 then show=0

non=n-o
npts=total(n)
ncat=n_elements(probcat)

;Compute value as a function of cost/loss ratio between 0 and 1
nb=50    ;number of increments between 0 and 1
CLratio=1./(nb-1.)*findgen(nb)
value=fltarr(nb)
maxvalue=make_array(nb,/float,value=-999.)  ;Envelope of value curves

if (show) then plot,[0,1],[0,0],title=title,yrange=[0,1], $
      xtitle='Cost/Loss Ratio',ytitle='Relative Value',linestyle=2

for icat=0,ncat-2 do begin
  p=probcat(icat+1)
  ;;;;p=pclim   ;Used for deterministic forecast value assessment
  k=max([0,where(probcat lt p)])
  Y=total(  o[0:k])/npts          ;misses
  Z=total(non[k+1:ncat-1])/npts   ;false alarms
  X=total(  o[k+1:ncat-1])/npts   ;hits
  for ib=0,nb-1 do begin
    if CLratio[ib] lt pclim then $
      value[ib]=(Clratio[ib]*(X+Z-1.)+Y) / (CLratio[ib]*(pclim-1.))
    if CLratio[ib] ge pclim then $
      value[ib]=(Clratio[ib]*(X+Z)+Y-pclim) / (pclim*(CLratio[ib]-1.))
    if value[ib] gt maxvalue[ib] then maxvalue[ib]=value[ib]
    endfor
  if (show) then oplot,CLratio,value,linestyle=0,thick=0.5
  endfor

if (show) then oplot,CLratio,maxvalue,thick=2

end