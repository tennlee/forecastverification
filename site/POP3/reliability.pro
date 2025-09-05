pro reliability, probcat,n,o,BS, $
  title=title,decomp=decomp,show=show,debug=debug

;Compute Brier score, and components: reliability, resolution, uncertainty
;Optional: make reliability diagram for probabilistic forecast of
;value at each gridpoint exceeding threshold

;Input:   probcat      probability categories
;         n            # forecasts in each fcst probability category
;         o            # occurrences in each fcst probability category
;Output:  BS           Brier Score computed as rel-res+unc (Wilks, 1995)
;Keywds:  title        optional title for plot
;         decomp       [obar, reliability, resolution, uncertainty]
;         show         show reliability diagram (default=don't)
;         debug        print,n,o, components

BS=-999.
decomp=[-999.,-999.,-999.,-999.]

;Get rid of empty bins
ss=where(n gt 0.,ng)
if ng eq 0 then goto, done
nn=n[ss]
oo=o[ss]
pcat=probcat[ss]

;Compute observed probability in each fcst prob. category
oo=oo/nn

ntot=total(nn)

;Compute reliability, resolution, and uncertainty
obar=total(nn*oo)/ntot
rel =total(nn*(pcat-oo)^2.)/ntot
res =total(nn*(oo-obar)^2.)/ntot
unc =obar*(1.-obar)
BS=rel-res+unc

if n_elements(debug) eq 0 then debug=0
if (debug) then begin
  print,'obar,rel,res,unc=',obar,rel,res,unc
  print,'BS=',BS
  endif

decomp=[obar,rel,res,unc]

if n_elements(show) eq 0 then show=0
if (show) then begin
  if n_elements(title) eq 0 then title=''
  ;Plot diagonal line
  plot,[0,1],[0,1],title=title, $
    xtitle='Forecast Probability',ytitle='Observed Frequency', $
    /xstyle,/ystyle,linestyle=0,thick=1
  ;Plot reliability
  oplot,pcat,oo,max_value=1.,thick=2
  ;Plot climatological frequency
  oplot,[0,1],[obar,obar],linestyle=2,thick=1
  ;Plot no-skill line
  oplot,[0,1],[obar/2,obar+(1.-obar)/2.],linestyle=1,thick=1
  plot,[0,1],[0,max(nn)],/nodata,/noerase,position=[.35,.7,.60,.90], $
    xtitle='Fcst Probability',ytitle='# Times Fcst', $
    charsize=0.8
  xbox=[-1,1,1,-1,-1] & ybox=[0,0,1,1,0]
  for i=0,ng-1 do polyfill,pcat[i]+0.02*xbox,nn[i]*ybox,color=255,/noclip
  endif

done:
end
