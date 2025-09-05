function rps,prob,obsocc

;Computes the Ranked Probability Score RPS, which is ~the mean squared
;difference between the observed and forecast CDFs.
;Reference: Stanski et al, 1989
;
;Input:  prob    forecast probability in each category
;                  dimensions (nc,nr,nbins) or (npts,nbins)
;        obsocc  observed occurrence in each category
;                  dimensions (nc,nr,nbins) or (npts,nbins)
;Output: rps     ranked probability score

sz=size(obsocc)
ndim=sz[0]
nbins=sz[ndim]

CDFo=fltarr(nbins)
CDFf=fltarr(nbins)

;Want to work with arrays of dimension (npts,nbins)
if ndim eq 2 then begin
  npts=sz[1]
  ob=obsocc
  prb=prob
endif else begin
  npts=sz[1]*sz[2]
  ob=reform(obsocc,npts,nbins)
  prb=reform(prob,npts,nbins)
endelse

;Loop on data
rps=0.
n=0.
CDFo=fltarr(nbins)
CDFf=fltarr(nbins)
for i=0,npts-1 do begin

  ;Define the CDF for this observation
  o=reform(ob[i,*],nbins)
  CDFo[0]=o[0]
  for j=1,nbins-1 do CDFo[j]=CDFo[j-1]+o[j]

  ;Invert the fcst prob array to create CDF
  p=reform(prb[i,*],nbins)
  CDFf[0]=p[0]
  for j=1,nbins-1 do CDFf[j]=CDFf[j-1]+p[j]

  ;Compute mean squared difference (nearly)
  msd=total((CDFf-CDFo)^2) / (nbins-1.)
;print,'o   =',o,'   CDFo=',CDFo,'   p   =',p,'   CDFf=',CDFf

  rps=rps+msd
  n=n+1

  continue:
  endfor

;Get mean RPS over all data
rps=rps/n

return,rps
end
