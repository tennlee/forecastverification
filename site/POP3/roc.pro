pro roc, n,o,h,f,area, $
  binormal=binormal,title=title,show=show,debug=debug

;Compute area under ROC curve.
;Optional: make relative operating characteristic (ROC) diagram for
;probabilistic forecast of value at each gridpoint exceeding threshold

;Input:   n            # forecasts in each fcst probability category
;         o            # occurrences in each fcst probability category
;         title        optional title for plot
;Output:  h            hit rate as a function of probability threshold
;         f            false alarm rate (POFD) "   "
;         area         area under ROC curve
;Keywds:  binormal     use binormal fit to estimate ROC area (default=not)
;         title        title of plot
;         show         show plot (default=don't)
;         debug        print values of n,o,n,f

if n_elements(binormal) eq 0 then binormal=0
if n_elements(title) eq 0 then title='Relative Operating Characteristic'
if n_elements(show) eq 0 then show=0
if n_elements(debug) eq 0 then debug=0

;Number of probability categories
ncat=n_elements(n)

;Compute hit rates and false alarm rates by probability category
non=n-o
h=make_array(ncat,/float,value=-999.)
f=make_array(ncat,/float,value=-999.)
for k=0,ncat-2 do begin
  W=total(non[0:k])
  Y=total(  o[0:k])
  Z=total(non[k+1:ncat-1])
  X=total(  o[k+1:ncat-1])
  if (X+Y) gt 0. then h[k]=X/(X+Y)
  if (Z+W) gt 0. then f[k]=Z/(Z+W)
  endfor
;Add endpoints at [0,0] and [1,1]
h[ncat-1]=0.
f[ncat-1]=0.
hs=[1.,h]
fs=[1.,f]

;Get rid of empty bins
ss=where(hs ge 0. and fs ge 0.,nok)
hs=hs[ss]
fs=fs[ss]

if (debug) then begin
  print,format='("      f=",16f8.4)',f
  print,format='("      h=",16f8.4)',h
  endif


;---------------------------------------------------------------------------

if (binormal) then begin
  print,'Using binormal estimation of ROC area'

  ;Compute area under the curve using binormal curve fitting (Ian Mason's notes)
  if nok ge 3 then begin
    ;Transform hit rate and false alarm rate to Z
    ;but don't transform values of 0 or fit won't work
    ss=where(h gt 0. and f gt 0.,npos)
    nz=npos-1
    zh=fltarr(nz)
    zf=fltarr(nz)
    for k=0,nz-1 do zh[k]=gauss_cvf(h[ss[k]])
    for k=0,nz-1 do zf[k]=gauss_cvf(f[ss[k]])
    ;Fit line to zh vs. zf
    if nz eq 1 then begin
      slope=1.              ;in 2-category case assume equal variances of PDFs
      intcpt=zh[0]-zf[0]
      endif
    if nz eq 2 then begin
      slope=(zh[1]-zh[0])/(zf[1]-zf[0])    ;2 points define line
      intcpt=zh[0]-slope*zf[0]
      endif
    if nz ge 3 then begin
      coeff= $               ;Regression fit
        regress(reform(zf,1,nz),zh,replicate(1.,nz),zhfit,const=intcpt)
      slope=coeff[0]
      print,'Binormal fit: slope, intercept=',slope,intcpt
      endif
    ;For integration take finer intervals in Z space
    zffit=-3.0+0.1*findgen(61)
    zhfit=intcpt+slope*zffit
    ;Now transform back to probability space before integrating
    hfit=fltarr(61)
    ffit=fltarr(61)
    for k=0,60 do hfit[k]=1.-gauss_pdf(zhfit[k])
    for k=0,60 do ffit[k]=1.-gauss_pdf(zffit[k])
    hfit=[1.,hfit,0.]
    ffit=[1.,ffit,0.]
    ;Now integrate under curve using trapezoid rule
    area=0.
    for k=0,61 do area=area+abs(ffit[k]-ffit[k+1])*(hfit[k]+hfit[k+1])/2.
  endif else begin
    area=0.5
  endelse
  print,'ROC area (binormal fit)=',area

endif else begin

  ;Compute area under the curve from line segments
  area=0.
  if nok ge 3 then for k=0,nok-2 do $
                   area=area+abs(fs[k]-fs[k+1])*(hs[k]+hs[k+1])/2. $
              else area=0.5
  print,'ROC area (line segments)=',area

endelse

;---------------------------------------------------------------------------

if (show) then begin
  plot,[0,1],[0,1],title=title, $
    xtitle='False Alarm Rate',ytitle='Hit Rate',/xstyle,/ystyle,linestyle=2
  if (binormal) then oplot,ffit,hfit,max_value=1., thick=2 $
                else oplot,fs,hs,max_value=1.,thick=2
  ;oplot,fs,hs,max_value=1.,psym=7
  endif

done:
end
