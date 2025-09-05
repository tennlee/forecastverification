pro pstats,obsocc,pop,refpop,ncat, $
   pod=pod,pofd=pofd, show=show

;Compute probabilistic verification statistics, including skill scores, on
;pooled data.
;
;Input:  obsocc      observed occurrence [npts]
;        pop         forecast probability of precipitation [npts]
;                    Note: negative values of obs or prob treated as bad data!
;        refpop      reference probability of precipitation (climatology) [npts]
;        ncat        number of probability categories (generally equal to
;                      number of ensemble members + 1 )
;Keywords:
;        pod, pofd   output probability of detection, probability of false
;                      detection [ncat,nthresh]
;        show        show='reliability', 'ROC', or 'value'

if n_elements(show) eq 0 then show=''

probcat=1./(ncat-1.)*findgen(ncat) & print,'probcat=',probcat
binsize=probcat[1]-probcat[0]
halfbinsize=binsize/2.

;Check compatibility of array sizes
npts=n_elements(obsocc)
if n_elements(pop) ne npts or n_elements(refpop) ne npts $
  then begin
  print,'pstats -- size mismatch between obsocc, pop, refpop, thresh -- quit'
  goto, done
  endif

;Keep only good observations (i.e., >=0; masked values=-1)
ss=where(obsocc ge 0,npts)
obs=obsocc[ss]
prob=pop[ss]
refprob=refpop[ss]

;Output POD and POFD
pod=fltarr(ncat+1)
pofd=fltarr(ncat+1)

;Working arrays for probabilistic verification
n   =make_array(ncat,/float,value=0.)
o   =make_array(ncat,/float,value=0.)
nref=make_array(ncat,/float,value=0.)
oref=make_array(ncat,/float,value=0.)

;Bin by forecast probability category
for k=0,ncat-1 do begin
  ;Which forecasts fall into probability category k
  ss=where(abs(prob-probcat[k]) lt halfbinsize,num)
  if num gt 0 then begin
    n[k]=num              ;# fcsts in category k
    o[k]=total(obs[ss])   ;# observed occurrences in category k
    endif
  sr=where(abs(refprob-probcat[k]) lt halfbinsize,num)
  if num gt 0 then begin
    nref[k]=num
    oref[k]=total(obs[sr])
    endif
  endfor


;Compute Brier score manually
BS=total((prob-obs)^2.)/npts

;Brier score decomposition into reliability, resolution, uncertainty
reliability,probcat,n,o,BSd,decomp=decomp,title='POP!Bhi',show=(show eq 'reliability')
rel=decomp[1]
res=decomp[2]
unc=decomp[3]
print,'BS-BSd=',BS-BSd

;Brier score for ref fcst (climatology), Brier skill score
BSref=total((refprob-obs)^2.)/npts
BSS=(BSref-BS)/BSref

;Area under ROC curve
roc,n,o,h,f,rocarea,/debug,show=(show eq 'ROC')
;roc,n,o,h,f,rocarea,/debug,/binormal,show=(show eq 'ROC')
pod=[1.,h]
pofd=[1.,f]

;Value plot
pclim=total(oref)/total(nref)
value,probcat,n,o,pclim,show=(show eq 'value')


print,'n=',n
print,'BS=',BS
print,'reliability=',rel
print,'resolution=',res
print,'uncertainty=',unc
print,'BSS=',BSS
print,'ROCarea=',rocarea

done:
end
