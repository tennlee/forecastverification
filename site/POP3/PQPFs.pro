pro PQPFs,models,forecast,sdate,edate,threshold, $
  btime=btime,ROCdiag=ROCdiag, random=random, best=best, minvol=minvol, $
  summer=summer, winter=winter, N=N, SE=SE

;Accumulate information on forecast and observed probabilities
;so that reliability diagram and Brier Skill Score can be computed
;
;Input:  model           array of model names
;                          ex: ['LAPS375','GASP','USAVM','UKGC','JMA','CMC','DWD']
;        forecast        '00_24','24_48', or ['00_24','24_48']
;        sdate,edate     (may be irrelavent -- see code)
;        threshold       single value or array of rain rates
;                              Expects thresholds of [1,2,5,10,25,50]
;Output: n               number of occurences of each probability category
;        o               number of observed events in each probability category
;Keywords:
;        btime           base time of forecast (00 or 12 (UTC), default=0)
;        /ROCdiag        plot a relative operating characteristic diagram
;        random          if used, # of models to combine randomly (10 samples/d)
;        best            if used, array of models to use (ex. ['UKGC','USAVM'])
;        minvol          Set to minimum volume of obs rain (isolate big storms)
;        /N or /SE       get results for northern (-20,115) to (-10,150)
;                          or southeastern (-45,135) to (-25,155) domain

;Computes statistics based on all data pooled, as well as daily average values

missng=1.e36

if n_params() lt 1 then models=['LAPS375','GASP','USAVM','UKGC','JMA','CMC','DWD']
if n_params() lt 2 then forecast='00_24'
if n_params() lt 3 then sdate=0
if n_params() lt 4 then edate=sdate
if n_params() lt 5 then threshold=[1,2,5,10,25,50]

if n_elements(btime) eq 0 then btime=0

model=models
;model=['LAPS','GASP','USAVM','UKGC','JMA','ECMWF','DWD']
;model=['LAPS375','GASP','USAVM','UKGC','JMA','CMC','DWD']
nmodels=n_elements(model)
focst=replicate(forecast[0],nmodels)

if n_elements(best) gt 0 then model=best

if n_elements(forecast) eq 2 then begin
  model=[model,model]
  focst=[focst,replicate(forecast[1],nmodels)]
  endif

if n_elements(random) eq 0 then random=0
if n_elements(minvol) eq 0 then minvol=0
season='annual'
if keyword_set(summer) then season='summer'
if keyword_set(winter) then season='winter'
region='AUST'
if keyword_set(N) then region='N'
if keyword_set(SE) then region='SE'

nthresh=n_elements(threshold)
nmodels=n_elements(model)
members=indgen(nmodels)

nuse=make_array(nmodels,/float,value=0.)
sumPS=make_array(nthresh,/float,value=0.)
sumPSref=make_array(nthresh,/float,value=0.)

nmem=nmodels
if random gt 0 then nmem=random
ncat=nmem+1
probcat=1./nmem*findgen(ncat)
print,format='("probcat=",16f8.4)',probcat

nc=45 & nr=35
pp=fltarr(nc,nr)
fcst=fltarr(nc,nr,nmodels)
obs=fltarr(nc,nr)
prob=fltarr(nc,nr,nthresh)
ref_fcst=fltarr(nc,nr,6)     ;Expects thresholds of [1,2,5,10,25,50]
ob=fltarr(nc,nr,nthresh)

n=make_array(ncat,nthresh,/float,value=0.)
o=make_array(ncat,nthresh,/float,value=0.)
ratio=make_array(ncat,nthresh,/float,value=-0.9999)
nref=make_array(ncat,nthresh,/float,value=0.)
oref=make_array(ncat,nthresh,/float,value=0.)
nd=make_array(ncat,nthresh,/float,value=0.)
od=make_array(ncat,nthresh,/float,value=0.)

mask=bytarr(nc,nr)
maskfile='mask_1.00.dat'
openr,unit7,maskfile,/get_lun
readu,unit7,mask
close,unit7 & free_lun,unit7
npts=total(mask)

;If /N or /SE then mask points outside those domains
if region eq 'N' then lims=[-20.,115.,-10.,150.] 
if region eq 'SE' then lims=[-45.,135.,-25.,155.] 
if region eq 'N' or region eq 'SE' then begin
  domainmask=mask*0
  lats=-45.+1.25*findgen(nr)
  lons=100.+1.25*findgen(nc)
  sslat=where(lats ge lims[0] and lats le lims[2],nla)
  sslon=where(lons ge lims[1] and lons le lims[3],nlo)
  domainmask[sslon[0]:sslon[nlo-1],sslat[0]:sslat[nla-1]]=1
  mask=mask*domainmask
  endif

maxdays=1000*10
dailySS  =make_array(maxdays,nthresh,/float,value=missng)
dailyarea=make_array(maxdays,nthresh,/float,value=missng)

factor=111.*0.8*111.*1.e-6     ;For getting rain volume

;For whole time series
;For MWR Poor Man's paper
;sdate=[19970901,19980901,19990901]
;edate=[19980831,19990831,19991231]
;For BMRC 2002 Modelling Workshop
;sdate=[20010901]
;edate=[20021031]

if season eq 'summer' then begin
;  sdate=[19971201,19981201,19991201]
;  edate=[19980228,19990228,20000229]
;For MWR Poor Man's paper
  sdate=[19971201,19981201,19991201]
  edate=[19980228,19990228,19991231]
endif
if season eq 'winter' then begin
;  sdate=[19980601,19990601,20000601]
;  edate=[19980831,19990831,20000831]
;For MWR Poor Man's paper
  sdate=[19980601,19990601]
  edate=[19980831,19990831]
endif
nyr=n_elements(sdate)

id=0L

for iyr=0,nyr-1 do begin

sday=jday(sdate[iyr])
ndys=ndays(sdate[iyr],edate[iyr])+1
rmonth=0

for idy=0,ndys-1 do begin
  day=sday+idy
  vdate=yyyymmdd(day)
  vmonth=fix((vdate mod 10000)/100)
  cmdate=strtrim(string(yyyymmdd(day)),2)
  cldate=strtrim(string(yyyymmdd(day-1)),2)        ;LAPS is 1 day behind @ 23Z
  c48date=strtrim(string(yyyymmdd(day-1)),2)        ;24-48 hour forecast
  cl48date=strtrim(string(yyyymmdd(day-2)),2)        ;LAPS is 1 day behind @ 23Z
  csdate=strtrim(string(yyyymmdd(day+1)),2)        ;analysis is 1 day ahead

print,'Getting forecast and observations for ',vdate

  for i=0,nmodels-1 do begin
    modelname=model[i]
    if model[i] eq 'GASP' and vdate ge 19990202 and vdate le 20001010 then $
      modelname='HRGASP'
;    if model[i] eq 'LAPS375' then modelname='LAPS'  ;not for BMRC2002ModWkshp

    read_qpf1deg,modelname,focst[i],vdate,pp,btime=btime,error=error

;****************************************************************************

    if error ne 0 then begin
      ;If problem is that it can't find HRGASP then look for GASP
      if modelname eq 'HRGASP' then begin
        read_qpf1deg,'GASP',focst[i],vdate,pp,btime=btime,error=error
        if error eq 0 then goto, read_data
        endif
      goto,continue    ;Don't use if fewer than nmodels models present!
      endif

;****************************************************************************

read_data:
    pp=pp*mask-(mask eq 0)
    fcst[*,*,i]=pp
  endfor
  fcst_save=fcst

  analfile='/bm/gdata/eee/rainval/1deg/qpf_anals/anal.'+csdate
  openr,unit,analfile,error=error,/get_lun
  if (error) then goto,continue    
  readu,unit,obs
  close,unit & free_lun,unit
  obs=obs*mask
  volume=total(obs)*factor
  obs=obs*mask-(mask eq 0)

  if volume lt minvol then goto, continue
  if (minvol) then print,vdate,'  Observed volume=',volume

  ;If using a random selection, choose models (10 times per day)
  ndone_random=0
  choose_again:
  if random gt 0 then begin
    members=randomx(nmodels,random)
    newfcst=fcst
    for m=0,random-1 do newfcst[*,*,m]=fcst_save[*,*,members[m]]
    fcst=newfcst
    nuse[members]=nuse[members]+1
    endif

  ;Get reference forecast (climatology)
  if vmonth ne rmonth then begin
    reffile='PQPF_ref.dat'
    case vmonth of
      1: reffile='PQPF_DJF_ref.dat'
      2: reffile='PQPF_JFM_ref.dat'
      3: reffile='PQPF_FMA_ref.dat'
      4: reffile='PQPF_MAM_ref.dat'
      5: reffile='PQPF_AMJ_ref.dat'
      6: reffile='PQPF_MJJ_ref.dat'
      7: reffile='PQPF_JJA_ref.dat'
      8: reffile='PQPF_JAS_ref.dat'
      9: reffile='PQPF_ASO_ref.dat'
     10: reffile='PQPF_SON_ref.dat'
     11: reffile='PQPF_OND_ref.dat'
     12: reffile='PQPF_NDJ_ref.dat'
    endcase
    print,'Getting reference (climatological frequency) from ',reffile
    openr,unit,reffile,/get_lun
    readu,unit,ref_fcst
    close,unit & free_lun,unit
    for ithr=0,5 do ref_fcst[*,*,ithr]=ref_fcst[*,*,ithr]*mask
    rmonth=vmonth
    endif

;Get probabilities of event occurring at each gridpoint
prob[*,*,*]=0.
for m=0,nmem-1 do begin
  for ithr=0,nthresh-1 do begin
    prob[*,*,ithr]=prob[*,*,ithr]+(fcst[*,*,m] ge threshold[ithr])
  endfor
endfor
prob=prob/nmem

;Get binary occurrence of observation and non-observation
ob[*,*,*]=0.
for ithr=0,nthresh-1 do begin
  ob[*,*,ithr]=(obs ge threshold[ithr])
endfor

nd[*,*]=0.
od[*,*]=0.

;Bin by forecast probability category 
binsize=probcat[1]-probcat[0]
for ithr=0,nthresh-1 do begin
obst=ob[*,*,ithr]
 for k=0,ncat-1 do begin
  ss=where(abs(prob[*,*,ithr]-probcat[k]) lt .01 and mask eq 1,num)
  if num gt 0 then begin
    nd[k,ithr]=num                           ;# fcsts in category k (daily)
    od[k,ithr]=total(obst[ss])               ;observed occurrences  (daily)
    n[k,ithr]=n[k,ithr]+nd[k,ithr]           ;# fcsts in category k (pooled)
    o[k,ithr]=o[k,ithr]+od[k,ithr]           ;observed occurrences  (pooled)
    endif
  sr=where(abs(ref_fcst[*,*,ithr]-probcat[k]) lt binsize/2 and mask eq 1,numref)
  if numref gt 0 then begin
    nref[k,ithr]=nref[k,ithr]+numref         ;# ref fcsts in category k (pooled)
    oref[k,ithr]=oref[k,ithr]+total(obst[sr])
    endif
 endfor
endfor

;Compute daily statistics
for ithr=0,nthresh-1 do begin
 odt=od[*,ithr]
 ndt=nd[*,ithr]
 non=ndt-odt
 h=make_array(ncat,/float,value=0.)
 f=make_array(ncat,/float,value=0.)
 for k=0,ncat-1 do begin
  ;Fill contingency table
  if k le ncat-2 then begin
    W=total(non[0:k])
    Y=total(odt[0:k])
    Z=total(non[k+1:ncat-1])
    X=total(odt[k+1:ncat-1])
    if (X+Y) gt 0. then h[k]=X/(X+Y)
    if (Z+W) gt 0. then f[k]=Z/(Z+W)
    endif
  endfor
  h[ncat-1]=0
  f[ncat-1]=0

  ;Get rid of empty bins
  ss=where(ndt gt 0.)
  ndt=ndt[ss]
  odt=odt[ss]
  h=h[ss]
  f=f[ss]
  probcatd=probcat[ss]
  h=[1.,h]
  f=[1.,f]
  ncatd=n_elements(h)

  ;Compute Brier skill score explicitly
  PS=total((prob[*,*,ithr]-ob[*,*,ithr])^2.)/npts
  PSref=total((ref_fcst[*,*,ithr]-ob[*,*,ithr])^2.)/npts
  SS=(PSref-PS)/PSref
  sumPS[ithr]   =sumPS[ithr]   +PS
  sumPSref[ithr]=sumPSref[ithr]+PSref

  ;Compute area under ROC curve
  if ncatd gt 2 then begin
    dx=abs(f[1:ncatd-1]-f[0:ncatd-2])
    area=0.
    for k=0,ncatd-2 do area=area+dx[k]*(h[k]+h[k+1])/2.
  endif else area=missng

dailySS[id,ithr]=SS
dailyarea[id,ithr]=area

endfor

print,format='(i8,4x," members=",3x,14a7)',vdate,model[members]
print,format='(12x," BSS=",14f7.3)',dailySS[id,*]
print,format='(12x,"area=",14f7.3)',dailyarea[id,*]

id=id+1
ndone_random=ndone_random+1
if random gt 0 and ndone_random lt 10 then goto, choose_again

continue:

endfor              ;End of loop on days
endfor              ;End of loop on years


;Get pooled results
print,format='("probcat=",16f8.4)',probcat
if n_elements(ROCdiag) eq 0 then ROCdiag=0
areapool=fltarr(nthresh)
BSSdecomp=fltarr(nthresh)
BSt=fltarr(nthresh)
for ithr=0,nthresh-1 do begin
  print,' '
  print,'Threshold=',threshold[ithr]
  print,format='("      n=",16i8  )',n[*,ithr]
  print,format='("      o=",16f8.4)',o[*,ithr]/n[*,ithr]
  title='Rain !Mb!N '+strtrim(string(threshold[ithr]),2)+' mm d!E-1!N'
  probcatt=probcat
  nt=n[*,ithr]
  ot=o[*,ithr]
  ntref=nref[*,ithr]
  otref=oref[*,ithr]
  reliability_diagram,probcatt,nt,ot,BS,title=title,decomp=decomp
  reliability_diagram,probcatt,ntref,otref,BSref,title=title
  ;BSSdecomp[ithr]=(BSref-BS)/BSref  ;Gives same result as direct calculation
  BSt[ithr]=BS
;  read,'Enter any number to continue, 0 to quit... ',number
;  if number eq 0 then goto, done
  if (ROCdiag) then begin
    probcatt=probcat
    nt=n[*,ithr]
    ot=o[*,ithr]
    ROC_diagram,probcatt,nt,ot,area,title=title
    areapool[ithr]=area
;    read,'Enter any number to continue, 0 to quit... ',number
;    if number eq 0 then goto, done
    endif
print,'obar, rel, res, unc =',decomp
endfor

BSSpool=(sumPSref-sumPS)/sumPSref

;Compute daily mean and S.D. values of Brier skill score and ROC area
meanSS  =make_array(nthresh,/float,value=missng)
meanarea=make_array(nthresh,/float,value=missng)
SDSS    =make_array(nthresh,/float,value=missng)
SDarea  =make_array(nthresh,/float,value=missng)
for ithr=0,nthresh-1 do begin
  SSd=dailySS[*,ithr]
  ss=where(SSd ne missng,nok)
  if nok ge 1 then begin
    result=moment(SSd[ss],sdev=sdev)
    meanSS[ithr]=result[0]
    SDSS[ithr]=sdev
    endif
  aread=dailyarea[*,ithr]
  ss=where(aread ne missng,nok)
  if nok ge 1 then begin
    result=moment(aread[ss],sdev=sdev)
    meanarea[ithr]=result[0]
    SDarea[ithr]=sdev
    endif
  endfor
print,' '
print,format='("Threshold      ",16i8  )',threshold
print,format='("Pooled BS      ",16f8.4)',BSt
print,format='("Pooled BSS     ",16f8.4)',BSSpool
;print,format='("BSS from decomp",16f8.4)',BSSdecomp
print,format='("Pooled ROC area",16f8.4)',areapool
print,format='("Daily mean BSS ",16f8.4)',meanSS
print,format='("Daily mean area",16f8.4)',meanarea
print,format='("S.D. daily BSS ",16f8.4)',SDSS
print,format='("S.D. daily area",16f8.4)',SDarea

if random gt 0 then print,'# times each model used:',nuse

;Write statistics to output file for later plotting
outputfile='POORMAN_'+forecast+'_'+season+'_'+region+'.pstats'
openw,unit,outputfile,/get_lun
printf,unit,nthresh,'  Number of thresholds'
printf,unit,ncat,'  Number of probability categories'
printf,unit,format='("  prob=",11f8.4)',probcat
for ithr=0,nthresh-1 do begin
  printf,unit,format='(f4.1," n=",11i8  )',threshold[ithr],n[*,ithr]
  for i=0,ncat-1 do if n[i,ithr] gt 0 then ratio[i,ithr]=o[i,ithr]/n[i,ithr]
  printf,unit,format='(f4.1," o=",11f8.4)',threshold[ithr],ratio[*,ithr]
endfor
printf,unit,' '
printf,unit,format='("Threshold      ",6i8  )',threshold
printf,unit,format='("Pooled BS      ",6f8.4)',BSt
printf,unit,format='("Pooled BSS     ",6f8.4)',BSSpool
printf,unit,format='("Pooled ROC area",6f8.4)',areapool
printf,unit,format='("Daily mean BSS ",6f8.4)',meanSS
printf,unit,format='("Daily mean area",6f8.4)',meanarea
printf,unit,format='("S.D. daily BSS ",6f8.4)',SDSS
printf,unit,format='("S.D. daily area",6f8.4)',SDarea
close,unit & free_lun,unit
print,'Wrote statistics to '+outputfile


done:

print,string(7b)
end
pro reliability_diagram, probcat,n,o,BS,title=title,decomp=decomp

;Make reliability diagram for probabilistic forecast of
;value at each gridpoint exceeding threshold

;Input:   probcat      Probability categories
;         n            # forecasts in each fcst probability category
;         o            # occurrences in each fcst probability category
;         title        optional title for plot
;Output:  BS           Brier Score computed as rel-res+unc (Wilks, 1995)
;         decomp       [obar, reliability,resolution, uncertainty]

sz=size(probcat)
ncat=sz[1]
nmem=ncat-1

;Compute observed probability in each fcst prob. category
for k=0,ncat-1 do o[k]=o[k]/n[k]        

;Get rid of empty bins
ss=where(o ge 0. and o le 1.,n1)
n=n[ss]
o=o[ss]
probcat=probcat[ss]

nn=total(n)

if n_elements(title) eq 0 then title='Reliability Diagram'
plot,[0,1],[0,1],title=title, $
  xtitle='Forecast Probability',ytitle='Observed Probability', $
  /xstyle,/ystyle,linestyle=2
oplot,probcat,o,max_value=1.,thick=2

;Compute reliability, resolution, and uncertainty
obar=total(n*o)/nn
rel =total(n*(probcat-o)^2.)/nn
res =total(n*(o-obar)^2.)/nn
unc =obar*(1.-obar)
BS=rel-res+unc

;print,'obar,rel,res,unc=',obar,rel,res,unc
;print,'BS=',BS
decomp=[obar,rel,res,unc]

end
pro ROC_diagram, probcat,n,o,area,title=title

;Make relative operating characteristic (ROC) diagram for probabilistic 
;forecast of value at each gridpoint exceeding threshold

;Input:   probcat      Probability categories
;         n            # forecasts in each fcst probability category
;         o            # occurrences in each fcst probability category
;         title        optional title for plot
;Output:  area         area under ROC curve

sz=size(probcat)
ncat=sz[1]
nmem=ncat-1

;Compute hit rates and false alarm rates by probability category
;(Follow Buizza et al. 1998 QJRMS)
non=n-o
h=fltarr(ncat)
f=fltarr(ncat)
for k=0,ncat-2 do begin
  W=total(non[0:k])
  Y=total(  o[0:k])
  Z=total(non[k+1:ncat-1])
  X=total(  o[k+1:ncat-1])
  h[k]=X/(X+Y)
  f[k]=Z/(Z+W)
  endfor
h[ncat-1]=0
f[ncat-1]=0

;Get rid of empty bins
for k=0,ncat-1 do o[k]=o[k]/n[k]        
ss=where(o ge 0. and o le 1.,nok)
n=n[ss]
o=o[ss]
h=h[ss]
f=f[ss]
probcat=probcat[ss]

h=[1.,h]
f=[1.,f]

print,format='("      f=",16f8.4)',f
print,format='("      h=",16f8.4)',h

if n_elements(title) eq 0 then title='Relative Operating Characteristic'
plot,[0,1],[0,1],title=title, $
  xtitle='False Alarm Rate',ytitle='Hit Rate', $
  /xstyle,/ystyle,linestyle=2
oplot,f,h,max_value=1.,thick=2
oplot,f,h,max_value=1.,psym=7

;Compute area under the curve (use trapezoid rule)
if n_elements(f) lt 3 then begin
  area=0.5
  goto, done
  endif
;dx=abs(f[1:ncat-1]-f[0:ncat-2])
dx=abs(f[1:nok-1]-f[0:nok-2])
area=0.
;for k=0,ncat-2 do area=area+dx[k]*(h[k]+h[k+1])/2.
for k=0,nok-2 do area=area+dx[k]*(h[k]+h[k+1])/2.

;print,'area=',area

done:
end
