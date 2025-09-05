pro prob_scores,show=show,h48=h48

;Read from Pertti's POP_3 file and compute a variety of probabilistic
;verification scores and diagrams:
;  Brier score
;  Brier skill score
;  Ranked probability score
;  Ranked probability skill score
;  Reliability diagram
;  ROC diagram
;  ROC area
;  Value score
;This is to put on the web page as an example

;To show a plot, set show='reliability', 'ROC', or 'value'
;To compute results for 48h forecasts set /h48 (default=24h fcsts)

;Notes from Pertti:
;------------------
;OBS is daily rainfall in mm
;P's are probabilities adding up to 1.0
;h0: No rain (precip le 0.2 mm)
;h1: precip between 0.3 and 4.4 mm
;h2: precip ge 4.5 mm

yyyy=0L & mm=0L & dd=0L
RR=0.
P1d=fltarr(3)
P2d=fltarr(3)
P24=fltarr(365,3)
P48=fltarr(365,3)
O24=fltarr(365,3)
O48=fltarr(365,3)
R24=fltarr(365)
R48=fltarr(365)
line=''

;Read the data from file
;Throw out bad data
n24=0
n48=0
openr,unit,'POP_3cat_2003.txt',/get_lun
readf,unit,line
for id=0,364 do begin
  readf,unit,yyyy,mm,dd,RR,P1d,P2d
  if RR lt 999 and P1d[0] ge 0. then begin
    P24[n24,*]=P1d
    O24[n24,0]=(RR le 0.2)
    O24[n24,1]=(RR ge 0.3 and RR le 4.4)
    O24[n24,2]=(RR ge 4.5)
    R24[n24]=RR
    n24=n24+1
    endif
  if RR lt 999 and P2d[0] ge 0. then begin
    P48[n48,*]=P2d
    O48[n48,0]=(RR le 0.2)
    O48[n48,1]=(RR ge 0.3 and RR le 4.4)
    O48[n48,2]=(RR ge 4.5)
    R48[n48]=RR
    n48=n48+1
    endif
  endfor
close,unit & free_lun,unit
print,'Number of good samples n24,n48=',n24,n48
n=n24
P24=P24[0:n-1,*]
P48=P48[0:n-1,*]
O24=O24[0:n-1,*]
O48=O48[0:n-1,*]
R24=R24[0:n-1,*]
R48=R48[0:n-1,*]

;Plot rain amount to see whether there is a seasonal signal
;If there is a strong seasonal cycle then the sample climatology
;for the whole year might not be an appropriate reference
plot,R24,ytitle='Rainfall (mm)',xtitle='Day',xrange=[0,n],/xstyle

;Get sample climatological probability
P24clim=fltarr(3)
P48clim=fltarr(3)
for icat=0,2 do begin
  P24clim[icat]=mean(O24[*,icat])
  P48clim[icat]=mean(O48[*,icat])
  endfor
print,'P24clim=',P24clim
print,'P48clim=',P48clim
P24clim=[[replicate(P24clim[0],n)],[replicate(P24clim[1],n)],[replicate(P24clim[2],n)]]
P48clim=[[replicate(P48clim[0],n)],[replicate(P48clim[1],n)],[replicate(P48clim[2],n)]]

;-------------------------------------------------------------------

if n_elements(show) eq 0 then show=''

;Compute statistics using pstats.pro
;Start with probability of rain (cats 1 + 2)
obs=O24[*,1]+O24[*,2]
pop=P24[*,1]+P24[*,2]
refpop=P24clim[*,1]+P24clim[*,2]
obshi=O24[*,2]
pophi=P24[*,2]
refpophi=P24clim[*,2]
RR=R24
RaPS=rps(P24,O24)
RPSref=rps(P24clim,O24)

if n_elements(h48) gt 0 then begin
 print,'Verifying 48h forecasts'
 obs=O48[*,1]+O48[*,2]
 pop=P48[*,1]+P48[*,2]
 refpop=P48clim[*,1]+P48clim[*,2]
 obshi=O48[*,2]
 pophi=P48[*,2]
 refpophi=P48clim[*,2]
 RR=R48
 RaPS=rps(P48,O48)
 RPSref=rps(P48clim,O48)
 endif

ncat=11
print,'' & print,'Statistics for POP (2nd and 3rd categories)'
pstats,obs,pop,refpop,ncat,show=show
print,'' & print,'Statistics for POP for highest rain category'
;If plotting results then make another window the same size as the first
if show ne '' then window,1,xsize=!d.x_size,ysize=!d.y_size
pstats,obshi,pophi,refpophi,ncat,show=show

RPSS=(RPSref-RaPS)/RPSref
print,'RPS=',RaPS
print,'RPSS=',RPSS

end