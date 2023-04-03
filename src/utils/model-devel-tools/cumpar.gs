'reinit'
"open gfsh-wstar.ctl"
"open gfsh-heateng.ctl"

ifiles=2
prefix=shallow
'set xlopts 1 7 0.13'
'set ylopts 1 7 0.12'
"set mproj off"
"set z 1 "
"d precip"
"d precip.2"
"d precip.4"
"d precip.3"
"d precip.5"
*"set ccolor 2"
"q pos"
"clear"
"d precip.3-precip.2"
return
"d precip.4"
"d precip.5"
return
************


"set y 1"
"set z 1 25"
"set axlim -20 25"
"d dtdt"
"d dtdt.2"
"d dtdt.3"
"q pos"
*"clear"
"d dqdt"
"d dqdt.2"
"d dqdt.3"
"draw title DTDT DQDT"
"q pos"
"clear"
"set axlim -10 10"
"d delh"
"d delh.2"
"d delh.3"
"draw title DELH"
*return
"clear"


*********** precipitation plot
is=1
loc=7.5
while(is<ifiles+1)
say is
"set dfile "is""
"set z 1"
"set y 1 161"
var=precip
*if(is=1); var=prliq;endif
'set grads off'
'set axlim 0 5' 
'set cthick 16'
'set cmark 0'
'set ccolor 'is+2''
'd 'var'.'is''
'define 'varm''is'=ave('var'.'is',y=1,y=161)'
'q define'
pre=sublin(result,is)
prem=subwrd(pre,2)
say prem
rc2=getinfo(is,loc,prem)
loc=loc-0.25 
is=is+1
endwhile
"draw title Precip (mm/h)"
'printim 'prefix'-precip.png white x3000 y2000'
'q pos'
'clear'

*********** vertical profiles 
"set t 1"
"set z 1 41"
"set y 1"
loc=7.5
nvar=11
is=1
while(is<ifiles+1)
say is
"set dfile "is""
ivar=1
while(ivar<nvar+1)
if(ivar=1); varx=dtdt;endif
if(ivar=2); varx=dqdt;endif
if(ivar=3); varx=delt;endif
if(ivar=4); varx=delq;endif
if(ivar=5); varx=mup;endif
if(ivar=6); varx=mdn;endif
if(ivar=7); varx=dqldt;endif
if(ivar=8); varx=upent;endif
if(ivar=9); varx=updet;endif
if(ivar=10); varx=dnent;endif
if(ivar=11); varx=dndet;endif

* average var.
'define 'varx'm'is'=ave('varx'.'is',y=1,y=161)'

ivar=ivar+1
endwhile
is=is+1
endwhile
*return


'set zlog on'
nvar=5
ivar=1
while(ivar<nvar+1)
if(ivar=1); varx=dTdt;unit="K/day";x1=-30;x2=30;xf=1;endif
if(ivar=2); varx=dQdt;unit="K/day";x1=-20;x2=20;xf=1;endif
if(ivar=3); varx=mup;unit="kg/s/m^2/100";x1=-0.5;x2=5;xf=100;endif
if(ivar=4); varx=mdn;unit="kg/s/m^2/100";x1=-2.5;x2=0.1;xf=100;endif
if(ivar=5); varx=dqldt;unit="K/day";x1=-1;x2=3;xf=1;endif
if(ivar=6); varx=delt;unit="m/s/day";x1=-30;x2=30;xf=1;endif
if(ivar=7); varx=delq;unit="m/s/day";x1=-30;x2=30;xf=1;endif

clear
loc=7.5
is=1
"set z 1 41"
while(is<ifiles+1)
say is
*"set dfile "is""
'set grads off'
'set axlim 'x1' 'x2'' 
'set cthick 10'
'set cmark 0'
'set ccolor 'is+2''
'd 'xf'*'varx'm'is''

rc2=getinfo2(is,loc)
loc=loc-0.25 

is=is+1
endwhile
"draw title "varx"  "unit""
"draw ylab pressure (hPa)"
file1=varx
'gxyat -y 2000 -x 3000 'prefix'-'file1'.png' 
*'printim 'prefix'-'file1'.png white x1500 y1000'
"q pos"
ivar=ivar+1
endwhile
"reset"



return
******************************************************
function getinfo(is,loc,prem)
"q file "is""
 name=sublin(result,1)
* name=subwrd(name,4) 
"set strsiz 0.16 0.17"
"set string "is+2" l 7 0"
"draw string 2.5 "loc" "name" mean="prem" mm/h"
return
******************************************************
function getinfo2(is,loc)
"q file "is""
 name=sublin(result,1)
* name=subwrd(name,4) 
"set strsiz 0.16 0.17"
"set string "is+2" l 7 0"
"draw string 3. "loc" "name""
return


*'gxyat -y 2000 -x 3000 fig.png'
