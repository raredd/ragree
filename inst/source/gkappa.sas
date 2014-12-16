*--------------------------------------------------------------------------
Macro for genenized kappa statistics
Missing score is allowed
Using longitudinal data format
Wanling Xie
September 2009
---------------------------------------------------------------------------;

%macro gkappa(data=, bygrp=, pnumber=, score=, outdata=);
proc sort data=&data; by &bygrp; 
proc freq data=&data noprint;
    tables &pnumber/out=patval;
    tables &score/out=outval2;
    tables &pnumber*&score/norow nocol nopercent out=outval;
    by &bygrp;
    title1 'Example from page 278-279 of Woolson and Clarke';

data outval; 
  set outval; 
  njk=COUNT; 
  njk2=njk*njk;

* merge number of scores for &pnumber k with the &pnumber*&score data; 
data patval; set patval; nk=COUNT; keep &bygrp &pnumber nk;
proc sort data=outval; by &bygrp &pnumber; 
proc sort data=patval; by &bygrp &pnumber; 
data outval;  
   merge outval patval; by &bygrp &pnumber;
   tempnk=njk*(nk-1); 

proc sort data=outval; by &bygrp &score; 
proc means data=outval noprint; 
var njk njk2 tempnk;
by &bygrp &score;  
output out=final sum=sumnjk sumnjk2 poss; 

data outval2; 
set outval2; 
pj=(PERCENT/100);
pj2=(PERCENT/100)**2; 

proc sort data=outval2; by &bygrp &score; 
data final; 
  merge final outval2; 
  by &bygrp &score; 
  actualj=sumnjk2-sumnjk; 
  pej=pj2; 
  agreej=actualj/poss; 
  kappaj=(agreej-pj)/(1-pj); 

/*
proc print data=final; var &score sumnjk sumnjk2 pj pj2 pej poss actualj agreej kappaj; 
by &bygrp;
title1 'Kappa statistic for each category'; 
run;
*/
proc means data=final noprint; 
var actualj pj2 poss;
by &bygrp; 
output out=meanval sum=actual pe poss;
data &outdata;
 set meanval;
 agree=actual/poss; 
 kappa=(agree-pe)/(1-pe); 
 DROP  _TYPE_    _FREQ_;
run;
%mend gkappa; 
