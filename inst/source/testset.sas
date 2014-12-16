options nocenter formdlim=' ' ls=80;

data one;
 input patient rater category $;
 cards;
 1 1 mild    
 1 2 mod
 1 3 mod
 1 4 mod
 2 1 mild    
 2 2 mild
 2 3 mild
 2 4 mild
 3 1 mild    
 3 2 sev
 3 3 sev
 3 4 sev
 4 1 mod    
 4 2 mod
 4 3 mod
 4 4 sev
 5 1 sev    
 5 2 sev
 5 3 sev
 5 4 sev
 6 1 mod    
 6 2 mod
 6 3 sev
 6 4 sev
 7 1 mild    
 7 2 mild
 7 3 sev
 7 4 sev
 8 1 mild    
 8 2 mild
 8 3 mild
 8 4 sev
 9 1 mild    
 9 2 mild
 9 3 mild
 9 4 mod
 10 1 mod    
 10 2 sev
 10 3 sev
 10 4 sev
 11 1 mild    
 11 2 sev
 11 3 sev
 11 4 sev
 12 1 mild    
 12 2 mild
 12 3 mild
 12 4 mild
 13 1 mod    
 13 2 sev
 13 3 sev
 13 4 sev
 14 1 mod    
 14 2 mod
 14 3 mod
 14 4 mod
 15 1 mild    
 15 2 mild
 15 3 mild
 15 4 sev
 16 1 mild
 16 2 mod
 16 3 sev
 16 4 sev
 17 1 mod
 17 2 sev
 17 3 sev
 17 4 sev
 18 1 mild
 18 2 mod
 18 3 mod
 18 4 mod
 19 1 sev
 19 2 sev
 19 3 sev
 19 4 sev
 20 1 mild
 20 2 mod
 20 3 mod
 20 4 sev
 21 1 mild
 21 2 mild
 21 3 mild
 21 4 mild
 22 1 mild
 22 2 mod
 22 3 mod
 22 4 mod
 23 1 mod
 23 2 mod
 23 3 mod
 23 4 sev
 24 1 sev
 24 2 sev
 24 3 sev
 24 4 sev
 25 1 mild
 25 2 mild
 25 3 mod
 25 4 mod
 26 1 mod
 26 2 mod
 26 3 mod
 26 4 mod
 27 1 mild
 27 2 mild
 27 3 mild
 27 4 sev
;
run;


proc freq data=one;
    tables patient/out=patval;
    tables category/out=outval2;
    tables patient*category /norow nocol nopercent out=outval;
    tables patient*category;
    title1 'Example from page 278-279 of Woolson and Clarke';
run ;

*proc print data=outval; 
proc print data=outval2; run ;
proc print data=patval;run;
data outval; 
  set outval; 
  njk=COUNT; 
  njk2=njk*njk;
run ;
proc print data=outval;
  var patient category njk njk2; 
run ;
* merge number of scores for patient k with the patient*category data; 
data patval; set patval; nk=COUNT; keep patient nk; run ;
proc sort data=outval; by patient; run;
proc sort data=patval; by patient; run;
data outval;  
   merge outval patval; by patient;
   tempnk=njk*(nk-1); 
run;
proc print data=outval;
  var patient category njk njk2 nk tempnk; 
run;
proc sort data=outval; by category;  run ;
proc means data=outval noprint; 
var njk njk2 tempnk;
by category;  
output out=final sum=sumnjk sumnjk2 poss; 
run ;
proc print data=final; var category sumnjk sumnjk2 poss; run;

data outval2; 
set outval2; 
pj=(PERCENT/100);
pj2=(PERCENT/100)**2;
pj3=(PERCENT/100)**3;
run;
proc print data=outval2; run;
proc sort data=outval2; by category;  run;
data final; 
  merge final outval2; 
  by category; 
  actualj=sumnjk2-sumnjk; 
  pej=pj2; 
  agreej=actualj/poss; 
  kappaj=(agreej-pj)/(1-pj); 
run ;

proc print data=final; var category sumnjk sumnjk2 pj pj2 poss actualj agreej kappaj; run ;
proc means data=final noprint; 
var actualj pj2 pj3 poss; 
output out=meanval sum=actual pe p3 poss; 
run;
data meanval;
 set meanval; 
 agree=actual/poss;
if pe ne 1 then kappa=(agree-pe)/(1-pe); 
  nr=4; 
  c1=2*nr-3; c2=2*(nr-2);
  pe2=pe**2;
  * variance of generalized kappa page 279 woolson;
 if pe ne 1 then vk=( 2*(pe-c1*pe2+ c2*p3))/(poss*(1-pe)**2);
 se=sqrt(vk);
run;

proc print data=meanval;
title1 'Agreement';
var  actual poss agree pe kappa c1 c2 vk se;
run;


%inc "Z:\macros\magree.sas";
%magree(data=one,items=patient,raters=rater, response=category,stat=kappa);


data one ; 
	set one ;
	marker = 1 ;
run ;

%inc "Z:\macros\gkappa.sas";
%gkappa(data=one, bygrp=marker, pnumber=patient, score=category, outdata=tmp);



*proportion of pairs that agree; 
*output from magree generalized kappa=0.40634; 
data meanval; 
  set meanval;
  gkappa= 0.40634; 
  agreeval=gkappa*(1-sumj)+sumj; 
run;
proc print data=meanval; 
  var sumj gkappa agreeval;
run;
