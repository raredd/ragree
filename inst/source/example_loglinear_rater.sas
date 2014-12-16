*example using loglinear models for interrater agreement;

options ls=130 noovp nocenter orientation=landscape formdlim=' ';

*example on how to run loglinear model analysis; 

*page 36 of VonEye book; 
data one; 
input x y count; 
if x=y then delta=1; 
if x ne y then delta=0; 
cards; 
1 1 4
1 2 3
1 3 1
1 4 0 
2 1 1
2 2 5
2 3 2
2 4 2
3 1 1
3 2 1
3 3 6
3 4 2
4 1 0
4 2 1
4 3 5
4 4 6
; 

run;
*page 65 of VonEye book;
data two; 
input r1 r2 r3  count; 
if r1=r2 then delta12=1; 
if r1 ne r2 then delta12=0;
if r1=r3 then delta13=1; 
if r1 ne r3 then delta13=0; 
if r2=r3 then delta23=1; 
if r2 ne r3 then delta23=0; 
if r1=r2=r3 then delta123=1;
else delta123=0; 
cards; 
1 1 1 4
1 1 2 3
1 1 3 6
1 2 1 2 
1 2 2 1
1 2 3 3
1 3 1 2
1 3 2 2
1 3 3 17
2 1 1 0
2 1 2 1
2 1 3 2
2 2 1 1 
2 2 2 1
2 2 3 1
2 3 1 0
2 3 2 0
2 3 3 4
3 1 1 0
3 1 2 1
3 1 3 3
3 2 1 0 
3 2 2 1
3 2 3 8
3 3 1 0
3 3 2 4
3 3 3 96
; 

run;
proc print data=two; 
  var r1 r2 r3 delta12 delta13 delta23 delta123;  
  title1;
  run;

data three; 
input a b count;
if a=b then delta=1;
else if a ne b then delta=0;
if a<b then trend=1; 
else if a>b then trend=-1;
else if a=b then trend=0; 
cards;
1 1 17
1 2 27
1 3 3
2 1 16
2 2 45
2 3 14
3 1 1
3 2 3
3 3 3
; 
run;
proc print data=three; 
  var a b count delta trend;   
  title1;
  run;
 


ods rtf file='C:\edie\Loglinear rater agreement\example_loglinear_rater.rtf'  style=SASDocPrinter;
ods noptitle; 
proc freq data=one; 
 table x*y/kappa; 
 weight count; 
 run;
proc genmod data=one; 
class x y; 
model count=x y/dist=poisson type3 type1;
*output out=residuals pred=pred stdreschi=stdreschi stdresdev=stdresdev reslik=reslik;
title1 'Base model main effects only'; 	
run;
proc genmod data=one; 
class x y; 
model count=x y delta/dist=poisson type3 type1;
output out=residuals pred=pred stdreschi=stdreschi stdresdev=stdresdev reslik=reslik;
title1 'Equal Weight Agreement Model';  	
run;
proc print data=residuals; 
run;

proc freq data=two; 
 table r1*r2*r3/kappa; 
 weight count; 
 run;
* genmod;
  * base model main effects only section 2.1 page 65; 
proc genmod data=two; 
class r1 r2 r3; 
model count=r1 r2 r3/dist=poisson type3 type1;
*output out=residuals pred=pred stdreschi=stdreschi stdresdev=stdresdev reslik=reslik;
title1 'Base model main effects only'; 	
run;
proc genmod data=two; 
class r1 r2 r3; 
model count=r1 r2 r3 delta12 delta13 delta23/dist=poisson type3 type1;
*output out=residuals pred=pred stdreschi=stdreschi stdresdev=stdresdev reslik=reslik;
title1 'Simultaneous agreement between pairs of raters. Section 2.4.1.1';  	
run;
proc genmod data=two; 
class r1 r2 r3; 
model count=r1 r2 r3 delta12 delta13 delta23 delta123/dist=poisson type3 type1;
*output out=residuals pred=pred stdreschi=stdreschi stdresdev=stdresdev reslik=reslik;
title1 'Agreement among all raters. Section 2.4.1.2';  	
run;
proc genmod data=two; 
class r1 r2 r3; 
model count=r1 r2 r3 delta123/dist=poisson type3 type1;
*output out=residuals pred=pred stdreschi=stdreschi stdresdev=stdresdev reslik=reslik;
title1 'Agreement among all raters-include only the delta for all raters. Section 2.4.1.2';  	
run;

proc freq data=three; 
 table a*b/kappa; 
 weight count;
 title;  
 run;
* genmod;
  * base model main effects only section 2.1 page 65; 

 *other models are in section 2.4.2 for rater specific trends;

proc genmod data=three; 
class a b; 
model count=a b/dist=poisson type3 type1;
*output out=residuals pred=pred stdreschi=stdreschi stdresdev=stdresdev reslik=reslik;
title1 'Base model main effects only'; 	
run;
proc genmod data=three; 
class a b; 
model count=a b delta trend/dist=poisson type3 type1;
*output out=residuals pred=pred stdreschi=stdreschi stdresdev=stdresdev reslik=reslik;
title1 'Equal Weight Agreement Model with trend Section 2.4.1.1';  	
run;
 
ods rtf close;
