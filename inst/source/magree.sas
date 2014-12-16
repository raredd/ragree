/*
* macro on web site http://support.sas.com/ctx/samples;
*/
%macro magree(version,
   data=_last_,
   items=,
   raters=,
   response=,
   stat=BOTH
   );

%if &version ne %then %put MAGREE macro Version 1.1;
options nonotes;

%let error=0;

/* Verify DATA= is specified and the data set exists */
%if &data ne %then %do;
  %if %sysfunc(exist(&data)) ne 1 %then %do;
    %put ERROR: DATA= data set, %upcase(&data), not found.;
    %goto exit;
  %end;
%end;
%else %do;
  %put ERROR: The DATA= parameter is required.;
  %goto exit;
%end;

/* Verify required parameters were specified */
%if &items= %then %do;
  %put ERROR: The ITEMS= parameter is required.;
  %goto exit;
%end;
%if &raters= %then %do;
  %put ERROR: The RATERS= parameter is required.;
  %goto exit;
%end;
%if &response= %then %do;
  %put ERROR: The RESPONSE= parameter is required.;
  %goto exit;
%end;

%if %upcase(&stat) ne KAPPA and 
    %upcase(&stat) ne KENDALL and
    %upcase(&stat) ne BOTH %then %do;
  %put ERROR: STAT= must be set to KAPPA, KENDALL, or BOTH.;
  %goto exit;
%end;

/* Check for existence of variable names */
%let dsid=%sysfunc(open(&data));
%if &dsid %then %do;
  %let chkitem=%sysfunc(varnum(&dsid,%upcase(&items)));
  %let chkrate=%sysfunc(varnum(&dsid,%upcase(&raters)));
  %let chkresp=%sysfunc(varnum(&dsid,%upcase(&response)));
  %let rc=%sysfunc(close(&dsid));
%end;
%else %do;
  %put ERROR: Could not open DATA= data set.;
  %goto exit;
%end;
%if &chkitem=0 %then %do;
    %put ERROR: Variable %upcase(&items) not found.;
    %goto exit;
%end;
%if &chkrate=0 %then %do;
    %put ERROR: Variable %upcase(&raters) not found.;
    %goto exit;
%end;
%if &chkresp=0 %then %do;
    %put ERROR: Variable %upcase(&response) not found.;
    %goto exit;
%end;

title2 "MAGREE macro";

%if &data=_last_ %then %let data=&syslast;

/* Get type of response variable so can format later.
 ======================================================================*/
%let dsid=%sysfunc(open(&data));
%if &dsid %then %do;
  %let ynum=%sysfunc(varnum(&dsid,&response));
  %let ytype=%sysfunc(vartype(&dsid,&ynum));
  %let rc=%sysfunc(close(&dsid));
%end;
%else %do;
  %put WARNING: Could not check type of &response variable.;
  %put %str(         ) Continuing assuming it is numeric.;
%end;

/* Verify that all subjects have the same number of ratings.  Note
 * that the raters variable must use the same values to identify the
 * raters in all subjects even if different raters rate the subjects.
 * But the number of raters must be the same across all subjects.
 * Create response by item summary table needed to compute kappas and 
 * get numbers of subjects, raters, and response categories.
 ======================================================================*/
data _nomiss; set &data; if &response ne " "; run;

proc freq data=_nomiss noprint;
  table &items*&raters   / sparse out=_balance;
  table &response*&items / out=_ycnts(drop=percent);
  table &items           / out=_n;
  table &raters          / out=_m;
  table &response        / out=_k(drop=count percent);
  run;

data _null_;
  set _balance nobs=_nobs;
  if count ne 1 then call symput('error',1);
  run;
%if &error %then %do;
  %put ERROR: Each rater must rate each subject exactly once.;
  %goto exit;
%end;

data _null_;
  set _n nobs=n;
  set _m nobs=m;
  set _k nobs=k;
  call symput('m',left(m));
  call symput('n',left(n));
  call symput('k',left(k));
  run;


/************************  Compute kappa  ******************************/
%if %upcase(&stat)=KAPPA or %upcase(&stat)=BOTH %then %do;

/* Create coded values (1,2,3,...) for response categories to use as
 * indices. 
 ======================================================================*/
data _ycnts;
  set _ycnts; 
  by &response;
  if first.&response then _code+1;
  run;

/* Compute kappa statistics for each category and overall
 ======================================================================*/
data _kappas;
  set _ycnts end=eof;
  array kapnum {&k} knum1-knum&k;
  array catsum {&k} sum1-sum&k;
  array kapj {&k} kappa1-kappa&k;  /* V6 limit: 999 categories */
  array pb {&k} pbar1-pbar&k;
  array zkapj {&k} zkj1-zkj&k;
  array prkapj {&k} prkj1-prkj&k;
  kapnum{_code} + count*(&m - count);
  catsum{_code} + count;
  if eof then do;
    knum=0; kden=0; pqqp=0;
    nmm = &n*&m*(&m-1);
    sekapj=sqrt(2/nmm);
    do j=1 to &k;
     put "catsum[j]=";
     put catsum{j};
      pb{j} = catsum{j}/(&m*&n);
     put "pbj[j]=";
     put pb{j};
      pq = pb{j}*(1-pb{j});
     put "pq=";
     put pq;
      kapj{j} = 1 - (kapnum{j}/(nmm*pq));
      zkapj{j} = kapj{j}/sekapj;
      prkapj{j} = 1-probnorm(zkapj{j});
      knum = knum + pq*kapj{j};
      kden = kden + pq; 
      pqqp = pqqp + pq*((1-pb{j})-pb{j});
     put "pqqp=";
     put pqqp;
    end;
    kappa=knum/kden;
     put "kden=";
     put kden;
     put "n*k*(k-1)=";
     put nmm;
     put "pqqp=";
     put pqqp; 
    sekap = (sqrt(2)/(kden*sqrt(nmm)))*sqrt(kden**2-pqqp);
    zkap=kappa/sekap; prkap=1-probnorm(zkap);
    keep &response kapp stderr z prob;
    do i=1 to &k;
       set _k;
       kapp=kapj{i}; stderr=sekapj; z=zkapj{i}; prob=prkapj{i}; output;
    end;
    &response=.; kapp=kappa; stderr=sekap; z=zkap; prob=prkap;
    output;
  end;
  run;

/* Print kappa statistics
 ======================================================================*/
%let fexist=0;
data _null_;
  if ("&ytype"="N" and cexist("work.formats._yfmt.format") ne 1) or
     ("&ytype"="C" and cexist("work.formats._yfmt.formatc") ne 1) then
     call symput("fexist","1");
  run;
%if &fexist %then %do;
  proc format;
     value  %if &ytype=C %then $_yfmt;
            %else               _yfmt;  .="Overall";
  run;
%end;

proc print noobs label;
  format prob pvalue. 
         &response %if &ytype=C %then $_yfmt.; %else _yfmt.; ;
  label kapp="Kappa" stderr="Standard Error" prob="Prob>Z";
  title3 "Kappa statistics for nominal response";
  run;
%end;


/*********** Compute Kendall's Coefficient of Concordance, W ***********/
%if %upcase(&stat)=KENDALL or %upcase(&stat)=BOTH %then %do;

%if &ytype=C %then %do;
  %put %str(ERROR: Kendall%'s Coefficient of Concordance requires a);
  %put %str(      ) numeric, ordinal response.;
  %goto exit;
%end;

/* Rank the data using average ties.  If data are already ranked, this
 * won't change anything.
 ======================================================================*/
proc sort data=&data out=_sortr;
  by &raters;
  run;
proc rank data=_sortr out=_ranked;
  by &raters; 
  var &response;
  run;

/* R-square from one-way ANOVA is Kendall's W.
 ======================================================================*/
proc anova data=_ranked outstat=_anova noprint; 
  class &items; 
  model &response = &items; 
  run;

/* Compute F statistic and p-value for testing W=0.
   Ties are handled correctly.
 ======================================================================*/
data _w; 
  set _anova;
  retain SSsubj;
  if _n_=2 then do;
     w=ss/(ss+SSsubj);
     f=(&m-1)*w/(1-w);
     numdf=&n-1-2/&m;
     dendf=(&m-1)*numdf;
     prob=1-probf(f,numdf,dendf);
     keep w f numdf dendf prob; 
     output;
  end;
  SSsubj=ss;
  run;    

/* Print Kendall's coefficient of concordance and test
 ======================================================================*/
proc print noobs label;
  var w f numdf dendf prob;
  format prob pvalue.;
  label w="Coeff of Concordance" numdf="Num DF" dendf="Denom DF"
        prob="Prob>F";
  title3 "Kendall's Coefficient of Concordance for ordinal response";
  run;
%end;

%exit:
options notes;
title;
%mend magree;


