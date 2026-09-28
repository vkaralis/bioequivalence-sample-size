% Estimate sample size using APPROXIMATE formulas (with 'z')

% The results from this analysis are very close (in accordance with) to those obtained from nQuery
%---------------------------------------------------------------------------------------------------------------------------------

% An 'ITERATION METHOD' can also be constructed                                 <-------------- !!!!
% see: CS_BE.doc pg. 46 (Mathematical formulas for Sample Size)
%---------------------------------------------------------------------------------------------------------------------------------


% !!! CORRECT FOR 'GMR=1' see ... Iterative                 <------------- !!!!


clc; clear;


a2s     = 0.10;               % 2-sided a => 1-sided: a/2
       a1s = a2s / 2 ;

power = 0.80;
       b = 1- power;



%--------------
za1s = norminv( 1- a1s  ); disp('za1s= '); disp(za1s);  % for a1/2
zb = norminv( 1- b  ); disp('zb= '); disp(zb);

%{
% Wang formulas
sd = 2.5;
diff = 1.0;
% N estimation
n = (  (za1s + zb)^2 * sd^2  ) / (2*diff^2);
disp(n)
%}

disp( 'Approximate solutions using z' ); disp('----------------------------------------'); disp(' ');




% --------------------------------------------------------------------------
% A1. Parallel / un-transformed (= Additive model) e.g., Clinical study
r = 1;                % allocatio ratio
mt = 5;             %mean value of T
mr = 4;             %mean value of R
d = 1.25;           % equivalence limit
s =0.1;              % residual variability (in Normal scale)

% ... CORRECT FOR mT <, > mR            <-----------------!!!

n = (r+1)/r * s^2 *    (za1s + zb)^2  / (   (  mt - mr) - d   )^2 ;
          nr = round(n);
          if nr > n
                    ntot = 2 * nr;
          else
                    ntot = 2 * ( round(n)+1 );
          end
disp('A1');   
disp('Parallel / un-transformed');
disp('N total='); disp(ntot)

          % OR more accurately: Using the a 'correction term' (from Julious 2004)
          n = (r+1)/r * s^2 *    (za1s + zb)^2  / (   (  mt - mr) - d   )^2     + za1s / 2 ;
                    nr = round(n);
                     if nr > n
                               ntot = 2 * nr;
                     else
                                ntot = 2 * ( round(n)+1 );
                     end 
          disp('N total (corr)='); disp(ntot)
          
          % check for 'odd' numbers
          if mod(ntot,2) == 1
                    ntot = ntot + 1;
                    disp('N total (mult. of 2)='); disp(ntot)
          end

          
          
          


%---
% A2. Parallel / LN--transformed (= multiplicative model) e.g., BE study
r = 1;                % allocatio ratio
gmr=1.05;             % (expected) point GMR of the study
d = 1.25;           % bioequivalence limit
s =0.25;              % residual variability (in LN scale)

% ... CORRECT FOR: GMR < 1              <-----------------!!!
n =  (r+1)/r * s^2 *  (za1s + zb)^2    / (  log(gmr) - log(d)  )^2 ;
          nr = round(n);
          if nr > n
                    ntot = 2 * nr;
          else
                    ntot = 2 * (round(n)+1);
          end
disp('A2');   
disp('Parallel / LN-transformed');
disp('N total='); disp(ntot)


          % OR more accurately: Using the a 'correction term' (from Julious 2004)
          n =  (r+1)/r * s^2 *  (za1s + zb)^2    / (  log(gmr) - log(d)  )^2    + za1s / 2 ;
          nr = round(n);
                    if nr > n
                              ntot = 2 * nr;
                    else
                              ntot = 2 * (round(n)+1);
                    end
          disp('N total (corr.)='); disp(ntot)

          % check for 'odd' numbers
          if mod(ntot,2) == 1
                    ntot = ntot + 1;
                    disp('N total (mult. of 2)='); disp(ntot)
          end

          
          
          



% ------------------------------------------------------------------------------------------


% B1. Cross-over / un-transformed (= Additive model)
mt = 5;             %mean value of T
mr = 4;             %mean value of R
s =0.5 ;              % residual variability (in Normal-scale)
d0 = 1.2 ;           % equivalence limit

          if mt > mr + d0
                    disp('Impossible. Provide other mT, mR, limits values');
          elseif mt > mr
                    d = d0 ;
          elseif mt < mr - d0
                    disp('Impossible. Provide other mT, mR, limits values');
          else
                    d = - d0  ;
          end

n =  2*  s^2 *  (za1s + zb)^2   /   (   (mt - mr) - d  )^2 ;
          nr = round(n);
          if nr > n
                    ntot = nr;
          else
                    ntot = (round(n)+1);
          end
disp('B1');   
disp('Cross-over / un-transformed');
disp('N total='); disp(ntot)

          % OR more accurately: Using the a 'correction term' (from Julious 2004)
          n =  2*  s^2 *  (za1s + zb)^2   /   (   (mt - mr) - d  )^2       + za1s / 2 ;
          nr = round(n);
          if nr > n
                    ntot = nr;
          else
                    ntot = (round(n)+1);
          end
          disp('N total (corr.)='); disp(ntot)
          
          % check for 'odd' numbers
          if mod(ntot,2) == 1
                    ntot = ntot + 1;
                    disp('N total (mult. of 2)='); disp(ntot)
          end
          
          
          
          

%---
% B2. Cross-over / LN--transformed (= multiplicative model)           <------------- !!
belup = 1.25 ;          % bioequivalence limit
bello  = 0.80 ;          % bioequivalence limit

s =0.5 ;                   % residual variability (in LN-scale)
gmr = 1.05 ;              % (expected) point GMR of the study

          if gmr >= belup
                    disp('Impossible. Provide different Input');
          elseif gmr <= bello
                    disp('Impossible. Provide different Input');
          end

          if gmr < 1
                    d = bello ;
          else
                    d = belup ;
          end

n =  2 *   s^2  * (za1s + zb)^2     / (   log(gmr) - log(d)   )^2  ;
          nr = round(n);
          if nr > n
                    ntot = nr;
          else
                    ntot = (round(n)+1);
          end
disp('B2');   
disp('Cross-over / LN-transformed (with Sw)');
disp('N total ='); disp(ntot)

          % OR more accurately: Using the a 'correction term' (from Julious 2004)
          n =  2 *   s^2  * (za1s + zb)^2     / (   log(gmr) - log(d)   )^2    + za1s / 2 ;
          nr = round(n);
          if nr > n
                    ntot = nr;
          else
                    ntot = (round(n)+1);
          end
          disp('N total (corr.)='); disp(ntot)

          
          % check for 'odd' numbers
          if mod(ntot,2) == 1
                    ntot = ntot + 1;
                    disp('N total (mult. of 2)='); disp(ntot)
          end
          
          
          
          
          

%---
% B3. Cross-over / LN--transformed (= multiplicative model) with CVw            <---------- !!
belup = 1.25 ;          % bioequivalence limit
bello  = 0.80 ;          % bioequivalence limit

cv =0.32;                  % CVw
          s = sqrt( log(1 + cv^2) ) ;
gmr = 1.06 ;              % (expected) point GMR of the study

          if gmr >= belup
                    disp('Impossible. Provide different Input');
          elseif gmr <= bello
                    disp('Impossible. Provide different Input');
          end

          if gmr < 1
                    d = bello ;
          else
                    d = belup ;
          end

% n =  2 * cv^2  * (za1s + zb)^2     / (   log(gmr) - log(d)   )^2 ;   % <--- !!! REQUIRES a CORRECTION CV=sqrt(exp(mse)-1)
n =  2 * s^2  * (za1s + zb)^2     / (   log(gmr) - log(d)   )^2    ;
          nr = round(n);
          if nr > n
                    ntot = nr;
          else
                    ntot = (round(n)+1);
          end
disp('B3');          
disp('Cross-over / LN-transformed (with CVw)');
disp('N total='); disp(ntot)


          % OR more accurately: Using the a 'correction term' (from Julious 2004)
          n =  2 * s^2  * (za1s + zb)^2     / (   log(gmr) - log(d)   )^2      + za1s / 2;
          nr = round(n);
          if nr > n
                    ntot = nr;
          else
                    ntot = (round(n)+1);
          end
          disp('N total (corr.)='); disp(ntot)
          
          
          % check for 'odd' numbers
          if mod(ntot,2) == 1
                    ntot = ntot + 1;
                    disp('N total (mult. of 2)='); disp(ntot)
          end
          
          


