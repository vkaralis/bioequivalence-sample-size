% Estimate sample size using ITERATIVE algorithm (initially: with 'z' ---> then with 't')
%---------------------------------------------------------------------------------------------------------------------------------

% !!! ~ 99% accordance with Diletti's tables ( the only small differences when: GMR=1,CV --> 0.0 )


%---------------------------------------------------------------------------------------------------------------------------------
% based on: Machin book pg.121          => see: CS_BE_Designs.doc pg 46

% Assuming: Cross over design (fully balanced: equal n1 & n2), ln-transform
% Input: a, b ,     CVw, GMR,           iterbase, itersize
%---------------------------------------------------------------------------------------------------------------------------------

% see: CS_BE.doc pg. 46 (Mathematical formulas for Sample Size)
%---------------------------------------------------------------------------------------------------------------------------------
% for 'BE3' method of: sample size estimator script (N_CrossOover_Parallel)
%---------------------------------------------------------------------------------------------------------------------------------


clc; clear;



% INPUT
% -----------

a2s     = 0.10;               % 2-sided a => 1-sided: a/2
       a1s = a2s / 2 ;

power = 0.90;
       b = 1- power;

belup = 1.25 ;          % upper BE limit
bello  = 0.80 ;          % lower BE limit

%---
gmr = 1.05 ;              % (expected) point GMR of the study              <--------------------

cvw = 0.20 ;                 % CVw                                                            <--------------------
          s = sqrt( log(1 + cvw^2) ) ;
%---          

maxiter = 100 ;     % maximum number of iterations
itersize = 0.5 ;      % difference in size limit for termination of the iterations


% Iterations based on:
          % 1. Approximate N (ntot1)=> 1
          % 2. Approximate N with correction term (ntot2) [from Julious 2004] => 2
                    iterbase = 2;                                                      %    <--------------------


                    
                    
                    

% CALCULATIONS
% -----------------------

za1s = norminv( 1- a1s  ); % disp('za1s= '); disp(za1s);  % for a1/2

if gmr ~= 1
          zb = norminv( 1- b  ); % disp('zb= '); disp(zb);
else
          zb = norminv( 1- b/2  ); 
end



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

% A. Approximate method (ntot1)

          n =  2 * s^2  * (za1s + zb)^2     / (   log(gmr) - log(d)   )^2    ;

          nr = round(n);
          if nr > n
                    ntot1 = nr;
          else
                    ntot1 = (round(n)+1);
          end
          %disp('N total='); disp(ntot1)


% B. Approximate method (ntot2): more accurate using a 'correction term' (from Julious 2004)

          n =  2 * s^2  * (za1s + zb)^2     / (   log(gmr) - log(d)   )^2      + za1s / 2;
          
          nr = round(n);
          if nr > n
                    ntot2 = nr;
          else
                    ntot2 = (round(n)+1);
          end
          %disp('N total (corr.)='); disp(ntot2)
          


          
          
% C. ITERATIONS (with 't' statistics)
%------------------------------------------

% set on which estimates (with OR without correction term) the iteration algorithm will be based)
          if iterbase == 1
          n0 = ntot1;
          else
          n0 = ntot2 ;
          end

          
% check if ntot2 <= 2              => No iterations
          if ntot2 <= 2
                    
                    disp ('!!!   N <= 2 and NO iterations were made'); disp(' '); disp(' ');
                    %nfinal = ntot2 ;
                    nfinal = ntot2 + 2 ;          % to be in accordance with Diletti's tables
                    count = 0 ;
                    
                    
                    
          else
          
                    status = 0;
                    count = 0;
                    for i=1 :maxiter

                              df = round( n0 ) - 2 ;
                              ta1s = tinv( 1-a1s, df ); % disp('t= '); disp( ta1s);
                              
                              if gmr ~= 1
                              tb = tinv( 1- b, df  ); % disp('zb= '); disp(zb);
                              else
                              tb = tinv( 1- b/2, df  ); 
                              end
                              %tb = tinv( 1- b, df  ); % disp('tb= '); disp(tb);
          
                              n(i) =  2 * s^2  * (ta1s + tb)^2     / (   log(gmr) - log(d)   )^2      + ta1s / 2;
                    
                                        if n(i) - n0 < itersize
                                                  nfinal = n(i); 
                                                  status = 1;
                                                  break
                                        else
                                                  n0 = n(i);
                                        end

                                        % check if the maximum # of iterations is exceeded
                                        if (  ( i == 100 )   &&   ( status == 0) )
                                                  disp('Maximum number of Iterations reached without convergence')
                                        end
                                         
                                        count = count + 1;
                                      
                    
                    end
          
                    
          end
          
          

% disp(n);          % make a plot !!!
niter = round( nfinal );
%disp(nfinal)



% ===========================================================


% OUTPUT
%--------------

disp( 'SAMPLE SIZE (total): '); disp('------------------------------'); disp(' ');

disp('A. Approximate method '); disp( ntot1 );
disp('B. Approximate method (with correction term) '); disp( ntot2 );
disp('C. Iterative method '); disp( niter );
          disp(' # of iterations= '); disp( count ); disp(' ');


% Round to an 'even' number (if odd number)
          % A
                    if mod( ntot1, 2 ) == 1
                              ntot1e = ntot1 + 1;
                    else
                              ntot1e = ntot1 ;
                    end
          % B
                    if mod( ntot2, 2 ) == 1
                              ntot2e = ntot2 + 1;
                    else
                              ntot2e = ntot2 ;
                    end    
          % C
                    if mod( niter, 2 ) == 1
                              nitere = niter + 1;
                    else
                              nitere = niter ;
                    end
                    
disp( 'SAMPLE SIZE (total / rounded): '); disp('-----------------------------------------'); disp(' ');

          disp('A. Approximate method '); disp( ntot1e );
          disp('B. Approximate method (with correction term) '); disp( ntot2e );
          disp('C. Iterative method '); disp( nitere );
         

          
