% Power estimate for: 2x2, Multiplicative
% ---------------------------------------------------------------------------


clc; clear;

% Input
n = [ 12 :2 :36 ]   ;                                 %   <----------
          df = n - 2;

gmr = 0.90 ;                                              %   <----------
cvw = 0.20 ;                                               %   <----------
          s = sqrt( log(1 + cvw^2) ) ;            % s = sqrt(MSE) from ANOVA

          
a = 0.05;           % significance level
hilim = 1.25;           % BE limit
lolim = 0.80 ; 


% -------------------------------------------------------------------------

% check whether GMR >= 1
        if gmr < 1
              d = lolim ;
        else
              d = hilim ;
        end


 ta = tinv( 1-a, df ); 
 
% n =  2 *   s^2  * (za1s + zb)^2     / (   log(gmr) - log(d)   )^2  =>

          tb =   - ta + sqrt(  n .* ( log(gmr) - log(d) )^2  ./  ( 2* s^2 )  ) ;
          tbc = - ta + sqrt(  (n - ta/2) .* ( log(gmr) - log(d) )^2  ./  ( 2* s^2 )  )  ;  % with the 'correction term' for better accuracy
                    % => Power
                    power   = 100* tcdf( tb, df ) 
                    powerc = 100* tcdf( tbc, df ) 



%----------------------------------------------------------------------------------------
% FIGURES

                    xmin =  min ( n )      - 0.01 ;
                    xmax = max  ( n )  + 0.01 ;
                    ymin =  0  -  1  ;
                    ymax = 100 + 1   ;    
                    
                    % Splines
                    xstep = xmin  :0.001  :xmax  ;  
                   
% Power (%) vs. N
figure(1);
hold on
grid on
plot1  = scatter(  n, powerc ,  110,  'o',  'filled',  'b'   );  sp1 = spline(  n, powerc,  xstep   );
plot2  = plot(    xstep, sp1, 'b' ,  'LineWidth',2 ) ;
xlabel( 'Sample size' , 'FontSize',16 );
ylabel( 'Power (%)' , 'FontSize',16 );
axis(  [ xmin xmax ymin ymax ]  ) ;
title('Power with the Correction term');
hold off

%{
figure(2);
hold on
grid on
plot3  = scatter(  n, power ,  120,  '^',  'r'   );  sp1 = spline(  n, power,  xstep   );
plot4  = plot(    xstep, sp1, 'r' ,  'LineWidth',2 ) ;
xlabel( 'Sample size' , 'FontSize',16 );
ylabel( 'Power (%)' , 'FontSize',16 );
axis(  [ xmin xmax ymin ymax ]  ) ;
title('Power (no Correction)');
hold off
%}
