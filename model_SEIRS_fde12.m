function [t, y] = model_SEIRS_fde12(N,alpha)

%   This code is in an appendix of the following paper:
%
%   [*] Rosa, S.; Torres, D.F.M. Numerical Fractional Optimal Control of 
%   Respiratory Syncytial Virus Infection in Octave/MATLAB. Mathematics
%   2023, 11, 1511. 
%
%   DOI: https://doi.org/10.3390/math11061511

%model_SEIRS_fde12   Octave/MATLAB code to solve a certain initial value 
%                    problem with the help of the fde12 function.

%   Revision: 1.1 - Date: March, 20 2023 

% initial conditions
y0=[0.426282; 0.0109566; 0.0275076; 0.535254];

% Values of parameters
miu = 0.0113; niu = 36; epsilon = 91; b0 = 85; b1 = 0.167; c1 = 0.167;
gama = 1.8; tfinal = 5; phi = pi/2;
ft = linspace(0,tfinal,N); h = tfinal/(N-1); 


% Correction of values of parameters
miu_ = miu^alpha; niu_ = niu^alpha; epsilon_ = epsilon^alpha;
gama_ = gama^alpha;

% time-dependent parameters
flambda = @(t) miu.^alpha.*(1 + c1.* cos( 2.* pi.* t + phi) );
fbeta = @(t) b0.^alpha.* (1 + b1.* cos( 2.* pi.* t + phi ) );


% Differential system of equations of the model 
fdefun = @(t,y,ft)[flambda(t)-miu_*y(1)-fbeta(t)*y(1)*y(3)+gama_*y(4); ...
    fbeta(t)*y(1)*y(3)-(miu_+epsilon_)*y(2);
    epsilon_*y(2)-(miu_+niu_)*y(3); ...
    niu_*y(3)-miu_*y(4)-gama_*y(4)]; 

% resolution of system with solver fde12 
[t,y] = fde12(alpha,fdefun,0,tfinal,y0,h,ft);


end
