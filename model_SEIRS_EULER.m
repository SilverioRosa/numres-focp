function [t,y] = model_SEIRS_EULER(N,alpha)

%   This code is in an appendix of the following paper:
%
%   [*] Rosa, S.; Torres, D.F.M. Numerical Fractional Optimal Control of 
%   Respiratory Syncytial Virus Infection in Octave/MATLAB. Mathematics
%   2023, 11, 1511. 
%
%   DOI: https://doi.org/10.3390/math11061511

%model_SEIRS_EULER   is an Octave/MATLAB code to solve a certain initial 
%                    value problem  by using the fractional forward Euler's
%                    method to approximate the variables of the 
%                    fractional system of equations.

%   Revision: 1.1 - Date: March, 20 2023 


% Values of parameters
miu = 0.0113; niu = 36; epsilon = 91;  gama = 1.8; tfinal = 5; b0 = 85;
b1 = 0.167; c1 = 0.167; phi = pi/2;

% initial conditions
S0 = 0.426282; E0 = 0.0109566; I0 = 0.0275076; R0 =  0.535254;

% Correction of values of parameters
miu_ = miu^alpha; niu_ = niu^alpha; epsilon_ = epsilon^alpha;
gama_ = gama^alpha;

% time-dependent parameters
flambda = @(t) miu_*(1 + c1 * cos( 2 * pi * t + phi) );
fbeta = @(t) b0^alpha.* (1 + b1 * cos( 2 * pi * t + phi ) );

% Initialization of variables
t = linspace(0,tfinal,N); h = tfinal/N; init = zeros(1,N);
S = init; E = init; I = init; R = init;
S(1) = S0; E(1) = E0; I(1) = I0; R(1) = R0;
beta = fbeta(t); lambda = flambda(t);


for j = 2:N
    aux_s = 0; aux_e = 0; aux_i = 0; aux_r = 0;
    for k = 1:j-1
        bk = (j-k+1)^alpha-(j-k)^alpha;

        % Differential system of equations of the model
        aux_s = aux_s+bk*(lambda(k)-miu_*S(k)-beta(k)*S(k)*I(k) ...
            +gama_*R(k));
        aux_e = aux_e+bk*(beta(k)*S(k)*I(k)-(miu_+epsilon_)*E(k));
        aux_i = aux_i+bk*(epsilon_*E(k)-(miu_+niu_)*I(k));
        aux_r = aux_r+bk*(niu_*I(k)-miu_*R(k)-gama_*R(k));
    end
    S(j) = S0+h^alpha/gamma(1+alpha)*aux_s;
    E(j) = E0+h^alpha/gamma(1+alpha)*aux_e;
    I(j) = I0+h^alpha/gamma(1+alpha)*aux_i;
    R(j) = R0+h^alpha/gamma(1+alpha)*aux_r;
end

y(1,:) = S; y(2,:) = E; y(3,:) = I; y(4,:) = R;
end

