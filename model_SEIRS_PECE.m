function [t,y] = model_SEIRS_PECE(N,alpha)

%   This code is in an appendix of the following paper:
%
%   [*] Rosa, S.; Torres, D.F.M. Numerical Fractional Optimal Control of 
%   Respiratory Syncytial Virus Infection in Octave/MATLAB. Mathematics 
%   2023, 11, 1511. 
%
%   DOI: https://doi.org/10.3390/math11061511

%model_SEIRS_PECE   Octave/MATLAB code to solve a certain initial value 
%                   problem by using the predict-evaluate-correct-evaluate 
%                  (PECE) method of Adams–Bashforth–Moulton.

%   Revision: 1.1 - Date: March, 20 2023 

% Values of parameters
miu = 0.0113; niu = 36; epsilon = 91;  gama = 1.8; tfinal = 5; 
b0 = 85; b1 = 0.167; c1=0.167; phi = pi/2;

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
beta = fbeta(t); lambda = flambda(t);
S = init; E = init; I = init; R = init; b = init; a = init;
S(1) = S0; E(1) = E0; I(1) = I0; R(1) = R0;
Sp = S; Ep = E; Ip = I; Rp = R;

% computation of coefficients a_k and b_k
for k = 1:N
    b(k) = k^alpha-(k-1)^alpha;
    a(k) = (k+1)^(alpha+1)-2*k^(alpha+1)+(k-1)^(alpha+1);
end

for j = 2:N

    % First part: prediction

    aux_s = 0; aux_e = 0; aux_i = 0; aux_r = 0;
    for k = 1:j

        % Differential system of equations of the model
        aux_s = aux_s+b(j-k+1)*(lambda(k)-miu_*S(k)...
            -beta(k)*S(k)*I(k)+gama_*R(k));
        aux_e = aux_e+b(j-k+1)*(beta(k)*S(k)*I(k)...
            -(miu_+epsilon_)*E(k));
        aux_i = aux_i+b(j-k+1)*(epsilon_*E(k)-(miu_+niu_)*I(k));
        aux_r = aux_r+b(j-k+1)*(niu_*I(k)-miu_*R(k)-gama_*R(k));
    end

    Sp(j) = S0+h^alpha/gamma(1+alpha)*aux_s;
    Ep(j) = E0+h^alpha/gamma(1+alpha)*aux_e;
    Ip(j) = I0+h^alpha/gamma(1+alpha)*aux_i;
    Rp(j) = R0+h^alpha/gamma(1+alpha)*aux_r;

    % Second part: correction

    aux_ss = lambda(j)-miu_*Sp(j)-beta(j)*Sp(j)*Ip(j)+gama_*Rp(j);
    aux_ee = beta(j)*Sp(j)*Ip(j)-(miu_+epsilon_)*Ep(j);
    aux_ii = epsilon_*Ep(j)-(miu_+niu_)*Ip(j);
    aux_rr = niu_*Ip(j)-miu_*Rp(j)-gama_*Rp(j);

    auxx = ((j-1)^(alpha+1)-(j-1-alpha)*j^alpha);
    aux_s0 = auxx*(lambda(1)-miu_*S(1)-beta(1)*S(1)*I(1)+gama_*R(1));
    aux_e0 = auxx* (beta(1)*S(1)*I(1)-(miu_+epsilon_)*E(1));
    aux_i0 = auxx*(epsilon_*E(1)-(miu_+niu_)*I(1));
    aux_r0 = auxx*(niu_*I(1)-miu_*R(1)-gama_*R(1));

    aux_s = 0; aux_e = 0; aux_i = 0; aux_r = 0;
    for k = 1:j-1

        % Differential system of equations of the model
        aux_s = aux_s+a(j-k)*(lambda(k)-miu_*S(k)-beta(k)*S(k)*I(k)...
            +gama_*R(k));
        aux_e = aux_e+a(j-k)*(beta(k)*S(k)*I(k)-(miu_+epsilon_)*E(k));
        aux_i = aux_i+a(j-k)*(epsilon_*E(k)-(miu_+niu_)*I(k));
        aux_r = aux_r+a(j-k)*(niu_*I(k)-miu_*R(k)-gama_*R(k));
    end

    S(j) = S0+h^alpha/gamma(2+alpha)*(aux_ss+aux_s0+aux_s);
    E(j) = E0+h^alpha/gamma(2+alpha)*(aux_ee+aux_e0+aux_e);
    I(j) = I0+h^alpha/gamma(2+alpha)*(aux_ii+aux_i0+aux_i);
    R(j) = R0+h^alpha/gamma(2+alpha)*(aux_rr+aux_r0+aux_r);
end

y(1,:) = S; y(2,:) = E; y(3,:) = I; y(4,:) = R;
end

