function  [t,y] = FOCP_PECE(N,alpha);

%    This code is in an appendix of the following paper:
%
%    [*] Rosa, S.; Torres, D.F.M. Numerical Fractional Optimal Control of 
%    Respiratory Syncytial Virus Infection in Octave/MATLAB. Mathematics 
%    2023, 11, 1511. 
%
%    DOI: https://doi.org/10.3390/math11061511

%FOCP_PECE   is our Octave/MATLAB code for the numerical solution of the 
%            fractional optimal control problem  with given initial 
%            conditions.

%   Revision: 1.1 - Date: March, 20 2023 

% values assumed as global
global tfinal miu niu epsilon gama  b0 b1 c1 phi k1 k2 S0 E0 I0 R0;


% Values of parameters
miu = 0.0113; niu = 36; epsilon = 91;  gama = 1.8; tfinal = 5;
b0 = 85; b1 = 0.167; phi = pi/2; c1 = .167;

% parameters of the algorithm
k1 = 1; k2 = 0.001; trmax = 1.0; tol = 0.001; test = 1;

% initial conditions
S0 = 0.426282; E0 = 0.0109566; I0 = 0.0275076; R0 =  0.535254;

% initialization of variables
t = linspace(0,tfinal,N);
init = zeros(1,N); S = init; E = init; I = init; R = init;
p1 = init; p2 = init; p3 = init; p4 = init; Ta = init;

% iterations of the numerical method

iter=0;

while test>tol,

    oldS = S; oldE = E; oldI = I; oldR = R;
    oldp1 = p1; oldp2 = p2; oldp3 = p3; oldp4 = p4; oldTa = Ta;

    % forward PECE iterations
    [y1] = system1_control(Ta,t,N,alpha);
    S = y1(1,:); E = y1(2,:); I = y1(3,:); R = y1(4,:);

    % backward PECE iterations
    [y2] = system2_adjoint(S,I,Ta,t,N,alpha);
    p1 = y2(1,:); p2 = y2(2,:); p3 = y2(3,:); p4 = y2(4,:);

    % new control
    Ta = projection((p3-p4).*I/(2*k2),trmax);
    Ta = ( Ta + oldTa ) / 2;

    % Relative error values for convergence
    vector = [max(abs(S-oldS))/(max(abs(S))),...
        max(abs(oldE-E))/(max(abs(E))),...
        max(abs(oldI-I))/(max(abs(I))),...
        max(abs(oldR-R))/(max(abs(R))),...
        max(abs(oldp1-p1))/(max(abs(p1))),...
        max(abs(oldp2-p2))/(max(abs(p2))),...
        max(abs(oldp3-p3))/(max(abs(p3))),...
        max(abs(oldp4-p4))/(max(abs(p4))), ...
        max(abs(oldTa-Ta))/(max(abs(Ta)))]*100;

    test = max(vector);
    iter=iter+1
end

y(1,:) = S; y(2,:) = E; y(3,:) = I; y(4,:) = R; y(5,:) = Ta;
y(6,:) = p1; y(7,:) = p2; y(8,:) = p3; y(9,:) = p4;

end


% function II: resolution of the fractional control system


function [y]= system1_control(Ta,t,N,alpha)

global  b0 b1 c1 phi miu gama epsilon niu tfinal S0 E0 I0 R0;

% time-dependent parameters
flambda = @(t) miu^alpha*(1 + c1 * cos( 2 * pi * t + phi) );
fbeta = @(t) b0^alpha.* (1 + b1 * cos( 2 * pi * t + phi ) );

% Correction of values of parameters
miu_ = miu^alpha; niu_ = niu^alpha; epsilon_ = epsilon^alpha;
gama_ = gama^alpha;

% initialization of variables
beta = fbeta(t); lambda = flambda(t);
h = tfinal/N; init = zeros(1,N);
S = init; E = init; I = init; R = init; a = init; b = init;
S(1) = S0; E(1) = E0; I(1) = I0; R(1) = R0;
Sp = init; Ep = init; Ip = init; Rp = init;

% computation of coefficients a_k and b_k
for k = 1:N
    b(k) = k^alpha-(k-1)^alpha;
    a(k) = (k+1)^(alpha+1)-2*k^(alpha+1)+(k-1)^(alpha+1);
end

for j = 2:N

    % First part: predict

    % differential equations of control system
    aux_s = 0; aux_e = 0; aux_i = 0; aux_r = 0;
    for k = 1:j
        aux_s = aux_s+b(j-k+1)*(lambda(k)-miu_*S(k)...
            -beta(k)*S(k)*I(k)+gama_*R(k));
        aux_e = aux_e+b(j-k+1)*(beta(k)*S(k)*I(k)...
            -(miu_+epsilon_)*E(k));
        aux_i = aux_i+b(j-k+1)*(epsilon_*E(k)-(miu_+niu_+Ta(k))*I(k));
        aux_r = aux_r+b(j-k+1)*(niu_*I(k)-miu_*R(k)-gama_*R(k)...
            +Ta(k)*I(k));
    end

    Sp(j) = S0+h^alpha/gamma(1+alpha)*aux_s;
    Ep(j) = E0+h^alpha/gamma(1+alpha)*aux_e;
    Ip(j) = I0+h^alpha/gamma(1+alpha)*aux_i;
    Rp(j) = R0+h^alpha/gamma(1+alpha)*aux_r;

    % Second part: correct

    aux_ss = lambda(j)-miu_*Sp(j)-beta(j)*Sp(j)*Ip(j)+gama_*Rp(j);
    aux_ee = beta(j)*Sp(j)*Ip(j)-(miu_+epsilon_)*Ep(j);
    aux_ii = epsilon_*Ep(j)-(miu_+niu_+Ta(j))*Ip(j);
    aux_rr = niu_*Ip(j)-miu_*Rp(j)-gama_*Rp(j)+Ta(j)*Ip(j);

    auxx = ((j-1)^(alpha+1)-(j-1-alpha)*j^alpha);
    aux_s0 = auxx*(lambda(1)-miu_*S(1)-beta(1)*S(1)*I(1)+gama_*R(1));
    aux_e0 = auxx* (beta(1)*S(1)*I(1)-(miu_+epsilon_)*E(1));
    aux_i0 = auxx*(epsilon_*E(1)-(miu_+niu_+Ta(1))*I(1));
    aux_r0 = auxx*(niu_*I(1)-miu_*R(1)-gama_*R(1)+Ta(1)*I(1));

    aux_s = 0; aux_e = 0; aux_i = 0; aux_r = 0;
    for k = 1:j-1
        aux_s = aux_s+a(j-k)*(lambda(k)-miu_*S(k)-beta(k)*S(k)*I(k)...
            +gama_*R(k));
        aux_e = aux_e+a(j-k)*(beta(k)*S(k)*I(k)-(miu_+epsilon_)*E(k));
        aux_i = aux_i+a(j-k)*(epsilon_*E(k)-(miu_+niu_+Ta(k))*I(k));
        aux_r = aux_r+a(j-k)*(niu_*I(k)-miu_*R(k)...
            -gama_*R(k)+Ta(k)*I(k));
    end

    S(j) = S0+h^alpha/gamma(2+alpha)*(aux_ss+aux_s0+aux_s);
    E(j) = E0+h^alpha/gamma(2+alpha)*(aux_ee+aux_e0+aux_e);
    I(j) = I0+h^alpha/gamma(2+alpha)*(aux_ii+aux_i0+aux_i);
    R(j) = R0+h^alpha/gamma(2+alpha)*(aux_rr+aux_r0+aux_r);
end

y(1,:) = S; y(2,:) = E; y(3,:) = I; y(4,:) = R;
end


% function III: resolution of the fractional adjoint system


function [y] = system2_adjoint(S,I,Ta,t,N,alpha)

global  miu gama epsilon niu tfinal k1 b0 b1 phi;

% time-dependent parameter
fbeta = @(t) b0^alpha.* (1 + b1 * cos( 2 * pi * t + phi ) );

% Correction of values of parameters
miu_=miu^alpha; niu_=niu^alpha; epsilon_=epsilon^alpha;
gama_=gama^alpha;

% initialization of variables
beta = fbeta(t);
h = tfinal/N; init = zeros(1,N); a = init; b = init;
p1 = init; p2 = init; p3 = init; p4 = init;
p1p = init; p2p = init; p3p = init; p4p = init;

% First part: predict

S = S(end:-1:1); I = I(end:-1:1);
Ta = Ta(end:-1:1); beta = beta(end:-1:1);

% computation of coefficients a_k and b_k
for k = 1:N
    b(k) = k^alpha-(k-1)^alpha;
    a(k) = (k+1)^(alpha+1)-2*k^(alpha+1)+(k-1)^(alpha+1);
end

for j = 2:N

    % differential equations of adjoint system
    aux_p1 = 0; aux_p2 = 0; aux_p3 = 0; aux_p4 = 0;
    for k = 1:j
        aux_p1 = aux_p1+b(j-k+1)*(-1)*(p1(k)*(miu_+beta(k)*I(k))- ...
            beta(k)*I(k)*p2(k));
        aux_p2 = aux_p2+b(j-k+1)*(-1)*(p2(k)*(miu_+epsilon_)...
            -epsilon_*p3(k));
        aux_p3 = aux_p3+b(j-k+1)*(-1)*(-k1+beta(k)*p1(k)*S(k)...
            -p2(k)*beta(k)*S(k)+p3(k)*(miu_+niu_+Ta(k))...
            -p4(k)*(niu_+Ta(k)));
        aux_p4 = aux_p4+b(j-k+1)*(-1)*(-gama_*p1(k)...
            +p4(k)*(miu_+gama_));
    end

    p1p(j) = h^alpha/gamma(1+alpha)*aux_p1;
    p2p(j) = h^alpha/gamma(1+alpha)*aux_p2;
    p3p(j) = h^alpha/gamma(1+alpha)*aux_p3;
    p4p(j) = h^alpha/gamma(1+alpha)*aux_p4;

    % Second part: correct

    aux_pp1 = (-1)*(p1p(j)*(miu_+beta(j)*I(j))-beta(j)*I(j)*p2p(j));
    aux_pp2 = (-1)*(p2p(j)*(miu_+epsilon_)-epsilon_*p3p(j));
    aux_pp3 = (-1)*(-k1+beta(j)*p1p(j)*S(j)-p2p(j)*beta(j)*S(j)...
        +p3p(j)*(miu_+niu_+Ta(j))-p4p(j)*(niu_+Ta(j)));
    aux_pp4 = (-1)*(-gama_*p1p(j)+p4p(j)*(miu_+gama_));

    auxx = (-1)*((j-1)^(alpha+1)-(j-1-alpha)*j^alpha);
    aux_p10 = auxx*(p1(1)*(miu_+beta(1)*I(1))-beta(1)*I(1)*p2(1));
    aux_p20 = auxx*(p2(1)*(miu_+epsilon_)-epsilon_*p3(1));
    aux_p30 = auxx*(-k1+beta(1)*p1(1)*S(1)-p2(1)*beta(1)*S(1)...
        +p3(1)*(miu_+niu_+Ta(1))-p4(1)*(niu_+Ta(1)));
    aux_p40 = auxx*(-gama_*p1(1)+p4(1)*(miu_+gama_));

    aux_p1 = 0; aux_p2 = 0; aux_p3 = 0; aux_p4 = 0;
    for k = 1:j-1
        aux_p1 = aux_p1+a(j-k)*(-1)*(p1(k)*(miu_+beta(k)*I(k))- ...
            beta(k)*I(k)*p2(k));
        aux_p2 = aux_p2+a(j-k)*(-1)*( p2(k)*(miu_+epsilon_)...
            -epsilon_*p3(k));
        aux_p3 = aux_p3+a(j-k)*(-1)*( -k1+beta(k)*p1(k)*S(k)...
            -p2(k)*beta(k)*S(k)+p3(k)*(miu_+niu_+Ta(k))...
            -p4(k)*(niu_+Ta(k)));
        aux_p4 = aux_p4+a(j-k)*(-1)*(-gama_*p1(k)+p4(k)*(miu_+gama_));
    end

    p1(j) = h^alpha/gamma(2+alpha)*(aux_pp1+aux_p10+aux_p1);
    p2(j) = h^alpha/gamma(2+alpha)*(aux_pp2+aux_p20+aux_p2);
    p3(j) = h^alpha/gamma(2+alpha)*(aux_pp3+aux_p30+aux_p3);
    p4(j) = h^alpha/gamma(2+alpha)*(aux_pp4+aux_p40+aux_p4);
end

y(1,:) = p1(end:-1:1); y(2,:) = p2(end:-1:1); y(3,:) = p3(end:-1:1);
y(4,:) = p4(end:-1:1);
end


% function IV: control projection over the set of admissible controls


function [v] = projection(vect,trmax)

isNeg = vect<0; vect(isNeg) = 0;
isHuge = vect>trmax; vect(isHuge) = trmax;
v = vect;

end