fclose all;
close all;
clear;
clc;

%% Load experimental data
input = readmatrix('../TestData/DCB-01_Toughness.csv');
nPts = size(input,1);

%% Unit system - Use Consistant Units
% mm-kg-s (force: N, stress: MPa)

%% Input parameters 
b = 25; % specimen width
h = 2.7; % half specimen thickness
L = 123; % Total specimen length
a0 = input(1,4); % Initial crack length
dFdd = input(20,3)/input(20,2); % Initial stiffness of DCB
C0 = 1/dFdd; % Initial compliance of DCB
ELong = 19.3e3; % Approximate longitudinal modulus
ETran = 6.11e3; % Approximate transverse modulus
GLT = 2.288e3; % Approximate shear modulus between long. and tran. directions

%% Computing Ef
% Procedure: 1) Assume a value of Ef based on E1, 1) compute Gamma, 2) Use 
% Gamma to Delta (root rotation correction), 3) Solve for Ef using Delta. 
% Iterate on new Ef. 

% Generate initial guess of Ef
EfOld = L^3/(4*b*(2*h)^3)*(dFdd);

error = 1;
iter = 1;
tol = 1e-5;

while (error > tol) && (iter < 100)
    Gamma = 1.18*sqrt(EfOld*ETran)/GLT;
    
    Delta = h*sqrt( (EfOld/11/GLT)*(3 - 2*(Gamma/(1+Gamma))^2) );
    
    EfNew = (C0 - 12*(a0+abs(Delta))/(5*b*h*GLT))^(-1) ...
                * 8*(a0+abs(Delta))^3/(b*h^3);
    
    error = (EfNew-EfOld)/EfOld;
    EfOld = EfNew;
    iter = iter+1;
end

Ef = EfOld;
    
%% Now compute effective crack length ae

% initialize ae and Cerr
ae = zeros(1, nPts);
Cerr = zeros(1,nPts);

for i=1:nPts
    delta = input(i,2);
    F = input(i, 3);
    
    alpha = 8/(b*h^3*Ef);
    beta = 12/(5*b*h*GLT);
    gamma = -delta/F;
    
    A = ( alpha^2 * (-108*gamma + 12*sqrt(3/alpha ...
            * (4*beta^3 + 27*gamma^2 * alpha))))^(1/3);
    
    ae(i) = A/(6*alpha) - 2*beta/A;
    
    Cerr(i) = (6*F^2)/(b^2*h) * ( (2*ae(i)^2)/(h^2*Ef) + 1/(5*GLT) );
end

figure();
subplot(2,2,1); hold on;
plot(ae(2:end),'.-', 'DisplayName', 'de Moura')
plot(input(:,4),'.-', 'DisplayName', 'Optical')
legend()

subplot(2,2,2); hold on;
plot(Cerr(2:end),'.-', 'DisplayName', 'de Moura')
plot(input(:,5), '.-', 'DisplayName', 'Optical')
legend()

subplot(2,2, 3:4); hold on;
plot(ae(2:end),Cerr(2:end),'.-', 'DisplayName', 'de Moura')
plot(input(:,4), input(:,5), '.-', 'DisplayName', 'Optical')
legend()
