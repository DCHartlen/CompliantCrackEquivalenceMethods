fclose all;
close all;
clear;
clc;

%% Load experimental data
input = readmatrix('../TestData/DCB-01_Toughness.csv');
nPts = size(input,1);

%% Unit system - Use Consistant Units
% mm-kg-s (force: N, stress: MPa)

% %% Input parameters 
% b = 25; % specimen width
% h = 2.7; % half specimen thickness
% L = 123; % Total specimen length
% a0 = input(1,4); % Initial crack length
% dFdd = input(20,3)/input(20,2); % Initial stiffness of DCB
% C0 = 1/dFdd; % Initial compliance of DCB
% ELong = 19.3e3; % Approximate longitudinal modulus
% ETran = 6.11e3; % Approximate transverse modulus
% GLT = 2.288e3; % Approximate shear modulus between long. and tran. directions
% I = b*h^3/12; % Second moment of inertia
% A = b*h; % x-sectional area
% 
% %% Estimate Ef
% % Procedure: 1) Assume a value of Ef based on E1, 1) compute Gamma, 2) Use 
% % Gamma to Delta (root rotation correction), 3) Solve for Ef using Delta. 
% % Iterate on new Ef. 
% 
% % Generate initial guess of Ef
% EfOld = L^3/(4*b*(2*h)^3)*ELong;
% 
% error = 1;
% iter = 1;
% tol = 1e-5;
% 
% while (error > tol) && (iter < 100)
%     Gamma = 1.18*sqrt(EfOld*ETran)/GLT;
%     
%     Delta = h*sqrt( (EfOld/11/GLT)*(3 - 2*(Gamma/(1+Gamma))^2) );
%     
%     EfNew = (C0 - 12*(a0+abs(Delta))/(5*b*h*GLT))^(-1) ...
%                 * 8*(a0+abs(Delta))^3/(b*h^3);
%     
%     error = (EfNew-EfOld)/EfOld;
%     EfOld = EfNew;
%     iter = iter+1;
% end
% 
% Ef = EfOld;
% % Ef = ELong

%% Iteratively solve for a0

h = 1.5;
b = 15;
Ef = 116e3;
GLT = 4e3;
I = b*h^3/12;
A = b*h;
C0 = (0.9375-.1)/4.9752;
ETran = 3e3;

aOld = 51;
iter = 1;
error = 1;
tol = 1e-6;
while (error > tol) && (iter < 100)
    aNew = ( (3/2*Ef*I) * (C0 / (1+(3/10)*(Ef/GLT)*(h^2/aOld^2))) )^(1/3);
    error = abs(aNew-aOld)/aOld;
    iter = iter+1;
    aOld = aNew;
end
a = aNew;

%% Solve for stress distribution factors x1, x2, x3
x1 = [...
    h/sqrt(6*GLT) * sqrt(5*Ef + 5*Ef*sqrt(1 - (36*GLT^2)/(5*Ef*ETran))),...
    h/sqrt(6*GLT) * sqrt(5*Ef - 5*Ef*sqrt(1 - (36*GLT^2)/(5*Ef*ETran)))]
x1 = max(x1)


x2 = (x1/2)*(-1 + sqrt(-1 + (10*Ef*h^2)/(3*GLT*x1^2)))

gam1 = x1 + 2*x2;
gam2 = x1^2 + 3*x1*x2 + 3*x2^2;
gam3 = x1^2 + 2*x1*x2 + 2*x2^2;

polyCoefs = [...
    7*(60*Ef*I),...
    (3*gam1 + 9*a)/(20*Ef*I),...
    (gam2 + 6*gam1*a)/(12*Ef*I),...
    (gam1*gam3 + 15*gam2*a)/(60*Ef*I) - (2*a)/(GLT*A),...
    (gam1*gam3*a)/(20*Ef*I) - (gam1*a)/(GLT*A) - (6*h)/(b*ETran),...
    -(2*h*gam1 + 6*h*a)/(b*ETran) ];

x3 = roots(polyCoefs);
x3 = max(x3(x3 == real(x3)))

