fclose all;
close all;
clear;
clc;

%% Load experimental data
input = readmatrix('../TestData/deGraciaTestDcb.csv');
nPts = size(input,1);

%% Iteratively solve for a0
h = 1.5;
b = 15;
Ef = 116e3;
GLT = 4.0e3;
ETran = 7e3;
I = b*h^3/12;
A = b*h;

input = [input, input(:,1)./input(:,2)];
input(1,3) = input(2,3);

C0 = mean(input(1,3));

aOld = 40;
iter = 1;
error = 1;
tol = 1e-6;
while (error > tol) && (iter < 100)
    aNew = ( (3/2*Ef*I) * (C0 / (1+(3/10)*(Ef/GLT)*(h^2/aOld^2))) )^(1/3);
    error = abs(aNew-aOld)/aOld;
    iter = iter+1;
    aOld = aNew;
end
a = aNew

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
    7/(60*Ef*I),...
    (3*gam1 + 9*a)/(20*Ef*I),...
    (gam2 + 6*gam1*a)/(12*Ef*I),...
    (gam1*gam3 + 15*gam2*a)/(60*Ef*I) - (2*a)/(GLT*A),...
    (gam1*gam3*a)/(20*Ef*I) - (gam1*a)/(GLT*A) - (6*h)/(b*ETran),...
    -(2*h*gam1 + 6*h*a)/(b*ETran) ];

x3 = roots(polyCoefs);
x3 = max(x3(x3 == real(x3)))


%% Determine crack length
beta1 = (x1^2 + 3*x1*x2 + 4*x1*x3 + 3*x2^2 + 8*x2*x3 + 5*x3^2) ...
    / (x1 + 2*x2 +2*x3);
beta2 = 3 / x3 / (x1 + 2*x2 + 2*x3);
beta3 = (x1^2*x2 + 3*x1*x3^2 + 3*x1*x2*x3 + 3*x2^2*x3 + 6*x2*x3^2 ...
+ 3*x3^3) / (x1 + 2*x2 +2*x3);
beta4 = (x1 + 2*x2 + 3*x3) / x3 / (x1 + 2*x2*x3);

% Determine current crack length using FPI at each point
crackLen = zeros(nPts,1);
crackLen(1) = a;
GI = zeros(nPts,1);
theta = zeros(nPts, 1);

for iPt = 2:nPts
    comp = input(iPt,3);
    P = input(iPt,2);
    
    goalSeek = @(a) ( (2*a^3)/(3*Ef*I) + (beta1*a^2)/(2*Ef*I) + (2*a)/(A*GLT) ...
    + (4*h*beta2*a)/(b*ETran) + (beta3*a)/(6*Ef*I) + (4*beta4*h)/(b*ETran) - comp);
    
    a = fzero(goalSeek, crackLen(iPt-1));
    
    crackLen(iPt) = a;
    
    GI(iPt) = (P^2*a^2)/(b*Ef*I) + (P^2*beta1*a)/(2*b*Ef*I) ...
        + (P^2)/(b*A*GLT) + (2*h*beta2*P^2)/(b^2*ETran) ...
        + (beta3*P^2)/(12*b*Ef*I);
    
    theta(iPt) = P/(12*Ef*I) * (beta3 + 3*beta1*a);
end
GI(1) = GI(2);

figure(); tiledlayout('flow');
nexttile; hold on;
plot(input(:,1),crackLen)
nexttile; hold on;
plot(crackLen-46,GI)


    