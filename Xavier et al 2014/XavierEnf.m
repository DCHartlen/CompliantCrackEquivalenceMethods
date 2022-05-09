fclose all;
close all;
clear;
clc;

% Load experimental data
input = readmatrix('../TestData/XavierEnfExpTestData.csv');
nPts = size(input,1);

% Mat properterties
E1 = 12.18e3;
E3 = 1.91e3;
G13 = 1.12e3;

% Geometry and specimen properties
b = 20;
h = 10;
L = 230;
a0 = 162;

A = 2*b*h;
I = 8*b*h^3/12;

C0 = polyfit(input(1:10,2), input(1:10,1), 1);
C0 = C0(1);
C0Corr = C0 - (3*L)/(5*G13*A);

Ef = (3*a0^3 + 2*L^3)/(12*I) * (C0 - (3*L)/(5*G13*A))^(-1);

input = [input,input(:,1)./input(:,2)];

GII = zeros(nPts,1);
for i = 1:nPts
    C = input(i,3);
    
    CCorr = C - (3*L)/(5*G13*A);

    ae(i,1) = (CCorr/C0Corr*a0^3 + (2/3)*(CCorr/C0Corr - 1)*L^3)^(1/3);

    GII(i,1) = (9*input(i,2)^2)/(16*b^2*Ef*h^3) * ...
        (CCorr/C0Corr*a0^3 + (2/3)*(CCorr/C0Corr-1)*L^3)^(2/3);
end

figure();
subplot(1,2,1); hold on;
plot(ae(5:end), '.-')
subplot(1,2,2); hold on;
plot(ae(5:end),GII(5:end), '.-')
