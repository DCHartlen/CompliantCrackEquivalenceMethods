%% MATLAB initialization
fclose all;
close all;
clear;
clc;

%% Load example force-displacement data
% Taken from de Gracia et all (2015)
forceDisp = readmatrix('deGraciaDcbTestData.csv');

%% Create the dcbCalculator object and specify necessary parameters
% Use consistant units. Unit system here is kg-mm-s (MPa)

% Method can be set to 'deGracia' for the method of de Gracia et al (2015)
% or 'deMoura' for the method of de Moura et al (2008)
dcbCalculator = dcbCrackEquivMethod('deGracia');
% Use name-value arguments to define specimen half height, h, and specimen
% width, b. 
dcbCalculator.setSpecimenGeom('h', 3/2, 'b', 15);
% Define material properties. 
dcbCalculator.setMaterialProps('E3', 7.0e3, 'G13', 4.0e3);

%% Extract resistance curve
% This portion of the code can be placed in a loop to process multiple
% specimens simultanously. 

% First, estimate initial compliance. Use linear portion of
% force-displacement curve, hence the indexing below. 
C0 = dcbCalculator.fitCompliance(forceDisp(1:20,1), forceDisp(1:20,2));
% Now fit specimen specific flexural modulus
Ef = dcbCalculator.fitFlexuralModulus(C0, 40);
% Compute energy release rate and effective crack length
[GI, aE] = dcbCalculator.computeCerr(forceDisp(:,1), forceDisp(:,2));

%% Plot force-displacement and resistance curves
figure();
subplot(1,2,1); hold on;
plot(forceDisp(:,1), forceDisp(:,2), '.-', 'displayName', 'Exp.')
title('Force-Displacement')
xlabel('Crosshead Displacement (mm)')
ylabel('Crosshead Force (N)')
grid on;

subplot(1,2,2); hold on;
plot(aE, GI, '.-', 'DisplayName', 'R-Curve')
title('R-Curve')
xlabel('Effective Crack Length (mm)')
ylabel('Energy Release Rate (J/mm^2)')
grid on