%% MATLAB initialization
fclose all;
close all;
clear;
clc;

%% Load example force-displacement data
% Taken from Xavier et al. (2014)
forceDisp = readmatrix('XavierEnfTestData.csv');

%% Create the dcbCalculator object and specify necessary parameters
% Use consistant units. Unit system here is kg-mm-s (MPa)

% Constructor for enfCrackEquivMethod takes no argument
enfCalculator = enfCrackEquivMethod();
% Use name-value arguments to define specimen half height (h), specimen
% width (b), and half length between sports (l)
enfCalculator.setSpecimenGeom('h', 20/2, 'b', 20, 'L', 230/2);
% Define material properties. Only shear modulus is required
enfCalculator.setMaterialProps('G13', 1.12e3);

%% Extract resistance curve
% This portion of the code can be placed in a loop to process multiple
% specimens simultanously. 

% Fit specimen-specific initial compliance. Use linear region of
% force-displacement curve. Hence the indexing below.
C0 = enfCalculator.fitCompliance(forceDisp(1:25,1), forceDisp(1:25,2));

% Evaluate energy release rate and effective crack length
[GII, aE] = enfCalculator.computeCerr(...
                forceDisp(:,1), forceDisp(:,2), 162.0, C0);
            
%% Plot force-displacement and resistance curves
figure();
subplot(1,2,1); hold on;
plot(forceDisp(:,1), forceDisp(:,2), '.-', 'displayName', 'Exp.')
title('Force-Displacement')
xlabel('Crosshead Displacement (mm)')
ylabel('Crosshead Force (N)')
grid on;

subplot(1,2,2); hold on;
plot(aE, GII, '.-', 'DisplayName', 'R-Curve')
title('R-Curve')
xlabel('Effective Crack Length (mm)')
ylabel('Energy Release Rate (J/mm^2)')
grid on
