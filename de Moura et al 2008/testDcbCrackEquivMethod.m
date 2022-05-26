fclose all;
% close all;
clear;
clc;

% load data
load('../TestData/Dcb_ForceDisp.mat')
nSpec = length(responseCurves);
% [responseCurves(:).('a0')] = ...
%     deal(37.64, 33.0, 25.4, 35.0, 37.8, 37.8, 37.0, 37.5);
[responseCurves(:).('a0')] = ...
    deal(37.64, 33.0, 30.0, 35.0, 40, 37.8, 33, 37.5);

% set up dcbCrackEquivMethod and fixed properties
dcbCalculator = dcbCrackEquivMethod();
dcbCalculator.setSpecimenGeom('h', 2.7, 'b', 25);
dcbCalculator.setMaterialProps('E3', 6.11e3, 'G13', 2.288e3)

% cycle through all DCB specimens
for iSpec = 1:4
    forceDisp = responseCurves(iSpec).data;
    
    C0 = dcbCalculator.fitCompliance(forceDisp(1:25,1), forceDisp(1:25,2));
    
    Ef = dcbCalculator.fitFlexuralModulus(C0, responseCurves(iSpec).a0);
    
    [responseCurves(iSpec).GI, responseCurves(iSpec).aE] =...
        dcbCalculator.computeCerr(forceDisp(:,1), forceDisp(:,2),...
            Ef);
end
% Sampling rate changes, need to update indices for C0 estimationg
for iSpec = 5:8
    forceDisp = responseCurves(iSpec).data;
    
    C0 = dcbCalculator.fitCompliance(forceDisp(1:50,1), forceDisp(1:50,2));
    
    Ef = dcbCalculator.fitFlexuralModulus(C0, responseCurves(iSpec).a0);
    
    [responseCurves(iSpec).GI, responseCurves(iSpec).aE] =...
        dcbCalculator.computeCerr(forceDisp(:,1), forceDisp(:,2),...
            Ef);
end

% % cycle through all DCB specimens (de Gracia Method)
% for iSpec = 1:4
%     forceDisp = responseCurves(iSpec).data;
%     
%     C0 = dcbCalculator.fitCompliance(forceDisp(1:25,1), forceDisp(1:25,2));
%     
%     Ef = dcbCalculator.fitFlexuralModulus(C0, responseCurves(iSpec).a0);
%     
%     [responseCurves(iSpec).GI, responseCurves(iSpec).aE] =...
%         dcbCalculator.computeCerrDeGracia(forceDisp(:,1), forceDisp(:,2),...
%             Ef, C0);
% end
% % Sampling rate changes, need to update indices for C0 estimationg (de
% % Gracia)
% for iSpec = 5:8
%     forceDisp = responseCurves(iSpec).data;
%     
%     C0 = dcbCalculator.fitCompliance(forceDisp(1:50,1), forceDisp(1:50,2));
%     
%     Ef = dcbCalculator.fitFlexuralModulus(C0, responseCurves(iSpec).a0);
%     
%     [responseCurves(iSpec).GI, responseCurves(iSpec).aE] =...
%         dcbCalculator.computeCerrDeGracia(forceDisp(:,1), forceDisp(:,2),...
%             Ef, C0);
% end

% Plot force disp and R-curve
figure();
subplot(2,1,1); hold on;
for iSpec = 1:nSpec
    plot(responseCurves(iSpec).data(:,1),...
        responseCurves(iSpec).data(:,2),...
        'DisplayName', responseCurves(iSpec).specId);
end
xlabel('Crosshead Disp (mm)')
ylabel('Crosshead Force (N)')
legend()
grid on;

subplot(2,1,2); hold on;
for iSpec = 1:nSpec
    plot(responseCurves(iSpec).aE, responseCurves(iSpec).GI,...
        'DisplayName', responseCurves(iSpec).specId);
end
xlabel('Effective Crack Length (mm)')
ylabel('Critical Energy Release Rate (J/mm^2)')
legend();
grid on;

% Plot force disp and corrected R-Curve
figure();
subplot(2,1,1); hold on;
for iSpec = 1:nSpec
    plot(responseCurves(iSpec).data(:,1),...
        responseCurves(iSpec).data(:,2),...
        'DisplayName', responseCurves(iSpec).specId);
end
xlabel('Crosshead Disp (mm)')
ylabel('Crosshead Force (N)')
legend()
grid on;

subplot(2,1,2); hold on;
for iSpec = 1:nSpec
    index = responseCurves(iSpec).GI > 0.1;
    rCurve = [responseCurves(iSpec).aE(index), responseCurves(iSpec).GI(index)];
    
    polyCoefs = polyfit(rCurve(1:25,2), rCurve(1:25,1), 1);
    offset = polyval(polyCoefs, 0.1);
    
    rCurve(:,1) = rCurve(:,1) - offset;
    
    plot(rCurve(:,1), rCurve(:,2),...
        'DisplayName', responseCurves(iSpec).specId);
end
xlabel('Effective Crack Length (mm)')
ylabel('Critical Energy Release Rate (J/mm^2)')
legend();
grid on;

