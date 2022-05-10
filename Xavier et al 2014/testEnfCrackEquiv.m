fclose all;
close all;
clear;
clc;

% load data
load('../TestData/ProcessedEnfData.mat')
nSpec = length(enfData);
[enfData(:).('a0')] = deal(34.7, 35.8, 33.5, 35.9, 36.2, 32.9);

% Setup enfCrackEquiv class and fixed properties
enfCalculator = enfCrackEquivMethod();
enfCalculator.setSpecimenGeom('h', 2.7, 'b', 25, 'L', 50);
enfCalculator.setMaterialProps('G13', 2.288e3)

% cycle through all ENF specimens
for iSpec = 1:nSpec
    forceDisp = enfData(iSpec).forceDisp;
    C0 = enfCalculator.fitCompliance(forceDisp(1:25,1), forceDisp(1:25,2));
    
    [enfData(iSpec).GII, enfData(iSpec).aE] =...
        enfCalculator.computeCerr(forceDisp(:,1), forceDisp(:,2),...
            enfData(iSpec).a0, C0);
end

% Plot force disp and R-curve
figure();
subplot(2,1,1); hold on;
for iSpec = 1:nSpec
    plot(enfData(iSpec).forceDisp(:,1), enfData(iSpec).forceDisp(:,2),...
        'DisplayName', enfData(iSpec).specId);
end
xlabel('Crosshead Disp (mm)')
ylabel('Crosshead Force (N)')
legend()
grid on;

subplot(2,1,2); hold on;
for iSpec = 1:nSpec
    plot(enfData(iSpec).aE, enfData(iSpec).GII,...
        'DisplayName', enfData(iSpec).specId);
end
xlabel('Effective Crack Length (mm)')
ylabel('Critical Energy Release Rate (J/mm^2)')
legend();
grid on;

% Plot force disp and corrected R-Curve
figure();
subplot(2,1,1); hold on;
for iSpec = 1:nSpec
    plot(enfData(iSpec).forceDisp(:,1), enfData(iSpec).forceDisp(:,2),...
        'DisplayName', enfData(iSpec).specId);
end
xlabel('Crosshead Disp (mm)')
ylabel('Crosshead Force (N)')
legend()
grid on;

subplot(2,1,2); hold on;
for iSpec = 1:nSpec
    index = enfData(iSpec).GII > 0.7;
    rCurve = [enfData(iSpec).aE(index), enfData(iSpec).GII(index)];
    
    polyCoefs = polyfit(rCurve(1:25,2), rCurve(1:25,1), 1);
    offset = polyval(polyCoefs, 0.7);
    
    rCurve(:,1) = rCurve(:,1) - offset;
    
    plot(rCurve(:,1), rCurve(:,2),...
        'DisplayName', enfData(iSpec).specId);
end
xlabel('Effective Crack Length (mm)')
ylabel('Critical Energy Release Rate (J/mm^2)')
legend();
grid on;





    

