fclose all;
close all;
clear;
clc;

% file names for crack-tracked method
fileNames = {...
    'DCB-01_Toughness.csv',...
    'DCB-02_Toughness.csv',...
    'DCB-03_Toughness.csv',...
    'DCB-05_Toughness.csv',...
    'DCB-06_Toughness.csv',...
    'DCB-07_Toughness.csv',...
    'DCB-08_Toughness.csv',...
    'DCB-09_Toughness.csv'};

% load crackEquivData
load('CrackEquivMethod.mat')
crackEquiv = responseCurves;

for i=1:8
    % load crack dta
    astmMethod = readmatrix(fileNames{i});
    figure('name', crackEquiv(i).specId);
    % plot uncorrected lengths
    subplot(1,2,1); hold on;
    plot(astmMethod(:,4), astmMethod(:,5),'.-', 'DisplayName','Astm')
    plot(crackEquiv(i).aE, crackEquiv(i).GI, '.-', 'DisplayName','Crack Equiv.')
    xlabel('Crack Length/ Effective Crack Length (mm)')
    ylabel('CERR (J/mm^2)')
    legend('Location','SouthEast');
    
    % plot crack extentions, not lengths
    index = crackEquiv(i).GI > 0.1;
    rCurve = [crackEquiv(i).aE(index), crackEquiv(i).GI(index)];
    
    polyCoefs = polyfit(rCurve(1:25,2), rCurve(1:25,1), 1);
    offset = polyval(polyCoefs, 0.1);

    rCurve(:,1) = rCurve(:,1) - offset;
    
    subplot(1,2,2); hold on;
        plot(astmMethod(:,4)-astmMethod(1,4), astmMethod(:,5),'.-', 'DisplayName','Astm')
    plot(rCurve(:,1), rCurve(:,2), '.-',...
        'DisplayName', 'CrackEquiv');
    xlabel('Crack Extension/Effective Crack Extension (mm)')
    ylabel('CERR (J/mm^2)')
    legend('Location','SouthEast');
end
    