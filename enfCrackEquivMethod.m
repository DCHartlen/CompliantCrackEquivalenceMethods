classdef enfCrackEquivMethod < handle
    % enfCrackEquiv
    % class-based approach for processing ENF data without crack tracking.
    % Based on the Xavier et al (2014) methodology
    
    properties
        E1
        E3
        G13
        h
        b
        L
        A
        I
        a0
        C0
    end
    
    methods
        function self = enfCrackEquivMethod(varargin)
            
        end
        
        function setMaterialProps(self, varargin)
            % setMaterialProps
            % Allows user to set approximate material properties
            nvArg = inputParser;
            addParameter(nvArg, 'E1', []);
            addParameter(nvArg, 'E3', []);
            addParameter(nvArg, 'G13', []);
            nvArg.KeepUnmatched = true;
            nvArg.CaseSensitive = false;
            parse(nvArg, varargin{:});
            nvArg = nvArg.Results;
            
            self.E1 = nvArg.E1;
            self.E3 = nvArg.E3;
            self.G13 = nvArg.G13; 
        end
        
        function setSpecimenGeom(self, varargin)
            % setSpecimenGeom
            % Allows user to set overall specimen geometries
            nvArg = inputParser;
            addParameter(nvArg, 'h', []);
            addParameter(nvArg, 'b', []);
            addParameter(nvArg, 'L', []);
            nvArg.KeepUnmatched = true;
            nvArg.CaseSensitive = false;
            parse(nvArg, varargin{:});
            nvArg = nvArg.Results;
            
            self.h = nvArg.h;
            self.b = nvArg.b;
            self.L = nvArg.L;
            
            self.A = 2*self.b*self.h;
            self.I = 8*self.b*self.h^3/12;
        end
        
        function C0 = fitCompliance(self, disp, force)
            % fitCompliance
            % Provided a subset of force-displacement curve, return
            % compliance
            
            % Use linear regression, does not assume a zero intercept
            C0 = polyfit(disp, force, 1);
            C0 = C0(1)^(-1);
            
        end
        
        function [GII, aEquiv] = computeCerr(self, disp, force, a0, C0)
            % computeCerr
            % Compute effective crack length and mode II CERR. Internally
            % calculates flexural modulus and compliance.
            
            % Compute correct compliances and flexural modulus
            Ef = (3*a0^3 + 2*self.L^3)/(12*self.I) * ...
                (C0 - (3*self.L)/(5*self.G13*self.A))^(-1);
            
            C0Corr = C0 - (3*self.L)/(5*self.G13*self.A);
            
            C = disp./force;
            CCorr = C - (3*self.L)/(5*self.G13*self.A);
            
            % For all force-displacement points, compute CERR and aE           
            aEquiv = (CCorr./C0Corr.*a0^3 + (2*self.L^3/3)...
                .*(CCorr./C0Corr - 1)).^(1/3);
            
            GII = (9.*force.^2) / (16*self.b^2*Ef*self.h^3) ...
                .* (CCorr./C0Corr.*a0.^3 + (2*self.L^3/3)...
                .*(CCorr./C0Corr - 1)).^(2/3);
            
        end
        
    end
    
end
            
            
            
            
            
            
            
            