classdef dcbCrackEquivMethod < handle
    % dcbCrackEquivMethod
    % A class-based approach to extracting the resistance curve of a DCB
    % specimen using the crack equivalence method of de Moura et al (2008)
    
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
        function self = dcbCrackEquivMethod(varargin)
            % dcbCrackEquivMethod
            % Constructor
            
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
            
            % Area and momemt of inerta for single arm
            self.A = self.b*self.h;
            self.I = self.b*self.h^3/12;
        end
        
        function C0 = fitCompliance(self, disp, force)
            % fitCompliance
            % Provided a subset of force-displacement curve, return
            % compliance
            
            % Use linear regression, does not assume a zero intercept
            C0 = polyfit(disp, force, 1);
            C0 = C0(1)^(-1);
            
        end
        
        function Ef = fitFlexuralModulus(self, C0, a0)
            % fitFlexuralModulus
            % Performs an iterative method to fit flexural modulus for
            % specific specimen crack length and initial compliance
            
            % set FPI properties
            error = 1;
            iter = 1;
            tol = 1e-5;
            
            % generate initial guess for flexural modulus
            EfOld = self.E3 .* 3;
            
            % Perform fixed point iteration to find subject-specific Ef
            while (error > tol) && (iter < 100)
                Gamma = 1.18*sqrt(EfOld*self.E3)/self.G13;
                
                Delta = self.h * sqrt( (EfOld/11/self.G13) ...
                    *(3 - 2*(Gamma/(1+Gamma))^2) );
                
                EfNew = (C0 - 12*(a0+abs(Delta))/....
                    (5*self.b*self.h*self.G13))^(-1) ...
                    * 8*(a0+abs(Delta))^3/(self.b*self.h^3);
                
                error = (EfNew-EfOld)/EfOld;
                EfOld = EfNew;
                iter = iter+1;
            end
            
            Ef = EfOld;
        end
        
        function [Gi, aEquiv] = computeCerr(self, disp, force, Ef)
            % computeCerr
            % Computes critical energy release rate (GI) and effective
            % crack length (aEquiv)
            
            alpha = 8 / (self.b*self.h^3*Ef);
            beta = 12 / (5*self.b*self.h*self.G13);
            gamma = -disp./force;
            
            Afact = ( alpha^2 .* (-108.*gamma + 12 .* sqrt(3/alpha ...
                .* (4*beta^3 + 27.*gamma.^2.*alpha))) ).^(1/3);
            
            aEquiv = Afact./(6*alpha) - (2*beta)./Afact;
            
            Gi = (6.*force.^2)./(self.b^2*self.h) .* ( (2 .* aEquiv.^2) ...
                ./ (self.h^2*Ef) + 1/(5*self.G13) );
            
        end
        
    end
    
end
        
            
            
                