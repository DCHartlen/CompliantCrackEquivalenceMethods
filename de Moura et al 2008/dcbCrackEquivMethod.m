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
        Ef
    end
    
    properties (Access = private)
        method
    end
    
    methods
        function self = dcbCrackEquivMethod(varargin)
            % dcbCrackEquivMethod
            % Constructor
            % If method type has been defined
            if length(varargin) == 1
                if strcmpi(varargin{1}, 'deMoura')
                    self.method = 'deMoura';
                    fprintf('Using method of de Moura et al. (2008)\n');
                elseif strcmpi(varargin{1}, 'deGracia')
                    self.method = 'deGracia';
                    fprintf('Using method of de Gracia et al. (2015)\n');
                else
                    error(['%s is not a defined option. Available options are "deMoura" and "deGracia"\n'],...
                        varargin{1})
                end
            else
                self.method = 'deMoura';
                fprintf('Defaulting to method of de Moura et al. (2008)\n')
            end        
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
            self.C0 = C0;
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
            self.Ef = Ef;
        end
        
        function [Gi, aEquiv] = computeCerr(self, disp, force)
            % computeCerr
            % Carries out calculation of critical energy release rate (Gi)
            % and equivalent crack length(aEquiv) using the method
            % specified during class initialization. 
            
            if strcmpi(self.method, 'deMoura')
                [Gi, aEquiv] = computeCerrDeMoura(self, disp, force, ...
                    self.Ef);
            elseif strcmpi(self.method, 'deGracia')
                [Gi, aEquiv] = computeCerrDeGracia(self, disp, force, ...
                    self.Ef, self.C0);
            else
                error('Method not defined. Aborting');
            end
        end
        
        function [Gi, aEquiv] = computeCerrDeMoura(self, disp, force, Ef)
            % computeCerrDeMoura
            % Computes critical energy release rate (GI) and effective
            % crack length (aEquiv) using the method of de Moura et al
            % (2008). 
            
            alpha = 8 / (self.b*self.h^3*Ef);
            beta = 12 / (5*self.b*self.h*self.G13);
            gamma = -disp./force;
            
            Afact = ( alpha^2 .* (-108.*gamma + 12 .* sqrt(3/alpha ...
                .* (4*beta^3 + 27.*gamma.^2.*alpha))) ).^(1/3);
            
            aEquiv = Afact./(6*alpha) - (2*beta)./Afact;
            
            Gi = (6.*force.^2)./(self.b^2*self.h) .* ( (2 .* aEquiv.^2) ...
                ./ (self.h^2*Ef) + 1/(5*self.G13) );
            
        end
        
        function [Gi, aEquiv] = computeCerrDeGracia(self, disp, force, Ef, C0)
            % computeCerrDeGracia
            % Computes critical energy release rate (GI) and effective
            % crack length (aEquiv) using the method of de Gracia et al
            % (2015).
            
            % Iteratively solve for a0
            aOld = 40;
            iter = 1;
            error = 1;
            tol = 1e-6;
            while (error > tol) && (iter < 100)
                aNew = ( (3/2*Ef*self.I) * (C0 ...
                    / (1+(3/10)*(Ef/self.G13)*(self.h^2/aOld^2))) )^(1/3);
                error = abs(aNew-aOld)/aOld;
                iter = iter+1;
                aOld = aNew;
            end
            a = aNew;
            
            % Solve for stress distribution factors x1, x2, x3
            x1 = [...
                self.h/sqrt(6*self.G13) * sqrt(5*Ef + 5*Ef*sqrt(1 ...
                - (36*self.G13^2)/(5*Ef*self.E3))), ...
                self.h/sqrt(6*self.G13) * sqrt(5*Ef ...
                - 5*Ef*sqrt(1 - (36*self.G13^2)/(5*Ef*self.E3)))];
            x1 = max(x1);
            
            
            x2 = (x1/2)*(-1 + sqrt(-1 + (10*Ef*self.h^2)/(3*self.G13*x1^2)));
            
            gam1 = x1 + 2*x2;
            gam2 = x1^2 + 3*x1*x2 + 3*x2^2;
            gam3 = x1^2 + 2*x1*x2 + 2*x2^2;
            
            polyCoefs = [...
                7/(60*Ef*self.I),...
                (3*gam1 + 9*a)/(20*Ef*self.I),...
                (gam2 + 6*gam1*a)/(12*Ef*self.I),...
                (gam1*gam3 + 15*gam2*a)/(60*Ef*self.I) - (2*a)/(self.G13*self.A),...
                (gam1*gam3*a)/(20*Ef*self.I) - (gam1*a)/(self.G13*self.A) - (6*self.h)/(self.b*self.E3),...
                -(2*self.h*gam1 + 6*self.h*a)/(self.b*self.E3) ];
            
            x3 = roots(polyCoefs);
            x3 = max(x3(x3 == real(x3)));
            
            % Determine Crack length
            beta1 = (x1^2 + 3*x1*x2 + 4*x1*x3 + 3*x2^2 + 8*x2*x3 + 5*x3^2) ...
                / (x1 + 2*x2 +2*x3);
            beta2 = 3 / x3 / (x1 + 2*x2 + 2*x3);
            beta3 = (x1^2*x2 + 3*x1*x3^2 + 3*x1*x2*x3 + 3*x2^2*x3 + 6*x2*x3^2 ...
                + 3*x3^3) / (x1 + 2*x2 +2*x3);
            beta4 = (x1 + 2*x2 + 3*x3) / x3 / (x1 + 2*x2*x3);
            
            % Determine current crack length using FPI at each point
            nPts = length(disp);
            crackLen = zeros(nPts,1);
            crackLen(1) = a;
            Gi = zeros(nPts,1);
            theta = zeros(nPts, 1);
            
            for iPt = 2:nPts
                comp = disp(iPt)/force(iPt);
                P = force(iPt);
                
                goalSeek = @(a) ( (2*a^3)/(3*Ef*self.I) + (beta1*a^2)/(2*Ef*self.I) + (2*a)/(self.A*self.G13) ...
                    + (4*self.h*beta2*a)/(self.b*self.E3) + (beta3*a)/(6*Ef*self.I) + (4*beta4*self.h)/(self.b*self.E3) - comp);
                
                a = fzero(goalSeek, crackLen(iPt-1));
                
                crackLen(iPt) = a;
                
                Gi(iPt) = (P^2*a^2)/(self.b*Ef*self.I) + (P^2*beta1*a)/(2*self.b*Ef*self.I) ...
                    + (P^2)/(self.b*self.A*self.G13) + (2*self.h*beta2*P^2)/(self.b^2*self.E3) ...
                    + (beta3*P^2)/(12*self.b*Ef*self.I);
                
                theta(iPt) = P/(12*Ef*self.I) * (beta3 + 3*beta1*a);
            end
            Gi(1) = Gi(2);
            aEquiv = crackLen;
            
        end
            
            
        end
        
end
        
            
            
                