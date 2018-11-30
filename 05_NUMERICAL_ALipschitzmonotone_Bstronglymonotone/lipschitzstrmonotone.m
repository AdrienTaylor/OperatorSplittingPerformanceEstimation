function [WC] = lipschitzstrmonotone(m,L,t,verbose)

%% Contraction factor for TOS.
% clear all; clc;

%% Solver/symbolic parameters

% set 'analysis' to (i)  symb to get the analytical LMI
%                   (ii) num to get the numerical solution minimizing the
%                   contraction factor
num      = 'numerical';
symb     = 'symbolical';
analysis = num;

% if numerical
% verbose  = 1;     % let the solver talk [0/1] ?
tol      = 1e-10; % accuracy for the solver

%% Algorithm' parameters (for numerics only)

% complete here if the analysis is to be performed numerically
alpha_value = 1;
theta_value = t;

%% Problem class parameters

% complete the options; the "_values" are not taken into account in the
% symbolical mode

AisNonZero  = 1;    % is A active ? [0/1]
AisSub      = 0;    % is A a subdifferential ? [0/1]
Acoco_value = 0;    % assumed beta>=0 (possibly 0)
Alips_value = Inf;  % assumed L>0 (possibly Inf)
Astrm_value = m;    % assumed m<min(L,1/beta) (possibly 0)

BisNonZero  = 1;    % is B active ? [0/1]
BisSub      = 0;    % is B a subdifferential ? [0/1]
Bcoco_value = 0;    % assumed beta>=0 (possibly 0)
Blips_value = L;    % assumed L>0 (possibly Inf)
Bstrm_value = 0;    % assumed m<min(L,1/beta) (possibly 0)

CisNonZero  = 0;    % is C active ? [0/1]
CisSub      = 0;    % is C a subdifferential ? [0/1]
Ccoco_value = 0;    % assumed beta>=0 (possibly 0)
Clips_value = Inf;  % assumed L>0 (possibly Inf)
Cstrm_value = 0;   % assumed m<min(L,1/beta) (possibly 0)

% EXAMPLES:
% vanilla gradient method: CisNonZero=1; BisNonZero=0; AisNonZero=0;
% vanilla prox- on A:      CisNonZero=0; BisNonZero=0; AisNonZero=1;
% vanilla prox- on B:      CisNonZero=0; BisNonZero=1; AisNonZero=0;
% FBS:                     CisNonZero=1; BisNonZero=0; AisNonZero=1;
% DRS:                     CisNonZero=0; BisNonZero=1; AisNonZero=1;

%% PEP settings (not to be modified)

switch analysis
    case num
        alpha = alpha_value;
        theta = theta_value;
        
        Acoco = Acoco_value;
        Alips = Alips_value;
        Astrm = Astrm_value;
        Bcoco = Bcoco_value;
        Blips = Blips_value;
        Bstrm = Bstrm_value;
        Ccoco = Ccoco_value;
        Clips = Clips_value;
        Cstrm = Cstrm_value;
        
    case symb
        alpha = sym('alpha');
        theta = sym('theta');
        
        Acoco = sym('Acoco');
        Alips = sym('ALips');
        Astrm = sym('Astrm');
        Bcoco = sym('Bcoco');
        Blips = sym('BLips');
        Bstrm = sym('Bstrm');
        Ccoco = sym('Ccoco');
        Clips = sym('CLips');
        Cstrm = sym('Cstrm');
end

% M matrices for the characteristics; corresponding to constraints of the
% form (DeltaX DeltaT) M (DeltaX DeltaT)^T >= 0

M_coco = @(beta)([0  1/2; 1/2   -beta]);
M_lips = @(L)([L^2    0;   0   -1]);
M_strm = @(m)([-m 1/2; 1/2      0]);
M_grad = @(m,L)([-m (1+m/L)/2; (1+m/L)/2 -1/L]);

% Algorithm' notations: (replace zi by xi and yi below)
%       zB = J_{alpha*B} z
%       zC = alpha C zB
%       zA = J_{alpha*A} (2 zB - z - zC)
%       z+ = z - theta ( zB - zA)
% we also use the following:
%       DA  = (2 zB - z - zC - zA)/alpha
%       DB  = (z  - zB) / alpha
%       DC  = zC / alpha

% if all operators are nonzero, the Gram matrix is defined as below:
%   P = [ z | zA | zB  | zC]
%   G = P^T * P,
% whereas if B, C and/or A is zero, we simply discard it in P
% (i.e., we remove the corresponding column(s)).

dimG = 1 + AisNonZero + BisNonZero + CisNonZero;
z    = zeros(1, dimG); z(1,1) = 1;
zC   = zeros(1, dimG); zC(1,1+AisNonZero+BisNonZero+CisNonZero) = CisNonZero;
if BisNonZero
    zB = zeros(1, dimG); zB(1,2+AisNonZero) = 1;
else
    zB = z;
end
if AisNonZero
    zA = zeros(1, dimG); zA(1,2) = 1;
else
    zA = 2*zB - z - zC;
end
zp   = z - theta * ( zB - zA);

DA   = (2 * zB - z - zC - zA) / alpha;
DB   = (z - zB)  / alpha;
DC   = zC / alpha;


%% PEP (dual form)

% VARIABLES (similar correspondences for B and C):
% tau     : contraction factor
% lamA(1) : cocoercivity of A (0 if not used)
% lamA(2) : Lipschitzness of A (0 if not used)
% lamA(3) : strong monotonicity of A (0 if not used)
% lamA(4) : smoothness and strong convexity of A (0 if not used)
switch analysis
    case num
        tau  = sdpvar(1);
        lamA = sdpvar(4,1);
        lamB = sdpvar(4,1);
        lamC = sdpvar(4,1);
    case symb
        tau  = sym('tau');
        lamA = sym('LA%d', [4 1]);
        lamB = sym('LB%d', [4 1]);
        lamC = sym('LC%d', [4 1]);
end

S = - tau * (z.'*z) + (zp.'*zp); % init. dual matrix

if AisNonZero
    if AisSub
        S = S + lamA(4) * [zA; DA].' * M_grad(Astrm,1/Acoco) * [zA; DA];
    else
        S = S + lamA(1) * [zA; DA].' * M_coco(Acoco) * [zA; DA] ...
            + lamA(3) * [zA; DA].' * M_strm(Astrm) * [zA; DA];
        if Alips ~= Inf
            S = S  + lamA(2) * [zA; DA].' * M_lips(Alips) * [zA; DA];
        end
        
    end
end
if BisNonZero
    if BisSub
        S = S + lamB(4) * [zB; DB].' * M_grad(Bstrm,1/Bcoco) * [zB; DB];
    else
        S = S + lamB(1) * [zB; DB].' * M_coco(Bcoco) * [zB; DB] ...
            + lamB(3) * [zB; DB].' * M_strm(Bstrm) * [zB; DB];
        if Blips ~= Inf
            S = S  + lamB(2) * [zB; DB].' * M_lips(Blips) * [zB; DB];
        end
    end
end

if CisNonZero
    if CisSub
        S = S + lamC(4) * [zB; DC].' * M_grad(Cstrm,1/Ccoco_value) * [zB; DC];
    else
        S = S + lamC(1) * [zB; DC].' * M_coco(Ccoco) * [zB; DC] ...
            + lamC(3) * [zB; DC].' * M_strm(Cstrm) * [zB; DC];
        if Clips ~= Inf
            S = S  + lamC(2) * [zB; DC].' * M_lips(Clips) * [zB; DC];
        end
    end
end

switch analysis
    case num
        cons = (S <= 0);
        cons = cons + (lamA >= 0);
        cons = cons + (lamB >= 0);
        cons = cons + (lamC >= 0);
        
        obj = tau;
        
        solver_opt = sdpsettings('solver','mosek','verbose',verbose,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',tol);
        solverDetails=optimize(cons,obj,solver_opt);
    case symb
        fprintf('This matrix should be negative semidefinite:\n');
        simplify(S)
        % Note that there are "trivial" simplifications to be done in S,
        % as for example removing the unnecessary parameters.
end
WC = double(obj);


































