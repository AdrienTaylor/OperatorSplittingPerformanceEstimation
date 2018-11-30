function [WC] = strmonotonelipschitz(m,L,t,verbose)

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
Alips_value = L;    % assumed L>0 (possibly Inf)
Astrm_value = 0;    % assumed m<min(L,1/beta) (possibly 0)

BisNonZero  = 1;    % is B active ? [0/1]
BisSub      = 0;    % is B a subdifferential ? [0/1]
Bcoco_value = 0;    % assumed beta>=0 (possibly 0)
Blips_value = Inf;  % assumed L>0 (possibly Inf)
Bstrm_value = m;    % assumed m<min(L,1/beta) (possibly 0)

CisNonZero  = 0;    % is C active ? [0/1]
CisSub      = 0;    % is C a subdifferential ? [0/1]
Ccoco_value = 0;    % assumed beta>=0 (possibly 0)
Clips_value = Inf;  % assumed L>0 (possibly Inf)
Cstrm_value = 0;    % assumed m<min(L,1/beta) (possibly 0)

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
M_lips = @(L)([1    0;   0   -1/L^2]);
M_strm = @(m)([-m 1/2; 1/2      0]);
M_grad = @(m,L)([-m (1+m/L)/2; (1+m/L)/2 -1/L]);

% Algorithm' notations: (replace zi by xi and yi below)
%       z1 = J_{alpha*B} z
%       z2 = J_{alpha*A} (2 z1 - z - alpha C z1)
%       z+ = z - theta z1 + theta z2
% we also use the following:
%       DX  =   x - y
%       DB  = Bx1 - By1
%       DC  = Cx1 - Cy1
%       DA  = Ax2 - Ay2
%       DX1 =  x1 - y1  = DX - alpha * DB
%       DX2 =  x2 - y2  = 2 * DX1 - DX - alpha * DC - alpha * DA
%       DXp =  x+ - y+  = DX - theta * DX1 + theta * DX2


% if all operators are nonzero, the Gram matrix is defined as below:
%   P = [ DX | DB | DC  | DA]
%   G = P^T * P,
% whereas if B, C and/or A is zero, we simply discard it in P
% (i.e., we remove the corresponding column(s)).

dimG = 1 + AisNonZero + BisNonZero + CisNonZero;
DX   = zeros(1, dimG); DX(1,1) = 1;
DB   = zeros(1, dimG); DB(1,1+BisNonZero) = BisNonZero;
DC   = zeros(1, dimG); DC(1,1+BisNonZero+CisNonZero) = CisNonZero;
DA   = zeros(1, dimG); DA(1,1+BisNonZero+CisNonZero+AisNonZero) = AisNonZero;

DX1  = DX - alpha * DB;
DX2  = 2 * DX1 - DX - alpha * DC - alpha * DA;
DXp  = DX - theta * DX1 + theta * DX2;

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

S = - tau * (DX.'*DX) + (DXp.'*DXp); % init. dual matrix

if AisNonZero
    if AisSub
        S = S + lamA(4) * [DX2; DA].' * M_grad(Astrm,Acoco) * [DX2; DA];
    else
        S = S + lamA(1) * [DX2; DA].' * M_coco(Acoco) * [DX2; DA] ...
                            + lamA(2) * [DX2; DA].' * M_lips(Alips) * [DX2; DA] ...
                            + lamA(3) * [DX2; DA].' * M_strm(Astrm) * [DX2; DA];
    end
end
if BisNonZero
    if BisSub
        S = S + lamB(4) * [DX1; DB].' * M_grad(Bstrm,Bcoco) * [DX1; DB];
    else
        S = S + lamB(1) * [DX1; DB].' * M_coco(Bcoco) * [DX1; DB] ...
                            + lamB(2) * [DX1; DB].' * M_lips(Blips) * [DX1; DB] ...
                            + lamB(3) * [DX1; DB].' * M_strm(Bstrm) * [DX1; DB];
    end
end

if CisNonZero
    if CisSub
        S = S + lamC(4) * [DX1; DC].' * M_grad(Cstrm,Ccoco_value) * [DX1; DC];
    else
        S = S + lamC(1) * [DX1; DC].' * M_coco(Ccoco) * [DX1; DC] ...
                            + lamC(2) * [DX1; DC].' * M_lips(Clips) * [DX1; DC] ...
                            + lamC(3) * [DX1; DC].' * M_strm(Cstrm) * [DX1; DC];
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



































