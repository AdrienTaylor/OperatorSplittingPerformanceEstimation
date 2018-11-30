function output = OSPEP_theta(alpha, AisNonZero, Acoco_value, Alips_value, Astrm_value, BisNonZero, Bcoco_value, Blips_value, Bstrm_value, CisNonZero, Ccoco_value, Clips_value, Cstrm_value, verbose)
% Problem class parameters

% complete the options; the "_values" are not taken into account in the
% symbolical mode

%AisNonZero   is A active ? [0/1]
%Acoco_value  assumed beta>=0 (possibly 0)
%Alips_value  assumed L>0 (possibly Inf)
%Astrm_value  assumed m<min(L,1/beta) (possibly 0)

%BisNonZero   is B active ? [0/1]
%Bcoco_value  assumed beta>=0 (possibly 0)
%Blips_value  assumed L>0 (possibly Inf)
%Bstrm_value  assumed m<min(L,1/beta) (possibly 0)

%CisNonZero   is C active ? [0/1]
%Ccoco_value  assumed beta>=0 (possibly 0)
%Clips_value  assumed L>0 (possibly Inf)
%Cstrm_value  assumed m<min(L,1/beta) (possibly 0)

% EXAMPLES:
% vanilla gradient method: CisNonZero=1; BisNonZero=0; AisNonZero=0;
% vanilla prox- on A:      CisNonZero=0; BisNonZero=0; AisNonZero=1;
% vanilla prox- on B:      CisNonZero=0; BisNonZero=1; AisNonZero=0;
% FBS:                     CisNonZero=1; BisNonZero=0; AisNonZero=1;
% DRS:                     CisNonZero=0; BisNonZero=1; AisNonZero=1;

% verbose  = 1;     % let the solver talk [0/1] ?
tol      = 1e-10; % accuracy for the solver




% PEP settings (not to be modified)

theta = sdpvar(1);

Acoco = Acoco_value;
Alips = Alips_value;
Astrm = Astrm_value;
Bcoco = Bcoco_value;
Blips = Blips_value;
Bstrm = Bstrm_value;
Ccoco = Ccoco_value;
Clips = Clips_value;
Cstrm = Cstrm_value;

% M matrices for the characteristics; corresponding to constraints of the
% form (DeltaX DeltaT) M (DeltaX DeltaT)^T >= 0

M_coco = @(beta)([0  1/2; 1/2   -beta]);
M_lips = @(L)([L^2    0;   0   -1]);
M_strm = @(m)([-m 1/2; 1/2      0]);

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

% PEP (dual form)

% VARIABLES (similar correspondences for B and C):
% tau     : contraction factor
% lamA(1) : cocoercivity of A (0 if not used)
% lamA(2) : Lipschitzness of A (0 if not used)
% lamA(3) : strong monotonicity of A (0 if not used)
% lamA(4) : smoothness and strong convexity of A (0 if not used)


tau  = sdpvar(1);
lamA = sdpvar(4,1);
lamB = sdpvar(4,1);
lamC = sdpvar(4,1);

S = - tau * (z.'*z); % init. dual matrix


if AisNonZero
    S = S + lamA(1) * [zA; DA].' * M_coco(Acoco) * [zA; DA] ...
        + lamA(3) * [zA; DA].' * M_strm(Astrm) * [zA; DA];
    if Alips ~= Inf
        S = S  + lamA(2) * [zA; DA].' * M_lips(Alips) * [zA; DA];
    end
end
if BisNonZero
    S = S + lamB(1) * [zB; DB].' * M_coco(Bcoco) * [zB; DB] ...
        + lamB(3) * [zB; DB].' * M_strm(Bstrm) * [zB; DB];
    if Blips ~= Inf
        S = S  + lamB(2) * [zB; DB].' * M_lips(Blips) * [zB; DB];
    end
end

if CisNonZero
    S = S + lamC(1) * [zB; DC].' * M_coco(Ccoco) * [zB; DC] ...
        + lamC(3) * [zB; DC].' * M_strm(Cstrm) * [zB; DC];
    if Clips ~= Inf
        S = S  + lamC(2) * [zB; DC].' * M_lips(Clips) * [zB; DC];
    end
end
Saug = -[-S zp.'; zp 1];% augmented S (Schur complement trick)
cons = (Saug <= 0);
cons = cons + (lamA >= 0);
cons = cons + (lamB >= 0);
cons = cons + (lamC >= 0);
if ~AisNonZero || ~BisNonZero
    cons = cons + (theta >= 0);
    cons = cons + (theta <= 2);
end
obj = tau;

solver_opt = sdpsettings('solver','mosek','verbose',verbose,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',tol);
optimize(cons,obj,solver_opt);



output = value(tau);
disp(['\theta value is ' , num2str(value(theta))]);





end




















