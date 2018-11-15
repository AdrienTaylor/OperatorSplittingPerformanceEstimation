function output = OSPEP_theta(alpha, AisNonZero, Acoco_value, Alips_value, Astrm_value, BisNonZero, Bcoco_value, Blips_value, Bstrm_value, CisNonZero, Ccoco_value, Clips_value, Cstrm_value) 
% Problem class parameters

% complete the options; the "_values" are not taken into account in the
% symbolical mode

%AisNonZero   is A active ? [0/1]
%Acoco_value  assumed beta>0 (possibly 0)
%Alips_value  assumed L>0 (possibly Inf)
%Astrm_value  assumed m<L (possibly 0)

%BisNonZero   is B active ? [0/1]
%Bcoco_value  assumed beta>0 (possibly 0)
%Blips_value  assumed L>0 (possibly Inf)
%Bstrm_value  assumed m<L (possibly 0)

%CisNonZero   is C active ? [0/1]
%Ccoco_value  assumed beta>0 (possibly 0)
%Clips_value  assumed L>0 (possibly Inf)
%Cstrm_value  assumed m<L (possibly 0)

% EXAMPLES:
% vanilla gradient method: CisNonZero=1; BisNonZero=0; AisNonZero=0;
% vanilla prox- on A:      CisNonZero=0; BisNonZero=0; AisNonZero=1;
% vanilla prox- on B:      CisNonZero=0; BisNonZero=1; AisNonZero=0;
% FBS:                     CisNonZero=1; BisNonZero=0; AisNonZero=1;
% DRS:                     CisNonZero=0; BisNonZero=1; AisNonZero=1;

verbose  = 1;     % let the solver talk [0/1] ?
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
M_lips = @(L)([1    0;   0   -1/L^2]);
M_strm = @(m)([-m 1/2; 1/2      0]);

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

% [            1,    -alpha*theta,    -alpha*theta,    -alpha*theta]
% [ -alpha*theta, alpha^2*theta^2, alpha^2*theta^2, alpha^2*theta^2]
% [ -alpha*theta, alpha^2*theta^2, alpha^2*theta^2, alpha^2*theta^2]
% [ -alpha*theta, alpha^2*theta^2, alpha^2*theta^2, alpha^2*theta^2]

sss = sdpvar(1);
mmm = [            1,    -alpha*theta,    -alpha*theta,    -alpha*theta;
    -alpha*theta,     alpha^2*sss,     alpha^2*sss,     alpha^2*sss;
    -alpha*theta,     alpha^2*sss,     alpha^2*sss,     alpha^2*sss;
    -alpha*theta,     alpha^2*sss,     alpha^2*sss,     alpha^2*sss];
mmm = mmm(1:dimG,1:dimG);
S = - tau * (DX.'*DX) + mmm; % init. dual matrix


if AisNonZero
    S = S + lamA(1) * [DX2; DA].' * M_coco(Acoco) * [DX2; DA] ...
        + lamA(2) * [DX2; DA].' * M_lips(Alips) * [DX2; DA] ...
        + lamA(3) * [DX2; DA].' * M_strm(Astrm) * [DX2; DA];
end
if BisNonZero
    S = S + lamB(1) * [DX1; DB].' * M_coco(Bcoco) * [DX1; DB] ...
        + lamB(2) * [DX1; DB].' * M_lips(Blips) * [DX1; DB] ...
        + lamB(3) * [DX1; DB].' * M_strm(Bstrm) * [DX1; DB];
end

if CisNonZero
    S = S + lamC(1) * [DX1; DC].' * M_coco(Ccoco) * [DX1; DC] ...
        + lamC(2) * [DX1; DC].' * M_lips(Clips) * [DX1; DC] ...
        + lamC(3) * [DX1; DC].' * M_strm(Cstrm) * [DX1; DC];
end

cons = (S <= 0);
cons = cons + (lamA >= 0);
cons = cons + (lamB >= 0);
cons = cons + (lamC >= 0);
cons = cons + (theta^2 <= sss);
cons = cons + (theta >= 0);
cons = cons + (theta <= 2);
obj = tau;

solver_opt = sdpsettings('solver','mosek','verbose',verbose,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',tol);
solverDetails=optimize(cons,obj,solver_opt);



output = value(tau);
disp(['\theta value is ' , num2str(value(theta))]);





end




















