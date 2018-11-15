%%
%Default values
AisNonZero = 1;
Astrm_value = 0;
Acoco_value = 0;
Alips_value = inf;

BisNonZero = 1;
Bstrm_value = 0;
Bcoco_value = 0;
Blips_value = inf;

CisNonZero = 1;
Cstrm_value = 0;
Ccoco_value = 1;
Clips_value = inf;
%%
%Figure 1
AisNonZero = 1;
Astrm_value = 1;
Acoco_value = 0;
Alips_value = inf;

BisNonZero = 1;
Bstrm_value = 0;
Bcoco_value = 0.01;
Blips_value = 5;

CisNonZero = 1;
Cstrm_value = 0;
Ccoco_value = 9;
Clips_value = inf;
%%
rho2 = @(alpha) OSPEP_theta(alpha, AisNonZero, Acoco_value, Alips_value, Astrm_value, BisNonZero, Bcoco_value, Blips_value, Bstrm_value, CisNonZero, Ccoco_value, Clips_value, Cstrm_value);

N = 100;
aa = 10.^linspace(-3,1,N);
vv = zeros(N,1);

for ii=1:N
    vv(ii) = rho2(aa(ii));
end

semilogx(aa,vv,'k', 'LineWidth',2);

xlabel('\alpha');
ylabel('\rho^2');

ylim([0.7,1])
xlim([1e-3,1e1])
yticks([0.7 0.8 0.9 1])
set(gcf,'Position',[100 100 600 200])

myprint('plot_fig1')

%Find minimizer with fminunc
fminunc(@(alpha) rho2(alpha),1)



%%
%Figure 2a
AisNonZero = 1;
Astrm_value = 1;
Acoco_value = .07;
Alips_value = inf;

BisNonZero = 1;
Bstrm_value = 4;
Bcoco_value = .02;
Blips_value = inf;

CisNonZero = 1;
Cstrm_value = 0;
Ccoco_value = 9;
Clips_value = inf;

%Plot range
yrange = [0.2,1];
xrange = [1e-3,1e1];
%%
%Figure 2b
AisNonZero = 1;
Astrm_value = 1;
Acoco_value = 0;
Alips_value = inf;

BisNonZero = 1;
Bstrm_value = 0;
Bcoco_value = .03;
Blips_value = inf;

CisNonZero = 0;
Cstrm_value = [];
Ccoco_value = [];
Clips_value = [];

%Plot range
yrange = [0.6,1];
xrange = [1e-3,1e1];
%%
%Figure 2c
AisNonZero = 1;
Astrm_value = 0;
Acoco_value = .03;
Alips_value = inf;

BisNonZero = 1;
Bstrm_value = 1;
Bcoco_value = 0;
Blips_value = inf;

CisNonZero = 0;
Cstrm_value = [];
Ccoco_value = [];
Clips_value = [];

%Plot range
yrange = [0.6,1];
xrange = [1e-3,1e1];
%%
%Figure 2d
AisNonZero = 1;
Astrm_value = 1;
Acoco_value = 0;
Alips_value = inf;

BisNonZero = 1;
Bstrm_value = 0;
Bcoco_value = 0;
Blips_value = 4;

CisNonZero = 0;
Cstrm_value = [];
Ccoco_value = [];
Clips_value = [];

%Plot range
yrange = [0.7,1];
xrange = [1e-3,1e1];
%%
%Figure 2e
AisNonZero = 1;
Astrm_value = 0;
Acoco_value = 0;
Alips_value = 4;

BisNonZero = 1;
Bstrm_value = 1;
Bcoco_value = 0;
Blips_value = inf;

CisNonZero = 0;
Cstrm_value = [];
Ccoco_value = [];
Clips_value = [];

%Plot range
yrange = [0.7,1];
xrange = [1e-3,1e1];
%%
%Figure 2f
AisNonZero = 0;
Astrm_value = [];
Acoco_value = [];
Alips_value = [];


BisNonZero = 1;
Bstrm_value = 1;
Bcoco_value = 0.1;
Blips_value = inf;

CisNonZero = 1;
Cstrm_value = 0;
Ccoco_value = 0;
Clips_value = 8;

%Plot range
yrange = [0.98,1];
xrange = [1e-4,1e-1];
%%
%Figure 2g
AisNonZero = 1;
Astrm_value = 1;
Acoco_value = .07;
Alips_value = 7;

BisNonZero = 1;
Bstrm_value = .3;
Bcoco_value = .02;
Blips_value = 2;

CisNonZero = 1;
Cstrm_value = 0.01;
Ccoco_value = 9;
Clips_value = .05;

%Plot range
yrange = [0.4,1];
xrange = [1e-3,1e2];
%%
%Figure 2h
AisNonZero = 0;
Astrm_value = [];
Acoco_value = [];
Alips_value = [];

BisNonZero = 1;
Bstrm_value = 1;
Bcoco_value = 0;
Blips_value = inf;

CisNonZero = 1;
Cstrm_value = 0;
Ccoco_value = .1;
Clips_value = inf;

%Plot range
yrange = [0.6,1];
xrange = [1e-3,1e1];
%%
%Figure 2i
AisNonZero = 1;
Astrm_value = 1;
Acoco_value = 0;
Alips_value = inf;

BisNonZero = 0;
Bstrm_value = [];
Bcoco_value = [];
Blips_value = [];

CisNonZero = 1;
Cstrm_value = 0;
Ccoco_value = .1;
Clips_value = inf;

%Plot range
yrange = [0.6,1];
xrange = [1e-3,1e1];
%%
rho2 = @(alpha) OSPEP_theta(alpha, AisNonZero, Acoco_value, Alips_value, Astrm_value, BisNonZero, Bcoco_value, Blips_value, Bstrm_value, CisNonZero, Ccoco_value, Clips_value, Cstrm_value);

N = 100;
aa = 10.^linspace(log10(xrange(1)),log10(xrange(2)),N);
vv = zeros(N,1);

for ii=1:N
    vv(ii) = rho2(aa(ii));
end

close all;
semilogx(aa,vv,'k', 'LineWidth',2);

ylim(yrange)
xlim(xrange)
set(gcf,'Position',[100 100 300 200])

myprint('plot_fig2g')

fminunc(@(alpha) rho2(alpha),1e-2)
