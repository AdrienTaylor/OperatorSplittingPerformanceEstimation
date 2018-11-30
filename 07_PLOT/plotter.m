%% Which case should we reproduce ? (type the name of the figure among:
% - fig1
% - fig2a, fig2b, fig2c, fig2d, fig2e, fig2f, fig2g, fig2h, fig2i
clc; clear all;
namefig    = 'fig2i';
optimalpha = 0; % do we optimize alpha with fminfun? [0/1]
verbose    = 0; % verbose solver ? [0/1]
switch namefig
    case 'fig1'%Figure 1
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
        yrange = [0.7,1];
        xrange = [1e-3,1e1];
    case 'fig2a'%Figure 2a
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
    case 'fig2b'%Figure 2b
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
    case 'fig2c'
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
    case 'fig2d' %Figure 2d
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
    case 'fig2e'%Figure 2e
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
    case 'fig2f'%Figure 2f
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
    case 'fig2g'%Figure 2g
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
    case 'fig2h'%Figure 2h
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
    case 'fig2i'%Figure 2i
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
end

%% Reproduce the figure and compute optimal alpha
rho2 = @(alpha) OSPEP_theta(alpha, AisNonZero, Acoco_value, Alips_value, Astrm_value, BisNonZero, Bcoco_value, Blips_value, Bstrm_value, CisNonZero, Ccoco_value, Clips_value, Cstrm_value, verbose);

N = 100;
aa = 10.^linspace(log10(xrange(1)),log10(xrange(2)),N);
vv = zeros(N,1);

for ii=1:N
    vv(ii) = rho2(aa(ii));
    if ~verbose
        fprintf('%d SDPs solved on a total of %d\n',ii,N);
    end
end

close all;
semilogx(aa,vv,'k', 'LineWidth',2);

ylim(yrange)
xlim(xrange)
set(gcf,'Position',[100 100 300 200])

if optimalpha
    fminunc(@(alpha) rho2(alpha),1e-2)
end