function [WC, feas, mode] = getUpperbound(mu,L,theta)

C    = sqrt(((2*(theta-1)*mu+theta-2)^2+L^2*(theta-2*(mu+1))^2)/(L^2+1));
if theta*(theta+C)/(mu+1)^2/C * (C+mu*((2*(theta-1)*mu+theta-2)-L^2*(theta-2*(mu+1)))/(L^2+1)) >=0
    mode = 1;
    WC   = ((theta+C)/2/(mu+1))^2;
elseif L<=1 && mu >= (L^2+1)/(L-1)^2 && theta<=-(2*(mu+1)*(L+1)*(mu+(mu-1)*L^2-2*mu*L-1))/(mu+L*(L^2+L+1)+2*mu^2*(L-1)+mu*L*(1-(L-3)*L)+1)
    mode = 2;
    WC   = (1-theta*(L+mu)/(L+1)/(mu+1))^2;
else
    mode = 3;
    WC   = (2-theta)/4/mu/(L^2+1) * ...
        (theta*(1-2*mu+L^2)-2*mu*(L^2-1))*...
        (theta*(1+2*mu+L^2)-2*(mu+1)*(L^2+1))/...
        (theta*(1+2*mu-L^2)-2*(mu+1)*(1-L^2));
end

switch mode
    case 1
        lambdaA  = theta*(theta+C)/(mu+1);
        lambdaBL = ((2-theta)*theta*mu*L^2*(theta+C)/(mu+1)/(L^2+1)/C)/L^2;
        lambdaBm = theta*(theta+C)/(mu+1)^2/C * ...
            (C+mu*((2*(theta-1)*mu+theta-2)-L^2*(theta-2*(mu+1)))/(L^2+1));
    case 2
        lambdaA  = 2*theta*(1+L)/(1-L)*(1-theta*(mu+L)/(L+1)/(mu+1));
        lambdaBL = (theta*L*(mu-1)/(mu+1)*(1-theta*(mu+L)/(L+1)/(mu+1)))/L^2;
        lambdaBm = 0;
    case 3
        lambdaA  = theta*(theta-2*mu*(theta+L^2-1)/(L^2+1));
        lambdaBL = ((2-theta)*theta*L^2/(L^2+1) * ...
            (theta * (L^2+1)-2*mu*(theta+L^2-1))/...
            ((2-theta)*(1-L^2)-2*mu*(theta+L^2-1)))/L^2;
        lambdaBm = 0;
end

islambdaAPos  = (lambdaA >=-10^-8);
islambdaBLPos = (lambdaBL>=-10^-8);
islambdaBmPos = (lambdaBm>=-10^-8);
[~ ,psd] = getDualMatrix(mu,L,theta,WC,lambdaBm,lambdaBL,lambdaA);

feas = min([islambdaAPos islambdaBLPos islambdaBmPos psd]);


end
function [S,psd] = getDualMatrix(mu,L,theta,tau,lambdaBm,lambdaBL,lambdaA)

S = [tau+lambdaBL-1             lambdaA/2-theta                 theta-lambdaBL-lambdaBm/2;
    lambdaA/2-theta             -theta^2+lambdaA+lambdaA*mu     theta^2-lambdaA;
    theta-lambdaBL-lambdaBm/2   theta^2-lambdaA                 -lambdaBL*L^2-theta^2+lambdaBL+lambdaBm];

psd = (min(eig(S))>=-10^-7);
end 