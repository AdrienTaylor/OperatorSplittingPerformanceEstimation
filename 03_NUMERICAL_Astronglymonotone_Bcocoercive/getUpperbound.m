function [WC, feas, mode] = getUpperbound(mu,beta,theta)

if mu*beta-mu+beta < 0 && theta<= 2*((beta+1)*(mu-beta-mu*beta))/(mu+mu*beta-beta-beta^2-2*mu*beta^ 2)
    mode = 1;
    WC   = (1-theta*beta/(beta+1))^2;
elseif mu*beta-mu-beta > 0 && theta <=2*(mu^2+beta^2+mu*beta+mu+beta-mu^2*beta^2)/(mu^2+beta^2+mu^2*beta+mu*beta^2+mu+beta-2*mu^2*beta^2)
    mode = 2;
    WC   = (1-theta*(1+mu*beta)/(mu+1)/(beta+1))^2;
elseif theta>=2*(mu*beta+mu+beta)/(2*mu*beta+mu+beta)
    mode = 3;
    WC   = (1-theta)^2;
elseif mu*beta+mu-beta<0 && theta<=2*(mu+1)*(beta-mu-mu*beta)/(beta+mu*beta-mu-mu^2-2*mu^2*beta)
    mode = 4;
    WC   = (1-theta*mu/(mu+1))^2;
else
    mode = 5;
    WC   = (2-theta)/4*((2-theta)*mu*(beta+1)-theta*beta*(mu-1))*((2-theta)*beta*(mu+1)-theta*mu*(beta-1))/((2-theta)*mu*beta*(mu+1)*(beta+1)-theta*mu^2*beta^2);
end

switch mode
    case 1
        lambdaA = 2*theta*(1+beta)/(1-beta)*(1-theta*beta/(beta+1));
        lambdaB = 2*theta*(1-theta*beta/(beta+1));
    case 2
        lambdaA = 2*theta*(beta+1)/(beta-1)*(1-theta*(1+mu*beta)/(mu+1)/beta+1);
        lambdaB = 2*theta*(mu-1)/(mu+1)*(1-theta*(1+mu*beta)/(mu+1)/(beta+1));
    case 3
        lambdaA = 2*theta*(theta-1);
        lambdaB = 2*theta*(theta-1);
    case 4
        lambdaA = 2*theta*(1-theta*mu/(mu+1));
        lambdaB = 2*theta*(1-mu)/(1+mu)*(1-theta*mu/(mu+1));
    case 5
        lambdaA = theta*((2-theta)*mu*(beta+1)+theta*beta*(1-mu))/beta;
        lambdaB = theta*(2-theta)/beta*((2-theta)*mu*(beta+1)+theta*beta*(1-mu))/(2*mu*beta*(1-theta)+(2-theta)*(mu+beta+1));
end

islambdaAPos = (lambdaA>=-10^-8);
islambdaBPos = (lambdaB>=-10^-8);
[~ ,psd] = getDualMatrix(mu,beta,theta,WC,lambdaB,lambdaA);

feas = min([islambdaAPos islambdaBPos psd]);


end
function [S,psd] = getDualMatrix(mu,beta,theta,tau,lambdaB,lambdaA)

S = [tau + beta*lambdaB-1    -theta+lambdaA/2        theta-(1/2+beta)*lambdaB;
    -theta+lambdaA/2         -theta^2+(1+mu)*lambdaA theta^2-lambdaA;
    theta-(1/2+beta)*lambdaB theta^2-lambdaA         -theta^2+(1+beta)*lambdaB];

psd = (min(eig(S))>=-10^-7);
end