function [WC, ind] = getLowerbound(mu,beta,theta)

LB1 = (1-theta*beta/(1+beta))^2;
LB2 = (1-theta*(1+mu*beta)/(1+mu)/(1+beta))^2;
LB3 = (1-theta)^2;
LB4 = (1-theta*mu/(1+mu))^2;

numeratorK   = beta^2*(mu-1)*(2*(theta-1)*mu+theta-2)^2+...
    2*beta*(theta-2)*mu^2*(2*(theta-1)*mu+theta-2)+...
    (theta-2)^2*(mu+1)*mu^2;
denominatorK = beta^2*(theta-2)*mu*(2*beta*(theta-2)-theta-2)+...
    beta^2*(2*beta+1)*(theta-2)^2+...
    (2*beta-1)*mu^3*(2*beta*(theta-1)+theta-2)^2-...
    (2*beta-1)*mu^2*(2*beta-theta+2)*(2*beta*(theta-1)+theta-2);

K = numeratorK/denominatorK;

if K>0 && K<=1/beta^2
    term1_a = 2*theta*mu+theta-2*beta*theta*K*mu-theta*K+2*K*mu+2*K-2*mu-2;
    term2_a = 4*(theta-2)^2*(mu+1)^2*(K-beta^2*K^2)+((theta-2)*(K-1)-2*mu*(theta-beta*theta*K+K-1))^2;
    denominatora = 2*(theta-2)*sqrt(K-beta^2*K^2);
    a = (term1_a + sqrt(term2_a))/denominatora;
    
    A = [mu -a; a mu];
    B =[beta*K -sqrt(K-K^2*beta^2); sqrt(K-K^2*beta^2) beta*K];
    JA = inv(A+eye(2));
    JB = inv(B+eye(2));
    T  = eye(2) - theta*JB + theta*JA*(2*JB-eye(2));
    TT = T*T.';
    LB5 = max(eig(TT));
    [WC, ind] = max([LB1 LB2 LB3 LB4 LB5]);
else
    [WC, ind] = max([LB1 LB2 LB3 LB4]);
end


end