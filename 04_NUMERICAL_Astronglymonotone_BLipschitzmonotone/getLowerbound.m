function [WC, ind] = getLowerbound(mu,L,theta)

LB1 = ((theta+sqrt(((2*(theta-1)*mu+theta-2)^2+L^2*(theta-2*(mu+1))^2)/(L^2+1)))/(mu+1)/2)^2;
LB2 = (1-theta*(L+mu)/(L+1)/(mu+1))^2;


numeratorK   = (1+L^2)*((-1+mu)*(-2+theta+2*(-1+theta)* mu)^2+L^4 *...
    (-1+mu)*(theta-2 *(1+mu))^2-2* L^2 *(4* (1+mu)* (1+mu^2)-...
    4*theta *(1+mu)* (1+mu^2)+theta^2 *(1+mu+2*mu^2)));

denominatorK = 2*L*((1 + L^2)^2 *(-2 + theta)^2 + (1 + L^2)^2 *...
    (-2 + theta)^2 * mu - 4 * (1 + L^2 * (-1 + theta)) *...
    (-1 + L^2 + theta) * mu^2 + 4* (-1 + L^2 + theta)^2 *mu^3);

K = numeratorK/denominatorK;

if  K>=0 && K<1
    
    JA = 1/(mu+1) * [1 0; 0 0];
    B = L * [K -sqrt(1-K^2); sqrt(1-K^2) K];
    JB = inversetwotimestwo(B+eye(2));
    T  = eye(2) - theta*JB + theta*JA*(2*JB-eye(2));
    TT = T*T.';
    LB3 = max(eig(TT));
    [WC, ind] = max([LB1 LB2 LB3]);
else
    [WC, ind] = max([LB1 LB2]);
end


end
function [inv] = inversetwotimestwo(A)
a = A(1,1); b=A(1,2); c=A(2,1); d=A(2,2);
inv = 1/(a*d-b*c) * [d -b; -c a];

end