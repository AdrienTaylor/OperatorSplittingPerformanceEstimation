%%
%Compare Upper and lower bounds of Section C

clear all;
clc;

tol    = 1e-8;
m_list = 1e-2:1e-2:2;
L_list = 1e-2:1e-2:2;
t_list = 1e-2:1e-2:2;
nb_comp = length(m_list)*length(L_list)*length(t_list);
cu_comp = 0;
nb_fail = 0;
nb_feas = 0;
UB_table = zeros(length(m_list), length(L_list), length(t_list));
LB_table = zeros(length(m_list), length(L_list), length(t_list));
fe_table = zeros(length(m_list), length(L_list), length(t_list));
Va_table = zeros(length(m_list), length(L_list), length(t_list));
mo_table = zeros(length(m_list), length(L_list), length(t_list));

toggle_every = 1e5;% toggle message every XXX verification;
for i = 1:length(m_list)
    for j = 1:length(L_list)
        for k = 1:length(t_list)
            mu    = m_list(i);
            L     = L_list(j);
            theta = t_list(k);
            
            [UB_table(i,j,k), fe_table(i,j,k), modeUB] = getUpperbound(mu,L,theta);
            mo_table(i,j,k)                 = modeUB;
            [LB_table(i,j,k)      , modeLB] = getLowerbound(mu,L,theta);
            Va_table(i,j,k)         = (abs(UB_table(i,j,k)-LB_table(i,j,k))<=tol);
            cu_comp                 = cu_comp + 1;
            if nb_fail == 0 && (abs(UB_table(i,j,k)-LB_table(i,j,k))>tol)
                first_L_failed = L; first_t_failed = theta; first_m_failed = mu;
            end
            nb_fail                 = nb_fail + (abs(UB_table(i,j,k)-LB_table(i,j,k))>tol);
            nb_feas                 = nb_feas + fe_table(i,j,k);
            if mod(cu_comp, toggle_every) == 0 
            fprintf('%d on a total of %d were done so far; SUMMARY: \n',cu_comp, nb_comp);
            fprintf('currently %d on %d trials with matching upper/lower bounds\n', cu_comp-nb_fail,cu_comp);
            fprintf('currently %d on %d upper bounds were feasible\n', nb_feas, cu_comp);
            end
        end
    end
end

%%
%Compare Bound with SDP

clear all;
clc;

verbose = 0; tol = 1e-8;
m_list = 1e-2:.1:1;
L_list = 1e-2:.1:1;
t_list = .1:.1:2;


nb_comp = length(m_list)*length(L_list)*length(t_list);
cu_comp = 0;
nb_fail = 0;
UB_table = zeros(length(m_list), length(L_list), length(t_list));
SL_table = zeros(length(m_list), length(L_list), length(t_list));
max_err = tol;

toggle_every = 1e1;% toggle message every XXX verification;
for i = 1:length(m_list)
    for j = 1:length(L_list)
        for k = 1:length(t_list)
            mu    = m_list(i);
            L     = L_list(j);
            theta = t_list(k);
            
            [UB_table(i,j,k), ~, ~] = getUpperbound(mu,L,theta);
            [LS_table(i,j,k)] = lipschitzstrmonotone(mu,L,theta,verbose);

            cu_comp                 = cu_comp + 1;
            if (abs(UB_table(i,j,k)-LS_table(i,j,k))>tol)
                nb_fail  = nb_fail+1;
                L_failed(nb_fail) = L;
                t_failed(nb_fail) = theta;
                m_failed(nb_fail) = mu;
            end
            max_err = max(max_err, abs(UB_table(i,j,k)-LS_table(i,j,k)));
            if mod(cu_comp, toggle_every) == 0 
            fprintf('%d on a total of %d were done so far; SUMMARY: \n',cu_comp, nb_comp);
            fprintf('currently %d out of %d trials with matching guarantees\n', cu_comp-nb_fail,cu_comp);
            fprintf('largest difference so far: %10.9f\n', max_err);
            end
        end
    end
end