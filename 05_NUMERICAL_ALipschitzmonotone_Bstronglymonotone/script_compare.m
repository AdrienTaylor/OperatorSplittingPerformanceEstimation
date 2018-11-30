clear all;
clc;

verbose = 0; tol = 1e-8;
m_list = 1e-2:.1:1;
L_list = 1e-2:.1:1;
t_list = .1:.1:2;



nb_comp = length(m_list)*length(L_list)*length(t_list);
cu_comp = 0;
nb_fail = 0;
LS_table = zeros(length(m_list), length(L_list), length(t_list));
SL_table = zeros(length(m_list), length(L_list), length(t_list));
max_err = tol;

toggle_every = 1e1;% toggle message every XXX verification;
for i = 1:length(m_list)
    for j = 1:length(L_list)
        for k = 1:length(t_list)
            mu    = m_list(i);
            L     = L_list(j);
            theta = t_list(k);
            [SL_table(i,j,k)] = strmonotonelipschitz(mu,L,theta,verbose);
            [LS_table(i,j,k)] = lipschitzstrmonotone(mu,L,theta,verbose);

            cu_comp                 = cu_comp + 1;
            if (abs(SL_table(i,j,k)-LS_table(i,j,k))>tol)
                nb_fail  = nb_fail+1;
                L_failed(nb_fail) = L;
                t_failed(nb_fail) = theta;
                m_failed(nb_fail) = mu;
            end
            max_err = max(max_err, abs(SL_table(i,j,k)-LS_table(i,j,k)));
            if mod(cu_comp, toggle_every) == 0 
            fprintf('%d on a total of %d were done so far; SUMMARY: \n',cu_comp, nb_comp);
            fprintf('currently %d out of %d trials with matching guarantees\n', cu_comp-nb_fail,cu_comp);
            fprintf('largest difference so far: %10.9f\n', max_err);
            end
        end
    end
end
