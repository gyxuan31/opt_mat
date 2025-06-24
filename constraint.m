function [c, ceq] = constraint(e, predicted_len, total_UE, num_RB, gamma, num_setreq, rb_counts)
    E = reshape(e, predicted_len, total_UE, num_RB);
    c = zeros(1, predicted_len); % <=0

    % total RB limit
    for t = 1:predicted_len
        e_sum = 0;
        for n = 1:total_UE
            for k = 1:num_RB
                e_sum = e_sum + E(t,n,k);
            end
        end
        c(t) = e_sum - num_RB*gamma;
    end

    % requirement
    for i = 1:num_setreq
        e_req = 0;
        for k = 1:num_RB
            e_req = e_req + E(t,n,k);
        end
        c(t+i) = rb_counts(i) - e_req;
    end

    ceq = [];
end
