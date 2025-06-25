function [c, ceq] = constraint(e, predicted_len, total_UE, num_RB, gamma, num_setreq, rb_counts, num_RU, RU_UE)
    E = reshape(e, predicted_len, total_UE, num_RB);
    c = zeros(1, predicted_len + predicted_len*num_RU*num_RB); % <=0


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
    
    % RB limit under every RU
    for t= 1:predicted_len
        for i = 1:num_RU
            for k = 1:num_RB
                count = 0;
                for u = RU_UE{i}
                    count = count +E(t,u,k);
                end
                c(predicted_len*t + num_RU*(i-1) + k) = count - 2;

            end
        end
    end

    % requirement
    % for i = 1:num_setreq
    %     e_req = 0;
    %     for k = 1:num_RB
    %         e_req = e_req + E(t,n,k);
    %     end
    %     c(t+i) = rb_counts(i) - e_req;
    % end

    ceq = [];
end
