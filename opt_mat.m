rng(1);

T = load('parameter.mat').T;
num_RU = load('parameter.mat').num_RU;
UERU = load('parameter.mat').UERU; % UE under every RU
total_UE = load('parameter.mat').total_UE;
num_RB = load('parameter.mat').num_RB;
num_ref = load('parameter.mat').num_ref;
gamma = load('parameter.mat').gamma;
num_setreq = load('parameter.mat').num_setreq;
B = double(load('parameter.mat').B);
P = load('parameter.mat').P;
sigmsqr = load('parameter.mat').sigmsqr;
eta = double(load('parameter.mat').eta);
predicted_len = load('parameter.mat').predicted_len;
rayleigh_gain = load('parameter.mat').rayleigh_gain;

user_RU_norm = zeros(1, total_UE); % RU index for every user
user_RU_op = zeros(1, total_UE);

load_dis = load('parameter.mat');
distance = load_dis.distance_true;
distance = reshape(distance, T, total_UE, num_RU);
load_pred = load('parameter.mat');
prediction = load_pred.prediction;
prediction = reshape(prediction, T-num_ref, predicted_len, total_UE, num_RU);


% ----- RANDOM - randomly generate e -----
for t = 1:T
    for i = 1:total_UE
        temp = zeros(1, num_RU);
        for j = 1:num_RU
            temp(j) = distance(t,i,j);
        end
        [~, user_RU_norm(i)] = min(temp);
    end
end

rb_counts = randi([0, num_RB/num_RU], 1, total_UE); % initial allocation
e_random = zeros(total_UE, num_RB);
record_random = [];
for i = 1:total_UE
    count = rb_counts(i);
    if count > 0
        selected_rbs = randperm(num_RB, count); 
        e_random(i, selected_rbs) = 1;
    end
end


for t = num_ref+1:T
    data_rate_random = zeros(1, total_UE);
    for n = 1:total_UE
        for k = 1:num_RB
            if e_random(n, k) == 1
                signal = P * distance(t, n, user_RU_norm(n)).^(-eta) * rayleigh_gain(n, k);
                interference = 0;
                for others = 1:total_UE
                    for i = 1:num_RU
                        if others ~= n && e_random(others, k) == 1 && user_RU_norm(others) ~= user_RU_norm(n)
                            interference = interference + ...
                                P * distance(t, n, user_RU_norm(i))^(-eta) * rayleigh_gain(n, k);
                        end
                    end
                end
                SINR = signal / (interference + sigmsqr);
                data_rate_random(n) = data_rate_random(n) + B * log(1 + SINR);
            end
        end
    end
    record_random = [record_random, sum(log(1+data_rate_random))];
end

% ----- STATIC -----


%  ----- OP -----
nvars = double(predicted_len * total_UE * num_RB);

lb = zeros(1, nvars);
ub = ones(1, nvars);

record_op = [];

for t = 1:T-num_ref
    disp(t);
    data_rate_op = zeros(1, total_UE);

    % OP
    pre_distance = reshape(prediction(t,:,:,:), predicted_len, total_UE, num_RU);
    for i = 1:predicted_len
        for u = 1:total_UE
            temp = zeros(1, num_RU);
            for j = 1:num_RU
                temp(j) = pre_distance(i,u,j);
            end
            [~, user_RU_op(u)] = min(temp);
        end
    end

    RU_UE_op = cell(1, num_RU); % save the index of UE under every RU
    for r = 1:num_RU
        idx = find(user_RU_op == r);
        RU_UE_op{r} = idx; % UE index of every RU
    end

    popSize = 50;
    initPop = zeros(popSize, nvars);
    initPop(popSize/2+1:end, :) = 1;

    objfun = @(e) compute_total_rate(round(e), predicted_len, total_UE, num_RB, pre_distance, rayleigh_gain, P, sigmsqr, eta, B, user_RU_op);
    confun = @(e) constraint(e, predicted_len, total_UE, num_RB, gamma, num_setreq, rb_counts, num_RU, RU_UE_op);
    options = optimoptions('ga', 'PopulationSize', popSize, 'MaxGenerations', 50, 'Display', 'iter',  'ConstraintTolerance', 1e-6,...
        'InitialPopulationMatrix', initPop);
    [e_opt, fval] = ga(objfun, nvars, [], [], [], [], lb, ub, confun, options);

    % fmin = @(e) compute_total_rate(e, predicted_len, total_UE, num_RB, pre_distance, rayleigh_gain, P, sigmsqr, eta, B, T, user_RU, num_RU);
    % options = optimoptions('fmincon', 'MaxIterations', 100, 'Display', 'iter-detailed', ...
    %     'MaxFunctionEvaluations', 1e5,'Algorithm','interior-point'); % , 'SpecifyObjectiveGradient',true
    % [e_opt, fval] = fmincon(fmin, repmat(e_norm, predicted_len, 1, 1), [], [], [], [], lb, ub, @constraints, options); % repmat(e_norm, predicted_len, 1, 1)
    e_opt;
    % for i = 1: nvars
    %     if e_opt(i) >= 0.5
    %         e_opt(i) = 1;
    %     else
    %         e_opt(i) = 0;
    %     end
    % end
    e_opt = reshape(e_opt, predicted_len, total_UE, num_RB);
    
    % check RB cannot allocated to 2 UEs under one RU
    RU_UE_norm = cell(1, num_RU); % save the index of UE under every RU
    for r = 1:num_RU
        idx = find(user_RU_norm == r);
        RU_UE_norm{r} = idx; % UE index of every RU
    end

    for i = 1:num_RU % set the repeat UE = 0
        for k = 1:num_RB
            ue_list = RU_UE_norm{i}; % UE under RU(i)
            allocated_UE = ue_list(e_opt(t, ue_list, k) > 0);
    
            if numel(allocated_UE) > 1
                dist_list = zeros(1, numel(allocated_UE));
                for idx = 1:numel(allocated_UE)
                    u = allocated_UE(idx);
                    dist_list(idx) = distance(t, u, i);
                end
                [~, min_idx] = min(dist_list);
                keep_ue = allocated_UE(min_idx); % keep the shortest UE
                for idx = 1:numel(allocated_UE) % set others = 0
                    u = allocated_UE(idx);
                    if u ~= keep_ue
                        e_opt(t, u, k) = 0;
                    else
                        e_opt(t, u, k) = 1;
                    end
                end
            end
        end
    end
    
    % calculate the real data rate
    for n = 1:total_UE
        for k = 1:num_RB
            if e_opt(1, n, k) >= 0.5
                signal = P * distance(t, n, user_RU_norm(n))^(-eta) * rayleigh_gain(n, k);
                interference = 0;
                for others = 1:total_UE
                    for i = 1:num_RU
                        if others ~= n && e_opt(1, others, k) >= 0.5 && user_RU_norm(others) ~= user_RU_norm(n)
                            interference = interference + ...
                                P * distance(t, n, user_RU_norm(i))^(-eta) * rayleigh_gain(n, k);
                        end
                    end
                end
                SINR = signal / (interference + sigmsqr);
                data_rate_op(n) = data_rate_op(n) + B * log(1 + SINR);
            end
        end
    end
    record_op = [record_op, sum(log(1+data_rate_op))];


    fprintf('Normal data rate: %.2f\n', record_random(t));
    fprintf('Optmed data rate: %.2f\n', record_op(t));
    
    subplot(1,2,1);
    hold on;
    imagesc(squeeze(e_random));
    xlabel('RB Index');
    ylabel('User Index');
    title('normal RB Allocation');
    xline(0.5 : 1 : double(num_RB), 'Color', [0.5 0.5 0.5]);
    yline(0.5 : 1 : double(total_UE), 'Color', [0.5 0.5 0.5]);
    hold off;

    subplot(1,2,2);
    hold on;
    imagesc(squeeze(e_opt(t, :, :)));
    xlabel('RB Index');
    ylabel('User Index');
    title('Optimized RB Allocation');
    xlim([0.5, double(num_RB)+0.5]); 
    ylim([0.5, double(total_UE)+0.5]); 
    xline(0.5 : 1 : double(num_RB), 'Color', [0.5 0.5 0.5]);
    yline(0.5 : 1 : double(total_UE), 'Color', [0.5 0.5 0.5]);
    colorbar;
    hold off;

end

fprintf('Normal data rate: %.2f\n', record_random);
fprintf('Optmed data rate: %.2f\n', record_op);

figure('Color','w')
hold on;
t_len = length(record_op);
plot(1:t_len, record_random(1:t_len), 'LineWidth', 2, 'Color', 'k');
plot(1:t_len, record_op, 'LineWidth', 2, 'Color', '#8fbc8f');

xlabel('Time Step');
ylabel('Total Data Rate (bps)');
title('Geometric Data Rate');
legend('Normal Allocation', 'Optimized Allocation','Location','southeast');

grid on;
box on;
% ylim([200, max(record_op)*1.05]);