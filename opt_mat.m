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

% calculate every UE connect which RU
for t = 1:T
    for i = 1:total_UE
        temp = zeros(1, num_RU);
        for j = 1:num_RU
            temp(j) = distance(t,i,j);
        end
        [~, user_RU_norm(i)] = min(temp);
    end
end

% init
rb_counts = randi([0, num_RB/num_RU], 1, total_UE); % initial allocation
e_random = zeros(total_UE, num_RB);
e_avg = zeros(total_UE, num_RB);
rec_dr_random = []; % record data rate
rec_dr_avg = [];
rec_dr_op = [];
rec_e_avg = zeros(T, total_UE, num_RB);
rec_e_op = zeros(T, total_UE, num_RB);

% RANDOM - randomly generate e - fix allocation
for i = 1:total_UE
    count = rb_counts(i);
    if count > 0
        selected_rbs = randperm(num_RB, count); 
        e_random(i, selected_rbs) = 1;
    end
end

% ---------- NORMAL BASELINE ----------
for t = num_ref+1:T
    RU_UE_norm = cell(1, num_RU); % save the index of UE under every RU
    for r = 1:num_RU
        idx = find(user_RU_norm == r);
        RU_UE_norm{r} = idx; % UE index of every RU
    end

    % ----- RANDOM -----
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
    rec_dr_random = [rec_dr_random, sum(log(1+data_rate_random))];

    % ----- AVERAGE -----
    for r = 1:num_RU
        ue_list = RU_UE_norm{r};
        num_ue = numel(ue_list);

        if num_ue > 0
            RB_per_UE = floor(num_RB / num_ue); % num of RB allocated to every UE
            remaining_RB = num_RB - RB_per_UE * num_ue;

            RB_pool = randperm(num_RB);

            for i = 1:num_ue
                u = ue_list(i);
                start_idx = (i - 1) * RB_per_UE + 1;
                end_idx = i * RB_per_UE;

                e_avg(u, RB_pool(start_idx:end_idx)) = 1;
            end

            for j = 1:remaining_RB % remaining RB randomly allocate
                u = ue_list(j);
                e_avg(u, RB_pool(RB_per_UE * num_ue + j)) = 1;
            end
        end
    end
    rec_e_avg(t,:,:) = e_avg(:,:);

    data_rate_static = zeros(1, total_UE);
    for n = 1:total_UE % calculate data rate
        for k = 1:num_RB
            if e_avg(n, k) == 1
                signal = P * distance(t, n, user_RU_norm(n))^(-eta) * rayleigh_gain(n, k);
                interference = 0;

                for others = 1:total_UE
                    for i = 1:num_RU
                        if others ~= n && e_avg(others, k) == 1 && user_RU_norm(others) ~= user_RU_norm(n)
                            interference = interference + ...
                                P * distance(t, n, user_RU_norm(i))^(-eta) * rayleigh_gain(n, k);
                        end
                    end
                end

                SINR = signal / (interference + sigmsqr);
                data_rate_static(n) = data_rate_static(n) + B * log(1 + SINR);
            end
        end
    end
    rec_dr_avg = [rec_dr_avg, sum(log(1 + data_rate_static))];
end

%  ---------- OP ----------
nvars = double(predicted_len * total_UE * num_RB);

lb = zeros(1, nvars);
ub = ones(1, nvars);



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

    popSize = 80;
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
    for i = 1: nvars
        if e_opt(i) >= 0.5
            e_opt(i) = 1;
        else
            e_opt(i) = 0;
        end
    end
    e_opt = reshape(e_opt, predicted_len, total_UE, num_RB);  % shape(predicted_len, total_UE, num_RB)

    % check RB cannot allocated to 2 UEs under one RU
    for i = 1:num_RU % set the repeat UE = 0
        for k = 1:num_RB
            ue_list = RU_UE_norm{i}; % UE under RU(i)
            allocated_UE = ue_list(e_opt(1, ue_list, k) > 0);

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
                        e_opt(1, u, k) = 0;
                    else
                        e_opt(1, u, k) = 1;
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
    rec_dr_op = [rec_dr_op, sum(log(1+data_rate_op))];
    rec_e_op(t,:,:) = e_opt(1,:,:)


    fprintf('Normal data rate: %.2f\n', rec_dr_random(t));
    fprintf('Normal data rate: %.2f\n', rec_dr_avg(t));
    fprintf('Optmed data rate: %.2f\n', rec_dr_op(t));

    subplot(1,2,1);
    hold on;
    imagesc(squeeze(e_random));
    xlabel('RB Index');
    ylabel('User Index');
    title('normal RB Allocation');
    xlim([0.5, double(num_RB)+0.5]); 
    ylim([0.5, double(total_UE)+0.5]); 
    xline(0.5 : 1 : double(num_RB), 'Color', [0.5 0.5 0.5]);
    yline(0.5 : 1 : double(total_UE), 'Color', [0.5 0.5 0.5]);
    hold off;

    subplot(1,2,2);
    hold on;
    imagesc(squeeze(e_opt(1, :, :)));
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

fprintf('Normal data rate: %.2f\n', rec_dr_random);
fprintf('Optmed data rate: %.2f\n', rec_dr_op);

% PLOT 1 - data rate geometric mean value
figure('Color','w')
hold on;
t_len = length(rec_dr_op);
plot(1:t_len, rec_dr_random(1:t_len), 'LineWidth', 2, 'Color', '#3480b8');
plot(1:t_len, rec_dr_avg(1:t_len), 'LineWidth', 2, 'Color', '#8fbc8f')
plot(1:t_len, rec_dr_op(1:t_len), 'LineWidth', 2, 'Color', '#c82423');

xlabel('Time Step');
ylabel('Geometric Mean of Data Rate');
legend('Static Allocation', 'Average Allocation', 'MPC-based Allocation', 'Location','southeast');

grid on;
box on;
hold off;
% ylim([200, max(record_op)*1.05]);

% PLOT 2 - resource utility
util_op = zeros(1, T);
util_random = zeros(1, T);
util_avg = zeros(1, T);
for t = 1:T
    % RANDOM
    util_random_list = any(e_random, 1); % (total_UE, num_RB)
    util_random(t) = sum(util_random_list) / num_RB;
    % AVG
    e_avag = squeeze(rec_e_avg(t,:,:));
    util_avg_list = any(e_avag, 1);
    util_avg(t) = sum(util_avg_list) / num_RB;
    % OP
    e_op = squeeze(rec_e_op(t,:,:)); % (T, total_UE, num_RB);
    util_op_list = any(e_op, 1);
    util_op(t) = sum(util_op_list) / num_RB;
end
figure('Color','w')
hold on;
t_len = length(util_op);
plot(1:t_len, util_random(1:t_len), 'LineWidth', 2, 'Color', '#3480b8');
plot(1:t_len, util_avg(1:t_len), 'LineWidth', 2, 'Color', '#8fbc8f')
plot(1:t_len, util_op(1:t_len), 'LineWidth', 2, 'Color', '#c82423');

xlabel('Time Step');
ylabel('RB Utility (%)');
legend('Static Allocation', 'Average Allocation', 'MPC-based Allocation', 'Location','southeast');

grid on;
box on;
hold off;

% PLOT 3 - data rate of every RU

        

% PLOT 4 - UE number increase
