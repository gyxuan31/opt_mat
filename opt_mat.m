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

for t = 1:T
    for i = 1:total_UE
        temp = zeros(1, num_RU);rand
        for j = 1:num_RU
            temp(j) = distance(t,i,j);
        end
        [~, user_RU_norm(i)] = min(temp);
    end
end


% ----- NORMAL - randomly generate e -----
rb_counts = randi([0, 5], 1, total_UE); % initial allocation
e_norm = zeros(total_UE, num_RB);
record_norm = [];
for i = 1:total_UE
    count = rb_counts(i);
    if count > 0
        selected_rbs = randperm(num_RB, count); 
        e_norm(i, selected_rbs) = 1;
    end
end

for t = num_ref+1:T
    data_rate_norm = zeros(1, total_UE);
    for n = 1:total_UE
        for k = 1:num_RB
            if e_norm(n, k) == 1
                signal = P * distance(t, n, user_RU_norm(n)).^(-eta) * rayleigh_gain(n, k);
                interference = 0;
                for others = 1:total_UE
                    for i = 1:num_RU
                        if others ~= n && e_norm(others, k) == 1 && user_RU_norm(others) ~= user_RU_norm(n)
                            interference = interference + ...
                                P * distance(t, n, user_RU_norm(i))^(-eta) * rayleigh_gain(n, k);
                        end
                    end
                end
                SINR = signal / (interference + sigmsqr);
                data_rate_norm(n) = data_rate_norm(n) + B * log(1 + SINR);
            end
        end
    end
    record_norm = [record_norm, sum(log(1+data_rate_norm))];
end

%  ----- OP -----
nvars = double(predicted_len * total_UE * num_RB)

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
    % user_RU_norm;
    % user_RU_op;

    popSize = 20;

    % initPop = rand(popSize, nvars);
    % initPop(1, :) = reshape(e_norm, 1, []);

    objfun = @(e) compute_total_rate(e, predicted_len, total_UE, num_RB, pre_distance, rayleigh_gain, P, sigmsqr, eta, B, T, user_RU_op, num_RU);
    confun = @(e) constraint(e, predicted_len, total_UE, num_RB, gamma, num_setreq, rb_counts);
    options = optimoptions('ga', 'PopulationSize', popSize, 'MaxGenerations', 50, 'Display', 'iter');
    [e_opt, fval] = ga(objfun, nvars, [], [], [], [], lb, ub, confun, options);

    % fmin = @(e) compute_total_rate(e, predicted_len, total_UE, num_RB, pre_distance, rayleigh_gain, P, sigmsqr, eta, B, T, user_RU, num_RU);
    % options = optimoptions('fmincon', 'MaxIterations', 100, 'Display', 'iter-detailed', ...
    %     'MaxFunctionEvaluations', 1e5,'Algorithm','interior-point'); % , 'SpecifyObjectiveGradient',true
    % [e_opt, fval] = fmincon(fmin, repmat(e_norm, predicted_len, 1, 1), [], [], [], [], lb, ub, @constraints, options); % repmat(e_norm, predicted_len, 1, 1)
    e_opt
    for i = 1: nvars
        if e_opt(i) >= 0.5
            e_opt(i) = 1;
        else
            e_opt(i) = 0;
        end
    end
    e_opt = reshape(e_opt, predicted_len, total_UE, num_RB);
    
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


    fprintf('Normal data rate: %.2f\n', record_norm(t));
    fprintf('Optmed data rate: %.2f\n', record_op(t));
    
    % subplot(1,2,1);
    % imagesc(squeeze(e_norm));
    % xlabel('RB Index');
    % ylabel('User Index');
    % title('Optimized RB Allocation (1 = Allocated)');
    % 
    % subplot(1,2,2);
    % imagesc(squeeze(e_opt(t, :, :)));
    % xlabel('RB Index');
    % ylabel('User Index');
    % title('Optimized RB Allocation (1 = Allocated)');
    % colorbar;

end

fprintf('Normal data rate: %.2f\n', record_norm);
fprintf('Optmed data rate: %.2f\n', record_op);

figure('Color','w')
hold on;
t_len = length(record_norm);
plot(1:t_len, record_norm, 'LineWidth', 2, 'Color', 'k');
plot(1:t_len, record_op, 'LineWidth', 2, 'Color', '#8fbc8f');

xlabel('Time Step');
ylabel('Total Data Rate (bps)');
title('Geometric Data Rate');
legend('Normal Allocation', 'Optimized Allocation','Location','southeast');

grid on;
ylim([200, max(record_op)*1.05]);



