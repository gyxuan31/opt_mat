clear all;
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
load_dis = load('parameter.mat');
distance = load_dis.distance_true;
distance = reshape(distance, T, total_UE, num_RU);
load_pred = load('parameter.mat');
prediction = load_pred.prediction;
prediction = reshape(prediction, T-num_ref, predicted_len, total_UE, num_RU);

rec_dr_random = load('output.mat').rec_dr_random;
rec_dr_avg = load('output.mat').rec_dr_avg;
rec_dr_op = load('output.mat').rec_dr_op;
rec_e_avg = load('output.mat').rec_e_avg;
rec_e_op = load('output.mat').rec_e_op;
rec_e_op = rec_e_op(1:T-num_ref,:,:);
e_random = load('output.mat').e_random;

% PLOT 1 - data rate geometric mean value
figure('Color','w')
hold on;
t_len = length(rec_dr_op);
plot(1:t_len, rec_dr_random(1:t_len), 'LineWidth', 2, 'Color', '#3480b8');
plot(1:t_len, rec_dr_avg(1:t_len), 'LineWidth', 2, 'Color', '#8fbc8f')
plot(1:t_len, rec_dr_op(1:t_len), 'LineWidth', 2, 'Color', '#c82423');

xlabel('Time Step');
ylabel('Geometric Mean of Data Rate');
xlim([1,T-num_ref])
legend('Static Allocation', 'Average Allocation', 'MPC-based Allocation', 'Location','southeast');

grid on;
box on;
hold off;
% ylim([200, max(record_op)*1.05]);

% PLOT 2 - total resource utilization
util_op = zeros(1, T);
util_random = zeros(1, T);
util_avg = zeros(1, T);
for t = 1:T-num_ref
    % RANDOM
    util_random_list = any(e_random, 1); % (total_UE, num_RB)
    temp = sum(util_random_list);
    util_random(t) = temp / double(num_RB);
    % AVG
    e_avag = squeeze(rec_e_avg(t,:,:))
    util_avg_list = any(e_avag, 1);
    util_avg(t) = sum(util_avg_list) / double(num_RB);
    % OP
    e_op = squeeze(rec_e_op(t,:,:)); % (T, total_UE, num_RB);
    util_op_list = any(e_op, 1);
    util_op(t) = sum(util_op_list) / double(num_RB);
end
figure('Color','w')
hold on;
t_len = length(util_op);
plot(1:t_len, util_random(1:t_len), 'LineWidth', 1.5, 'Color', '#3480b8');
plot(1:t_len, util_avg(1:t_len), 'LineWidth', 1.5, 'Color', '#8fbc8f')
plot(1:t_len, util_op(1:t_len), 'LineWidth', 1.5, 'Color', '#c82423');
xlim([1,100]); 
ylim([0,1]);
xlabel('Time Step');
ylabel('RB Utilization (%)');
legend('Static Allocation', 'Average Allocation', 'MPC-based Allocation', 'Location','southeast');

grid on;
box on;
hold off;

% PLOT - every RU utilization
util_ru_op = zeros(1, T);
util_ru_random = zeros(1, T);
util_ru_avg = zeros(1, T);

% calculate user_RU
for t = num_ref+1:T
    for i = 1:total_UE
        temp = zeros(1, num_RU);
        for j = 1:num_RU
            temp(j) = distance(t,i,j);
        end
        [~, user_RU_norm(i)] = min(temp);
    end
end

figure('Color','w')
for a = 1:num_RU
    util_op = zeros(1, T);
    util_random = zeros(1, T);
    util_avg = zeros(1, T);
    for t = 1:T-num_ref
        RU_UE_norm = cell(1, num_RU); % save the index of UE under every RU
        for r = 1:num_RU
            idx = find(user_RU_norm == r);
            RU_UE_norm{r} = idx; % UE index of every RU
        end
        
        % RANDOM
        e_ran = e_random(RU_UE_norm{a},:);
        util_random_list = any(e_ran, 1); % (total_UE, num_RB)
        temp = sum(util_random_list);
        util_random(t) = temp / double(num_RB);
        % AVG
        e_avag = squeeze(rec_e_avg(t,:,:));
        e_avag = e_avag(RU_UE_norm{a},:);
        util_avg_list = any(e_avag, 1);
        util_avg(t) = sum(util_avg_list) / double(num_RB);
        % OP
        e_op = squeeze(rec_e_op(t,:,:)); % (T, total_UE, num_RB);
        e_op = e_op(RU_UE_norm{a}, :);
        util_op_list = any(e_op, 1);
        util_op(t) = sum(util_op_list) / double(num_RB);
    end
    subplot(1,double(num_RU), double(a));
    hold on;
    t_len = length(util_op);
    plot(1:t_len, util_random(1:t_len), 'LineWidth', 1.5, 'Color', '#3480b8');
    plot(1:t_len, util_avg(1:t_len), 'LineWidth', 1.5, 'Color', '#8fbc8f')
    plot(1:t_len, util_op(1:t_len), 'LineWidth', 1.5, 'Color', '#c82423');
    xlim([1,T-num_ref]); 
    ylim([0,1]);
    xlabel('Time Step');
    ylabel('RB Utilization (%)');
  
    grid on;
    box on;

end
legend('Static Allocation', 'Average Allocation', 'MPC-based Allocation', 'Location','northeast');
hold off;


% PLOT 4 - UE number increase
multi_rec_dr_random = load('multi_output.mat').multi_rec_dr_random;
multi_rec_dr_avg = load('multi_output.mat').multi_rec_dr_avg;
multi_rec_dr_op = load('multi_output.mat').multi_rec_dr_op;
multi_num_UE = load('multi_UE.mat').multi_num_UE;

figure('Color','w')
hold on;
t_len = length(multi_rec_dr_op) - 1;
plot(1:t_len, multi_rec_dr_random(1:t_len), 'LineWidth', 1.5, 'Color', '#3480b8');
plot(1:t_len, multi_rec_dr_avg(1:t_len), 'LineWidth', 1.5, 'Color', '#8fbc8f')
plot(1:t_len, multi_rec_dr_op(1:t_len), 'LineWidth', 1.5, 'Color', '#c82423');

xlabel('UE number');
ylabel('Geometric Mean of Data Rate');
legend('Static Allocation', 'Average Allocation', 'MPC-based Allocation', 'Location','southeast');

grid on;
box on;
hold off;
