rng(1);
m = 5; n = 24; 
SCALE = 10000;
B = [25000 12000 12000 11000 17000];
P =[0.281	0.175	0.234	0.059	0.083	0.336	0.406	0.131	0.290	0.367	0.375	0.036	0.016	0.071	0.368	0.041	0.177	0.402	0.223	0.290	0.132	0.288	0.350	0.008;
0.459	0.286	0.383	0.096	0.136	0.549	0.663	0.215	0.474	0.601	0.613	0.058	0.027	0.116	0.602	0.067	0.289	0.656	0.365	0.474	0.216	0.470	0.572	0.013;
0.137	0.085	0.114	0.029	0.041	0.164	0.198	0.064	0.142	0.179	0.183	0.017	0.008	0.035	0.180	0.020	0.086	0.196	0.109	0.141	0.065	0.140	0.171	0.004;
0.589	0.366	0.491	0.123	0.174	0.703	0.850	0.275	0.608	0.770	0.786	0.075	0.034	0.149	0.771	0.086	0.370	0.841	0.468	0.608	0.277	0.603	0.733	0.016;
0.018	0.011	0.015	0.004	0.005	0.022	0.027	0.009	0.019	0.024	0.025	0.002	0.001	0.005	0.024	0.003	0.012	0.026	0.015	0.019	0.009	0.019	0.023	0.001];
T = sin(linspace(-pi, pi, n)) * SCALE;
T = T - min(T) + SCALE;

c = [61000 80000 61000 23000 64000];
R = [0.15 1.18 0.57 2.08 2.43];

nvars = m * n;
lb = zeros(1, nvars);
ub = repmat(T, 1, m); 

fitnessFcn = @(x) ad_objfun(x, m, n, P, R, B);

nonlcon = @(x) ad_confun(x, m, n, P, R, T, c, B);

% options = optimoptions('ga', ...
%     'PopulationSize', 150, ...
%     'MaxGenerations', 100, ...
%     'Display', 'iter');
% x_opt = ga(fitnessFcn, nvars, [], [], [], [], lb, ub, nonlcon, options);
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','MaxFunctionEvaluations',1e7, 'MaxIterations', 1e7);
x_opt = fmincon(fitnessFcn, ones(1,nvars)*10000,[],[], [],[],lb,ub,nonlcon, options);
D = reshape(x_opt, [m, n]);
si = R(:) .* sum(P .* D, 2);

disp('Total allocation per user:'); disp(sum(D,2));
disp('s_i values:'); disp(si);

figure;
imagesc(D);
xlabel('Hour'); ylabel('Service Provider');
title('Allocation Matrix D');
colorbar;
