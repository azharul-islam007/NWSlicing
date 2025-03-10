close all;
ul_enable_nw_characteristics = [...
0.393	0.373	0.474	0.476
0.679	0.881	0.855	0.809
0.877	0.963	0.993	0.996
0.991	0.995	0.999	0.992
0.996	0.998	0.988	0.988
0.984	0.994	0.983	0.986
0.988	0.995	0.999	0.993
0.957	0.991	0.987	0.999
0.98	0.995	0.998	0.991
0.999	0.999	0.991	1
];

ul_disable_nw_characteristics = [...
0.393	0.373	0.474	0.476
0.945	0.91	0.959	0.928
0.913	0.914	0.918	0.918
0.762	0.797	0.814	0.82
0.828	0.897	0.904	0.904
0.867	0.874	0.877	0.867
0.933	0.958	0.963	0.967
0.82	0.819	0.819	0.822
0.861	0.867	0.865	0.867
0.945	0.967	0.968	0.971
];

figure(1);
hold on;
sim_round = 1: 10;
lambda_e_lst = 300 : 100 : 600;
Markers = {'o', 's', 'v', '^'};
for i = [1  4]
    lambda_e = lambda_e_lst(i);
    plot(ul_enable_nw_characteristics(:, i) * 100, 'Marker', Markers{i}, 'Color', 'blue', ...
        'DisplayName', sprintf('ignore..=false,lambda=%d', lambda_e));
    plot(ul_disable_nw_characteristics(:, i) * 100, 'Marker', Markers{i}, 'Color', 'red', ...
        'DisplayName', sprintf('ignore..=true,lambda=%d', lambda_e));    
end
grid;
legend('Location', 'southeast');
xlabel('simulation round');
ylabel('UL PRB utilization (%)');
hold off;