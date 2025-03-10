close all;
dl_enable_nw_characteristics = [...
0.119	0.104	0.211	0.284
0.246	0.217	0.383	0.409
0.514	0.432	0.687	0.707
0.71	0.566	0.747	0.869
0.837	0.711	0.869	0.92
0.888	0.905	0.962	0.918
0.941	0.883	0.971	0.95
0.913	0.949	0.989	0.976
0.978	0.951	0.98	0.993
0.993	0.938	0.977	1  
];

dl_disable_nw_characteristics = [...
0.119	0.104	0.211	0.284
0.232	0.24	0.362	0.359
0.418	0.445	0.483	0.492
0.362	0.412	0.64	0.717
0.513	0.613	0.778	0.839
0.6	0.608	0.619	0.743
0.758	0.829	0.903	0.924
0.797	0.806	0.816	0.815
0.793	0.823	0.841	0.858
0.908	0.926	0.947	0.96
];

figure(1);
hold on;
sim_round = 1: 10;
lambda_e_lst = 300 : 100 : 600;
Markers = {'o', 's', 'v', '^'};
for i = [1  4]
    lambda_e = lambda_e_lst(i);
    plot(dl_enable_nw_characteristics(:, i) * 100, 'Marker', Markers{i}, 'Color', 'blue', ...
        'DisplayName', sprintf('ignore..=false,lambda=%d', lambda_e));
    plot(dl_disable_nw_characteristics(:, i) * 100, 'Marker', Markers{i}, 'Color', 'red', ...
        'DisplayName', sprintf('ignore..=true,lambda=%d', lambda_e));    
end
grid;
legend('Location', 'southeast');
xlabel('simulation round');
ylabel('DL PRB utilization (%)');
hold off;