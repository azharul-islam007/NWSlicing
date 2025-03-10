close all;
ul_prb_utilization = [...
0.884401	0.826914
0.918494	0.837709
0.926829	0.855933
0.923048	0.85396
];

figure(1);
hold on;
lambda_e_lst = 300 : 100 : 600;
plot(lambda_e_lst', ul_prb_utilization(:, 1) * 100, 'Color', 'blue', 'Marker', 'o', ...
    'DisplayName', 'ignore..=false');
plot(lambda_e_lst', ul_prb_utilization(:, 2) * 100, 'Color', 'red', 'Marker', '^', ...
    'DisplayName', 'ignore..=true');
xlabel('lambda');
ylabel('UL PRB utilization (%)');
legend('Location', 'southeast');
grid;
hold off;