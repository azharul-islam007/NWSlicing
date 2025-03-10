close all;
dl_prb_utilization = [...
0.713895	0.549983
0.665822	0.580695
0.777606	0.660124
0.802623	0.698949
];

figure(1);
hold on;
lambda_e_lst = 300 : 100 : 600;
plot(lambda_e_lst', dl_prb_utilization(:, 1) * 100, 'Color', 'blue', 'Marker', 'o', ...
    'DisplayName', 'ignore..=false');
plot(lambda_e_lst', dl_prb_utilization(:, 2) * 100, 'Color', 'red', 'Marker', '^', ...
    'DisplayName', 'ignore..=true');
xlabel('lambda');
ylabel('DL PRB utilization (%)');
legend('Location', 'southeast');
grid;
hold off;