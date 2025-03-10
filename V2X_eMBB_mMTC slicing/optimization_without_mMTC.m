data = [
    1.36, 0.06; % V2X and eMBB for episode 1
    1.44, 0.00; % V2X and eMBB for episode 2
    0.33, 2.03; % V2X and eMBB for episode 3
    1.39, 2.64; % V2X and eMBB for episode 4
    1.38, 0.15; % V2X and eMBB for episode 5
    1.80, 0.36; % V2X and eMBB for episode 6
    2.42, 1.89; % V2X and eMBB for episode 7
    1.41, 2.05; % V2X and eMBB for episode 8
    1.37, 0.00; % V2X and eMBB for episode 9
    0.03, 1.96; % V2X and eMBB for episode 10
    1.20, 0.01; % V2X and eMBB for episode 11
    0.10, 1.61; % V2X and eMBB for episode 12
    2.03, 2.44; % V2X and eMBB for episode 13
    1.46, 0.01; % V2X and eMBB for episode 14
    1.48, 0.37  % V2X and eMBB for episode 15
];

% The q_values are the episode numbers in this case
q_values = 1:15;

% Create a new figure
figure;

% Plot the V2X values against the Q values
plot(q_values, data(:, 1), 'LineWidth', 2, 'Color', [0 0.5 0.5]);
hold on;

% Plot the eMBB values against the Q values
plot(q_values, data(:, 2), 'LineWidth', 2, 'Color', [0.5 0 0]);

% Add labels and title to the plot
xlabel('actions');
ylabel('V2X and eMBB values');
%title('Line Plot of Q values vs V2X and eMBB values');

% Add a legend
legend('V2X', 'eMBB');

% Display the grid
grid on;

hold off;