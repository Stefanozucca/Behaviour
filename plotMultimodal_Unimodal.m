function plotMultimodal_Unimodal(Unimodal, Multimodal)
% Generate synthetic data for demonstration
N = size(Multimodal,1); % Number of subjects
% Multimodal = AllMultimodal_Zone_cumulative(:,5); % Random time spent in left arm
% Unimodal = AllUnimodal_Zone_cumulative(:,5); % Random time spent in right arm

% Compute mean and SEM
mean_Multimodal = mean(Multimodal);
mean_Unimodal = mean(Unimodal);
sem_Multimodal = std(Multimodal) / sqrt(N);
sem_Unimodal = std(Unimodal) / sqrt(N);

% Create the figure
figure;
hold on;

% Plot data points and connect them with lines
for i = 1:N
    plot([1, 2], [Unimodal(i), Multimodal(i)], '-bo', 'MarkerSize', 8, 'LineWidth', 1.5,'Color', [0.5, 0.5, 0.5],'MarkerFaceColor', [0.5, 0.5, 0.5]);
end
hold on

% Plot mean and SEM
errorbar(1, mean_Unimodal, sem_Unimodal, 'k', 'LineWidth', 2, 'CapSize', 10, 'Marker', 'o', 'MarkerFaceColor', 'k');
errorbar(2, mean_Multimodal, sem_Multimodal, 'k', 'LineWidth', 2, 'CapSize', 10, 'Marker', 'o', 'MarkerFaceColor', 'k');
plot([1, 2], [mean_Unimodal, mean_Multimodal], 'k-', 'LineWidth', 2);


% Customize plot
xlim([0.5, 2.5]);
ylim([0 200]);
xticks([1 2]);
xticklabels({'Unimodal', 'Multimodal'});
ylabel('Time Spent (s)');
% grid on;
% title('Time Spent in Left and Right Arm for Each Subject');

hold off;

end
