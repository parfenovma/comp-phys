clear; clc; close all;

blue_color = [0.12, 0.47, 0.71];
red_color  = [0.84, 0.15, 0.16];
set(0, 'DefaultFigureColor', 'w');
set(0, 'DefaultAxesColor', 'w');
set(0, 'DefaultTextFontName', 'DejaVu Sans');
set(0, 'DefaultAxesFontName', 'DejaVu Sans');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontSize', 12);
set(0, 'DefaultAxesLineWidth', 1.5);
set(0, 'DefaultLineLineWidth', 1.7);
set(0, 'DefaultTextFontWeight', 'normal');
set(0, 'DefaultAxesTitleFontWeight', 'normal'); 
set(0, 'DefaultAxesTickDir', 'in');
set(0, 'DefaultAxesTickLength', [0.015 0.025]);

set(0, 'DefaultAxesXGrid', 'on');
set(0, 'DefaultAxesYGrid', 'on');
set(0, 'DefaultAxesGridLineStyle', '-');
set(0, 'DefaultAxesGridColor', [0.6 0.6 0.6]);
set(0, 'DefaultAxesGridAlpha', 0.5);


a0 = 0.1;
f0 = 2.0;
w0 = 2 * pi * f0;
T = 0.5;
N = 16;
h = T / N;

l_indices = 0 : (N - 1);
t_l = l_indices * h;
p_l = 2 * a0 * sin(3 * w0 * t_l) + a0 * cos(5 * w0 * t_l);

p_T_discrete = zeros(1, N);
for n = 0 : (N - 1)
    sum_val = 0;
    for l = 0 : (N - 1)
        exponent = -1i * 2 * pi * n * l / N;
        sum_val = sum_val + p_l(l + 1) * exp(exponent);
    end
    p_T_discrete(n + 1) = sum_val / N;
end

p_restored = zeros(1, N);
for l = 0 : (N - 1)
    sum_val = 0;
    for n = 0 : (N - 1)
        exponent = 1i * 2 * pi * n * l / N;
        sum_val = sum_val + p_T_discrete(n + 1) * exp(exponent);
    end
    p_restored(l + 1) = sum_val;
end

p_restored_real = real(p_restored);

h_fine = h / 10;

t_fine = 0 : h_fine : (T - h_fine); 

p_interpolated = zeros(1, length(t_fine));

for idx = 1 : length(t_fine)
    t = t_fine(idx);
    sum_sinc = 0;
    for l = 0 : (N - 1)
        sum_sinc = sum_sinc + p_l(l + 1) * sinc((t - l * h) / h);
    end
    p_interpolated(idx) = sum_sinc;
end

p_analytical_fine = 2 * a0 * sin(3 * w0 * t_fine) + a0 * cos(5 * w0 * t_fine);


figure('Position', [110, 55, 1500, 1050]);

axes('Position', [0.0781, 0.5651, 0.9069, 0.3890]);
plot(t_l, p_l, 'bo-', 'MarkerSize', 8, 'DisplayName', 'Исходный дискретный p(l)');
hold on;
plot(t_l, p_restored_real, 'rx', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Восстановленный ОДПФ');
hold off;

title('\bf a) \rm  Проверка обратимости прямого и обратного ДПФ');
xlabel('Время t, мкс');
ylabel('Амплитуда, МПа');
grid on;
legend('Location', 'best'); 


axes('Position', [0.0781, 0.0747, 0.9069, 0.3890]);
plot(t_l, p_l, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'DisplayName', 'Узлы дискретизации p(l)');
hold on;
plot(t_fine, p_interpolated, 'r-', 'LineWidth', 2, 'DisplayName', 'Интерполяция (sinc)');
plot(t_fine, p_analytical_fine, '--', 'Color', [0.5, 0.6, 1.0], 'LineWidth', 4, 'DisplayName', 'Аналитический сигнал');
hold off;
title('\bf b) \rm  Восстановление непрерывного сигнала (Теорема Котельникова)');
xlabel('Время t, мкс');
ylabel('Амплитуда, МПа');
grid on;
legend('Location', 'best');


print(gcf, 'fig_4_m.png', '-dpng', '-r300');

% evaluate absolute error
error_manual = p_restored_real - p_l;
max_err_manual = max(abs(error_manual));
figure('Position', [150, 150, 900, 800]);
stem(t_l, error_manual, 'filled', 'Color', red_color, 'MarkerFaceColor', red_color, 'LineWidth', 1.5);
hold on;
plot(xlim, [0 0], 'k-', 'LineWidth', 1);
hold off;
title('Абсолютная ошибка восстановления точек для собственной реализации');
xlabel('Время \itt\rm, мкс');
ylabel('\Delta\itp\rm, МПа');
grid on;
print(gcf, 'fig_4_error.png', '-dpng', '-r300');
