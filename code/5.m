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
P_fft_raw = fft(p_l);
P_fft_normalized = P_fft_raw / N;
amplitude_fft = abs(P_fft_normalized);
fs = 1 / h;
freqs_fft = (0 : (N - 1)) * (fs / N);
p_ifft = ifft(P_fft_raw);
p_ifft_real = real(p_ifft);


figure('Position', [110, 55, 1600, 1300]);
axes('Position', [0.0650, 0.5513, 0.9225, 0.4114]);
half_N = (N / 2) + 1;
stem(freqs_fft(1:half_N), amplitude_fft(1:half_N), 'filled', 'b', 'DisplayName', 'Спектр ДПФ');
hold on;
plot(6.0, 0.1, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Аналитическая амплитуда (3f_0)');
plot(10.0, 0.05, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Аналитическая амплитуда (5f_0)');
hold off;
title('\bf a) \rm  Амплитудный спектр, полученный с помощью БПФ (нормированный на N)');
xlabel('Частота, МГц');
ylabel('Амплитуда, МПа');
grid on;
legend('Location', 'northeast');

axes('Position', [0.0650, 0.0598, 0.9225, 0.4114]);
plot(t_l, p_l, 'k-', 'LineWidth', 2, 'DisplayName', 'Исходный сигнал p(l)');
hold on;
plot(t_l, p_ifft_real, 'r--', 'LineWidth', 2, 'DisplayName', 'Восстановленный ifft');
hold off;
title('\bf b) \rm  Проверка обратимости прямого и обратного БПФ (FFT / IFFT)');
xlabel('Время \itt\rm, мкс');
ylabel('Амплитуда, МПа');
grid on;
legend('Location', 'northeast');
print(gcf, 'fig_5_m.png', '-dpng', '-r300');

% evaluate absolute error
error_fft    = p_ifft_real - p_l;
max_err_fft    = max(abs(error_fft));
figure('Position', [150, 150, 900, 800]);
stem(t_l, error_fft, 'filled', 'Color', red_color, 'MarkerFaceColor', red_color, 'LineWidth', 1.5);
hold on;
plot(xlim, [0 0], 'k-', 'LineWidth', 1);
hold off;
title('Абсолютная ошибка восстановления точек для встроенной реализации');
xlabel('Время \itt\rm, мкс');
ylabel('\Delta\itp\rm, МПа');
grid on;
text_str2 = sprintf('Max Error = %.2e', max_err_fft);
legend(text_str2, 'Location', 'best');
print(gcf, 'fig_5_error.png', '-dpng', '-r300');
