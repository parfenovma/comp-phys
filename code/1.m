clear; clc; close all;

blue_color = [0.12, 0.47, 0.71];
red_color  = [0.84, 0.15, 0.16];
set(0, 'DefaultFigureColor', 'w');
set(0, 'DefaultAxesColor', 'w');
set(0, 'DefaultTextFontName', 'DejaVu Sans');
set(0, 'DefaultAxesFontName', 'DejaVu Sans');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontSize', 12);
set(0, 'DefaultTextFontWeight', 'normal');
set(0, 'DefaultAxesTitleFontWeight', 'normal'); 
set(0, 'DefaultAxesTickDir', 'in');
set(0, 'DefaultAxesLineWidth', 1.5);
set(0, 'DefaultLineLineWidth', 1.7);
set(0, 'DefaultAxesTickLength', [0.015 0.025]);

set(0, 'DefaultAxesXGrid', 'on');
set(0, 'DefaultAxesYGrid', 'on');
set(0, 'DefaultAxesGridLineStyle', '-');
set(0, 'DefaultAxesGridColor', [0.6 0.6 0.6]);
set(0, 'DefaultAxesGridAlpha', 0.5);

a_0 = 0.1;
f_0 = 2.0;
w_0 = 2 * 3.14 * f_0;
T = 1 / f_0;
t = 0: 2*T: 1000;
p_t = 2 * a_0 * sin(3*w_0 * t) + a_0 * cos(5 * w_0 * t);
freqs_f0_multipliers = -6: 1: 6;
freqs = freqs_f0_multipliers * f_0;

coefficients = zeros(length(freqs), length(t));
coefficients(freqs_f0_multipliers == 3)  = -0.1i;
coefficients(freqs_f0_multipliers == -3) =  0.1i;
coefficients(freqs_f0_multipliers == 5)  =  0.05;
coefficients(freqs_f0_multipliers == -5) =  0.05;

c_amp = abs(coefficients);
c_phase = angle(coefficients);
c_real = real(coefficients);
c_im = imag(coefficients);

figure('Position', [110, 55, 1200, 1200]); 

axes('Position', [0.0667, 0.7066, 0.9208, 0.2606]); 
plot(t, p_t, 'b', 'LineWidth', 1.7);
title('Исходный акустический сигнал p(t)');
xlabel('Время t, мкс');
ylabel('p(t), МПа');
grid on;


axes('Position', [0.08, 0.38, 0.39, 0.22]); 
stem(freqs, c_amp, 'filled', 'Color', blue_color, 'MarkerFaceColor', blue_color);
hold on; plot(xlim, [0 0], 'Color', red_color, 'LineWidth', 1.7); hold off;
title('\bf a) \rm Амплитудный спектр |p_T(f_n)|');
xlabel('Частота f_n, МГц');
grid on;

axes('Position', [0.57, 0.38, 0.39, 0.22]); 
stem(freqs, c_phase, 'filled', 'Color', blue_color, 'MarkerFaceColor', blue_color);
hold on; plot(xlim, [0 0], 'Color', red_color, 'LineWidth', 1.7); hold off;
title('\bf b) \rm Фазовый спектр arg(p_T(f_n))')
xlabel('Частота f_n, МГц');
ylabel('Фаза p_T(f_n)');
grid on;

axes('Position', [0.08, 0.07, 0.39, 0.22]); 
stem(freqs, c_real, 'filled', 'Color', blue_color, 'MarkerFaceColor', blue_color);
hold on; plot(xlim, [0 0], 'Color', red_color, 'LineWidth', 1.7); hold off;
title('\bf c) \rm  Действительная часть спектра Re(p_T(f_n))')
xlabel('Частота f_n, МГц');
grid on;

axes('Position', [0.57, 0.07, 0.39, 0.22]); 
stem(freqs, c_im, 'filled', 'Color', blue_color, 'MarkerFaceColor', blue_color);
hold on; plot(xlim, [0 0], 'Color', red_color, 'LineWidth', 1.7); hold off;
title('\bf d) \rm  Мнимая часть спектра Im(p_T(f_n))')
xlabel('Частота f_n, МГц');
grid on;


print(gcf, 'fig_1_m.png', '-dpng', '-r300');
