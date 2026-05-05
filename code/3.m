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
fs = 1 / h;

l_indices = 0:1:N;
t_l = l_indices * h;
p_l = 2 * a0 * sin(3*w0 * t_l) + a0 * cos(5 * w0 * t_l);
p_T_discrete = zeros(N);

half_N = N / 2;
for n = 0 : half_N
    sum_val = 0;
    for l = 0:(N-1)
        exponent = -1i * 2 * pi * n * l / N;
        sum_val = sum_val + p_l(l+1) * exp(exponent);
    end
    p_T_discrete(n+1) = sum_val / N;
end

for n = half_N+1:(N-1)
    p_T_discrete(n) = conj(p_T_discrete((N - n) + 1));
end

n_extended = -N:1:N;
freqs_extended = n_extended .* (fs / N);
p_T_extended = p_T_discrete(mod(n_extended, N) + 1);
p_T_extended(abs(p_T_extended) < 1e-10) = 0;

ampliude_ext = abs(p_T_extended);
phase_ext = angle(p_T_extended);
real_part_ext = real(p_T_extended);
imag_part_ext = imag(p_T_extended);


figure('Position', [110, 55, 1600, 1300]);

axes('Position', [0.0628, 0.5508, 0.4298, 0.4129]); 
stem(freqs_extended, ampliude_ext, 'filled', 'Color', blue_color, 'MarkerFaceColor', blue_color);
hold on;
y_bounds = ylim;
plot([fs/2, fs/2], y_bounds, '--r', 'LineWidth', 1.7, 'Color', red_color);
plot([-fs/2, -fs/2], y_bounds, '--r', 'LineWidth', 1.7, 'Color', red_color);
hold off;
title('\bf a) \rm  Амплитудный спектр ДПФ');
xlabel('Частота, МГц');
ylabel('Амплитуда, МПа');
grid on;

axes('Position', [0.5577, 0.5508, 0.4298, 0.4129]); 
stem(freqs_extended, phase_ext, 'filled', 'Color', blue_color, 'MarkerFaceColor', blue_color)
title('\bf b) \rm  Фазовый спектр ДПФ');
xlabel('Частота, МГц');
ylabel('Фаза, рад');
grid on;

axes('Position', [0.0628, 0.0583, 0.4298, 0.4129]); 
stem(freqs_extended, real_part_ext, 'filled', 'Color', blue_color, 'MarkerFaceColor', blue_color)
title('\bf c) \rm  Действительная часть ДПФ');
xlabel('Частота, МГц');
grid on;

axes('Position', [0.5577, 0.0583, 0.4298, 0.4129]); 
stem(freqs_extended, imag_part_ext, 'filled', 'Color', blue_color, 'MarkerFaceColor', blue_color)
title('\bf d) \rm  Мнимая часть ДПФ');
xlabel('Частота, МГц');
grid on;

print(gcf, 'fig_3_m.png', '-dpng', '-r300');
