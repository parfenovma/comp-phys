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
l_indices = 0 : (N - 1);
t_l = l_indices * h;
p_l = 2 * a0 * sin(3 * w0 * t_l) + a0 * cos(5 * w0 * t_l);
freqs_f0_multipliers = -6 : 6;
freqs = freqs_f0_multipliers * f0;
coeffs_deriv = zeros(1, length(freqs));
coeffs_deriv(freqs_f0_multipliers == 3)  =  0.6 * pi * f0;
coeffs_deriv(freqs_f0_multipliers == -3) =  0.6 * pi * f0;
coeffs_deriv(freqs_f0_multipliers == 5)  =  1i * 0.5 * pi * f0;
coeffs_deriv(freqs_f0_multipliers == -5) = -1i * 0.5 * pi * f0;


figure('Position', [110, 55, 1600, 1300]);

axes('Position', [0.0391, 0.5390, 0.4498, 0.4156]);
stem(freqs, abs(coeffs_deriv), 'filled', 'Color', blue_color, 'MarkerFaceColor', blue_color);
title('Амплитудный спектр dp/dt');

axes('Position', [0.5377, 0.5390, 0.4498, 0.4156]);
stem(freqs, angle(coeffs_deriv), 'filled', 'Color', blue_color, 'MarkerFaceColor', blue_color);
title('Фазовый спектр dp/dt');

axes('Position', [0.0391, 0.0484, 0.4498, 0.4156]);
stem(freqs, real(coeffs_deriv), 'filled', 'Color', blue_color, 'MarkerFaceColor', blue_color);
title('Действительная часть dp/dt');

axes('Position', [0.5377, 0.0484, 0.4498, 0.4156]);
stem(freqs, imag(coeffs_deriv), 'filled', 'Color', blue_color, 'MarkerFaceColor', blue_color);
title('Мнимая часть dp/dt');

print(gcf, 'fig_add_1_m.png', '-dpng', '-r300');

p_T_discrete = zeros(1, N);
for n = 0 : (N - 1)
    sum_val = 0;
    for l = 0 : (N - 1)
        exponent = -1i * 2 * pi * n * l / N;
        sum_val = sum_val + p_l(l + 1) * exp(exponent);
    end
    p_T_discrete(n + 1) = sum_val / N;
end

freqs_mapped = zeros(1, N);
for n = 0 : (N - 1)
    if n <= N/2
        freqs_mapped(n + 1) = n * (fs / N);
    else
        freqs_mapped(n + 1) = (n - N) * (fs / N);
    end
end

omega_mapped = 2 * pi * freqs_mapped;

p_T_deriv_discrete = p_T_discrete .* (1i * omega_mapped); 

n_extended = -N : N;
freqs_extended = n_extended * (fs / N);
p_T_deriv_extended = p_T_deriv_discrete(mod(n_extended, N) + 1);
p_T_deriv_extended(abs(p_T_deriv_extended) < 1e-10) = 0;

figure('Position', [110, 55, 1500, 1050]);

axes('Position', [0.0391, 0.5390, 0.4498, 0.4141]);
stem(freqs_extended, abs(p_T_deriv_extended), 'filled', 'Color', blue_color, 'MarkerFaceColor', blue_color);
hold on;
y_bounds = ylim;
plot([fs/2, fs/2], y_bounds, '--', 'Color', red_color, 'LineWidth', 1.5);
plot([-fs/2, -fs/2], y_bounds, '--', 'Color', red_color, 'LineWidth', 1.5);
hold off;
title('Амплитуда P_discrete в интервале [-1/h, 1/h]');


axes('Position', [0.5377, 0.5390, 0.4498, 0.4141]);
stem(freqs_extended, angle(p_T_deriv_extended), 'filled', 'Color', blue_color, 'MarkerFaceColor', blue_color);
title('Фаза P_{discrete}');


axes('Position', [0.0391, 0.0484, 0.4498, 0.4141]);
stem(freqs_extended, real(p_T_deriv_extended), 'filled', 'Color', blue_color, 'MarkerFaceColor', blue_color);
title('Re(P_{discrete})');


axes('Position', [0.5377, 0.0484, 0.4498, 0.4141]);
stem(freqs_extended, imag(p_T_deriv_extended), 'filled', 'Color', blue_color, 'MarkerFaceColor', blue_color);
title('Im(P_{discrete})');

print(gcf, 'fig_add_2_m.png', '-dpng', '-r300');

P_fft = fft(p_l);

% there is no fftfreq in matlab, so these freqs came out from python script
freqs_fft = [0 : (N/2 - 1), -N/2 : -1] * (fs / N); 
omega_fft = 2 * pi * freqs_fft;
P_deriv_fft = P_fft .* (1i * omega_fft);
p_deriv_ifft = ifft(P_deriv_fft);
p_deriv_ifft_real = real(p_deriv_ifft);
p_deriv_analytical = 0.6 * w0 * cos(3 * w0 * t_l) - 0.5 * w0 * sin(5 * w0 * t_l);

figure('Position', [110, 55, 1500, 750]);
plot(t_l, p_deriv_analytical, 'k-', 'LineWidth', 3, 'DisplayName', 'Аналитическая dp/dt');
hold on;
plot(t_l, p_deriv_ifft_real, 'r--o', 'MarkerSize', 6, 'LineWidth', 2, 'DisplayName', 'Восстановленная ifft');
hold off;

title('Восстановление производной через обратное БПФ');
xlabel('Время t, мкс');
ylabel('Амплитуда, МПа/мкс');
legend('Location', 'northeast');

print(gcf, 'fig_add_3_m.png', '-dpng', '-r300');
