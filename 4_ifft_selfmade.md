## Обратное ДПФ и интерполяция по теореме Котельникова

В этом разделе решается обратная задача: восстановление временного сигнала из его дискретного спектра (ОДПФ) и реконструкция непрерывной формы сигнала между узлами дискретизации.

### Программная реализация обратного ДПФ (ОДПФ)
Вычисление сеточной функции $p(l)$ по найденным амплитудам гармоник проводится по формуле (4):
$$ p(l) = \sum_{n=0}^{N-1} p_T(n) \exp\left( i \frac{2\pi nl}{N} \right) $$

Теоретически, прямое и обратное преобразования формируют строго обратимую пару. Это можно видеть на графиках в конце данного раздела. 

### Восстановление непрерывного сигнала

Для восстановления непрерывной функции $p(t)$ используется ряд Котельникова (формула 5):
$$ p(t) = \sum_{l=0}^{N-1} p(l) \text{sinc} \left( \frac{\pi}{h} (t - lh) \right) $$

Функция $\text{sinc}(x) = \frac{\sin(x)}{x}$ в частотной области представляет собой идеальный фильтр нижних частот (ФНЧ) прямоугольной формы. Суммирование сдвинутых функций $\text{sinc}$, умноженных на значения выборки $p(l)$, позволяет "сгладить" ступенчатый дискретный сигнал.

В библиотеке `NumPy` и в реализации sinc в `Matlab` используется нормированная функция `sinc(x)`, которая вычисляется как $\frac{\sin(\pi x)}{\pi x}$. Поэтому при программировании аргумент функции умножать на $\pi$ не нужно, формула принимает вид `numpy.sinc((t - l*h) / h)`.

### Графики результатов
Вычисления производились на сетке с шагом $h/10$. Видно, что интерполированная кривая проходит через исходные узлы дискретизации и совпадает с исходной аналитической функцией, что доказывает теорему Котельникова.

![ОДПФ и интерполяция](pic/fig_4_m.png)

\newpage

### Код № 4

```matlab
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

title('Проверка обратимости прямого и обратного ДПФ');
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
title('Восстановление непрерывного сигнала (Теорема Котельникова)');
xlabel('Время t, мкс');
ylabel('Амплитуда, МПа');
grid on;
legend('Location', 'best');

print(gcf, 'fig_4_m.png', '-dpng', '-r300');
```
