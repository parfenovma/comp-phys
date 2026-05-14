## Программная реализация ДПФ

Реализуем алгоритм ДПФ для вычисления спектра сеточной функции $p(l)$ по формуле из задания:

$$
\overline{p_T(n)} = \frac{1}{N} \sum_{l=0}^{N-1} p(l)W_N^{-nl}, \quad W_N^{nl} = \exp\left( i\frac{2\pi nl}{N} \right) \qquad
$$


### Использование свойств симметрии
Поскольку исходный сигнал $p(t)$ является действительной функцией, его спектр обладает свойством эрмитовой симметрией. Для $N$-точечного ДПФ:
$$ p_T(N - n) = p_T^*(n) $$
где $*$ обозначает комплексное сопряжение. 

Будем сначала вычислять спектральные коэффициенты $p_T(n)$ по формуле суммы только для первой половины спектра (от $n = 0$ до $n = N/2$), а вторую половину получим операцией комплексного сопряжения.


### Построение графиков
Графики построены в расширенном интервале от $-1/h$ до $1/h$ (т.е. от $-f_s$ до $f_s$). Видно, что амплитуды получились в два раза меньше рассчитанных аналитически (рис. 3.a). Это следствие того, что мы рассматриваем интервал, включающий отрицательные частоты.
Например: $$A_{analytical}⋅cos⁡(5ω_0t)=\frac{A_{analytical}}{2}⋅e^{i5ω0t}+\frac{A_{analytical}}{2}⋅e^{−i5ω0t}$$

![Собственное ДПФ](pic/fig_3_m.png)

\newpage

### Код №3

```matlab
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
stem(freqs_extended, ampliude_ext);
hold on;
y_bounds = ylim;
plot([fs/2, fs/2], y_bounds);
plot([-fs/2, -fs/2], y_bounds);
hold off;
title('Амплитудный спектр ДПФ');
xlabel('Частота, МГц');
ylabel('Амплитуда, МПа');
grid on;

axes('Position', [0.5577, 0.5508, 0.4298, 0.4129]); 
stem(freqs_extended, phase_ext)
title('Фазовый спектр ДПФ');
xlabel('Частота, МГц');
ylabel('Фаза, рад');
grid on;

axes('Position', [0.0628, 0.0583, 0.4298, 0.4129]); 
stem(freqs_extended, real_part_ext)
title('Действительная часть ДПФ');
xlabel('Частота, МГц');
grid on;

axes('Position', [0.5577, 0.0583, 0.4298, 0.4129]); 
stem(freqs_extended, imag_part_ext)
title('Мнимая часть ДПФ');
xlabel('Частота, МГц');
grid on;

print(gcf, 'fig_3_m.png', '-dpng', '-r300');
```