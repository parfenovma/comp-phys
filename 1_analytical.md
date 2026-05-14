# Самостоятельная работа: Дискретное преобразование Фурье. Вариант 1

## Аналитический расчет коэффициентов ряда Фурье


$$
\left\{
\begin{aligned}
a_0 &= 0.1 \text{ МПа} \\
f_0 &= 2 \text{ МГц} \\
\omega_0 &= 2\pi f_0 \\
p(t) &= 2a_0 \sin(3\omega_0 t) + a_0 \cos(5\omega_0 t)
\end{aligned}
\right.
$$

### Определение периода $T$
Сигнал состоит из двух гармоник с частотами $f_1 = 3f_0 = 6$ МГц и $f_2 = 5f_0 = 10$ МГц. 

Определим базовую частоту:
$$
f_0 = \text{НОД(6 МГц, 10 МГц)} = \text{2МГц}
$$


Следовательно, период сигнала $T$ равен:
$$ 
T = \frac{1}{f_0} = \frac{1}{2 \cdot 10^6} = 0.5 \text{ мкс}
$$

### Расчет аналитических коэффициентов спектра $p_T(f_n)$
Используем формулы Эйлера для представления тригонометрических функций через комплексные экспоненты:

$$
p(t) = 2a_0 \sin(3\omega_0 t) + a_0 \cos(5\omega_0 t)
$$


$$
\implies 
\left[
\begin{aligned}
 \sin(x) = \frac{e^{ix} - e^{-ix}}{2i} = -\frac{i}{2}e^{ix} + \frac{i}{2}e^{-ix} \\
 \cos(x) = \frac{e^{ix} + e^{-ix}}{2} = \frac{1}{2}e^{ix} + \frac{1}{2}e^{-ix}
\end{aligned}
\right]
$$
$$ \implies p(t) = 2(0.1) \left[ -\frac{i}{2}e^{i 3\omega_0 t} + \frac{i}{2}e^{-i 3\omega_0 t} \right] + 0.1 \left[ \frac{1}{2}e^{i 5\omega_0 t} + \frac{1}{2}e^{-i 5\omega_0 t} \right] $$
$$ \implies p(t) = -0.1i \cdot e^{i 3\omega_0 t} + 0.1i \cdot e^{-i 3\omega_0 t} + 0.05 \cdot e^{i 5\omega_0 t} + 0.05 \cdot e^{-i 5\omega_0 t} $$

Из полученного выражения напрямую следуют аналитические значения комплексных коэффициентов спектра $p_T(f_n)$ на частотах $f_n$:


$$
\begin{aligned}
& 3f_0 (6 МГц): p_T(3f_0) = -0.1i  & \text{Фаза}& = -\pi/2 \\
& -3f_0 (-6 МГц): p_T(-3f_0) = 0.1i; &  \text{Фаза}&= \pi/2 \\
&5f_0 (10 МГц): p_T(5f_0) = 0.05; & \text{Фаза}& = 0 \\
&-5f_0 (-10 МГц): p_T(-5f_0) = 0.05; &  \text{Фаза}&= 0\\
\end{aligned}
$$

### Графики функции и аналитического спектра

Амплитудный и действительный спектры четные (рис. 1.a, 1.c), а фазовый и мнимый спектры нечетные (рис. 1.b, 1.d).

![Аналитический спектр и сигнал](pic/fig_1_m.png)

\newpage

### Код №1
```matlab
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
plot(t, p_t);
title('Исходный акустический сигнал p(t)');
xlabel('Время t, мкс');
ylabel('p(t), МПа');
grid on;

axes('Position', [0.08, 0.38, 0.39, 0.22]); 
stem(freqs, c_amp);
hold on; plot(xlim, [0 0]); hold off;
title('Амплитудный спектр |p_T(f_n)|');
xlabel('Частота f_n, МГц');
grid on;

axes('Position', [0.57, 0.38, 0.39, 0.22]); 
stem(freqs, c_phase);
hold on; plot(xlim, [0 0]); hold off;
title('Фазовый спектр arg(p_T(f_n))')
xlabel('Частота f_n, МГц');
ylabel('Фаза p_T(f_n)');
grid on;

axes('Position', [0.08, 0.07, 0.39, 0.22]); 
stem(freqs, c_real);
hold on; plot(xlim, [0 0]); hold off;
title('Действительная часть спектра Re(p_T(f_n))')
xlabel('Частота f_n, МГц');
grid on;

axes('Position', [0.57, 0.07, 0.39, 0.22]); 
stem(freqs, c_im);
hold on; plot(xlim, [0 0]); hold off;
title('Мнимая часть спектра Im(p_T(f_n))')
xlabel('Частота f_n, МГц');
grid on;


print(gcf, 'fig_1_m.png', '-dpng', '-r300');
```
