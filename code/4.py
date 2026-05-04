import numpy as np
import matplotlib.pyplot as plt


a0 = 0.1
f0 = 2.0
w0 = 2 * np.pi * f0
T = 0.5
N = 16
h = T / N

l_indices = np.arange(N)
t_l = l_indices * h
p_l = 2 * a0 * np.sin(3 * w0 * t_l) + a0 * np.cos(5 * w0 * t_l)

p_T_discrete = np.zeros(N, dtype=complex)
for n in range(N):
    sum_val = 0j
    for l in range(N):
        exponent = -1j * 2 * np.pi * n * l / N
        sum_val += p_l[l] * np.exp(exponent)
    p_T_discrete[n] = sum_val / N


p_restored = np.zeros(N, dtype=complex)
for l in range(N):
    sum_val = 0j
    for n in range(N):
        exponent = 1j * 2 * np.pi * n * l / N
        sum_val += p_T_discrete[n] * np.exp(exponent)
    p_restored[l] = sum_val


p_restored_real = np.real(p_restored)


h_fine = h / 10
t_fine = np.arange(0, T, h_fine)
p_interpolated = np.zeros(len(t_fine))

for idx, t in enumerate(t_fine):
    sum_sinc = 0
    for l in range(N):
        sum_sinc += p_l[l] * np.sinc((t - l * h) / h)
    p_interpolated[idx] = sum_sinc


p_analytical_fine = 2 * a0 * np.sin(3 * w0 * t_fine) + a0 * np.cos(5 * w0 * t_fine)


fig = plt.figure(figsize=(10, 8))

plt.subplot(2, 1, 1)
plt.plot(t_l, p_l, 'bo-', markersize=8, label='Исходный дискретный $p(l)$')
plt.plot(t_l, p_restored_real, 'rX', markersize=6, label='Восстановленный ОДПФ')
plt.title('Проверка обратимости прямого и обратного ДПФ')
plt.xlabel('Время $t$, мкс')
plt.ylabel('Амплитуда, МПа')
plt.grid(True)
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(t_l, p_l, 'ko', markersize=6, label='Узлы дискретизации $p(l)$')
plt.plot(t_fine, p_interpolated, 'r-', linewidth=2, label='Интерполяция (sinc)')
plt.plot(t_fine, p_analytical_fine, 'b--', alpha=0.5, linewidth=4, label='Аналитический сигнал')

plt.title('Восстановление непрерывного сигнала (Теорема Котельникова)')
plt.xlabel('Время $t$, мкс')
plt.ylabel('Амплитуда, МПа')
plt.grid(True)
plt.legend()

plt.tight_layout()
print(fig.properties)
for i, ax in enumerate(fig.axes):
    bbox = ax.get_position()
    print(f"График {i+1}: [{bbox.x0:.4f}, {bbox.y0:.4f}, {bbox.width:.4f}, {bbox.height:.4f}]")
plt.savefig('fig_4.png', dpi=300)