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


P_fft_raw = np.fft.fft(p_l)


P_fft_normalized = P_fft_raw / N
amplitude_fft = np.abs(P_fft_normalized)


fs = 1 / h
freqs_fft = np.arange(N) * (fs / N)


p_ifft = np.fft.ifft(P_fft_raw)
p_ifft_real = np.real(p_ifft)


fig = plt.figure(figsize=(12, 10))

plt.subplot(2, 1, 1)

half_N = N // 2 + 1
plt.stem(freqs_fft[:half_N], amplitude_fft[:half_N], basefmt=" ",)

plt.plot(6.0, 0.1, 'ro', label='Аналитическая амплитуда (3f0)')
plt.plot(10.0, 0.05, 'ro', label='Аналитическая амплитуда (5f0)')

plt.title('Амплитудный спектр, полученный с помощью БПФ (нормированный на $N$)')
plt.xlabel('Частота, МГц')
plt.ylabel('Амплитуда, МПа')
plt.grid(True)
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(t_l, p_l, 'k-', linewidth=2, label='Исходный сигнал $p(l)$')
plt.plot(t_l, p_ifft_real, 'r--', linewidth=2, label='Восстановленный $ifft$')
plt.title('Проверка обратимости прямого и обратного БПФ (FFT / IFFT)')
plt.xlabel('Время $t$, мкс')
plt.ylabel('Амплитуда, МПа')
plt.grid(True)
plt.legend()

plt.tight_layout()
print(fig.properties)
for i, ax in enumerate(fig.axes):
    bbox = ax.get_position()
    print(f"График {i+1}: [{bbox.x0:.4f}, {bbox.y0:.4f}, {bbox.width:.4f}, {bbox.height:.4f}]")
plt.savefig('fig_5.png', dpi=300)
