import numpy as np
import matplotlib.pyplot as plt


a0 = 0.1  # MPa
f0 = 2.0  # MHz
w0 = 2 * np.pi * f0
T = 0.5  # мкс
N = 16   # points
h = T / N  # 0.03125 microseconds
fs = 1 / h  # 32 MHz


l_indices = np.arange(N)
t_l = l_indices * h
p_l = 2 * a0 * np.sin(3 * w0 * t_l) + a0 * np.cos(5 * w0 * t_l)

p_T_discrete = np.zeros(N, dtype=complex)

half_N = N // 2
for n in range(half_N + 1):
    sum_val = 0j
    for l in range(N):
        exponent = -1j * 2 * np.pi * n * l / N
        sum_val += p_l[l] * np.exp(exponent)
    p_T_discrete[n] = sum_val / N


for n in range(half_N + 1, N):
    p_T_discrete[n] = np.conj(p_T_discrete[N - n])

n_extended = np.arange(-N, N + 1)
freqs_extended = n_extended * (fs / N)
p_T_extended = p_T_discrete[n_extended % N]

p_T_extended[np.abs(p_T_extended) < 1e-10] = 0

amplitude_ext = np.abs(p_T_extended)
phase_ext = np.angle(p_T_extended)
real_part_ext = np.real(p_T_extended)
imag_part_ext = np.imag(p_T_extended)


fig = plt.figure(figsize=(12, 10))

plt.subplot(2, 2, 1)
plt.stem(freqs_extended, amplitude_ext, basefmt=" ")
plt.title('Амплитудный спектр ДПФ')
plt.xlabel('Частота, МГц')
plt.ylabel('Амплитуда, МПа')
plt.axvline(x=fs/2, color='r', linestyle='--', alpha=0.5)
plt.axvline(x=-fs/2, color='r', linestyle='--', alpha=0.5)
plt.grid(True)

plt.subplot(2, 2, 2)
plt.stem(freqs_extended, phase_ext, basefmt=" ")
plt.title('Фазовый спектр ДПФ')
plt.xlabel('Частота, МГц')
plt.ylabel('Фаза, рад')
plt.grid(True)

plt.subplot(2, 2, 3)
plt.stem(freqs_extended, real_part_ext, basefmt=" ")
plt.title('Действительная часть ДПФ')
plt.xlabel('Частота, МГц')
plt.grid(True)

plt.subplot(2, 2, 4)
plt.stem(freqs_extended, imag_part_ext, basefmt=" ")
plt.title('Мнимая часть ДПФ')
plt.xlabel('Частота, МГц')
plt.grid(True)

plt.tight_layout()
plt.savefig('fig_3.png', dpi=300)
