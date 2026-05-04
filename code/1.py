import numpy as np
import matplotlib.pyplot as plt

a0 = 0.1 # MPa
f0 = 2.0 # MHz
w0 = 2 * np.pi * f0
T = 1 / f0 # microseconds

t = np.linspace(0, 2*T, 1000)
p_t = 2 * a0 * np.sin(3 * w0 * t) + a0 * np.cos(5 * w0 * t)

freqs_f0_multipliers = np.arange(-6, 7) 
freqs = freqs_f0_multipliers * f0

coeffs = np.zeros(len(freqs), dtype=complex)

coeffs[freqs_f0_multipliers == 3] = -0.1j
coeffs[freqs_f0_multipliers == -3] = 0.1j
coeffs[freqs_f0_multipliers == 5] = 0.05
coeffs[freqs_f0_multipliers == -5] = 0.05

amplitude = np.abs(coeffs)
phase = np.angle(coeffs)
real_part = np.real(coeffs)
imag_part = np.imag(coeffs)

fig = plt.figure(figsize=(12, 12))

plt.subplot(3, 1, 1)
plt.plot(t, p_t, 'b-', linewidth=2)
plt.title('Исходный акустический сигнал $p(t)$', fontsize=14)
plt.xlabel('Время $t$, мкс')
plt.ylabel('$p(t)$, МПа')
plt.grid(True)


plt.subplot(3, 2, 3)
plt.stem(freqs, amplitude)
plt.title('Амплитудный спектр $|p_T(f_n)|$')
plt.xlabel('Частота, МГц')
plt.grid(True)

plt.subplot(3, 2, 4)
plt.stem(freqs, phase)
plt.title('Фазовый спектр $arg(p_T(f_n))$')
plt.xlabel('Частота, МГц')
plt.ylabel('Радианы')
plt.grid(True)

# Действительная и мнимая части
plt.subplot(3, 2, 5)
plt.stem(freqs, real_part)
plt.title('Действительная часть $Re(p_T(f_n))$')
plt.xlabel('Частота, МГц')
plt.grid(True)

plt.subplot(3, 2, 6)
plt.stem(freqs, imag_part)
plt.title('Мнимая часть $Im(p_T(f_n))$')
plt.xlabel('Частота, МГц')
plt.grid(True)

plt.tight_layout()
for i, ax in enumerate(fig.axes):
    bbox = ax.get_position()
    print(f"График {i+1}: [{bbox.x0:.4f}, {bbox.y0:.4f}, {bbox.width:.4f}, {bbox.height:.4f}]")
plt.savefig('fig_1.png', dpi=300)
