import matplotlib
import math as m
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')
plt.rcParams['axes.grid'] = True


def fft(s_t, fs, n=8192, shift=True):
    s_w = np.abs(sp.fft.fft(s_t, n=n))  # FFT
    s_w = 20 * np.log10(s_w / max(s_w))  # logarithmic scale
    if shift:
        s_w = np.concatenate((s_w[m.ceil(s_w.size / 2):], s_w[0:m.ceil(s_w.size / 2)]))  # Matlab FFT shift
        df = np.arange(-fs / 2, fs / 2, fs / s_w.size)  # X axis: frequency scale
    else:
        df = np.arange(0, fs, fs / s_w.size)
    return df, s_w

#input
N = 4096                         # period length
A_i = np.ones([4, 4]) * 100      # amplitude

#const
M = 4                                       # count of battery ACC
K = 4                                       # count of signals DDS
B = 12                                      # DDS phase width
f_s = 100e6                                 # frequency
w_dds = 16                                  # DDS output width
w_ampl = 8
w_mult = 14
w_acc = 15
A_max = 100                                 # max amplitude

f_i = np.ones([M, K]) * 10e6            # frequency
fi_i = np.zeros([M, K])                 # phase offset

t = np.arange(N)
DDS_out = np.zeros((M, K, N), dtype = complex)
DDS_mult = np.zeros((M, K, N), dtype = complex)
DDS_acc = np.zeros((M, N), dtype = complex)
DDS_res = np.zeros((M, N), dtype = complex)

A_i_code = np.round((A_i / A_max) * (2**(w_ampl) - 1))

FCW = np.round(f_i / f_s * 2**B)
PCW = np.round(fi_i / 360 * 2**B)

#
I_max = (2**w_dds / 2 - 1) * (2**w_ampl)
W_P = int(np.ceil(np.log2(I_max)))
# Q = 2 ** (W_P - w_mult + 1)
Q = 2 ** (w_dds - 1 + w_ampl - w_mult + 1)
A_max_o = np.round(I_max / Q)


for i in range(M):
    for j in range(K):
        for k in range(N):
            DDS_out[i, j, k] = np.round((2 ** (w_dds - 1) - 2) * np.exp((1j * 2 * np.pi * FCW[i, j] * t[k] + PCW[i, j]) / 2 ** B))
            DDS_mult[i, j, k] = np.round(DDS_out[i, j, k] * A_i_code[i, j] / Q)
            # print(np.log2(np.max(np.abs(DDS_mult[j, :]))))
            # print(DDS_mult[j, :])
    DDS_acc[i] = DDS_mult[i, :].sum(axis=0)
    # x, y = fft(DDS_mult[i, :], f_s, n=N)
# print(DDS_acc[0])

DDS_res = DDS_acc.transpose().reshape([1, M * N]).transpose()
# print(DDS_res)

# real_part = DDS_acc[0].real
# imaginary_part = DDS_acc[0].imag

# x, y = fft(DDS_mult, f_s, n=N)
# plt.plot(x, y)
# plt.show()


#output
print(DDS_res)


# save as file

# out = np.stack([DDS_res.real, DDS_res.imag]).transpose()[0, :, :]
# out[:, :] = out
# print(out, out.shape)
# print(DDS_acc.real.shape)
# np.savetxt('freal.dat', DDS_acc.real.transpose(), fmt='%d %d %d %d', delimiter='\n')
# np.savetxt('imag.dat', DDS_acc.imag.transpose(), fmt='%d %d %d %d', delimiter='\n')


# with open('DDS_res_real.dat', 'w') as file:                                       /////////////////////////////////
#     pass
# with open('DDS_res_imag.dat', 'w') as file:
#     pass
#
# # file_r = open('DDS_res_real.dat', 'a')                                        ////////////////////////////////
# file_i = open('DDS_res_imag.dat', 'a')
#
# data_to_save = np.array([DDS_acc[0].real, DDS_acc[0].imag]) #[real_part, imaginary_part])
# print(data_to_save)
# for i in range(103):
#     # np.savetxt(file_i, imaginary_part, fmt='%d')                                  //////////////////////
#     # np.savetxt(file_r, real_part, fmt='%d')
#     np.savetxt(file_i, data_to_save, fmt='%d')
# file_i.close()
# file_r.close()

np.savetxt("DDS_res.txt", DDS_res, fmt='%d')
# np.savetxt("DDS_res_imag.dat", imaginary_part, fmt='%d', delimiter='\n')


# plt.plot(t, all_DDS.real)
# plt.plot(t, all_DDS.imag)
# plt.grid(True)
# plt.savefig("DDS_out.png")

