import matplotlib
import math as m
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

matplotlib.use('TkAgg')
plt.rcParams['axes.grid'] = True


def fft(s_t, fs, n=8192, shift=False):
    s_w = np.abs(sp.fft.fft(s_t, n=n))  # FFT
    # s_w = 20 * np.log10(s_w / max(s_w))  # logarithmic scale
    if shift:
        s_w = np.concatenate((s_w[m.ceil(s_w.size / 2):], s_w[0:m.ceil(s_w.size / 2)]))  # Matlab FFT shift
        df = np.arange(-fs / 2, fs / 2, fs / s_w.size)  # X axis: frequency scale
    else:
        df = np.arange(0, fs, fs / s_w.size)
    return df, s_w


N = 1024                           # period length
M = 4                                       # count of battery ACC
K = 4                                       # count of signals DDS
B = 12                                      # DDS phase width
f_s = 50                                   # frequency
w_dds = 16                                  # DDS output width
w_ampl = 8                                  # bit depth ampl
w_mult = 14                                 # bit depth mult
w_acc = 15                                  # bit depth acc
A_max = 100                                 # max amplitude


A_i = np.asarray([[10, 15, 20, 25],
                  [30, 35, 40, 45],
                  [50, 55, 60, 65],
                  [70, 75, 80, 85]])       # amplitude
# A_i = np.ones([M, K]) * 100


f_i = np.asarray([[1, 2, 3, 4],
                  [5, 6, 7, 8],
                  [9, 10, 11, 12],
                  [13, 14, 15, 16]])    # frequency
# f_i = np.ones([M, K]) * 5
fi_i = np.zeros([M, K])                 # phase offset


t = np.arange(N)
DDS_out = np.zeros((M, K, N), dtype = complex)
DDS_mult = np.zeros((M, K, N), dtype = complex)
DDS_acc = np.zeros((M, N), dtype = complex)
DDS_res = np.zeros((M, N), dtype = complex)

A_i_code = np.round((A_i / A_max) * (2**(w_ampl) - 1))
# print(A_i_code)

FCW = np.round(f_i / f_s * 2**B)
PCW = np.round(fi_i / 360 * 2**B)


#for Sema

# print(FCW)
# print(A_i_code)

for m in range(M):
    for k in range(K):
        print(f"sudo ./mem 0x83c03{(K * m + k) << 2:03x} 0x{int(A_i_code[m, k]):02x}{int(PCW[m, k]):03x}{int(FCW[m, k]):03x}")

              # f".DATA({{{w_ampl}'d{int(A_i_code[m, k]):d},"
              # f"{B}'d{int(PCW[m, k]):d},{B}'d{int(FCW[m, k]):d}}}));"
              # f"#10; // ampl {A_i[m, k]}%, phi0 {fi_i[m, k]} deg, freq {f_i[m, k]} MHz")
        # f"sudo ./mem 0x83c03{(K * m + k) << 2:03x} 0x{int(A_i_code[m, k]):02x}{int(PCW[m, k]):03x}{int(FCW[m, k]):03x}")


I_max = (2**w_dds / 2 - 1) * (2**w_ampl)
W_P = int(np.ceil(np.log2(I_max)))
Q = 2 ** (w_dds - 1 + w_ampl - w_mult + 1)
A_max_o = np.round(I_max / Q)



#reading

read_m = []
with open("write.txt") as f:
    for line in f:
        read_m.append([int(x, 16) - 2 ** 16 if int(x, 16) >= 2 ** 15 else int(x, 16) for x in line.split()])
read_m = np.array(read_m).transpose()#[:,::4]
# print(read_m)

o_ACC = np.zeros((read_m.shape[0] >> 1, read_m.shape[1]), dtype=complex)
# print(read_m.shape, o_ACC.shape)
for i in range(4):
        o_ACC[i, :] = (read_m[i*2, :] + 1j*read_m[i*2+1, :])# * np.hamming(read_m.shape[1])
# print(o_ACC)

# sys.exit()

output = []
tmp = []
k = 0
with open('output_data.dat', 'r') as file:
    for line in file:
        number_list = [int(num, 16) - 2 ** 32 if int(num, 16) >= 2 ** 26 else int(num, 16) for num in line.split()]
        # print(number_list)
        tmp.append(number_list)
output = np.array(tmp).transpose()


# y = abs(x0 + 1i * x1)
y = abs(output[0, :] + 1j * output[1, :])
y = y[:4096].reshape(256, 4).transpose()
for i in range(len(y)):
    y[i] = y[i] / np.max(y[i])
output = 20 * np.log10(y)

for i in range(M):
    for j in range(K):
        for k in range(N):
            DDS_out[i, j, k] = np.round((2 ** (w_dds - 1) - 2) * np.exp((1j * 2 * np.pi * FCW[i, j] * t[k] + PCW[i, j]) / 2 ** B))
            DDS_mult[i, j, k] = np.round(DDS_out[i, j, k] * A_i_code[i, j] / Q)
        # plt.figure(2)
        # plt.title("MULT")
        # plt.xlabel("time")
        # plt.ylabel("MULT")
        # plt.plot(np.arange(N), DDS_mult[0, j, :].real)
        # # plt.plot(np.arange(N), DDS_mult[0, j, :].imag)
        # plt.legend(["1", "2", "3", "4"])

        # plt.figure(9)
        # plt.title("MULT")
        # plt.xlabel("time")
        # plt.ylabel("MULT")
        # # plt.plot(np.arange(N), DDS_mult[0, j, :].real)
        # plt.plot(np.arange(N), DDS_mult[0, j, :].imag)
        # plt.legend(["1", "2", "3", "4"])
    DDS_acc[i] = DDS_mult[i, :].sum(axis=0)

    # ACC signal

    plt.figure(3)
    plt.title("ACC")
    plt.xlabel("time")
    plt.ylabel("ACC")
    plt.plot(np.arange(N), DDS_acc[0].real)# + DDS_acc[1].real + DDS_acc[2].real + DDS_acc[3].real)
    plt.plot(np.arange(N), DDS_acc[0].imag)# + DDS_acc[1].imag + DDS_acc[2].imag + DDS_acc[3].imag)
    plt.legend(["1", "2", "3", "4"])

    #FFT ACC

    plt.figure(4)
    plt.title("ACC_FFT")
    plt.xlabel("Frequency")
    plt.ylabel("Power")
    x, y = fft(DDS_acc[i], f_s, n=1024)
    plt.plot(x, y.flatten())
    plt.legend(np.arange(M) + 1)

    #output from vivado FFT

    # plt.figure(5)
    # plt.title("output_FFT")
    # plt.xlabel("Frequency")
    # plt.ylabel("Power")
    # plt.plot(np.arange(output[i].size), output[i])
    # plt.legend(np.arange(M) + 1)

    #read signal

    plt.figure(6)
    plt.title("read")
    plt.xlabel("Frequency")
    plt.ylabel("Power")
    plt.plot(np.arange(o_ACC[i].size), o_ACC[0].real)
    plt.plot(np.arange(o_ACC[i].size), o_ACC[0].imag)
    plt.legend(np.arange(M) + 1)

    #FFT read signal
    plt.figure(7)
    plt.title("read_FFT")
    plt.xlabel("Frequency")
    plt.ylabel("Power")
    x, y = fft(o_ACC[i], f_s, n=1024)
    plt.plot(x, y.flatten())
    plt.legend(np.arange(M) + 1)
    # plt.legend(["1"])


plt.show()

DDS_res = DDS_acc.transpose().reshape([1, M * N]).transpose()



