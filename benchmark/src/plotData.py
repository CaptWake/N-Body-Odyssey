import numpy as np
import matplotlib.pyplot as plt
import re

PATH = '/home/dadi/Documents/Scuola/Informatica/Magistrale/Anno 1/Semestre 2/Sistemi paralleli e distribuiti/Progetto/resDim.txt'

N_TEST = 14
N_CONFIG = 7
N_PROG = 1

TITLE = "Graph"
LEGEND = ['seq AP', 'seq AP avx', 'AP OMP stat', 'AP OMP dyn', 'BH seq', 'BH stat', 'BH dyn']

def read_file(file_path):
    with open(file_path, 'r') as file:
        data = file.read().split('\n')
    return data

def filter_data(data):
    data = read_file(file_path)
    filtered_data = list(filter(lambda x: "#" in x, data))
    for i in range(len (filtered_data)):
        match = re.search("\d", filtered_data[i])
        if match:
            start_index = match.start()
            end_index = len(filtered_data[i]) - 2
            filtered_data[i] = filtered_data[i][start_index:end_index]
    return filtered_data

file_path = PATH
result = []
for i in range(N_PROG):
    result.append(np.zeros((N_CONFIG, N_TEST)))
data = read_file(file_path)
data = filter_data(data)
x = np.arange(0, N_TEST)
print(len(data))

for i in range(len(data)):
    result[i // (N_CONFIG * N_TEST)][i % N_CONFIG][(i - (i // (N_CONFIG * N_TEST) * N_CONFIG * N_TEST)) // N_CONFIG] = data[i]

fig, axs = plt.subplots(N_PROG, 1)
fig.suptitle(TITLE)

for i in range(N_PROG):
    for j in range(N_CONFIG):
        axs.plot(x, result[i][j], marker='*')
        axs.grid()
fig.legend(LEGEND, shadow=True)
plt.show()