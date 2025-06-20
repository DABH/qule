import multiprocessing
import os

# NOTE: Make sure log.csv exists with the correct headers before running this file

matrix_types = ["rand", "Had", "QFT", "grover"]
alpha = [0.1]
methods = ["GDP*", "GDLM*", "DLR*"]
# n = [4,8,16]
n = [8]
num_tests = 5
noise_mags = [0.01, 0.01, 0.001, 0.0001, 0.00001]
for matrix_type in matrix_types:
	for a in alpha:
		for method in methods:
			for i in n:
					for noise_mag in noise_mags:
						command = "./qule --matrix_type " + matrix_type + " --alpha " + str(a) + " --method " + method + " --n " + str(i) + " --num_tests " + str(num_tests) + " --add_noise --noise_mag " + str(noise_mag)
						print("running: ", command)
						p = multiprocessing.Process(target=os.system, args=(command,))
						p.start()
						# p.join()