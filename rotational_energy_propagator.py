import sys
import os
from subprocess import call
import numpy as np
from scipy.special import legendre
import argparse

parser = argparse.ArgumentParser(
	description="It computes rotational energy propagator for a linear rotor.",
	epilog="Enjoy the program! :)")
parser.add_argument(
	"--temperature",
	type=float,
	help="The temperature of the system in Kelvin",
	required=True)
parser.add_argument(
	"--trotter_number",
	type=int,
	metavar='P',
	help="The Trotter number.",
	required=True)
parser.add_argument(
	"--rotational_constant",
	type=float,
	metavar='B Constant',
	help="The rotational B constant in wave-number.",
	required=True)
parser.add_argument(
	"--size",
	type=int,
	metavar="NUMBER",
	help="The size of the grid",
	required=True)
parser.add_argument(
	"--spin_isomer",
	type=str,
	default="spinless",
	help="Coupling between the nuclear spin and the angular motion.")
args = parser.parse_args()

if __name__ == '__main__':
	temperature=args.temperature
	trotter_number=args.trotter_number
	rotational_constant=args.rotational_constant
	size_theta=args.size
	spin_isomer=args.spin_isomer

	if (spin_isomer=="spinless"):
		isomer="-"
		basis_type=""
	if (spin_isomer=="para"):
		isomer = "-p-"
		basis_type="even"
	if (spin_isomer=="ortho"):
		isomer = "-o-"
		basis_type="odd"

	CMRECIP2KL=1.4387672	   	# cm^-1 to Kelvin conversion factor
	beta=1.0/temperature
	tau=beta/trotter_number
	rotational_constant=rotational_constant*CMRECIP2KL

	tol = 10e-16
	maxj = 1000
	for count, j in enumerate(range(maxj)):
		boltzmann_term = np.exp(-tau * rotational_constant * j * (j + 1.0))
		if (boltzmann_term < tol):
			print(j, boltzmann_term, -tau * rotational_constant * j * (j + 1.0))
			jmax = j
			print('jmax=', jmax)
			break
	if (count == (maxj - 1)):
		jmax = maxj
		print('!!! Warning: maxj is reached')

	dens = np.zeros(size_theta)
	erot = np.zeros(size_theta)

	cost = np.linspace(-1.0, 1.0, size_theta, endpoint=True)

	for ic in range(size_theta):
		sum = 0.0
		sum_eng = 0.0
		for j in range(jmax + 1):
			Nj = (2.0 * j + 1.0)
			Pn = legendre(j)
			tmp = Nj * Pn(cost[ic]) * np.exp(-tau * rotational_constant * j * (j + 1.0))
			sum += tmp
			sum_eng += tmp * j * (j + 1.0)
		dens[ic] = sum
		erot[ic] = sum_eng

	dens = dens / (4.0 * np.pi)
	erot = rotational_constant * erot / (4.0 * np.pi * dens * trotter_number)
	dens_comb = np.array([cost, dens, erot])

	home = os.path.expanduser("~")
	plot_dir_path = "academic-project/output/"
	final_result_path = home + "/" + plot_dir_path + "final-pigs-outputs-for-plotting/"
	output_file = final_result_path + "rotational-energy-propagator-linear-rotor-rotational-constant" + str(args.rotational_constant) + "wavenumber-temperatute" + str(args.temperature) + "kelvin-trotter-number" + str(trotter_number) + ".txt" 
	print(output_file)

	np.savetxt(
		output_file,
		dens_comb.T,
		delimiter=' ',
		header='First col. --> ei.ej; 2nd and 3rd are the density and energy estimator, respectively. ')
	print("Done!")
