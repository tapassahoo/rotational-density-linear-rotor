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
	print(rotational_constant)
	exit()

	tol = 10e-16
	maxj = 1000
	for id, j in enumerate(range(maxj)):
		exp1 = math.exp(-tau * bconst * j * (j + 1.0))
		if (exp1 < tol):
			print(j, exp1, -tau * bconst * j * (j + 1.0))
			jmax = j
			print('jmax=', jmax)
			break
	if (id == (maxj - 1)):
		jmax = maxj
		print('!!! Warning: maxj is reached')

	dens = np.zeros(size_theta)
	erot = np.zeros(size_theta)
	erotsq = np.zeros(size_theta)

	cstep = 2.0 / float(size_theta - 1)
	cost = np.array([ic * cstep - 1.0 for ic in range(size_theta)])

	for ic in range(size_theta):
		sum = 0.0
		sum_eng = 0.0
		sum_engsq = 0.0
		for j in range(jmax + 1):
			Nj = (2.0 * j + 1.0)
			Pn = legendre(j)
			tmp = Nj * Pn(cost[ic]) * math.exp(-tau * bconst * j * (j + 1.0))
			sum += tmp
			sum_eng += tmp * j * (j + 1.0)
			sum_engsq += tmp * j * (j + 1.0) * j * (j + 1.0)
		dens[ic] = sum
		erot[ic] = sum_eng
		erotsq[ic] = sum_engsq

	dens = dens / (4.0 * math.pi)
	erot = bconst * erot / (4.0 * math.pi * dens * nslice)
	erotsq = bconst * bconst * erotsq / \
		(4.0 * math.pi * nslice * nslice * dens)
	dens_comb = np.array([cost, dens, erot, erotsq])
	np.savetxt(
		'linden1.dat',
		dens_comb.T,
		delimiter=' ',
		header='First col. --> ei.ej; 2nd and 3rd and 4th cols are the density and energy estimator and heat capacity estimator, respectively. ')
	exit()

	cstep1 = 2.0 / float(size_theta1 - 1)
	pstep1 = 2.0 * math.pi / float(size_phi1 - 1)
	cost1 = np.array([ic * cstep1 - 1.0 for ic in range(size_theta1)])
	sint1 = np.sqrt(1.0 - cost1 * cost1)
	phi1 = np.array([ip * pstep1 for ip in range(size_phi1)])

	cosg = np.zeros((size_theta1 * size_phi1 * size_theta1 * size_phi1), float)
	dens1 = np.zeros(
		(size_theta1 *
		 size_phi1 *
		 size_theta1 *
		 size_phi1),
		float)
	erot1 = np.zeros(
		(size_theta1 *
		 size_phi1 *
		 size_theta1 *
		 size_phi1),
		float)
	erotsq1 = np.zeros(
		(size_theta1 *
		 size_phi1 *
		 size_theta1 *
		 size_phi1),
		float)
	indices = np.zeros(
		(size_theta1 *
		 size_phi1 *
		 size_theta1 *
		 size_phi1),
		int)

	np.savetxt(
		'linden_indices.dat',
		(size_theta1,
		 size_phi1),
		fmt='%6d',
		delimiter=' ')
	ii = 0
	for ic1 in range(size_theta1):
		for ip1 in range(size_phi1):
			for ic2 in range(size_theta1):
				for ip2 in range(size_phi1):
					phi12 = phi1[ip1] - phi1[ip2]

					cosg[ii] = cost1[ic1] * cost1[ic2] + \
						sint1[ic1] * sint1[ic2] * math.cos(phi12)
					sum = 0.0
					sum_eng = 0.0
					sum_engsq = 0.0
					for j in range(jmax + 1):
						Nj = (2.0 * j + 1.0)
						Pn = legendre(j)
						tmp = Nj * Pn(cosg[ii]) * \
							math.exp(-tau * bconst * j * (j + 1.0))
						sum += tmp
						sum_eng += tmp * j * (j + 1.0)
						sum_engsq += tmp * j * (j + 1.0) * j * (j + 1.0)
					dens1[ii] = sum
					erot1[ii] = sum_eng
					erotsq1[ii] = sum_engsq
					indices[ii] = ii
					ii = ii + 1
	dens1 = dens1 / (4.0 * math.pi)
	#erot1 = bconst*erot1/(4.0*math.pi*dens1*nslice)
	erot1 = bconst * erot1 / (4.0 * math.pi * nslice)
	#erotsq1 = bconst*bconst*erotsq1/(4.0*math.pi*nslice*nslice*dens1)
	erotsq1 = bconst * bconst * erotsq1 / (4.0 * math.pi * nslice * nslice)
	dens_comb1 = np.array([cosg, dens1, erot1, erotsq1])
	#np.savetxt('linden.dat',dens_comb1.T, fmt='%10.8e', delimiter=' ',header='First col. --> ei.ej; 2nd and 3rd and 4th cols are the density and energy estimator and heat capacity estimator, respectively. ')
	np.savetxt('linden.dat', dens_comb1.T, fmt='%10.8e', delimiter=' ')
	cmd = 'cat linden.dat >> linden_indices.dat'
	os.system(cmd)
	call(['mv', 'linden_indices.dat', 'linden.out'])
