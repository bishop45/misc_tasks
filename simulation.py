# -*- coding: utf-8 -*-
import pylab as pl
import math
import numpy as np

#解析解
def exact(x,pe):
	return	(-math.exp(pe)+math.exp(pe*x))/(1-math.exp(pe))
#中心差分
def central_difference(idx,pe,dx):
	if (2.0-pe*dx == 0):
		return 0
	else:
		c2 = 1.0/(1.0-(((2.0+pe*dx)/(2.0-pe*dx))**(1.0/dx)))
		c1 = 1.0 - c2
		return c1 + c2 *(((2.0+pe*dx)/(2.0-pe*dx))**idx)
#風上差分
def upwind_difference(idx,pe,dx):
	c2 = 1.0/(1.0-(1.0+pe*dx)**(1.0/dx))
	c1 = 1.0 - c2
	return c1 + c2*((1.0+pe*dx)**idx)
#数値粘性
def numerical_viscosity(idx,pe,dx):
	return upwind_difference(idx+1,pe,dx) - 2*upwind_difference(idx,pe,dx) + upwind_difference(idx-1,pe,dx) / (2.0 * dx)

def main():
	l_dn = [10,20]
	l_pe = [1,5,10,20,50]
	for dn in l_dn:
		for pe in l_pe:
			phi = np.zeros((5,dn+1))
			for i in range(0,dn+1):
		 		phi[0][i] = 1.0/dn*i
		 		phi[1][i] = exact(phi[0][i],pe)
		 		phi[2][i] = central_difference(i,pe,1.0/dn)
				phi[3][i] = upwind_difference(i,pe,1.0/dn)
				phi[4][i] = numerical_viscosity(i,pe,1.0/dn)

			pl.scatter(phi[0], phi[1], color="r", label="exact")
			pl.scatter(phi[0], phi[2], color="g", label="central difference")
			pl.scatter(phi[0], phi[3], color="b", label="upwind difference")
			pl.title("dn:{0},pe:{1}.".format(dn,pe))
			pl.xlabel('x')
			pl.ylabel('phi')
			pl.legend (bbox_to_anchor=(0.5, -0.3), loc='center')
			pl.subplots_adjust(bottom=0.3)
			pl.savefig("dn{0}pe{1}.png".format(dn,pe))
			pl.clf()

			pl.scatter(phi[0], phi[3], color="b", label="upwind difference")
			pl.scatter(phi[0], phi[4], color="y",label = "numerical_viscosity")
			pl.xlim([-0.2,1.2])
			pl.ylim([-0.2,10])
			pl.title("dn:{0},pe:{1}.".format(dn,pe))
			pl.xlabel('x')
			pl.ylabel('Numerical_viscosity')
			pl.legend (bbox_to_anchor=(0.5, -0.3), loc='center')
			pl.subplots_adjust(bottom=0.3)
			pl.savefig("Numerical_viscosity_dn{0}pe{1}.png".format(dn,pe))
			pl.clf()

if __name__ == '__main__':
	main()
