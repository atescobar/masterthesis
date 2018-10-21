import numpy as np 


def falsePosition(r, x1, x2, tol=1e-3):
	print("trying to find root")
	f1 = r(x1)
	f2 = r(x2)
	print(np.sign(f1))
	print(np.sign(f2))
	if (True):#((f1 < 0 and f2 > p) or (f2 < 0 and f1 > p)):
		print("root bracketed")
		it =500
		for i in range(it):
			x3 = x2-f2*(x2-x1)/(f2-f1)
			f3 = f(x3)

			if(f2*f3 > 0):
				x2 = x3

			if (f1*f3> 0):
				x1 = x3

			if (abs(f3)<tol):
				return x3
		print("convergence problem: too many iterations, "+str(it))

	else:
		print("root is not bracketed")
	return None

