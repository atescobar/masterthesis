import numpy as np 


def falsePosition(f, x1, x2, tol=1):
	f1 = f(x1)
	f2 = f(x2)
	if ((f1 < 0 and f2 > 0) or (f2 < 0 and f1 > 0)):
		print("root bracketed")
		it =1000
		for i in range(it):
			x3 = x2-f2*(x2-x1)/(f2-f1)
			f3 = f(x3)

			s1 = np.sign(f1)
			s2 = np.sign(f2)
			s3 = np.sign(f2)

			if(s2*s3 > 0):
				f2 = f(x3)
				print(f2)
			else:
				f1 = f(x3)
				print(f1)

			if (s1*s3 > 0):
				f1 = f(x3)
				print(f1)
			else:
				f2 = f(x3)
				print(f2)

			if((i%10)==1):
				print(f3)

			if (abs(f3)<tol):
				return x3
		

		print("convergence problem: too many iterations, "+str(it))

	else:
		print("root is not bracketed")
	return None

