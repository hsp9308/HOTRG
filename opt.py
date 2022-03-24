import time
import numpy as np
from scipy.optimize import minimize_scalar
import subprocess
import sys


def func(x):
	params = x
	s = subprocess.check_output(['./TRG_it'] + ['%d' % int(sys.argv[1])] + [str(params)])
	print(s)
	s1 = complex(s)
	sr = np.real(s1)
	si = np.imag(s1)
	print(sr*sr+si*si)
	print(params)
	return (sr*sr+si*si)



start_time=time.time()
res = minimize_scalar(func, bounds=(0.04,0.05),method='bounded',options={'disp':True,'xtol':1e-8})

print(res.x)
print("--------  %s seconds elapsed  --------" % (time.time() - start_time))
