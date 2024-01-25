"""
Examples
Below, we apply our algorithms to compute
* the modular polynomial phi19 (over QQ)
* the modular polynomial phi103 (over a prime field for a suitable prime)
"""




load("main.sage")

#allow sage to do parallel computations
os.environ['SAGE_NUM_THREADS'] = '10'


#modular polynomial phi_{19} over ZZ
phi19 = modular_polynomial(19,verbose=False,ncpu=10); 


#modular polynomial phi_{103} over F_p with p = 79103
ell = 103
[n,u,v] = [7,5,0]
p = 79103
phi103 = modular_polynomial_modp(p,ell,n,u,v)