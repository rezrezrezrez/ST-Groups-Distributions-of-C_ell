from sympy import Matrix, symbols, BlockDiagMatrix
import numpy as np
import sympy as sp
from sympy.utilities.iterables import multiset_combinations
from math import factorial
from sage.all import *, binomial as sg_binomial


def ID(p):
    deg = p+1

    F.<z> = CyclotomicField(p^2)
    genus = euler_phi(p^2)
    P = PolynomialRing(F,genus/2,'u')

    blocks = []
    for i in [0..(genus/2)-1]:
        blocks.append(diagonal_matrix([P.gens()[i],(P.gens()[i])^(-1)]))

    return block_diagonal_matrix(blocks)


def IDsympy(p):
    genus = euler_phi(p**2)

    u = symbols(f"u1:{genus//2 + 1}")  
# creates a list containing [u1, u2, u3, ...], adjust range if you want to start the ordering differently

    blocks = []
    for i in range(genus // 2):
        zeta = u[i]
        block = sp.Matrix(
            [
                [zeta, 0],
                [0, 1/zeta]
            ]
        )
        blocks.append(block)

    identityMatrix = BlockDiagMatrix(
        *blocks
    ).as_explicit() 

    return identityMatrix


def charPoly_CoeffEig_diag(matrix, n, term):
# "n" is dimension of the matrix, "term" is the coefficient you're after
    k = n - term

    eigenvalues = matrix.diagonal() 
    combs = multiset_combinations(eigenvalues, k)

    total = sum(np.prod(comb) for comb in combs)

    return total * ((-1) ** (k))


def multinomial(*args):
    factorials = [math.factorial(i) for i in range(n+1)]

    total = factorials[sum(args)] 

    for arg in args:
        total //= factorials[arg] 

    return total


def combinations_sum(n, num_elements):
    if num_elements == 1:
        yield (n,)
    else:
        for i in range(n+1):
            for combination in combinations_sum(n - i, num_elements - 1):
                yield (i,) + combination


def nth_moment(n):
    Mn = 0;

    for combination in combinations_sum(n, ...): 
#be sure to change the number (not "n") to match the scenario
        ... = combination 
# "unpack" the combination variables (which are ...-tuples)

        uProd = np.prod(...) 
# the full nth moment of coefficient expression MUST go inside the parenthesis of np.prod()

        Mn += multinomial(*combination) * uProd

    return Mn


def normalizedtracea1(p,bound): 
	tracelist=[]

	for q in prime_range(p+1,2^bound):
		qlist=[q]
		P.<x,y,z> = GF(q)[]
		C = Curve(y^p*z-x^(p+1)+x*z^p)

		if Mod(q-1,p^2)==0: #only count points if q=1 mod p^2
			numpoints=C.count_points(1)
			N1=numpoints[0]; S1=N1-(q+1)
			normtrace1=S1/sqrt(q)
			qlist.append(normtrace1)
		else: 
			qlist.append(0)
			tracelist.append(qlist)

	return tracelist
