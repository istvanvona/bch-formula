
def Dynkin_BCH(X, Y, maxN = None, commutator = lambda x,y: x*y - y*x, addition = sum):
	"""
	The explicit series expansion of Baker-Campbell-Haussdorf-formula by Dynkin for Lie-algebras.

	Ref.: https://en.wikipedia.org/wiki/Baker%E2%80%93Campbell%E2%80%93Hausdorff_formula#An_explicit_Baker.E2.80.93Campbell.E2.80.93Hausdorff_formula
	"""	
	def rs_partitions(N, n):
		for p in Partitions(N, max_length = 2*n, min_length = n):	
			for q in Permutations(p + [0]*(2*n-len(p))):
				res = [(q[:n][i], q[n:][i]) for i in range(n)]
				if not(0 in map(lambda x:x[0] or x[1], res)):
					yield res
	
	multicommutator = 	lambda *args: (1 if len(args)%2 else -1)*reduce(commutator, args[::-1]) 

	N = 0
	while True if maxN is None else (N <= maxN):
		N += 1
		subterms_n = []
		for n in range(1, N+1):	
			subterms_rs = []
			for rs in rs_partitions(N, n):
				alg_coeff = multicommutator(*sum([[X]*r + [Y]*s for r,s in rs], []))
				fact = reduce(lambda x,y:x*y, [factorial(r)*factorial(s) for r,s in rs], 1)
				subterms_rs.append(alg_coeff/SR(fact))
			subterms_n.append((1 if n%2 else -1)/SR(n)*addition(subterms_rs))
		yield addition(subterms_n)/SR(N)
				





"""
def poincare1d_BCH():
	from sympy import Symbol, Matrix, series, exp, var

	ncvar = lambda s: Symbol(s, commutative = False)

	A,B,C = map(ncvar, ['A','B','C'])
	chi,t0,x0 = var('chi,t0,x0')

	P = LieAlgebra([A,B,C], {1:[(A,B,C),(A,C,B)]})


	BCH_series = sum([term.simplify() for term in Dynkin_BCH(chi*A, (t0*B+ x0*C), 6, commutator = P.commutator)])
	match_series = series(-chi/(1-exp(chi)), chi, 0, 8)
	
	expr = lambda A,B,C: exp(chi*(A + (t0-x0)/2*(B-C) - 1/(1-exp(chi))*(t0*B+x0*C)))
		
	sexpr = expr(A,B,C)

	A = Matrix([[0,1,0],
				[1,0,0],
				[0,0,0]])
	
	B = Matrix([[0,0,1],
    			[0,0,0],
    			[0,0,0]])

	C = Matrix([[0,0,0],
    			[0,0,1],
    			[0,0,0]])
	
	
	mexpr = expr(A,B,C).simplify()

	return BCH_series, match_series, sexpr, (A, B, C), mexpr

"""
