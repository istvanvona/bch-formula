from sympy import tanh



class LieAlgebra:
	
	
	def __init__(self, basis, structure, zero = 0, one = 1):
		for e in basis:
			if basis.count(e) > 1:
				raise RuntimeError('Basis elements must be different.')
		self.zero = zero
		self.one = one
		self.basis = basis

		comfort = lambda x: [x] if type(x) is tuple else x

		self.structure = dict()

		def safe_put(f,idx):
			orig = self.structure.get(idx)
			if orig is None:
				self.structure[idx] = f
			elif orig != f:
				raise RuntimeError('Structure constants are inconsistent.')

	
		for f, idxs in structure.items():
			for a,b,c in comfort(idxs):
				if a == b and f != self.zero:
					raise RuntimeError('Structure constants must be antisymmetric in first two indices.')
				ia, ib, ic = list(map(self.basis.index, [a,b,c]))
				
				if ia < ib:
					safe_put(f, (ia,ib,ic))
				else:
					safe_put(-f, (ib,ia,ic))
	
	def structure_constant(self, a, b, c):
		ia, ib, ic = list(map(self.basis.index, [a,b,c]))
		if ia < ib:
			return self.structure.get((ia,ib,ic), self.zero)
		else:
			return -self.structure.get((ib,ia,ic), self.zero)
		
	def _base_commutator(self, a, b):
		return sum([self.structure_constant(a,b,c)*c for c in self.basis], self.zero)
		
	

	def _list_commutator(self, a, b):
		from sympy import prod	
		for i, A in enumerate(a):
			for j, B in enumerate(b):
				yield prod(a[:i], self.one) * prod(b[:j], self.one) * self._base_commutator(A, B) * prod(b[j+1:], self.one) * prod(a[i+1:], self.one)
				
			
		

	def commutator(self, a, b):
		
		from sympy import prod

		def split(iterable, condition):
			true, false = [], []	
			for x in iterable:
				(false, true)[condition(x)].append(x)
			return true, false

		def expand_powers(term):
			from sympy import S
			for f in term.as_ordered_factors():
				for k,v in f.as_powers_dict().items():
					if S(v).is_integer:
						for i in range(v):
							yield k
					else:
						yield f
						

		def factorize(term):
			elems, coeffs = split(expand_powers(term), lambda x: x in self.basis)
			return prod(coeffs), elems
					

		res = []
		for i in a.as_ordered_terms():
			C_i, X_i  = factorize(i)
			for j in b.as_ordered_terms():
				C_j, X_j = factorize(j)
				res.append(C_i*C_j*sum(self._list_commutator(X_i, X_j), self.zero))

		return sum(res, self.zero)			
		

def rs_partitions(N, n):
	from sympy.utilities.iterables import multiset_permutations, partitions
	k = 2*n
	for m, p in partitions(N, m = k, size = True):	
		if m >= n:	
			for q in multiset_permutations(sum([[k]*v for k,v  in p.items()], []) + [0]*(k-m)):
				res = [(q[:n][i], q[n:][i]) for i in range(n)]
				if not(0 in map(lambda x:x[0] or x[1], res)):
					yield res
	
	
def multicommutator(*args, commutator = lambda x,y: x*y-y*x):
	"""
	Applying commutator multiple times:
		multicommutator(X1, X2, X3) ->  [X1, [X2, X3]]
	"""	
	from sympy.core.compatibility import reduce
	return (1 if len(args)%2 else -1)*reduce(commutator, args[::-1]) 



def Dynkin_BCH(X, Y, maxN = None, commutator = lambda x,y: x*y - y*x, addition = sum):
	"""
	The explicit series expansion of Baker-Campbell-Haussdorf-formula by Dynkin for Lie-algebras.

	Ref.: https://en.wikipedia.org/wiki/Baker%E2%80%93Campbell%E2%80%93Hausdorff_formula#An_explicit_Baker.E2.80.93Campbell.E2.80.93Hausdorff_formula
	"""	
	
	N = 0
	while True if maxN is None else (N <= maxN):
		N += 1
		subterms_n = []
		for n in range(1, N+1):
			from sympy import S
			subterms_rs = []
			for rs in rs_partitions(N, n):
				alg_coeff = multicommutator(*sum([[X]*r + [Y]*s for r,s in rs], []), commutator = commutator)
				from sympy.core.compatibility import reduce
				from sympy import factorial
				fact = reduce(lambda x,y:x*y, [factorial(r)*factorial(s) for r,s in rs], 1)
				subterms_rs.append(alg_coeff/S(fact))
			subterms_n.append((1 if n%2 else -1)/S(n)*addition(subterms_rs))
		yield addition(subterms_n)/S(N)
				



def poincare1d_BCH():
	"""
		Calculations are not quite right, see Poincare1D.ipynb.
	"""
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


