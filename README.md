<!--
author:   Jutho Haegeman

email:    jutho.haegeman@ugent.be

language: en

-->

[![LiaScript](https://raw.githubusercontent.com/LiaScript/LiaScript/master/badges/course.svg)](https://LiaScript.github.io/course/?https://raw.githubusercontent.com/Jutho/VectorFunctionSpacesSummary/main/README.md)

# Vector and Function Spaces - Summary

For every chapter of the lecture notes "Vector- and function spaces", this document lists the following items:

1. Important concepts: these mostly correspond to definitions, sometimes the definition is hidden inside a proposition or theorem. You do not need to know the literal definition, but you should be able to understand and use this term or concept correctly

2. Important theorems and propositions: it is indicated which proofs need to be known (actively) or understood

3. Additional topics for applications / exercises: what you need to know for the exercises, in particular if it is additional material that you do not need to know for theory

## 1. Elementary algebraic structures

**Important concepts:**

* Set, subset, intersection, union
* Map, domain, codomain, argument, image, injective, surjective, bijective
* Set cardinality: finite, countably infinite, uncountable $=$ uncountably infinite
* Equivalence relation (and partial order relation)
* Group, ring, field
* Structure preserving map $=$ homomorphism, isomorphism, endomorphism, automorphism group
* Permutation, parity or signature of permutation
* Vector space, linear combination, linear span, complete set, linear independence, basis, dimension, coordinates (coordinate isomorphism)
* Subspace, disjoint subspaces, direct sum, complement, codimension, quotient space, 
* Algebra, commutator

**Lemmas, propositions and theorems:**

No proofs from Chapter 1

**For applications / exercises**

## 2. Linear maps and matrices

**Important concepts:**

* Linear map, property of linearity ($=$ additivity + homogeneity), composition, identity
* Kernel ($=$ null space), image ($=$ range), rank, nullity
* Matrix, matrix representation of a linear map
* Matrix multiplication, transpose, hermitian conjugate, (anti)symmetric and (anti)Hermitian matrix
* determinant (Leibniz formula) and trace
* inverse, singular matrix ($=$ zero determinant) and nonsingular matrix
* General linear group ($=$ invertible matrices), basis transform $=$ similarity transform
* Linear functional, dual space ($=$ linear functionals), basis transform of linear functional, dual linear map ($\cong$ matrix transpose)
* Antilinear map
* System of linear equations, homogeneous and inhomogeneous, over- and underdetermined, upper and lower triangular matrix, LU decomposition, block matrix

**Lemmas, propositions and theorems:**

No proofs from Chapter 2. Nonetheless important:

* Theorem 2.11: Rank-nullity theorem $=$ dimension theorem

**For applications / exercises**

* Computing determinant of jacobian for integration measures
* Using Guassian elemination, Schur complements, Sherman-Morrison-Woodbury formula

## 3. Linear operators and eigenvalues

**Important concepts:**

* Projector and its relation to direct sum
* Polynomial of a matrix, annihilating polynomial, minimal annihilating polynomial
* Eigenvalue, eigenvector, eigenspace, spectrum, geometric multiplicity, algebraic multiplicity, characteristic polynomial, spectral decomposition ($=$ eigenvalue decomposition $=$ diagonalisation), spectral projector, diagonalisable versus defective matrix
* Invariant subspace, generalised eigenspace, Jordan normal form
* Companion matrix, diagonalised by Vandermonde matrix (not the generalisation for the case where eigenvalues coincide and the companion matrix is defective)
* Left eigenvector
* Function of an operator/matrix, matrix exponential

**Lemmas, propositions and theorems:**

Proofs to most lemmas and propositions are rather short and/or constructive and need to be known. Excluded from this is the construction of the Jordan normal form (Subsection 3.2.5).

Particularly important are:

* Proposition 3.3: projector versus direct sum
* Proposition 3.4
* Proposition 3.7: linear independence of different eigenvectors
* Theorem 3.12: Cayley Hamilton: passive knowledge of proof (understand but not reproduce)
* Proposition 3.15: invariant subspaces
* Proposition 3.19: Jordan decomposition of companion matrix: only passive knowledge

**For applications / exercises**

* Using Cayley-Hamilton theorem
* Computing eigenvalues and eigenvectors
* Computing a matrix function via eigenvalue or Jordan decomposition
* Solving an (autonomous linear) initial value problem or recurrence relation (in particular: Remark 3.62 and 3.69)

## 4. Norms and distances

**Important concepts:**

* Norm, Hölder p-norm
* Metric ($=$ distance), distance in a normed vector space, isometric map, limit of a sequence, continuity of a map, open subset, closed subset, closure, dense subset, separable metric space
* Equivalence of norms
* Cauchy sequence, metric complete vector space (=Banach space), closed subspace, dense subspace, complete set, absolutely converging series, Schauder basis
* Norm for linear maps, Frobenius norm, bounded map ($\equiv$ continuous map), subordinate norm, submultiplicative norm for operators, induced norm ($=$ operator norm), spectral radius
* Condition number

**Lemmas, propositions, theorems:**

No proofs from Chapter 4. Nonetheless important

* In a metric complete vector space (=Banach space)
  
  - Cauchy sequence $\iff$ convergent sequence
  - Absolutely convergent series $\implies$ convergent series

* Finite-dimensional vector spaces:
  
  - all norms are equivalent
  - always metric complete
  - finite-dimensional subspaces are always closed

* Linear maps are continuous if and only if bounded
* Bounded linear maps with operator norm are metric complete
* Operator norm is subordinate and submultiplicative
* Gelfand formula for spectral radius

## 5. Inner products and orthogonality

**Important concepts:**

* Bilinear form, sesquilinear form, quadratic form, symmetric and Hermitian, degenerate, positive (semi)definite versus indefinite, matrix congruence (= basis transform for bilinear forms)
* Inner product, standard/Euclidean inner product (in $\mathbb{C}^n$, in $\ell^2$, in $L^2$), metric/Gram matrix, inner product norm, Hilbert space
* Orthogonality, orthonormal set, orthogonal complement, orthogonal projection, orthogonal direct sum decomposition, orthogonal projector
* Orthonormal basis, Plancherel/Parseval identity, Gram-Schmidt orthonormalisation, QR decomposition
* (Anti)-isomorphism between dual space (= bounded linear functionals) and original Hilbert space
* Bounded linear maps, operator norm expressed using inner product, bounded linear maps have closed kernels
* Adjoint of a linear map, self-adjoint operators, isometric and unitary maps, normal operator
* Least squares solution, Moore-Penrose pseudoinverse

**Lemmas, propositions, theorems:**

Chapter 5 contains some of the most important theorems and propositions (Cauchy-Schwarz inequality, orthogonal projection theorem, expansion theorem, …). The corresponding proofs are often short and constructive. Shorter proofs need to be known actively, longer proofs passively.

Important active proofs:
* Theorem 5.2: Cauchy-Schwarz inequality
* Corollary 5.3: inner product norm
* Proposition 5.5: parallellogram law
* Proposition 5.7 + corollary 5.8: linear independence of orthogonal vectors
* Theorem 5.9: Pythagoras
* Lemma 5.21 and resulting from that Lemma 5.22: constructing orthogonal projection, Bessel's inequality
* Proposition 5.28: construction of the adjoint
* Proposition 5.29
* Proposition 5.30: relation between kernels and images of a map and its adjoint via orthogonal complements
* Proposition 5.31: characterisation of an linear isometry
* Proposition 5.34: characterisation of a normal operator
* Corollary 5.35 and 5.36: structure of eigenvalues and eigenvectors of normal operators
* Proposition 5.37, corollary 5.38 and 5.39: more properties of normal operators

Important passive proofs:
* Theorem 5.6: polarisation identity
* Theorem 5.13: orthogonal projection
* Theorem 5.16: orthogonal projectors and direct sum
* Theorem 5.23: expansion theorem
* Theorem 5.27: Riesz representation theorem

