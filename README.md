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
