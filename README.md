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

**For applications / exercises**

* Computing inner products, applying Gram-Schmidt (e.g. with custom inner products or between functions in a function space)
* Knowing and using the properties of Hermitian and normal operators

## 6. Unitary similarity and unitary equivalence

Chapter 6 discusses some of the common matrix decompositions in finite-dimensional linear algebra. A significant fraction of the material deals with algorithmic details and remarks, and serves as background for or addition to e.g. Python for Scientists. It is not of high importance for either the theory or exercises of this course. However, some of the main theorems regarding Schur decomposition and singular value decomposition are important.

**Important concepts:**

* Unitary and orthogonal group
* Discrete Fourier transform as unitary transformation and circulant matrices
* Schur decomposition and its relation to eigenvalue decomposition for normal matrices
* Canonical form of a bilinear map, intertia or signature
* Singular value decomposition: full, thin and compact; relation to rank, operator norm, condition number; applicability in the context of least squares solution (pseudo-inverse) and low rank approximations

**Lemmas, propositions, theorems:**

Important active proofs:

* Theorem 6.3: Schur decomposition
* Proposition 6.4: Normal matrices and Schur decomposition
* Proposition 6.6: Canonical form for congruence
* Proposition 6.8: Singular value decomposition (also remark 6.9 for its relation to eigenvalue / Schur decomposition)
* Proposition 6.9: SVD and rank
* Proposition 6.10: SVD and operator norm
* Proposition 6.11: SVD and Frobenius norm
* Proposition 6.12: SVD and condition number
* Theorem 6.14: SVD and low rank approximation in operator norm (not Theorem 6.15 in Frobenius norm)

Important passive proofs:

* Proposition 6.13: SVD and minimum norm least squares solution

## 7. Function spaces

Chapter 7 discusses some of the practicalities and issues related to working with functions spaces. It is mostly there to raise some awareness about these issues. The practical aspects are situated in Section 7.3 (properties of orthogonal polynomials) and parts of Section 7.4 (elementary properties of Fourier coefficients / Fourier series in 7.4.1, also Proposition 7.24 and 7.27 in section 7.4.3; relation with Discrete Fourier Transform in 7.4.5). All other (sub)sections are there as background material to raise some awareness about potential issues, and do not require an active understanding.

**Summary / Important concepts:**

* Function spaces can be given a proper norm and, for L^2, an inner product. The non-trivial step involves `identifying' functions that are equal almost everywhere, as one and the same (technically, working with equivalence classes of functions that are equal almost everywhere).
* The function space L^2 has interesting dense subspaces such as smooth or continuous functions, which are the ones we typically deal with
* Orthogonal polynomials and their general properties (recurrence relation and structure of roots=zeros)
* Fourier series to represent periodic functions; trigonometric polynomials are dense = Fourier modes are complete and thus the expansion theorem (Chapter 5 applies); unitary transformation between square integrable functions on an interval,and square summable sequence of Fourier coefficients (Parseval); slow convergence for discontinuities (Gibbs phenomenon); relation with discrete Fourier transform
* Unbounded operators are only defined on a subspace of the full Hilbert space = domain; the interesting class of operators are those for which that domain is still dense (=> no vector is orthogonal to the domain).
* Also the adjoint of an unbounded operator needs a domain, namely all $w$ for which we can make $\langle w, \hat{A} v \rangle = \langle \hat{A}^\dagger w, v\rangle$ work for all $v$ in the domain of $\hat{A}$. For differential operators, this relates to choosing boundary conditions for $v$ and $w$.
* Understanding the difference between being Hermitian/symmetric ($\langle w, \hat{A} v \rangle = \langle A w,v\rangle$ for all $v$ and $w$ in domain of $A$) and being self adjoint, again in the case of differential operators.

* The spectrum of an operator in an infinite-dimensional $\hat{A}$ Hilbert space consists of three parts:

  - The point spectrum: actual eigenvalues $\lambda$ with normalizable eigenvectors $v$: $\hat{A} v = \lambda v$
  - The continuous spectrum: values $\lambda$ for which we can find approximate eigenvectors, but no exact eigenvectors that we can properly normalize; we can find $v_\epsilon$ such that $\lVert \hat{A} v_\epsilon - \lambda v_\epsilon \rVert < \epsilon$ for all $\epsilon>0$, but the limit $\epsilon \to 0$ of $v_\epsilon$ is not well defined
  - The residual spectrum: very unintuitive and related to the fact that, on infinite-dimensional Hilbert spaces $\nu(\hat{A})$ (dimension of the kernel) and $\nu(\hat{A}^\dagger)$ do not need to be the same; the residual spectrum consists of values $\lambda$ for which no eigenvectors or approximate eigenvectors exist, but for which $\overline{\lambda}$ is in the point spectrum or continuous spectrum of $\hat{A}^\dagger$.

* For a self adjoint operator, the residual spectrum is empty, and the point spectrum and continuous spectrum only contain real numbers.

**Lemmas, propositions, theorems:**

Important active proofs:

* Proposition 7.8
* Proposition 7.9
* Proposition 7.10 (Cristoffel-Darboux formule zeker niet vanbuiten kennen)
* Proposition 7.11

For the next proposition on properties of the Fourier coefficients: proving the relation that the Fourier coefficients satisfy is an easy calculation, which you should be able to actively do. You do not need to know the technical conditions under which these manipulations are allowed, and which requirements they impose on the function $f$ or the Fourier coefficients $(\widehat{f}_k)$.

* Proposition 7.17
* Proposition 7.19
* Proposition 7.24
* Proposition 7.27

**For applications / exercises**

* Computing inner products and applying Gram-Schmidt to a small set of functions.
* Deriving relations of specific families orthogonal polynomials, e.g. deriving orthonormalization relation or recurrence relation from generating function.
* Computing simple Fourier coefficients using the elementary properties

## 8. Linear differential operators

This chapter provides an in-depth study of differential operators, and their role in the study of linear differential equations with boundary conditions.

**Important concepts:**

* Homogeneous and inhomogeneous (linear) differential equation, homogeneous and inhomogeneous boundary conditions, formal adjoint of a linear differential operator, Sturm-Liouville operator, separated boundary conditions.
* A $p$th order differential equation is well balanced (likely to have a solution that exists and is unique for any right hand side) if it has exactly $p$ boundary conditions; this is a necessary but not sufficient condition. Having $p$ boundary conditions is also necessary (but not sufficient) to be a self-adjoint operator.
* Initial value problem, fundamental matrix solution, time-ordered exponential, Wronskian, Floquet theorem
* Boundary value problem, Dirichlet and Neumann conditions for second order problems, Green's function, Green's operator as inverse of differential operator
* Sturm-Liouville eigenvalue problems: a regular Sturm-Liouville operator admits a spectral decompositoin where the eigenvectors provide a complete orthonormal basis for the Hilbert space, and can thus be used to compute e.g. the exponential or inverse (= Green's function) (Remark 8.31 and 8.32)

**Lemmas, propositions, theorems:**

Important active proofs:

* Proposition 8.3 and Corollary 8.4 (it is sufficient if you can prove this for the case $p=2$)
* Proposition 8.5 and Corollary 8.6
* Constructing the bilinear concomitant of the Sturm-Liouville operator (eq 8.30)
* Verifying that it is self-adjoint with respect to separated or periodic boundary conditions.
* Lemma 8.12
* Theorem 8.13
* Proposition 8.14
* Proposition 8.15
* Proposition 8.16
* Proposition 8.17
* Theorem 8.18

Important passive proofs:

* Proposition 8.9 and its generalisation to the inhomogeneous case, Proposition 8.19

No theorems beyond subsection 8.2.4; you need to understand examples and use the concepts (see above) in exercises.

**For applications / exercises**

* Be able to apply the Frobenius method, in particular for second order problems (Remarks 8.25 and 8.26). In particular, this requires that you can recognize/identify a regular singular point, and that you can derive the indicial equation. Also applying the simpler version (a Taylor series ansatz) when there is no singular point (coefficient of highest derivative is strictly positive)

* Understand how to use a Green's function, and how to construct it in the case of a second order problem with separated boundary conditions (p305). Actually, no exercises were made on this in class due to lack of time, so there is only the example in the theory section, and the exercise exam will not contain questions about this.

## 9. Fourier transforms and distributions

This chapter introduces the Fourier transform, first in the classical sense (as a unitary operator on $L^2(\mathbb{R})$). Then, the main concepts from the theory of distributions is introduced, which provides the mathematical framework for working with 'generalised functions' such as the Dirac delta 'function'. Within the setting of distributions, we can significantly extend the concept of derivatives, limits, series and Fourier transforms beyond their classical meaning, and we provide several examples of this. Finally, we revisit the different types of Fourier transforms and introduce a fourth type that combines very naturally with the three types that we have already seen. Then, we find various relations between these different types, in which the use of distributions plays a prominent role.

Sections 9.2.11 and 9.3.4 have only been briefly covered in class, 9.3.5 and 9.4 have not been covered at all. None of those need to be known for the exam.

**Important concepts:**

* Fourier transform, convolution, Fourier transform as unitary operator on $L^2(\mathbb{R})$, Parseval and Plancherel relation for Fourier transform
* Fourier transform of Gaussian distribution, characteristic function
* Test function, compact support, distribution, regular versus singular distribution, Dirac-delta distribution (and its derivatives), Heaviside function/distribution, Cauchy principal value, distributional derivative, distributional limit, distributional Fourier series and Fourier transform
* Fourier transforms on different domains: discrete Fourier transform, Fourier series, discrete-time Fourier transform, (continuous-time) Fourier transform.
* Sampling, Nyquist rate, reconstruction via sinc (Whittaker–Shannon interpolation formula). 

**Lemmas, propositions, theorems:**

For most of this chapter, no technical aspects of the theorems or proofs need to be known.

What you need to know is:

* Properties of Fourier transform: Proposition 9.1, Theorem 9.3 (convolution), Proposition 9.4 (derivative) (no technical requirements on $f$ or $\widehat{f}$)
* Computing Fourier transform of Gaussian distribution, passive understanding of proving central limit theorem
* Using the definition of translation, scaling, derivative, coordinate transform, limit of distributions on examples:
  Examples 9.9, 9.10, 9.11, 9.12, 9.13, 9.14, 9.15, 9.16, 9.17, 9.19
* Theorem 9.15 (given the structure of the complex logarithm and the distributional derivative of $\log |x|$)
* Using Proposition 9.16 (Dirac sequence) on examples such as Example 9.20, 9.21, 9.22
* Passive understanding of 9.2.5 (Cauchy Principal Value) and 9.2.9 (Dirac Comb distribution)
* Fourier transform of distributions on examples: Example 9.23, 9.24
* Theorem 9.18 (Poisson summation formula) given Dirac Comb distribution
* Sampling: proving proposition 9.20, corollary 9.21 and proposition 9.22.

**For applications / exercises**

Unfortunately, we did not have time for further exercises on this chapter. 