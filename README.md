<!--
author:   Jutho Haegeman

email:    jutho.haegeman@ugent.be

language: en

-->

[![LiaScript](https://raw.githubusercontent.com/LiaScript/LiaScript/master/badges/course.svg)](https://LiaScript.github.io/course/?https://raw.githubusercontent.com/Jutho/VectorFunctionSpacesSummary/main/README.md)

# Vector and Function Spaces - Summary

For every chapter of the lecture notes "Vector- and function spaces", this document lists the following items:

1. A brief summary

2. Not covered in class: sections which we have not or only briefly covered, and are not part of the exam (neither theory nor exercise). This corresponds to (subsections) that (should) have a ($\star$) indicator in the lecture notes.

3. Important concepts: these mostly correspond to definitions, sometimes the definition is hidden inside a proposition or theorem. You do not need to know the literal definition, but you should be able to understand and use this term or concept correctly, both for theory and exercise.

4. Important theorems and propositions: list of proofs that need to be known (actively) or understood (passively) for the theory exam

5. Additional topics for applications / exercises: what you need to know for the exercises, in particular if it is additional material that you do not need to know for theory

## 1. Elementary algebraic structures

**Summary:**

This chapter lists the basic mathematical structures and terminology that we will need and use, and doing so, also specifies the convention and notations that we will use for those. This should be almost completely repition; there will be no direct theory questions about this, but you should of course be able to use and understand this terminology

**Not convered in class:**

Section 1.2.4 was barely covered in class. The definition of the kernel of a group homomorphism generalised to the kernel of a vector space homomorphism (=linear map). The concept of a normal subgroup and a quotient group also generalise to that of a subspace and a quotient space, because a vector space is an abelian group with respect to vector addition. We have not at all discussed exact sequences.

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

**Summary:**

This chapter discusses general properties of vector space homomorphisms = linear maps, and in particular their representation as matrices in the case of finite-dimensional vector spaces. Most of this chapter should be repition, some parts are probably new (linear functionals and dual space, antilinear maps, determinants and inverses of block matrices). You need to be able to use the concepts from this chapter (for exercises and further theory), but there will not be any direct questions about it on the theory exam.

**Not covered in class:**

Section 2.5.4 (Double dual space) and Section 2.7.3 (Conjugate vector space) have not been covered in class.

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

**Summary:**

This chapter discusses general properties of vector space endomorphisms = linear operators = linear maps from a vector space to itself. Operators can be composed with themselves, giving rise to powers and polymials. This is important to then also introduce eigenvalues and eigenvectors, and finally, to generalise arbitrary scalar functions ($f:\mathbb{C}\to \mathbb{C}$) to operators/square matrices. Parts of this will be repition (projectors, eigenvalues, ...), other parts will be new (generalised eigenspaces, Jordan form, functions of operators).

**Not covered in class:**

Section 3.3.3 (Derivatives of matrix functions) was not covered in class. Section 3.4 was covered, but rather quickly, and mainly with how it should be used in applications and exercises (see below). There will be no theory questions from Section 3.4.

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

**Summary:**

This chapter introduces the concept of a norm, and then discusses at length how this can be used to make sense of limits of sequences of vectors, mostly in infinite-dimensional vector spaces, where there are some surprises and not everything is very intuitive. Then, we discuss the implications for linear maps on such normed vector spaces. Most of the proofs are very technical; none of them need to be known for the exam. However, some important concepts are important for the remainder of the course (especially the theory-part).

**Not covered in class:**

Section 4.4.3 was not covered in class. If you want some theoretical background about Markov chains and the Google PageRank algorithm, feel free to read ;-).

**Important concepts:**

* Norm, Hölder p-norm
* Metric ($=$ distance), distance in a normed vector space, isometric map, limit of a sequence, continuity of a map, open subset, closed subset, closure, dense subset, separable metric space
* Equivalence of norms
* Cauchy sequence, metric complete vector space (=Banach space), closed subspace, dense subspace, complete set, absolutely converging series, Schauder basis
* Norm for linear maps, Frobenius norm, bounded map ($\equiv$ continuous map), subordinate norm, submultiplicative norm for operators, induced norm ($=$ operator norm), spectral radius
* Condition number

**Lemmas, propositions, theorems:**

No proofs from Chapter 4. The following insights are nonetheless important:

* In a metric complete vector space (=Banach space)

  - Cauchy sequence $\iff$ convergent sequence
  - Absolutely convergent series $\implies$ convergent series

* Finite-dimensional vector spaces:

  - all norms are equivalent
  - always metric complete
  - finite-dimensional subspaces are always closed

* Infinite-dimensional proper subspaces of a Banach space can be dense: there is no intuitive or visual way to interpret this concept, there are vectors which are not in this subspace, but nonetheless they are arbitrary close to it (measured using the norm of the surrounding vector space). In a Hilbert space (next chapter), this means in particular that no vector can be made orthogonal to a dense subspace, since a unit vector that is orthogonal to a subspace has at least distance 1 to any vector in that subspace.

* Linear maps are continuous if and only if bounded
* Bounded linear maps with operator norm are metric complete
* Operator norm is subordinate and submultiplicative
* Gelfand formula for spectral radius

**For applications / exercises**
* Computing a norm. But this chapter is mostly there for theoretical purposes.

## 5. Inner products and orthogonality

**Summary:**

This important chapter introduces the concept of an inner product and the structures that follows from it, notabily, the concept of orthogonality and orthogonal projections. In particular, working with a basis in an infinite-dimensional vector space becomes more easy with an inner product (and associated norm) instead of "just a norm", when using orthogonality. In a Hilbert space (=metric complete inner product space), a set of vectors that is complete (the linear span defines a dense subspace) can be turned into an orthonormal set (using Gram-Schmidt) which then defines a basis (expansion theorem). Linear maps and linear operators between Hilbert spaces can also have more structure. Bounded linear functionals (elements from the dual space) are one-to-one associated with vectors in the primal space via the inner product. Bounded linear maps have adjoints. Linear operators can satisfy relations with their adjoint, e.g. they can be equal (self-adjoint), the adjoint can be the inverse (unitary) or they can commute with the adjoint (normal), which then imposes particular constraints on the spectrum and the eigenvectors. On the practical side, the orthogonal projection allows to construct least square solutions to overdetermined systems.

**Not covered in class:**

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

**Summary:**

Chapter 6 discusses some of the common matrix decompositions in finite-dimensional linear algebra. A significant fraction of the material deals with algorithmic details and remarks, and serves as background for or addition to e.g. Python for Scientists. It is not of high importance for either the theory or exercises of this course. However, some of the main theorems regarding Schur decomposition and singular value decomposition are important.

**Not covered in class:**

Subsections 6.4.4 (practical considerations regarding Schur decomposition) and 6.6.6 (polar decomposition) as well as Section 6.7 (Krylov methods) were not covered in class.

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

**Summary / Important concepts:**

Chapter 7 discusses some of the practicalities and issues related to working with functions spaces. It is mostly there to raise some awareness about these issues. The practical aspects are situated in Section 7.3 (properties of orthogonal polynomials) and parts of Section 7.4 (elementary properties of Fourier coefficients / Fourier series in 7.4.1, also Proposition 7.24 and 7.27 in section 7.4.3; relation with Discrete Fourier Transform in 7.4.5). All other (sub)sections are there as background material to raise some awareness about potential issues, and do not require an active understanding.

* Function spaces can be given a proper norm and, for $L^2$, an inner product. The non-trivial step involves `identifying' functions that are equal almost everywhere, as one and the same (technically, working with equivalence classes of functions that are equal almost everywhere).
* The function space $L^2$ has interesting dense subspaces such as smooth or continuous functions, which are the ones we typically deal with
* Orthogonal polynomials and their general properties (recurrence relation and structure of roots=zeros)
* Fourier series to represent periodic functions; trigonometric polynomials are dense = Fourier modes are complete and thus the expansion theorem (Chapter 5 applies); unitary transformation between square integrable functions on an interval,and square summable sequence of Fourier coefficients (Parseval); slow convergence for discontinuities (Gibbs phenomenon); relation with discrete Fourier transform
* Unbounded operators are only defined on a subspace of the full Hilbert space = domain; the interesting class of operators are those for which that domain is still dense (=> no vector is orthogonal to the domain).
* Also the adjoint of an unbounded operator needs a domain, namely all $w$ for which we can make $\langle w, \hat{A} v \rangle = \langle \hat{A}^\dagger w, v\rangle$ work for all $v$ in the domain of $\hat{A}$. For differential operators, this relates to choosing boundary conditions for $v$ and $w$.
* Understanding the difference between being Hermitian/symmetric ($\langle w, \hat{A} v \rangle = \langle A w,v\rangle$ for all $v$ and $w$ in domain of $A$) and being self adjoint, again in the case of differential operators.

* The spectrum of an operator $\hat{A}$ in an infinite-dimensionalHilbert space consists of three parts:

  - The point spectrum: actual eigenvalues $\lambda$ with normalizable eigenvectors $v$: $\hat{A} v = \lambda v$
  - The continuous spectrum: values $\lambda$ for which we can find approximate eigenvectors, but no exact eigenvectors that we can properly normalize; we can find $v_\epsilon$ such that $\lVert \hat{A} v_\epsilon - \lambda v_\epsilon \rVert < \epsilon$ for all $\epsilon>0$, but the limit $\epsilon \to 0$ of $v_\epsilon$ is not well defined
  - The residual spectrum: very unintuitive and related to the fact that, on infinite-dimensional Hilbert spaces $\nu(\hat{A})$ (dimension of the kernel) and $\nu(\hat{A}^\dagger)$ do not need to be the same; the residual spectrum consists of values $\lambda$ for which no eigenvectors or approximate eigenvectors exist, but for which $\overline{\lambda}$ is in the point spectrum or continuous spectrum of $\hat{A}^\dagger$.

* For a self adjoint operator, the residual spectrum is empty, and the point spectrum and continuous spectrum only contain real numbers.

**Not covered in class:**

Subsection 7.3.6 (Gaussian quadrature) was only briefly covered; there will be no questions about this on the exam.

**Lemmas, propositions, theorems:**

Important active proofs:

* Proposition 7.8
* Proposition 7.9
* Proposition 7.10 (Cristoffel-Darboux formula does not need to be known by heart)
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

**Summary:**

This chapter provides an in-depth study of differential operators, and their role in the study of linear differential equations with boundary conditions. We discuss how to construct the "formal adjoint" of a differential operator using partial integration, and how it relates to the boundary condition to actually construct the adjoint. We discuss how to decompose the solution of a differential equation in different parts, and how to study the existence and uniqueness of these different contributions. The adjoint plays a role via the Fredholm alternative theorem. Also, we discuss when a second order differential operator is self-adjoint (known as a Sturm-Liouville operator).

To better understand the role of boundary conditions, we take an extended detour via initial value problems, for which we can formally construct the solution (as a path ordered exponential), and we also discuss practical recipes (via a Taylor expansion or a generalisation thereof, known as Frobenius method). We find that we need $p$ boundary conditions for a $p$th order differential equation to be well balanced. Finally then, we can construct the solution to the inhomogeneous differential equation (with homogeneous boundary conditions) using the Green's function.

We then move on using differential operators in eigenvalue problems. The purpose thereof is for solving higher-dimensional partial differential equations, for which the spectral decomposion of a differential operator (if it exists) turns out to be very useful.

**Not covered in class:**

Subsections 8.2.5, 8.2.6 and 8.2.7 were not really covered, only the resulting method of Frobenius (at the end of 8.2.7), which is useful for the exercises. The example application in 8.2.8 (Bessel function) was also not covered.

Subsections 8.4.3, 8.4.4 and section 8.5 were only briefly discussed in class.

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

**Summary:**

This chapter introduces the Fourier transform, first in the classical sense (as a unitary operator on $L^2(\mathbb{R})$). Then, the main concepts from the theory of distributions is introduced, which provides the mathematical framework for working with 'generalised functions' such as the Dirac delta 'function'. Within the setting of distributions, we can significantly extend the concept of derivatives, limits, series and Fourier transforms beyond their classical meaning, and we provide several examples of this. Finally, we revisit the different types of Fourier transforms and introduce a fourth type that combines very naturally with the three types that we have already seen. Then, we find various relations between these different types, in which the use of distributions plays a prominent role.

**Not covered in class:**

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

**For applications / exercises:**

Unfortunately, we did not have time for further exercises on this chapter.

## Glossary

| English                                                              | Dutch                                                                                       | Typical notation                                                     | Wikipedia                                                                                    |
|----------------------------------------------------------------------|---------------------------------------------------------------------------------------------|----------------------------------------------------------------------|----------------------------------------------------------------------------------------------|
| Set                                                                  | Verzameling                                                                                 | $A, B, \dots$                                                        | <a href="https://en.wikipedia.org/wiki/Set_(mathematics)">`Set_(mathematics)`</a>                                              |
| Subset                                                               | Deelverzameling                                                                             | $B \subseteq A$                                                      | <a href="https://en.wikipedia.org/wiki/Subset">`Subset`</a>                                                         |
| Proper subset                                                        | Eigenlijke deelverzameling                                                                  |                                                                      |                                                                                              |
| Disjoint subsets                                                     | Disjuncte verzamelingen                                                                     | $A \cap B = \{\}$                                                    | <a href="https://en.wikipedia.org/wiki/Disjoint_sets">`Disjoint_sets`</a>                                                  |
| Cartesian product                                                    | Cartesisch product                                                                          | $A \times B$                                                         | <a href="https://en.wikipedia.org/wiki/Cartesian_product">`Cartesian_product`</a>                                              |
| Cardinality                                                          | Kardinaliteit                                                                               | $\lvert A \rvert$                                                    | <a href="https://en.wikipedia.org/wiki/Cardinality">`Cardinality`</a>                                                    |
| Countable                                                            | Aftelbaar                                                                                   |                                                                      | <a href="https://en.wikipedia.org/wiki/Countable_set">`Countable_set`</a>                                                  |
| Uncountable                                                          | Onaftelbaar                                                                                 |                                                                      | <a href="https://en.wikipedia.org/wiki/Uncountable_set">`Uncountable_set`</a>                                                |
| Map                                                                  | Afbeelding                                                                                  | $\varphi:A \to B$                                                    | <a href="https://en.wikipedia.org/wiki/Map_(mathematics)">`Map_(mathematics)`</a>                                              |
| Identity map                                                         | Identiteitsafbeelding                                                                       |                                                                      | <a href="https://en.wikipedia.org/wiki/Identity_function">`Identity_function`</a>                                              |
| Equivalence relation                                                 | Equivalentierelatie                                                                         | $a \sim b$                                                           | <a href="https://en.wikipedia.org/wiki/Equivalence_relation">`Equivalence_relation`</a>                                           |
| Equivalence class                                                    | Equivalentieklasse                                                                          | $[a]$                                                                | <a href="https://en.wikipedia.org/wiki/Equivalence_class">`Equivalence_class`</a>                                              |
| Partial order                                                        | Partiële ordening                                                                           | $a \preccurlyeq b$                                                   | <a href="https://en.wikipedia.org/wiki/Partially_ordered_set">`Partially_ordered_set`</a>                                          |
| Binary operation                                                     | Binaire bewerking                                                                           | $\cdot, \circ, +, \times, \dots$                                     | <a href="https://en.wikipedia.org/wiki/Binary_operation">`Binary_operation`</a>                                               |
| Group                                                                | Groep                                                                                       | $G, H, \dots$                                                        | <a href="https://en.wikipedia.org/wiki/Group_(mathematics)">`Group_(mathematics)`</a>                                            |
| Abelian group                                                        | Abelse groep                                                                                |                                                                      | <a href="https://en.wikipedia.org/wiki/Abelian_group">`Abelian_group`</a>                                                  |
| Subgroup                                                             | Deelgroep                                                                                   | $H \subseteq G$                                                      | <a href="https://en.wikipedia.org/wiki/Subgroup">`Subgroup`</a>                                                       |
| Homomorphism                                                         | Homomorfisme                                                                                | $\varphi, \chi, \dots$                                               | <a href="https://en.wikipedia.org/wiki/Homomorphism">`Homomorphism`</a>                                                   |
| Isomorphism                                                          | Isomorphisme                                                                                | $\phi$                                                               | <a href="https://en.wikipedia.org/wiki/Isomorphism">`Isomorphism`</a>                                                    |
| Isomorphic                                                           | Isomorf                                                                                     | $A \cong B$                                                          |                                                                                              |
| Automorphism group                                                   | Automorfismegroup                                                                           | $\mathrm{Aut}(A)$                                                    | <a href="https://en.wikipedia.org/wiki/Automorphism_group">`Automorphism_group`</a>                                             |
| Involution                                                           | Involutie                                                                                   |                                                                      | <a href="https://en.wikipedia.org/wiki/Involution_(mathematics)">`Involution_(mathematics)`</a>                                       |
| Permutation group                                                    | Permutatiegroep                                                                             | $S_n$                                                                | <a href="https://en.wikipedia.org/wiki/Permutation_group">`Permutation_group`</a>                                              |
| Group action                                                         | Groepswerking, groepsactie                                                                  | $g \triangleright a$, $a \triangleleft g$                            | <a href="https://en.wikipedia.org/wiki/Group_action">`Group_action`</a>                                                   |
| Ring                                                                 | Ring                                                                                        | $R$, $\mathbb{Z}$                                                    | <a href="https://en.wikipedia.org/wiki/Ring_(mathematics)">`Ring_(mathematics)`</a>                                             |
| Field                                                                | Veld                                                                                        | $\mathbb{F}$, $\mathbb{R}$, $\mathbb{C}$                             | <a href="https://en.wikipedia.org/wiki/Field_(mathematics)">`Field_(mathematics)`</a>                                            |
| Vector space, linear space                                           | Vectorruimte, lineaire ruimte                                                               | $V, W, \dots$                                                        | <a href="https://en.wikipedia.org/wiki/Vector_space">`Vector_space`</a>                                                   |
| (Linear) subspace                                                    | (Lineaire) deelruimte                                                                       | $W \preccurlyeq V$                                                   | <a href="https://en.wikipedia.org/wiki/Linear_subspace">`Linear_subspace`</a>                                                |
| Linear map                                                           | Lineaire afbeelding                                                                         | $\hat{A}, \hat{B}, \dots$                                            | <a href="https://en.wikipedia.org/wiki/Linear_map">`Linear_map`</a>                                                     |
| (Linear) operator                                                    | (Lineaire) operator                                                                         |                                                                      |                                                                                              |
| Linear transformation                                                | Lineaire transformatie                                                                      | $\hat{T}, \dots$                                                     |                                                                                              |
| Linear combination                                                   | Lineaire combinatie                                                                         |                                                                      | <a href="https://en.wikipedia.org/wiki/Linear_combination">`Linear_combination`</a>                                             |
| (Linear) span                                                        | Lineair omhulsel, lineair opspansel                                                         | $\mathrm{span}(S), \mathbb{F} S$                                     | <a href="https://en.wikipedia.org/wiki/Linear_span">`Linear_span`</a>                                                    |
| Basis                                                                | Basis                                                                                       | $B_V$                                                                | <a href="https://en.wikipedia.org/wiki/Basis_(linear_algebra)">`Basis_(linear_algebra)`</a>                                         |
| Dimension                                                            | Dimensie                                                                                    | $\dim(V)$                                                            | <a href="https://en.wikipedia.org/wiki/Dimension_(vector_space)">`Dimension_(vector_space)`</a>                                       |
| Direct sum                                                           | Directe som                                                                                 | $V \oplus W$                                                         | <a href="https://en.wikipedia.org/wiki/Direct_sum">`Direct_sum`</a>                                                     |
| Quotient space                                                       | Quotiëntruimte                                                                              |                                                                      | <a href="https://en.wikipedia.org/wiki/Quotient_space_(linear_algebra)">`Quotient_space_(linear_algebra)`</a>                                |
| Algebra                                                              | Algebra                                                                                     |                                                                      | <a href="https://en.wikipedia.org/wiki/Algebra_over_a_field">`Algebra_over_a_field`</a>                                           |
| Identity operator                                                    | Identiteitsoperator                                                                         | $\hat{1}, \hat{1}_V$                                                 | <a href="https://en.wikipedia.org/wiki/Identity_function">`Identity_function`</a>                                              |
| Kernel, null space                                                   | Kern, nulruimte                                                                             | $\ker(\hat{A})$                                                      | <a href="https://en.wikipedia.org/wiki/Kernel_(linear_algebra)">`Kernel_(linear_algebra)`</a>                                        |
| Image, range                                                         | Beeldruimte, bereik                                                                         | $\mathrm{im}(\hat{A})$                                               | <a href="https://en.wikipedia.org/wiki/Image_(mathematics)">`Image_(mathematics)`</a>                                            |
| Rank                                                                 | Rang                                                                                        | $\rho(\hat{A})$                                                      | <a href="https://en.wikipedia.org/wiki/Rank_(linear_algebra)">`Rank_(linear_algebra)`</a>                                          |
| Nullity                                                              | Nulliteit                                                                                   | $\nu(\hat{A})$                                                       | <a href="https://en.wikipedia.org/wiki/Kernel_(linear_algebra)#nullity">`Kernel_(linear_algebra)#nullity`</a>                                |
| Matrix                                                               | Matrix                                                                                      | $\mathsf{A}, \mathsf{B}, \dots$                                      | <a href="https://en.wikipedia.org/wiki/Matrix_(mathematics)">`Matrix_(mathematics)`</a>                                           |
| Identity matrix                                                      | Identiteitsmatrix                                                                           | $\mathsf{I}, \mathsf{I}_n$                                           | <a href="https://en.wikipedia.org/wiki/Identity_matrix">`Identity_matrix`</a>                                                |
| Transpose                                                            | Getransponeerde                                                                             | $\mathsf{A}^{\mathsf{T}}$                                            | <a href="https://en.wikipedia.org/wiki/Transpose">`Transpose`</a>                                                      |
| Hermitian conjugate, conjugate transpose, Hermitian adjoint, adjoint | Hermitisch geconjugeerde, geconjugeerd getransponeerde, Hermitisch toegevoegde, toegevoegde | $\mathsf{A}^{\mathsf{H}}$                                            | <a href="https://en.wikipedia.org/wiki/Conjugate_transpose">`Conjugate_transpose`</a>                                            |
| Symmetric                                                            | Symmetrisch                                                                                 |                                                                      | <a href="https://en.wikipedia.org/wiki/Symmetric_matrix">`Symmetric_matrix`</a>                                               |
| Skew-symmetric, antisymmetric                                        | Scheefsymmetrisch, antisymmetrisch                                                          |                                                                      | <a href="https://en.wikipedia.org/wiki/Skew-symmetric_matrix">`Skew-symmetric_matrix`</a>                                          |
| Hermitian                                                            | Hermitisch                                                                                  |                                                                      | <a href="https://en.wikipedia.org/wiki/Hermitian_matrix">`Hermitian_matrix`</a>                                               |
| Skew-Hermitian, anti-Hermitian                                       | Scheefhermitisch, antihermitisch                                                            |                                                                      | <a href="https://en.wikipedia.org/wiki/Skew-Hermitian_matrix">`Skew-Hermitian_matrix`</a>                                          |
| Row space, column space                                              | Rijruimte, kolomruimte                                                                      |                                                                      | <a href="https://en.wikipedia.org/wiki/Row_and_column_spaces">`Row_and_column_spaces`</a>                                          |
| Trace                                                                | Spoor                                                                                       | $\mathrm{tr}(\mathsf{A})$                                            | <a href="https://en.wikipedia.org/wiki/Trace_(linear_algebra)">`Trace_(linear_algebra)`</a>                                         |
| Determinant                                                          | Determinant                                                                                 | $\det(\mathsf{A})$                                                   | <a href="https://en.wikipedia.org/wiki/Determinant">`Determinant`</a>                                                    |
| Adjugate                                                             | Geadjugeerde                                                                                | $\mathrm{adj}(\mathsf{A})$                                           | <a href="https://en.wikipedia.org/wiki/Adjugate_matrix">`Adjugate_matrix`</a>                                                |
| Matrix inverse, invertible matrix                                    | Matrixinverse, inverteerbare matrix                                                         | $\mathsf{A}^{-1}$                                                    | <a href="https://en.wikipedia.org/wiki/Invertible_matrix">`Invertible_matrix`</a>                                              |
| Singular, degenerate                                                 | Singulier, ontaard                                                                          |                                                                      |                                                                                              |
| General linear group                                                 | Algemene lineaire group                                                                     | $\mathrm{GL}(V), \mathrm{GL}(n; \mathbb{F})$                         | <a href="https://en.wikipedia.org/wiki/General_linear_group">`General_linear_group`</a>                                           |
| Special linear group                                                 | Speciale lineaire group                                                                     | $\mathrm{SL}(n; bbF)$                                                | <a href="https://en.wikipedia.org/wiki/Special_linear_group">`Special_linear_group`</a>                                           |
| Similar matrices                                                     | Gelijkvormige matrices                                                                      | $\tilde{\mathsf{A}} = \mathsf{T}^{-1} \mathsf{A} \mathsf{T}$         | <a href="https://en.wikipedia.org/wiki/Matrix_similarity">`Matrix_similarity`</a>                                              |
| Linear functional, linear form, dual vector, covector                | Lineaire functionaal, lineaire vorm, duale vector, covector                                 | $\xi, \chi, \theta,\dots$                                            | <a href="https://en.wikipedia.org/wiki/Linear_form">`Linear_form`</a>                                                    |
| Dual (vector) space                                                  | Duale (vector)ruimte                                                                        | $V^\ast$                                                             | <a href="https://en.wikipedia.org/wiki/Dual_space">`Dual_space`</a>                                                     |
| Dual linear map                                                      | Duale lineaire afbeelding                                                                   | $\hat{A}^\ast$                                                       | <a href="https://en.wikipedia.org/wiki/Transpose_of_a_linear_map">`Transpose_of_a_linear_map`</a>                                      |
| Antilinear map                                                       | Antilineaire afbeelding                                                                     |                                                                      | <a href="https://en.wikipedia.org/wiki/Antilinear_map">`Antilinear_map`</a>                                                 |
| System of linear equations, linear system                            | Stelsel van lineaire vergelijkingen, linear systeem                                         |                                                                      | <a href="https://en.wikipedia.org/wiki/System_of_linear_equations">`System_of_linear_equations`</a>                                     |
| Overdetermined                                                       | Overgedetermineerd                                                                          |                                                                      | <a href="https://en.wikipedia.org/wiki/Overdetermined_system">`Overdetermined_system`</a>                                          |
| Underdetermined                                                      | Ondergedetermineerd                                                                         |                                                                      | <a href="https://en.wikipedia.org/wiki/Underdetermined_system">`Underdetermined_system`</a>                                         |
| Homogeneous, inhomogeneous                                           | Homogeen, inhomogeen                                                                        |                                                                      |                                                                                              |
| Triangular matrix, upper or lower triangular                         | Driehoeksmatrix, bovendriehoeksmatrix of benedendriehoeksmatrix                             | $\mathsf{T}$                                                         | <a href="https://en.wikipedia.org/wiki/Triangular_matrix">`Triangular_matrix`</a>                                              |
| Polynomial                                                           | polynoom, veelterm                                                                          | $p$                                                                  | <a href="https://en.wikipedia.org/wiki/Polynomial">`Polynomial`</a>                                                     |
| Monic polynomial                                                     | Monische veelterm                                                                           |                                                                      | <a href="https://en.wikipedia.org/wiki/Monic_polynomial">`Monic_polynomial`</a>                                               |
| Idempotent                                                           | Idempotent                                                                                  |                                                                      | <a href="https://en.wikipedia.org/wiki/Idempotence">`Idempotence`</a>                                                    |
| Nilpotent                                                            | Nilpotent                                                                                   |                                                                      | <a href="https://en.wikipedia.org/wiki/Nilpotent">`Nilpotent`</a>                                                      |
| Projection operator                                                  | Projectie-operator                                                                          | $\hat{P}, \hat{Q}$                                                   | <a href="https://en.wikipedia.org/wiki/Projection_(linear_algebra)">`Projection_(linear_algebra)`</a>                                    |
| Minimal (annihilating) polynomial                                    | Minimale (annihilerende) polynoom, minimaalveelterm                                         | $m_{\hat{A}}(z)$                                                     | <a href="https://en.wikipedia.org/wiki/Minimal_polynomial_(linear_algebra)">`Minimal_polynomial_(linear_algebra)`</a>                            |
| Characteristic polynomial                                            | Karakteristieke polynoom, karakteristieke veelterm                                          | $k_{\hat{A}}(z)$                                                     | <a href="https://en.wikipedia.org/wiki/Characteristic_polynomial">`Characteristic_polynomial`</a>                                      |
| Eigenvalue, eigenvector, eigenspace                                  | Eigenwaarde, eigenvector, eigenruimte                                                       | $\lambda$, $v_\lambda$, $V_\lambda$                                  | <a href="https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors">`Eigenvalues_and_eigenvectors`</a>                                   |
| Algebraic multiplicity                                               | Algebraische multipliciteit                                                                 | $q_\lambda$                                                          |                                                                                              |
| Geometric multiplicity                                               | Meetkundige multipliciteit                                                                  | $r_\lambda$                                                          |                                                                                              |
| Spectrum                                                             | Spectrum                                                                                    | $\sigma_{\hat{A}}$                                                   | <a href="https://en.wikipedia.org/wiki/Spectral_theory">`Spectral_theory`</a>                                                |
| Spectral decomposition, eigendecomposition                           | Spectrale decompositie, eigenwaarde-ontbidinding                                            | $\mathsf{A} \mathsf{V} = \mathsf{V} \mathsf{D}$                      | <a href="https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix">`Eigendecomposition_of_a_matrix`</a>                                 |
| Diagonalizable matrix                                                | Diagonaliseerbare matrix                                                                    |                                                                      | <a href="https://en.wikipedia.org/wiki/Diagonalizable_matrix">`Diagonalizable_matrix`</a>                                          |
| Defective matrix                                                     | Defecte matrix (?), niet-diagonaliseerbare matrix                                           |                                                                      | <a href="https://en.wikipedia.org/wiki/Defective_matrix">`Defective_matrix`</a>                                               |
| Invariant subspace                                                   | Invariante deelruimte                                                                       |                                                                      | <a href="https://en.wikipedia.org/wiki/Invariant_subspace">`Invariant_subspace`</a>                                             |
| Generalized eigenvector, generalized eigenspace                      | Veralgemeende eigenvector, veralgemeende eigenruimte (?)                                    | $u_{\lambda,k}$, $U_\lambda$                                         | <a href="https://en.wikipedia.org/wiki/Generalized_eigenvector">`Generalized_eigenvector`</a>                                        |
| Jordan normal form, Jordan canonical form                            | Jordan-normaalvorm                                                                          | $\mathsf{A} \mathsf{V} = \mathsf{V} \mathsf{J}$                      | <a href="https://en.wikipedia.org/wiki/Jordan_normal_form">`Jordan_normal_form`</a>                                             |
| Companion matrix                                                     | Begeleidende matrix (?)                                                                     | $\mathsf{C}_p$                                                       | <a href="https://en.wikipedia.org/wiki/Companion_matrix">`Companion_matrix`</a>                                               |
| Norm                                                                 | Norm                                                                                        | $\lVert v \rVert$                                                    | <a href="https://en.wikipedia.org/wiki/Norm_(mathematics)">`Norm_(mathematics)`</a>                                             |
| Metric                                                               | Metriek                                                                                     | $d(v,w)$                                                             | <a href="https://en.wikipedia.org/wiki/Metric_space">`Metric_space`</a>                                                   |
| (Linear) isometry                                                    | (Lineaire) isometrie                                                                        | $\hat{W}$                                                            | <a href="https://en.wikipedia.org/wiki/Isometry#Linear_isometry">`Isometry#Linear_isometry`</a>                                       |
| Limit of a sequence                                                  | Limiet van een rij                                                                          | $\lim_{n\to +\infty} a_n$                                            | <a href="https://en.wikipedia.org/wiki/Limit_of_a_sequence">`Limit_of_a_sequence`</a>                                            |
| Continuity of a map                                                  | Continuïteit van een afbeelding                                                             |                                                                      | <a href="https://en.wikipedia.org/wiki/Continuous_function#Continuous_functions_between_metric_spaces">`Continuous_function#Continuous_functions_between_metric_spaces`</a> |
| Open set                                                             | Open verzameling                                                                            |                                                                      | <a href="https://en.wikipedia.org/wiki/Open_set">`Open_set`</a>                                                       |
| Closed set                                                           | Gesloten verzameling                                                                        |                                                                      | <a href="https://en.wikipedia.org/wiki/Closed_set">`Closed_set`</a>                                                     |
| Compact set                                                          | Compacte verzameling                                                                        |                                                                      | <a href="https://en.wikipedia.org/wiki/Compact_space">`Compact_space`</a>                                                  |
| Closure                                                              | Sluiting                                                                                    |                                                                      | <a href="https://en.wikipedia.org/wiki/Closure_(topology)">`Closure_(topology)`</a>                                             |
| Dense subset                                                         | Dichte deelverzameling                                                                      |                                                                      | <a href="https://en.wikipedia.org/wiki/Dense_set">`Dense_set`</a>                                                      |
| Separable (metric) space                                             | Separabele (metrische) ruimte                                                               |                                                                      | <a href="https://en.wikipedia.org/wiki/Separable_space">`Separable_space`</a>                                                |
| Cauchy sequence                                                      | Cauchyrij                                                                                   |                                                                      | <a href="https://en.wikipedia.org/wiki/Cauchy_sequence">`Cauchy_sequence`</a>                                                |
| Complete metric space                                                | Complete metrische ruimte                                                                   |                                                                      | <a href="https://en.wikipedia.org/wiki/Complete_metric_space">`Complete_metric_space`</a>                                          |
| Banach space                                                         | Banach ruimnte                                                                              |                                                                      | <a href="https://en.wikipedia.org/wiki/Banach_space">`Banach_space`</a>                                                   |
| Absolute convergence of a series                                     | Absolute convergentie van een reeks                                                         |                                                                      | <a href="https://en.wikipedia.org/wiki/Absolute_convergence">`Absolute_convergence`</a>                                           |
| Bounded                                                              | Begrensd                                                                                    |                                                                      | <a href="https://en.wikipedia.org/wiki/Bounded_operator">`Bounded_operator`</a>                                               |
| Operator norm                                                        | Operatornorm                                                                                |                                                                      | <a href="https://en.wikipedia.org/wiki/Operator_norm">`Operator_norm`</a>                                                  |
| Spectral radius                                                      | Spectraalstraal                                                                             | $\rho_{\hat{A}}$                                                     | <a href="https://en.wikipedia.org/wiki/Spectral_radius">`Spectral_radius`</a>                                                |
| Condition number                                                     | Conditiegetal                                                                               | $\kappa(\hat{A})$                                                    | <a href="https://en.wikipedia.org/wiki/Condition_number">`Condition_number`</a>                                               |
| Bilinear form                                                        | Bilineaire vorm                                                                             | $B(v,w)$                                                             | <a href="https://en.wikipedia.org/wiki/Bilinear_form">`Bilinear_form`</a>                                                  |
| Quadratic form                                                       | Kwadratische vorm                                                                           | $q(v)$                                                               | <a href="https://en.wikipedia.org/wiki/Quadratic_form">`Quadratic_form`</a>                                                 |
| Sesquilinear form                                                    | Sesquilineaire vorm                                                                         | $C(v,w)$                                                             | <a href="https://en.wikipedia.org/wiki/Sesquilinear_form">`Sesquilinear_form`</a>                                              |
| (Positive) definitene, indefinite                                    | (Positief) definiet, indefiniet                                                             |                                                                      | <a href="https://en.wikipedia.org/wiki/Definite_quadratic_form">`Definite_quadratic_form`</a>                                        |
| Congruent matrices                                                   | Congruente matrices                                                                         | $\tilde{\mathsf{A}} = \mathsf{V}^{\mathsf{T}} \mathsf{A} \mathsf{V}$ | <a href="https://en.wikipedia.org/wiki/Matrix_congruence">`Matrix_congruence`</a>                                              |
| Inner product, scalar product                                        | Inwendig product, inproduct, scalair product                                                | $\langle v, w \rangle$                                               | <a href="https://en.wikipedia.org/wiki/Inner_product_space">`Inner_product_space`</a>                                            |
| Standard inner product, dot product                                  | Standaardinproduct                                                                          | $\mathbf{v} \cdot \mathbf{w}$                                        | <a href="https://en.wikipedia.org/wiki/Dot_product">`Dot_product`</a>                                                    |
| Gram matrix                                                          | Grammatrix                                                                                  | $g_{ij}$                                                             | <a href="https://en.wikipedia.org/wiki/Gram_matrix">`Gram_matrix`</a>                                                    |
| Hilbert space                                                        | Hilbertruimte                                                                               |                                                                      | <a href="https://en.wikipedia.org/wiki/Hilbert_space">`Hilbert_space`</a>                                                  |
| Orthogonality                                                        | Orthogonaliteit                                                                             | $v \perp w$                                                          | <a href="https://en.wikipedia.org/wiki/Orthogonality_(mathematics)">`Orthogonality_(mathematics)`</a>                                    |
| Orthonormality                                                       | Orthonormaliteit                                                                            |                                                                      | <a href="https://en.wikipedia.org/wiki/Orthonormality">`Orthonormality`</a>                                                 |
| Orthogonal complement                                                | Orthogonaal complement                                                                      | $S^\perp$                                                            | <a href="https://en.wikipedia.org/wiki/Orthogonal_complement">`Orthogonal_complement`</a>                                          |
| Orthogonal projection                                                | Orthogonale projectie                                                                       |                                                                      |                                                                                              |
| QR decomposition                                                     | QR-decompositie                                                                             | $\mathsf{A} = \mathsf{Q}\mathsf{R}$                                  | <a href="https://en.wikipedia.org/wiki/QR_decomposition">`QR_decomposition`</a>                                               |
| (Hermitian) adjoint                                                  | (Hermitisch) toegevoegde                                                                    | $\hat{A}^\dagger$                                                    | <a href="https://en.wikipedia.org/wiki/Hermitian_adjoint">`Hermitian_adjoint`</a>                                              |
| Self-adjoint operator                                                | Zelftoegevoegde operator                                                                    |                                                                      | <a href="https://en.wikipedia.org/wiki/Self-adjoint_operator">`Self-adjoint_operator`</a>                                          |
| Unitary                                                              | Unitair                                                                                     |                                                                      | <a href="https://en.wikipedia.org/wiki/Unitary_operator">`Unitary_operator`</a>                                               |
| (Linear) least squares method/solution                               | Kleinste-kwadratenmethode/-oplossing                                                        |                                                                      | <a href="https://en.wikipedia.org/wiki/Linear_least_squares">`Linear_least_squares`</a>                                           |
| Discrete Fourier transform                                           | Discrete Fouriertransformatie                                                               |                                                                      | <a href="https://en.wikipedia.org/wiki/Discrete_Fourier_transform">`Discrete_Fourier_transform`</a>                                     |
| Schur decomposition                                                  | Schurdecompositie                                                                           | $\mathsf{A}\mathsf{U} = \mathsf{U} \mathsf{T}$                       | <a href="https://en.wikipedia.org/wiki/Schur_decomposition">`Schur_decomposition`</a>                                            |
| Singular value decomposition                                         | Singulierewaardenontbinding                                                                 | $\mathsf{A} = \mathsf{U} \mathsf{S} \mathsf{V}^{\mathsf{H}}$         | <a href="https://en.wikipedia.org/wiki/Singular_value_decomposition">`Singular_value_decomposition`</a>                                   |
| Sylvester’s law of inertia                                           | Sylvester’s traagheidswet                                                                   |                                                                      | <a href="https://en.wikipedia.org/wiki/Sylvester%27s_law_of_inertia">`Sylvester%27s_law_of_inertia`</a>                                   |
| Low-rank approximation                                               | Lage-rangbenadering (?)                                                                     |                                                                      | <a href="https://en.wikipedia.org/wiki/Low-rank_approximation">`Low-rank_approximation`</a>                                         |
| Trigonometric polynomial                                             | Trigonometrische polynoom                                                                   |                                                                      | <a href="https://en.wikipedia.org/wiki/Trigonometric_polynomial">`Trigonometric_polynomial`</a>                                       |
| Fourier series                                                       | Fourierreeks                                                                                |                                                                      | <a href="https://en.wikipedia.org/wiki/Fourier_series">`Fourier_series`</a>                                                 |
| Point spectrum, continuous spectrum, residual spectrum               | Puntspectrum, continu spectrum, residueel spectrum                                          |                                                                      | <a href="https://en.wikipedia.org/wiki/Decomposition_of_spectrum_(functional_analysis)">`Decomposition_of_spectrum_(functional_analysis)`</a>                |
| Differential equation                                                | Differentiaalvergelijking                                                                   |                                                                      | <a href="https://en.wikipedia.org/wiki/Differential_equation">`Differential_equation`</a>                                          |
| Differential operator                                                | Differentiaaloperator                                                                       | $\hat{L}$                                                            | <a href="https://en.wikipedia.org/wiki/Differential_operator">`Differential_operator`</a>                                          |
| Formal adjoint                                                       | Formeel toegevoegde                                                                         |                                                                      |                                                                                              |
| Boundary conditions                                                  | Randvoorwaarden                                                                             |                                                                      |                                                                                              |
| Boundary value problem                                               | Randvoorwaardeprobleem                                                                      |                                                                      | <a href="https://en.wikipedia.org/wiki/Boundary_value_problem">`Boundary_value_problem`</a>                                         |
| Initial value problem                                                | Beginvoorwaardeprobleem                                                                     |                                                                      | <a href="https://en.wikipedia.org/wiki/Initial_value_problem">`Initial_value_problem`</a>                                          |
| Green’s function                                                     | Greense functie                                                                             |                                                                      | <a href="https://en.wikipedia.org/wiki/Green%27s_function">`Green%27s_function`</a>                                             |
| Sturm-Liouville operator                                             | Sturm-Liouville operator                                                                    |                                                                      | <a href="https://en.wikipedia.org/wiki/Sturm–Liouville_theory">`Sturm–Liouville_theory`</a>                                         |
| (Continuous) Fourier transform                                       | (Continue) Fouriertransformatie                                                             | $\widehat{f}$                                                        | <a href="https://en.wikipedia.org/wiki/Fourier_transform">`Fourier_transform`</a>                                              |
| Convolution                                                          | Convolutie                                                                                  |                                                                      | <a href="https://en.wikipedia.org/wiki/Convolution">`Convolution`</a>                                                    |
| Gaussian / normal distribution                                       | Gaussische / normale verdeling                                                              |                                                                      | <a href="https://en.wikipedia.org/wiki/Normal_distribution">`Normal_distribution`</a>                                            |
| Characteristic function                                              | Karateristieke functie                                                                      |                                                                      | <a href="https://en.wikipedia.org/wiki/Characteristic_function_(probability_theory)">`Characteristic_function_(probability_theory)`</a>                   |
| Distribution                                                         | Distributie                                                                                 |                                                                      | <a href="https://en.wikipedia.org/wiki/Distribution_(mathematics)">`Distribution_(mathematics)`</a>                                     |
| Test function                                                        | Testfunctie                                                                                 |                                                                      |                                                                                              |
| Cauchy principal value                                               | Cauchy hoofdwaarde                                                                          |                                                                      | <a href="https://en.wikipedia.org/wiki/Cauchy_principal_value">`Cauchy_principal_value`</a>                                         |
| Discrete-time Fourier transform                                      | Discrete-tijd-Fouriertransformatie                                                          |                                                                      | <a href="https://en.wikipedia.org/wiki/Discrete-time_Fourier_transform">`Discrete-time_Fourier_transform`</a>                                |
| Sampling                                                             | Bemonsteren                                                                                 |                                                                      | <a href="https://en.wikipedia.org/wiki/Sampling_(signal_processing)">`Sampling_(signal_processing)`</a>                                   |
| Nyquist rate / frequency                                             | Nyquist frequentie                                                                          |                                                                      | <a href="https://en.wikipedia.org/wiki/Nyquist_rate">`Nyquist_rate`</a>                                                   |
| Aliasing                                                             | Aliasing, vouwvervorming                                                                    |                                                                      | <a href="https://en.wikipedia.org/wiki/Aliasing">`Aliasing`</a>                                                       |
