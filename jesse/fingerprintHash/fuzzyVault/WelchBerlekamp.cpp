
#include "WelchBerlekamp.h"

//static bool LinearSolve (vec_ZZ_p & solution, mat_ZZ_p & A);


/*
 * WelchBerlekamp.cpp
 * Authors: Soren Johnson and Leonid Reyzin (reyzin@cs.bu.edu)
 *
 * This code and explanatory notes are hosted at
 * http://www.cs.bu.edu/~reyzin/code/fuzzy.html
 *
 * Uses Victor Shoup's NTL, see http://www.shoup.net
 *
 * This function decodes [n,n-2t,2t+1] Reed-Solomon codes over GF(p) using
 * the Welch-Berlekamp algorithm, as described in Section 3 of Notes for
 * Lecture 10 of Madhu Sudan's 2001 course "Algorithmic Introduction to
 * Coding Theory" http://people.csail.mit.edu/madhu/FT01/ (that description
 * is based, in turn, on a description in "Highly Resilient Correctors for
 * Multivariate Polynomials," by Gemmel and Sudan, Information Processing
 * Letters, 43(4):169-174, 1992).
 *
 * Namely, given vectors X and Y of length n each, it finds a polynomial f
 * over GF(p) of degree at most n-2t-1 such that f(X[i])=Y[i] for at least
 * n-t values of i.  Such a polynomial f, if it exists, is unique by the
 * funamental theorem of algebra (if there were two such polynomials f_1,
 * f_2, then f_1-f_2 would have at least n-2t roots, but degree at most
 * n-2t-1).
 *
 * Assumes (does not check!) that no value appears in X more than once,
 * that X and Y are of the same length, and that ZZ_p::init has been run
 * and that ZZ_p::modulus() is a prime.
 *
 * The polynomial f is placed into the "result" argument. Returns true if
 * successful, false if such a polynomial does not exist.
 *
 * This algorithm takes Theta(n^3) operations over GF(p).  It is not the
 * best known algorithm for this problem; faster algorithms are referenced
 * in the lecture notes mentioned above.
 *
 * NOTE for those using this code as part of Improved Juels-Sudan Secure
 * Sketch: the meaning of n and t in this function is different from the
 * meaning of n and t in ijs.cpp and ijsio.cpp, in order to be consistent
 * with standard coding notation.
 */
bool WelchBerlekamp (ZZ_pX & result, const vec_ZZ_p & X, const vec_ZZ_p & Y,  int t)
{
  int n = X.length();

  if (t<0) 
    return false;
  if (t==0) {                   // special case t=0: then no errors are allowed,
    interpolate(result, X, Y);  // and result is just an interpolation of X, Y.
    return true;                // We special-case it here to avoid
  }                             // creating a special case later in the code
  if (n<2*t) // degree n-2t-1<-1 doesn't make sense if n<2t.
    return false;
            
  // Note that the case of n=2t is possible: the result should be the 0
  // polynomial if at least half of the Y values are 0, and false otherwise
  // The decoding algorithm handles that correctly below

  int i, j;
  mat_ZZ_p A;  // the matrix to hold the linear system we'll be solving
  vec_ZZ_p NE; // coefficients of N and E as one vector

  A.SetDims(n, n + 1);

  /* We know that N has degree at most n-t-1,
     E has degree exactly t and is monic,
     and N(X[i])=Y[i]E(X[i]) for all i.  Thus, we need to solve
     the following matrix equation to find the coefficients of N and E:

    [1 x0 x0^2 ... x0^(n-t-1) -y0 -y0x0 -y0x0^2 ... -y0x0^(t-1)] [N0      ]   [y0x0^t]
    [1 x1 x1^2 ... x1^(n-t-1) -y1 -y1x1 -y1x1^2 ... -y1x1^(t-1)] [N1      ]   [y1x1^t]
    ...                                                          ...          ...
    ...                                                          [N(n-t-1)]   ...
    ...                                                          [E0      ] = ...
    ...                                                          [E1      ]   ...
    ...                                                          ...          ...
    [1 xm xm^2 ... xm^(n-t-1) -ym -ymxm -ymxm^2 ... -ymxm^(t-1)] [E(t-1)  ]   [ymxm^t]

    where m = n - 1 (it just got too messy to write (n-1) at the bottom).

  */


  int numNCoeffs=n-t; 
  // It is important that numNCoeffs<n because t<=0 was already dealt with above

  // populate the augmented matrix that includes the right-hand side as the last column
  for (i = 0; i < n; i++) {
    A[i][0] = 1;
    for (j = 1; j < numNCoeffs; j++) // set A[i][j] to have x_i^j for 0<=j<numNCoeffs
      mul(A[i][j], A[i][j - 1], X[i]);
    A[i][j]=-Y[i]; // this is allowed because j=numNCoeffs<n
    j++;
    for (; j<n; j++) // set A[i][j] to have -y_i x_i^{j-numNCoeffs} for numNCoeffs < j < n
      mul(A[i][j], A[i][j-1], X[i]);
    mul(A[i][n], -A[i][n - 1], X[i]); // set A[i][n] to have y_i x_i^t
  }

  // Solve the linear system -- this can fail if it has no solutions
  if (!LinearSolve(NE, A))
    return false;

  ZZ_pX N; 
  ZZ_pX E; // the error locating polynomial, degree t
  N.SetMaxLength(n-t);
  E.SetMaxLength(t+1);
  // set the coefficients of N and E that we solved for
  for (i = 0; i < numNCoeffs; i++)
    SetCoeff(N, i, NE[i]);
  for (i = 0; i < t; i++)
    SetCoeff(E, i, NE[numNCoeffs+i]);
  // set the leading coefficient of E to 1
  SetCoeff(E, t);

  ZZ_pX remainder;
  DivRem(result, remainder, N, E);
  // if remainder is 0, result is the polynomial we want.  Proof: N=pE
  // (where p is the result); since N(x_i)= y_i E(x_i) for n different x_i,
  // and E can have at most t roots, we have p(x_i) E(x_i) = y_i E(x_i) for
  // at least n-t nonzero values of E(x_i), and hence p(x_i) = y_i for at
  // least n-t values

  // If the remainder is nonzero, then the polynomial we want does not
  // exist, by the proof of Welch-Berlekamp (simply because any solution
  // (N, E) must have the ratio p if p with the desired properties exists)
  if (IsZero(remainder))
    return true;
  else 
    return false;
}


/*
 * Finds a solution to a system of linear equations represtented by an
 * n-by-n+1 matrix A: namely, denoting by B the left n-by-n submatrix of A
 * and by C the last column of A, finds a column vector x such that Bx=C.
 * If more than one solution exists, chooses one arbitrarily by setting some
 * values to 0.  If no solutions exists, returns false.  Otherwise, places
 * a solution into the first agrument and returns true.  Assumes (does not
 * check!) that A is an n by n+1 matrix and that ZZ_p::init has been run
 * and that ZZ_p::modulus() is a prime.
 *
 * NOTE that the matrix A changes (gets converted to row echelon form).
 */
static
bool LinearSolve (vec_ZZ_p & solution, mat_ZZ_p & A) {
  int n = A.NumRows();
  solution.SetLength(n);

  ZZ_p t;


  int firstDeterminedValue = n; 
  // we will be determining values of the solution
  // from n-1 down to 0.  At any given time,
  // values from firstDeterminedValue to n-1 have been
  // found. Initializing to n means
  // no values have been found yet.
  // To put it another way, the variabe firstDeterminedValue
  // stores the position of first nonzero entry in the row just examined
  // (except at initialization)

  int rank=gauss (A); // can start at rank-1, because below that are all zeroes
  for (int row = rank-1; row >=0; row--) {
    // remove all the known variables from the equation
    t = A[row][n];
    int col;
    for (int col = n-1; col>=firstDeterminedValue; col--)
      t-=(A[row][col]*solution[col]);
    
    // now we need to find the first nonzero coefficient in this row
    // if it exists before firstDeterminedValue
    // because the matrix is in row echelon form, the first nonzero
    // coefficient cannot be before the diagonal
    for (col=row; col<firstDeterminedValue; col++)
      if (!IsZero(A[row][col]))
	break;

    if (col<firstDeterminedValue) { // this means we found a nonzero coefficient
    // we can determine the variables in position from col to firstDeterminedValue 
      // if this loop executes even once, then the system is undertermined
      // we arbitrarily set the undetermined variables to 0, because it make math easier
      for (int j = col+1; j<firstDeterminedValue; j++)
	clear(solution[j]);
      // Now determine the variable at the nonzero coefficient
      div(solution[col], t, A[row][col]);
      firstDeterminedValue = col;
    }

    else // this means there are no nonzero coefficients before firstDeterminedValue.
      // Because we skip all the zero rows at the bottom, the matrix is in
      // row echelon form, and firstDeterminedValue is equal to the
      // position of first nonzero entry in row+1 (unless it is equal to n),
      // this means we are at a row with all zeroes except in column n
      // The system has no solution.
      return false;
    
  }

  // set the remaining undetermined values, if any, to 0
  for (int col=0; col<firstDeterminedValue; col++) 
    clear (solution[col]);

  return true;
}

