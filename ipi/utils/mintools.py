"""Contains simple helper algorithms for minimization.

Copyright (C) 2013, Joshua More and Michele Ceriotti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http.//www.gnu.org/licenses/>.

Algorithms implemented by Michele Ceriotti and Benjamin Helfrecht, 2015

Functions:
        bracket: Determines the 3 points that bracket the function minimum
        min_brent:  Does one-D minimization (line search) based on bisection
            method with derivatives. Uses 'bracket' function.
        min_approx: Does approximate n-D minimization (line search) based
            on sufficient function decrease in the search direction
        min_trm: Does approximate n-D minimization inside a trust-region

        BFGS: Constructs an approximate inverse Hessian to determine
            new search directions. Minimizes the function using 'min_approx' function.
        BFGS-TRM: Constructs an approximate inverse Hessian to determine
            new search directions. Minimizes the function using 'min_trm' function.
        L-BFGS: Uses the limited memory BFGS algorithm (L-BFGS) to
            compute new search directions. Minimizes using 'min_approx'
        L-BFGS_nls: L-BFGS algorithm without line search
            *** This function is less stable than L-BFGS and not any more efficient ***
        bracket_neb: Modified 'bracket' routine to make
            compatible with functions with unknown gradient
        min_brent_neb: Modified 'min_brent' routine to make
            compatible with functions with unknown gradient

        bracket, bracket_neb, min_brent, min_brent_neb,and BFGS subroutines adapted from:
            Press, W. H., Teukolsky, S. A., Vetterling, W. T., and Flannery, B. P. (1992).
            Numerical Recipes in C: The Art of Scientific Computing.
            Cambridge: Cambridge University Press

        LBFGS subroutine adapted from:
            Nocedal, J. (1980). Updating Quasi-Newton Matrices with
            Limited Storage. Mathematics of Computation, 35, 773-782.
            DOI: http://dx.doi.org/10.1090/S0025-5718-1980-0572855-7
        powell: powell formula to update the hessian
             (R. Fletcher. Practical Methods of Optimization. 2nd ed.(1987)
        nichols: nichols algorithm for optimization (minimum or transition state)
        Simons, J. and Nichols, J. (1990), Int. J. Quantum Chem., 38: 263-276.
"""

# TODO: CLEAN UP BFGS, L-BFGS, L-BFGS_nls TO NOT EXIT WITHIN MINTOOLS.PY
#       BUT USE UNIVERSAL SOFTEXIT

__all__ = ["min_brent"]

import numpy as np
import math

# import time
from ipi.utils.softexit import softexit
from ipi.utils.messages import verbosity, info, warning


# Bracketing function
def bracket(fdf, fdf0=None, x0=0.0, init_step=1.0e-3):
    """Given an initial point, determines the initial bracket for the minimum
    Arguments:
           fdf: function to minimize, derivative of function to minimize
           x0: initial point
           fdf0: value of function and its derivative at x0
           init_step: intial step size
    """

    # Constants
    gold = 1.618034  # Golden ratio
    glimit = 100.0  # Limit for magnification of parabolic fit step
    tiny = 1.0e-20  # Prevent division by zero

    # x0: initial point of evaluation (e.g. initial atomic position)
    # ax, bx, cx: bracketing points with ax < bx < cx
    # bracketing finished if an ax, bx, cx with ax < bx < cx is found
    # fa, fb, fc: value of function at ax, bx, cx

    if fdf0 is None:
        fdf0 = fdf(x0)
    ax = x0
    fa, dfa = fdf0
    bx = x0 + init_step
    fb, dfb = fdf(bx)
    info(" @BRACKET: Started bracketing", verbosity.debug)
    info(" @BRACKET: Evaluated first step", verbosity.debug)

    # Switch direction to move downhill, if necessary, and rearrange
    if fb > fa:
        tmp = ax
        ax = bx
        bx = tmp
        tmp = fb
        fb = fa
        fa = tmp
        tmp = dfb
        dfb = dfa
        dfa = tmp

    # Initial guess for third bracketing point
    cx = bx + gold * (bx - ax)
    fc, dfc = fdf(cx)
    info(
        " @BRACKET: Evaluated initial bracket: (%f:%f, %f:%f, %f:%f)"
        % (ax, fa, bx, fb, cx, fc),
        verbosity.debug,
    )

    # Loop until acceptable bracketing condition is achieved
    # u is a point between two of the bracketing points,
    # Use parabolic extrapolation to find u. "tiny" prevents possible division by zero
    while fb > fc:
        r = (bx - ax) * (fb - fc)
        q = (bx - cx) * (fb - fa)
        u = bx - ((bx - cx) * q - (bx - ax) * r) / (
            2.0 * math.copysign(max(abs(q - r), tiny), (q - r))
        )  # Point from parabolic fit
        ulim = bx + glimit * (
            cx - bx
        )  # Limit for parabolic fit point; *Can test various possibilities*

        # Find minimums between b and c or a and u
        # If parabolic fit unsuccessful, use default step magnification
        # Otherwise:
        # - Parabolic fit between c and its allowed limit
        # - Limit u to maximum allowed value
        # - Use default magnification
        if ((bx - u) * (u - cx)) > 0.0:
            fu, dfu = fdf(u)
            info(" @BRACKET: Evaluated new bracket point", verbosity.debug)
            # Found minimum between b and c?
            # -b------u-----c shift:
            # -a------b-----c
            if fu < fc:
                ax = bx
                bx = u
                fa = fb
                fb = fu
                dfa = dfb
                dfb = dfu
                info(
                    " @BRACKET: Bracketing completed: (%f:%f, %f:%f, %f:%f)"
                    % (ax, fa, bx, fb, cx, fc),
                    verbosity.debug,
                )
                return (ax, bx, cx, fb, dfb)
                # minimum between a and u?
                # -a-----b-----u-----c shift:
                # -a-----b-----c
            elif fu > fb:
                cx = u
                fc = fu
                dfc = dfu
                info(" @BRACKET: Bracketing completed", verbosity.debug)
                return (ax, bx, cx, fb, dfb)
                # parabolic extrapolation was not successful.
                # Use golden value (initial guess, default magnification).
            u = cx + gold * (cx - bx)
            fu, dfu = fdf(u)
            info(" @BRACKET: Evaluated new bracket point", verbosity.debug)
        # c < u (parabolic fit) < ulim?
        elif ((cx - u) * (u - ulim)) > 0.0:
            fu, dfu = fdf(u)
            info(" @BRACKET: Evaluated new bracket point", verbosity.debug)
            # minimum between c and u+gold(u-cx)?
            # -c----u----u+gold(u-cx) shift:
            # -b----c----u
            if fu < fc:
                bx = cx
                cx = u
                u = cx + gold * (cx - bx)
                fb = fc
                fc = fu
                dfb = dfc
                dfc = dfu
                fu, dfu = fdf(u)
                info(" @BRACKET: Evaluated new bracket point", verbosity.debug)
        # u >= ulim? limit u.
        elif ((u - ulim) * (ulim - cx)) >= 0.0:
            u = ulim
            fu, dfu = fdf(u)
            info(" @BRACKET: Evaluated new bracket point", verbosity.debug)
        # reject parabolic u, use golden value (default magnification)
        else:
            u = cx + gold * (cx - bx)
            fu, dfu = fdf(u)
            info(" @BRACKET: Evaluated new bracket point", verbosity.debug)

        # Shift points, so that points to continue are ax, bx, cx
        ax = bx
        bx = cx
        cx = u
        fa = fb
        fb = fc
        fc = fu
        dfa = dfb
        dfb = dfc
        dfc = dfu

    info(
        " @BRACKET: Bracketing completed: (%f:%f, %f:%f, %f:%f)"
        % (ax, fa, bx, fb, cx, fc),
        verbosity.debug,
    )
    return (ax, bx, cx, fb, dfb)


# One dimensional minimization function using function derivatives
# and Brent's method


def min_brent(fdf, fdf0, x0, tol, itmax, init_step):
    """Given a maximum number of iterations and a convergence tolerance,
    minimizes the specified function
    Arguments:
           x0: initial x-value
           fdf: function to minimize
           fdf0: initial function value
           tol: convergence tolerance
           itmax: maximum allowed iterations
           init_step: initial step size
    """

    # Initializations and constants
    # gold = 0.3819660  # Golden ratio
    # Safeguard against trying to find fractional precision for min that is exactly zero
    zeps = 1.0e-10
    e = 0.0  # Size of step before last

    # Call initial bracketing routine
    (ax, bx, cx, fb, dfb) = bracket(fdf, fdf0, x0, init_step)

    # Set bracket points
    if ax < cx:
        a = ax
    else:
        a = cx
    if ax > cx:
        b = ax
    else:
        b = cx

    # Initial points to evaluate
    # f* is evaluation of function
    # df* is the evaluation of the derivative
    # x = point with least function value so far
    # w = point with 2nd least function value
    # v = previous value of w
    # u = point at which function was evaluated most recently
    x = w = v = bx
    fw = fv = fx = fb  # Function
    dfw = dfv = dfx = dfb  # Function derivative

    # Main loop
    j = 1
    info(" @MINIMIZE: Started 1D minimization", verbosity.debug)
    while j <= itmax:
        # Determine tolerance
        xm = 0.5 * (a + b)
        tol1 = tol * abs(x) + zeps
        tol2 = 2.0 * tol1

        # Test for satisfactory completion: |b-a|<=tol*abs(x)
        if abs(x - xm) <= (tol2 - 0.5 * (b - a)):
            info(" @MINIMIZE: Finished minimization, energy = %f" % fx, verbosity.debug)
            # return (x, fx, dfx)
            fx, dfx = fdf(x)  # Evaluate again to update lm.dforces object
            return

        # Initialize d values (used to determine step size) to outside of bracket
        d = 0.0
        if abs(e) > tol1:
            d1 = 2.0 * (b - a)
            d2 = d1

            # Secant method with both d points (find zero point of derivative)
            if dfw != dfx:
                d1 = (w - x) * dfx / (dfx - dfw)
            if dfv != dfx:
                d2 = (v - x) * dfx / (dfx - dfv)

            # Choose estimate based on derivative at x and move distance on step
            # before last
            # estimates should be within the bracket
            u1 = x + d1
            u2 = x + d2
            ok1 = ((a - u1) * (u1 - b) > 0.0) and (dfx * d1 <= 0.0)
            ok2 = ((a - u2) * (u2 - b) > 0.0) and (dfx * d2 <= 0.0)
            olde = e
            e = d

            # Take an acceptable d
            # (x+d1 and x+d2 within bracket and d1,d2 have different sign from dfx);
            # if both are acceptable, choose smallest
            if ok1 or ok2:
                if ok1 and ok2:
                    if abs(d1) < abs(d2):
                        d = d1
                    else:
                        d = d2
                elif ok1:
                    d = d1
                else:
                    d = d2
                # movement must be less than half of the movement of the step before last
                # (better not to punish algorithm for one bad step)
                # (e: last step, olde: step before last)
                if abs(d) <= abs(0.5 * olde):
                    u = x + d
                    if ((u - a) < tol2) or ((b - u) < tol2):
                        d = math.copysign(tol1, (xm - x))
                    else:
                        # new d with d = 0.5*e
                        # which segment is decided by sign of derivative
                        if dfx >= 0.0:
                            e = a - x
                        else:
                            e = b - x
                        d = 0.5 * e
            # conditions for d1,d2 not fulfilled, new d = 0.5*e
            else:
                if dfx >= 0.0:
                    e = a - x
                else:
                    e = b - x
                d = 0.5 * e
        else:
            if dfx >= 0.0:
                e = a - x
            else:
                e = b - x
            d = 0.5 * e
        # one function evaluation per iteration
        if abs(d) >= tol1:
            u = x + d
            fu, dfu = fdf(u)
        else:
            u = x + math.copysign(tol1, d)
            fu, dfu = fdf(u)

            # If minimum step in downhill direction goes uphill, minimum has been found
            if fu > fx:
                info(
                    " @MINIMIZE: Finished minimization, energy = %f" % fx,
                    verbosity.debug,
                )
                # return (x, fx,dfx)
                fx, dfx = fdf(x)  # Evaluate again to update lm.dforces object
                return
        # order for next step: a < u (later x) < b
        if fu <= fx:
            if u >= x:
                a = x
            else:
                b = x
            # shift for next step
            v = w
            fv = fw
            dfv = dfw
            w = x
            fw = fx
            dfw = dfx
            x = u
            fx = fu
            dfx = dfu
        else:
            if u < x:
                a = u
            else:
                b = u
            if (fu <= fw) or (w == x):
                v = w
                fv = fw
                dfv = dfw
                w = u
                fw = fu
                dfw = dfu
            elif (fu < fv) or (v == x) or (v == w):
                v = u
                fv = fu
                dfv = dfu
        j += 1

    # Exit if maximum number of iterations exceeded
    info(
        " @MINIMIZE: Error -- maximum iterations for minimization (%d) exceeded, \
        exiting minimization"
        % itmax,
        verbosity.low,
    )
    info(" @MINIMIZE: Finished minimization, energy = %f" % fx, verbosity.debug)
    # return (x, fx,dfx)
    return


# Approximate line search


def min_approx(fdf, x0, fdf0, d0, big_step, tol, itmax):
    """Given an n-dimensional function and its gradient, and an
    initial point and a direction, finds a new point where the function
    is thought to be 'sufficiently' minimized, i.e. carries out an
    approximate minimization.
        Arguments:
            fdf: function and its gradient
            fdf0: initial function and gradient value
            d0: n-dimensional initial direction
            x0: n-dimensional initial point
            big_step: maximum step size
            tol: tolerance for exiting line search
            itmax: maximum number of iterations for the line search
    """

    # Initializations and constants
    info(" @MINIMIZE: Started approx. line search", verbosity.debug)
    n = len(x0.flatten())
    if fdf0 is None:
        fdf0 = fdf(x0)
    f0, df0 = fdf0
    if d0 is None:
        d0 = -df0 / np.linalg.norm(df0.flatten())
    x = np.zeros(n)
    alf = 1.0e-4

    # Step size
    stepsum = np.linalg.norm(d0.flatten())

    # Scale if attempted step is too large
    if stepsum > big_step:
        info(" @MINIMIZE: Scaled step size for line search", verbosity.debug)
        d0 *= big_step / stepsum

    slope = np.dot(df0.flatten(), d0.flatten())

    if slope >= 0.0:
        info(" @MINIMIZE: Warning -- gradient is >= 0 (%f)" % slope, verbosity.low)

    test = np.amax(
        np.divide(
            np.absolute(d0.flatten()), np.maximum(np.absolute(x0.flatten()), np.ones(n))
        )
    )

    # Setup to try Newton step first
    alamin = tol / test
    alam = 1.0

    # Minimization Loop
    i = 1
    while i < itmax:
        x = np.add(x0, (alam * d0))
        fx, dfx = fdf(x)
        info(" @MINIMIZE: Calculated energy", verbosity.debug)

        # Check for convergence on change in x
        if alam < alamin:
            x = x0
            info(
                " @MINIMIZE: Convergence in position, exited line search",
                verbosity.debug,
            )
            return (x, fx, dfx)

        # Sufficient function decrease
        elif fx <= (f0 + alf * alam * slope):
            info(
                " @MINIMIZE: Sufficient function decrease, exited line search",
                verbosity.debug,
            )
            return (x, fx, dfx)

        # No convergence; backtrack
        else:
            info(
                " @MINIMIZE: No convergence on step; backtrack to find point",
                verbosity.debug,
            )

            # First backtrack
            if alam == 1.0:
                f2 = None
                alam2 = None
                tmplam = -slope / (2.0 * (fx - f0 - slope))

            # Subsequent backtracks
            # coefficient should lie between:
            # 0.1*alam and 0.5*alam (= 0.1*lambda_1 and 0.5*lambda_1),
            # otherwise step lengths are too small
            else:
                rhs1 = fx - f0 - alam * slope
                rhs2 = f2 - f0 - alam2 * slope
                a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2)
                b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (
                    alam - alam2
                )
                if a == 0.0:
                    tmplam = -slope / (2.0 * b)

                else:
                    disc = b * b - 3.0 * a * slope
                    if disc < 0.0:
                        tmplam = 0.5 * alam

                    elif b <= 0.0:
                        tmplam = (-b + np.sqrt(disc)) / (3.0 * a)

                    else:
                        tmplam = -slope / (b + np.sqrt(disc))

                    # Coefficient less than 0.5 * lambda_1
                    if tmplam > (0.5 * alam):
                        tmplam = 0.5 * alam

        alam2 = alam
        f2 = fx

        # Coefficient greater than 0.1 * lambda_1
        alam = max(tmplam, 0.1 * alam)

        i += 1

    info(
        " @MINIMIZE: Error - maximum iterations for line search (%d) exceeded, \
        exiting search"
        % itmax,
        verbosity.low,
    )
    info(" @MINIMIZE: Finished minimization, energy = %f" % fx, verbosity.debug)
    return (x, fx, dfx)


# BFGS algorithm with approximate line search


def BFGS(x0, d0, fdf, fdf0, invhessian, big_step, tol, itmax):
    """BFGS minimization. Uses approximate line minimizations.
    Does one step.
        Arguments:
            x0: initial point
            d0: initial direction for line minimization
            fdf: function and gradient (mapper)
            fdf0: initial function and gradient value
            big_step: limit on step length
            tol: convergence tolerance
            itmax: maximum number of allowed iterations
    """

    info(" @MINIMIZE: Started BFGS", verbosity.debug)
    zeps = 1.0e-13
    u0, g0 = fdf0

    # Maximum step size
    n = len(x0.flatten())
    linesum = np.dot(x0.flatten(), x0.flatten())
    big_step = big_step * max(np.sqrt(linesum), n)

    # Perform approximate line minimization in direction d0
    x, u, g = min_approx(fdf, x0, fdf0, d0, big_step, tol, itmax)
    d_x = np.subtract(x, x0)

    # Update invhessian.
    d_g = np.subtract(g, g0)
    hdg = np.dot(invhessian, d_g.flatten())

    fac = np.vdot(d_g, d_x)
    fae = np.dot(d_g.flatten(), hdg)
    sumdg = np.vdot(d_g, d_g)
    sumxi = np.vdot(d_x, d_x)

    # Skip update if not 'fac' sufficiently positive
    if fac > np.sqrt(zeps * sumdg * sumxi):
        fac = 1.0 / fac
        fad = 1.0 / fae

        # Compute BFGS term
        dg = np.subtract((fac * d_x).flatten(), fad * hdg)
        invhessian += (
            np.outer(d_x, d_x) * fac - np.outer(hdg, hdg) * fad + np.outer(dg, dg) * fae
        )
        info(" @MINIMIZE: Updated invhessian", verbosity.debug)
    else:
        info(
            " @MINIMIZE: Skipped invhessian update; direction x gradient insufficient",
            verbosity.debug,
        )

    # Update direction
    d = np.dot(invhessian, -g.flatten())
    d0[:] = d.reshape(d_x.shape)
    info(" @MINIMIZE: Updated search direction", verbosity.debug)


# BFGS algorithm trust radius method
def BFGSTRM(x0, u0, f0, h0, tr, mapper, big_step):
    """Input: x0 = previous accepted positions
          u0 = previous accepted energy
          f0 = previous accepted forces
          h0 = previous accepted hessian
          tr = trust radius
      mapper = function to evaluate energy and forces
    big_step = limit on step length"""

    # Make one movement, evaluate if it has to be accepted or not.
    # If accepted, update tr and Hessian.
    # If not only update the tr and restart the loop
    accept = False
    while not accept:
        # Find new movement direction candidate
        d_x = min_trm(f0, h0, tr)
        # Make movement  and get new energy (u)  and forces(f) using mapper
        x = x0 + d_x
        u, g = mapper(x)
        f = -g

        # Compute energy gain

        true_gain = u - u0
        expected_gain = -np.dot(f0.flatten(), d_x.flatten())
        expected_gain += 0.5 * np.dot(
            d_x.reshape((1, d_x.size)), np.dot(h0, d_x.reshape((d_x.size, 1)))
        )
        harmonic_gain = -0.5 * np.dot(d_x.flatten(), (f0 + f).flatten())

        # Compute quality:
        d_x_norm = np.linalg.norm(d_x)

        if d_x_norm > 0.05:
            quality = true_gain / expected_gain
        else:
            quality = harmonic_gain / expected_gain
        accept = quality > 0.1

        # Update TrustRadius (tr)
        if quality < 0.25:
            tr[0] = 0.5 * d_x_norm
        elif quality > 0.75 and d_x_norm > 0.9 * tr:
            tr[0] = 2.0 * tr
            if tr > big_step:
                tr[0] = big_step

    # After accept, update the Hessian
    d_f = np.subtract(f, f0)
    TRM_UPDATE(d_x.flatten(), d_f.flatten(), h0)


# TRM functions (TRM_UPDATE and min_trm)


def TRM_UPDATE(dx, df, h):
    """Input: DX = X -X_old
           DF = F -F_old
           DG = -DF
           H  = hessian
    Task: updated hessian"""

    dx = dx[:, np.newaxis]  # dimension nx1
    dx_t = dx.T  # dimension 1xn
    dg = -df[:, np.newaxis]
    dg_t = dg.T

    # Bakken and Helgaker, JCP, 117,9160. Eq 44
    h1 = np.dot(dg, dg_t)
    h1 = h1 / (np.dot(dg_t, dx))
    h2a = np.dot(h, dx)
    h2b = np.dot(dx_t, h)
    h2 = np.dot(h2a, h2b)
    h2 = h2 / np.dot(dx_t, h2a)

    h += h1 - h2


def min_trm(f, h, tr):
    """Return the minimum of
    E(dx) = -(F * dx + 0.5 * ( dx * H * dx ),
    whithin dx**2 <tr

    IN    f  = forces        (n,)
          h  = hessian       (nxn)
          tr = trust-radius

    OUT   DX = displacement in cartesian basis

    INTERNAL
             ndim = dimension
             d    = hessian eigenvalues
             w    = hessian eigenvector (in columns)
             g    = gradient in cartesian basis
             gE   = gradient in eigenvector basis
             DX   = displacement in cartesian basis
             DXE  = displacement in eigenvector basis
    """

    # Resize
    ndim = f.size
    shape = f.shape
    f = f.reshape((1, ndim))

    # Diagonalize
    d, w = np.linalg.eigh(h)
    d = d[:, np.newaxis]  # dimension nx1

    gEt = np.dot(f, w)  # Change of basis  ##
    gE = gEt.T  # dimension nx1

    # Count negative, zero and positive eigenvalues
    neg = (d < -1e-7).sum()
    zero = (d < 1e-7).sum() - neg
    # pos = d.size - neg - zero

    # Pull out zero-mode gE
    if zero > 0:
        gE[neg : neg + zero] = np.zeros((zero, 1))

    # Real work starts here
    DXE = np.zeros((ndim, 1))

    for i in range(0, ndim):
        if np.absolute(d[i]) > 1e-5:
            DXE[i] = gE[i] / d[i]

    min_d = np.amin(d)

    # Check if h is positive definite and use trivial result if within trust radius
    if min_d > 0.0:
        if neg != 0:
            print("problem in 'find'!!!")
        if np.linalg.norm(DXE) < tr:
            DX = np.dot(w, DXE)
            DX = DX.reshape(shape)
            return DX

    # If we don't have luck, let's start with the iteration
    lamb_min = max(0.0, -min_d)
    lamb_max = 1e30
    lamb = min(lamb_min + 0.5, 0.5 * (lamb_min + lamb_max))

    for i in range(0, 100):
        DXE = gE / (d + lamb)
        y = np.sum(DXE**2) - tr**2
        dy = -2.0 * np.sum((DXE**2) / (d + lamb))

        if np.absolute(y / dy) < 1e-5 or np.absolute(y) < 1e-13:
            break

        if y < 0.0:
            lamb_max = min(lamb, lamb_max)
        else:
            lamb_min = max(lamb, lamb_min)

        if dy > 0.0 or lamb_min > lamb_max:
            print("Problem in find. II")

        lamb = lamb - y / dy
        if lamb <= lamb_min or lamb >= lamb_max:
            lamb = 0.5 * (lamb_min + lamb_max)
    #  print 'iter',i,lamb, lamb_max,lamb_min,y,dy

    DX = np.dot(w, DXE)
    DX = DX.reshape(shape)
    return DX


# L-BFGS algorithm with approximate line search


def L_BFGS(x0, d0, fdf, qlist, glist, fdf0, big_step, tol, itmax, m, scale, k):
    """L-BFGS minimization. Uses approximate line minimizations.
    Does one step.
        Arguments:
            fdf = function and gradient
            fdf0 = initial function and gradient value
            d0 = initial direction for line minimization
            x0 = initial point
            qlist = list of previous positions used for reduced Hessian^-1 construction
            glist = list of previous gradients used for reduced Hessian^-1 construction
            m = number of corrections to store and use
            k = iteration (MD step) number
            big_step = limit on step length
            tol = convergence tolerance
            itmax = maximum number of allowed iterations
    """

    zeps = 1.0e-10
    n = len(x0.flatten())
    alpha = np.zeros(m)
    beta = np.zeros(m)
    rho = np.zeros(m)

    u0, g0 = fdf0

    # Maximum step size
    linesum = np.dot(x0.flatten(), x0.flatten())
    big_step = big_step * max(np.sqrt(linesum), n)

    # MC try to resolve the stuck BFGS bug
    if np.dot(g0.flatten(), d0.flatten()) > 0.0:
        # reset search direction if we are moving uphill!
        info(" @MINIMIZE: moving uphill, resetting search direction! ", verbosity.debug)
        d0 = g0 / np.sqrt(np.dot(g0.flatten(), g0.flatten()))

    print("@ GEOP step ", big_step)
    # Perform approximate line minimization in direction d0
    x, u, g = min_approx(fdf, x0, fdf0, d0, big_step, tol, itmax)

    # Compute difference of positions (gradients)
    # Build list of previous 'd_positions (d_gradients)'

    d_x = np.subtract(x, x0)
    if k < m:
        qlist[k] = d_x.flatten()
    else:
        qlist_aux = np.roll(qlist, -1, axis=0)
        qlist[:] = qlist_aux
        qlist[m - 1] = d_x.flatten()

    d_g = np.subtract(g, g0)
    if k < m:
        glist[k] = d_g.flatten()
    else:
        glist_aux = np.roll(glist, -1, axis=0)
        glist[:] = glist_aux
        glist[m - 1] = d_g.flatten()

    # Update direction.
    # 1_Determine bounds for L-BFGS 'two loop recursion'
    if k < (m - 1):
        bound1 = k
        bound2 = k + 1
    else:
        bound1 = m - 1
        bound2 = m

    # 2
    q = g.flatten()

    # 3_Loops
    fac = np.dot(d_g.flatten(), d_x.flatten())
    sumdg = np.dot(d_g.flatten(), d_g.flatten())
    sumdx = np.dot(d_x.flatten(), d_x.flatten())

    # Skip update if not 'fac' sufficiently positive
    if fac > np.sqrt(zeps * sumdg * sumdx):
        # Begin two loop recursion:
        # First loop
        for j in range(bound1, -1, -1):
            rho[j] = 1.0 / np.dot(glist[j], qlist[j])
            alpha[j] = rho[j] * np.dot(qlist[j], q)
            q = q - alpha[j] * glist[j]

        info(" @MINIMIZE: First L-BFGS loop recursion completed", verbosity.debug)

        if scale == 0:
            hk = 1.0
        elif scale == 1:
            hk = np.dot(glist[0], qlist[0]) / np.dot(glist[0], glist[0])
        elif scale == 2:
            hk = np.dot(glist[bound1], qlist[bound1]) / np.dot(
                glist[bound1], glist[bound1]
            )

        d = hk * q

        # Second loop
        for j in range(0, bound2, 1):
            beta[j] = rho[j] * np.dot(glist[j], d)
            d = d + qlist[j] * (alpha[j] - beta[j])

        info(" @MINIMIZE: Second L-BFGS loop recursion completed", verbosity.debug)
        d = -1.0 * d.reshape(d0.shape)

    else:
        info(
            " @MINIMIZE: Skipped direction update; direction * gradient insufficient",
            verbosity.debug,
        )
        # d = d0
        d = -1.0 * d_x

    d0[:] = d
    info(" @MINIMIZE: Updated search direction", verbosity.debug)


# Damped BFGS to use in NEB. Has no line search and no TRM.
def Damped_BFGS(x0, fdf, fdf0, hessian, big_step):
    """BFGS, damped as described in Nocedal, Wright (2nd ed.) Procedure 18.2
    The purpose is mostly using it for NEB optimization, but it is capable of
    plain geometry optimization also.
    Written for a DIRECT Hessian B, not for the inverse H.

    Currently it doesn't use min_approx, TRM or any other step determination,
    just the simplest (invhessian dot gradient) step, as in aimsChain.
    The reason is that both LS and TRM require energy, and the energy
    of NEB springs is ill-defined because of all projections that we do.
    This may be improved later, but first we need to have NEB working.

    Inside this function I use flattened vectors, restoring shape only when needed.
    I always keep x0 in the original shape.

    Does one step.

    Arguments:
      x0: initial point
      fdf: function and gradient (mapper)
      fdf0: initial function and gradient values
      hessian: approximate Hessian for the BFGS algorithm
      big_step: limit on step length. It is defined differently
                  compared to other optimization algorithms, take care.

    Returns:
      quality: minus cosine of the (gradient, dx) angle.
               Needed for the step length adjustment.
    """

    info(" @DampedBFGS: Started.", verbosity.debug)
    _, g0 = fdf0
    g0 = g0.flatten()

    # Nocedal's notation
    B = hessian

    # Calculate direction
    # When the inverse itself is not needed, people recommend solve(), not inv().
    info(" @DampedBFGS: sk = np.linalg.solve(B, -g0) ...", verbosity.debug)
    info(
        "              The code randomly crashes here with some versions of Numpy "
        "based on OpenBLAS.\n"
        "              If this happens, use Numpy based on MKL, e.g. from Anaconda.",
        verbosity.debug,
    )
    info("Operands:", verbosity.debug)
    info("%s,  %s" % (type(B), str(B.shape)), verbosity.debug)
    info("%s,  %s" % (type(g0), str(g0.shape)), verbosity.debug)
    sk = np.linalg.solve(B, -g0)
    info(" @DampedBFGS: Calculated direction.", verbosity.debug)

    # Cosine of the (f, dx) angle
    quality = -np.dot(sk / np.linalg.norm(sk), g0 / np.linalg.norm(g0))
    info(" @DampedBFGS: Direction quality: %.4f." % quality, verbosity.debug)

    # I use maximal cartesian atomic displacement as a measure of step length
    maxdispl = np.amax(np.linalg.norm(sk.reshape(-1, 3), axis=1))
    info(" @DampedBFGS: big_step = %.6f" % big_step, verbosity.debug)
    if maxdispl > big_step:
        info(
            " @DampedBFGS: maxdispl before scaling: %.6f bohr" % maxdispl,
            verbosity.debug,
        )
        sk *= big_step / maxdispl

    info(
        " @DampedBFGS: maxdispl:                %.6f bohr"
        % (np.amax(np.linalg.norm(sk.reshape(-1, 3), axis=1))),
        verbosity.debug,
    )

    # Force call
    _, g = fdf(x0 + sk.reshape(x0.shape))
    g = g.flatten()
    # coordinates CHECKED

    # Update hessian
    yk = g - g0
    skyk = np.dot(sk, yk)

    # Equation 18.15 in Nocedal
    theta = 1.0
    Bsk = np.dot(B, sk)
    sBs = np.dot(sk, Bsk)
    # Damped update if rhok isn't sufficiently positive
    if skyk < 0.2 * sBs:
        theta = (0.8 * sBs) / (sBs - skyk)
        info(
            " @DampedBFGS: damping update of the Hessian; "
            "(direction dot d_gradient) is small. "
            "theta := %.6f" % theta,
            verbosity.debug,
        )
        yk = theta * yk + (1 - theta) * Bsk
        skyk = np.dot(sk, yk)
    else:
        info(" @DampedBFGS: Update of the Hessian, no damping.", verbosity.debug)

    info(" @DampedBFGS: (s_k dot y_k) before reciprocating: %e" % skyk, verbosity.debug)
    try:
        rhok = 1.0 / skyk
    except:
        warning(" @DampedBFGS: caught ZeroDivisionError in 1/skyk.", verbosity.high)
        rhok = 1e5

    # Compute BFGS term (eq. 18.16 in Nocedal)
    B += np.outer(yk, yk) * rhok - np.outer(Bsk, Bsk) / sBs

    # If small numbers are found on the diagonal of the Hessian,
    # add small positive numbers. Somewhat dirty solution,
    # but it increased stability in some tests.
    # 1 Ha/Bohr^2 is ~97.2 eV/ang^2.
    eigvals = np.real(np.linalg.eigvals(B))
    if np.any(eigvals < 1e-1):
        info(" @DampedBFGS: stabilizing the diagonal of the Hessian.", verbosity.debug)
        B += 1e-2 * np.eye(len(B))

    return quality


def FIRE(
    x0,
    fdf,
    fdf0,
    v=None,
    a=0.1,
    N_dn=0,
    N_up=0,
    dt=0.1,
    maxstep=0.5,
    dtmax=1.0,
    dtmin=1e-5,
    Ndelay=5,
    Nmax=2000,
    finc=1.1,
    fdec=0.5,
    astart=0.1,
    fa=0.99,
):
    """FIRE algorithm based on
    Bitzek et al, Phys. Rev. Lett. 97, 170201 (2006) and
    Guénolé, J. et al.  Comp. Mat. Sci. 175, 109584 (2020).
    Semi-implicit Euler integration used.
    Done by Guoyuan Liu <liuthepro@outlook.com>, May 2021.

    FIRE does not rely on energy, therefore it is suitable for NEB calculation, where
    the energy is not conservative. Basic principle: accelerate towards force gradient
    (downhill direction) and stop immediately when going uphill.
    Try adjusting dt, dtmax, dtmin for optimal performance.

    Arguments:
        x0: initial beads positions
        fdf: energy and function mapper. call fdf(x) to update beads position and froces
        fdf0: initial value of energy and gradient
        v: current velocity
        a: velocity mixing factor, in the paper it is called alpha
        fa: a decrement factor
        astart: initial a value
        N_dn: number of steps since last downhill direction
        N_up: number of steps since last uphill direction
        dt: time interval
        dtmax: max dt (increase when uphill)
        dtmin: min dt (decrease when downhill)
        finc: dt increment factor
        fdec: dt decrement factor
        Ndelay: min steps required to be in one direction before adjust dt and a
        Nmax: max consecutive steps in uphill direction before trigger exit

    Returns:
        v, a, N, dt since they are dynamically adjusted
    """
    info(" @FIRE being called", verbosity.debug)
    _, g0 = fdf0
    force = -g0

    p = np.vdot(force, v)
    # downhill
    if p > 0.0:
        N_dn += 1
        N_up = 0
        if N_dn > Ndelay:
            dt = min(dt * finc, dtmax)
            a = a * fa
    # uphill
    else:
        N_dn = 0
        N_up += 1
        if N_up > Nmax:
            softexit.trigger("@FIRE is stuck for %d steps. We stop here." % N_up)
        dt = max(dt * fdec, dtmin)
        a = astart
        # correct uphill motion
        x0 -= 0.5 * dt * v
        # stop moving in uphill direction
        v = np.zeros(v.shape)

    # accelerate
    v += dt * force
    # change velocity direction with inertia
    if p > 0.0:
        f_unit = force / np.linalg.norm(force)
        v = (1 - a) * v + a * np.linalg.norm(v) * f_unit
    # update posistion
    dx = dt * v
    # check max dx
    normdx = np.linalg.norm(dx)
    if normdx > maxstep:
        dx = maxstep * dx / normdx
    x0 += dx

    info(" @FIRE: calling a gradient mapper to update position", verbosity.debug)
    fdf(x0)

    return v, a, N_dn, N_up, dt


# Bracketing for NEB, TODO: DEBUG THIS IF USING SD OR CG OPTIONS FOR NEB
def bracket_neb(fdf, fdf0=None, x0=0.0, init_step=1.0e-3):
    """Given an initial point, determines the initial bracket for the minimum
    Arguments:
           fdf: function to minimize
           x0: initial point
           fdf0: value of function at x0
           init_step: initial step size
    """

    # Constants
    gold = 1.618034  # Golden ratio
    glimit = 100.0  # Limit for magnification of parabolic fit step
    tiny = 1.0e-20  # Prevent division by zero

    # x0: initial point of evaluation (e.g. initial atomic position)
    # ax, bx, cx: bracketing points with ax < bx < cx
    # fa, fb, fc: value of function at ax, bx, cx
    if fdf0 is None:
        fdf0 = fdf(x0)
    ax = x0
    fa = fdf0
    bx = x0 + init_step
    fb = fdf(bx)[1]
    info(" @BRACKET: Started bracketing", verbosity.debug)
    info(" @BRACKET: Evaluated first step", verbosity.debug)

    # Switch direction to move downhill, if necessary, and rearrange
    if fb > fa:
        tmp = ax
        ax = bx
        bx = tmp
        tmp = fb
        fb = fa
        fa = tmp

    # Initial guess for third bracketing point
    cx = bx + gold * (bx - ax)
    fc = fdf(cx)[1]
    info(
        " @BRACKET: Evaluated initial bracket: (%f:%f, %f:%f, %f:%f)"
        % (ax, fa, bx, fb, cx, fc),
        verbosity.debug,
    )

    # Loop until acceptable bracketing condition is achieved
    # u is a point between two of the bracketing points
    # Use parabolic extrapolation to find u. "tiny" prevents possible division by zero
    while fb > fc:
        r = (bx - ax) * (fb - fc)
        q = (bx - cx) * (fb - fa)
        u = bx - ((bx - cx) * q - (bx - ax) * r) / (
            2.0 * math.copysign(max(abs(q - r), tiny), (q - r))
        )  # Point from parabolic fit
        ulim = bx + glimit * (
            cx - bx
        )  # Limit for parabolic fit point; *Can test various possibilities*

        # Find minimums between b and c or a and u
        # If parabolic fit unsuccessful, use default step magnification
        # Otherwise:
        # - Parabolic fit between c and its allowed limit
        # - Limit u to maximum allowed value
        # - Use default magnification
        if ((bx - u) * (u - cx)) > 0.0:
            fu = fdf(u)[1]
            info(" @BRACKET: Evaluated new bracket point", verbosity.debug)
            # Found minimum between b and c?
            # -b------u-----c shift:
            # -a------b-----c
            if fu < fc:
                ax = bx
                bx = u
                fa = fb
                fb = fu
                info(
                    " @BRACKET: Bracketing completed: (%f:%f, %f:%f, %f:%f)"
                    % (ax, fa, bx, fb, cx, fc),
                    verbosity.debug,
                )
                return (ax, bx, cx, fb)
                # minimum between a and u?
            # -a-----b-----u-----c shift:
            # -a-----b-----c
            elif fu > fb:
                cx = u
                fc = fu
                info(" @BRACKET: Bracketing completed", verbosity.debug)
                return (ax, bx, cx, fb)
                # parabolic extrapolation was not successful.
                # Use golden value (initial guess, default magnification).
            u = cx + gold * (cx - bx)
            fu = fdf(u)[1]
            info(" @BRACKET: Evaluated new bracket point", verbosity.debug)
        # c < u (parabolic fit) < ulim?
        elif ((cx - u) * (u - ulim)) > 0.0:
            fu = fdf(u)[1]
            info(" @BRACKET: Evaluated new bracket point", verbosity.debug)
            # minimum between c and u+gold(u-cx)?
            # -c----u----u+gold(u-cx) shift:
            # -b----c----u
            if fu < fc:
                bx = cx
                cx = u
                u = cx + gold * (cx - bx)
                fb = fc
                fc = fu
                fu = fdf(u)[1]
                info(" @BRACKET: Evaluated new bracket point", verbosity.debug)
        # u >= ulim? limit u.
        elif ((u - ulim) * (ulim - cx)) >= 0.0:
            u = ulim
            fu = fdf(u)[1]
            info(" @BRACKET: Evaluated new bracket point", verbosity.debug)
        # reject parabolic u, use golden value (default magnification)
        else:
            u = cx + gold * (cx - bx)
            fu = fdf(u)[1]
            info(" @BRACKET: Evaluated new bracket point", verbosity.debug)

        # Shift points, so that points to continue are ax, bx, cx
        ax = bx
        bx = cx
        cx = u
        fa = fb
        fb = fc
        fc = fu

    info(
        " @BRACKET: Bracketing completed: (%f:%f, %f:%f, %f:%f)"
        % (ax, fa, bx, fb, cx, fc),
        verbosity.debug,
    )
    return (ax, bx, cx, fb)


# Minimize using only forces; for NEB
def min_brent_neb(fdf, fdf0=None, x0=0.0, tol=1.0e-6, itmax=100, init_step=1.0e-3):
    """Given a maximum number of iterations and a convergence tolerance,
    minimizes the specified function
    Arguments:
           x0: initial x-value
           fdf: function to minimize
           fdf0: initial function value
           tol: convergence tolerance
           itmax: maximum allowed iterations
           init_step: initial step size
    """

    # Initializations and constants
    gold = 0.3819660
    zeps = 1e-10
    e = 0.0  # Step size for step before last
    d = 0.0

    (ax, bx, cx, fb) = bracket_neb(fdf, fdf0, x0, init_step)

    # Set bracket points
    if ax < cx:
        a = ax
    else:
        a = cx
    if ax > cx:
        b = ax
    else:
        b = cx

    # Initial points to evaluate
    # f* is evaluation of arbitrary function
    # x = point with least function value so far
    # w = point with 2nd least function value
    # v = previous value of w
    # u = point at which function was evaluated most recently
    # d = used to determine step size
    x = w = v = bx
    fw = fv = fx = fb

    # Main loop
    j = 1
    while j <= itmax:
        # Determine tolerance
        xm = 0.5 * (a + b)
        tol1 = tol * abs(x) + zeps
        tol2 = 2.0 * tol1

        # Test for satisfactory completion: |b-a|<=tol*abs(x)
        if abs(x - xm) <= (tol2 - 0.5 * (b - a)):
            xmin = x
            return xmin, fx

        # Complete an iteration if error is greater than tolerance
        # and construct parabolic fit from parameters
        # assumption:
        # function given at points x,w,v is approximately parabolic near the minimum

        # x-p/q is point, where derivative of fitted parabola is zero
        if abs(e) > tol1:
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = 2.0 * (q - r)

            if q > 0.0:
                p = -p
            q = abs(q)
            etmp = e
            e = d

            # Determine acceptability of parabolic fit
            # movement must be less than half of the movement of the step before last:
            # |p/q| <= |0.5*etmp| (old e)
            if (
                (abs(p) >= abs(0.5 * q * etmp))
                or (p <= (q * (a - x)))
                or (p >= (q * (b - x)))
            ):
                # step into larger of the two segments
                if x >= xm:
                    e = a - x
                else:
                    e = b - x
                # conditions for parabolic fit not fulfilled, new d: golden section
                d = gold * e

            # Take parabolic step (conditions fulfilled)
            # minimum of fitted parabola at point u
            else:
                d = p / q
                u = x + d
                if ((u - a) < tol2) or ((b - u) < tol2):
                    d = abs(tol1) * (xm - x) / abs(xm - x)
        # if abs(e) <= tol1 (first step in any case)
        else:
            # step into larger of the two segments
            if x < xm:
                e = a - x
            else:
                e = b - x
            d = gold * e
        if abs(d) >= tol1:
            u = x + d
        else:
            u = x + abs(tol1) * d / abs(d)
            # one function evaluation per iteration, derivative is computed, too..
            # include count...
        fu = fdf(u)[1]
        # order for next step: a < u (later x) < b
        if fu <= fx:
            if u >= x:
                a = x
            else:
                b = x

            # Reassign bracketing points
            v = w
            w = x
            x = u
            fv = fw
            fw = fx
            fx = fu
        else:
            if u < x:
                a = u
            else:
                b = u
            if fu <= fw or w == x:
                v = w
                w = u
                fv = fw
                fw = fu
            elif (fu <= fv) or (v == x) or (v == w):
                v = u
                fv = fu
        j += 1

    # Exit if maximum number of iterations exceeded
    xmin = x
    return xmin, fx


# L-BFGS without line search; WARNING: UNSTABLE


def L_BFGS_nls(
    x0,
    d0,
    fdf,
    qlist,
    glist,
    fdf0=None,
    big_step=100,
    tol=1.0e-6,
    itmax=100,
    init_step=1.0e-3,
    m=0,
    k=0,
):
    """L-BFGS minimization without line search
    Does one step.
        Arguments:
            fdf: function and gradient
            fdf0: initial function and gradient value
            d0: initial direction for line minimization
            x0: initial point
            qlist: list of previous positions used for reduced  Hessian^-1 construction
            glist: list of previous gradients used for reduced  Hessian^-1 construction
            m: number of corrections to store and use
            k: iteration (MD step) number
            big_step: limit on step length
            tol: convergence tolerance
            itmax: maximum number of allowed iterations
            init_step: initial step size
    """

    # Original function value, gradient, other initializations
    # zeps = 1.0e-10
    if fdf0 is None:
        fdf0 = fdf(x0)
    f0, df0 = fdf0
    n = len(x0.flatten())
    dg = np.zeros(n)
    g = df0
    x = np.zeros(n)
    # linesum = np.dot(x0.flatten(), x0.flatten())
    alpha = np.zeros(m)
    beta = np.zeros(m)
    rho = np.zeros(m)
    q = np.zeros(n)

    # Initial line direction
    xi = d0
    dg = df0

    # Step size
    stepsize = np.linalg.norm(d0.flatten())

    # First iteration; use initial step
    if k == 0:
        scale = 1.0
        while (
            np.linalg.norm(g.flatten()) >= np.linalg.norm(df0.flatten())
            or np.isnan(np.linalg.norm(g.flatten()))
            or np.isinf(np.linalg.norm(g.flatten()))
        ):
            x = np.add(
                x0,
                (scale * init_step * d0 / np.linalg.norm(d0.flatten())),
            )
            scale *= 0.1
            fx, g = fdf(x)
    else:
        # Scale if attempted step is too large
        if stepsize > big_step:
            d0 = big_step * d0 / np.linalg.norm(d0.flatten())
            info(" @MINIMIZE: Scaled step size", verbosity.debug)

        x = np.add(x0, d0)
        print("step size:", np.linalg.norm(d0.flatten()))
        fx, g = fdf(x)

    info(" @MINIMIZE: Started L-BFGS", verbosity.debug)
    info(" @MINIMIZE: Updated gradient", verbosity.debug)

    # Update line direction (xi) and current point (x0)
    xi = np.subtract(x, x0)

    # Build list of previous positions
    if k < m:
        qlist[k] = xi.flatten()
    else:
        qlist = np.roll(qlist, -1, axis=0)
        qlist[m - 1] = xi.flatten()

    # Update current point
    x0 = x

    # Compute difference of gradients
    q = g.flatten()
    dg = np.subtract(g, dg)

    # Build list of previous gradients
    if k < m:
        glist[k] = dg.flatten()
    else:
        glist = np.roll(glist, -1, axis=0)
        glist[m - 1] = dg.flatten()

    # Determine bounds for L-BFGS 'two loop recursion'
    if k < (m - 1):
        bound1 = k
        bound2 = k + 1
    else:
        bound1 = m - 1
        bound2 = m

    # Begin two loop recursion:
    # First loop
    for j in range(bound1, -1, -1):
        rho[j] = 1.0 / np.dot(glist[j], qlist[j])
        alpha[j] = rho[j] * np.dot(qlist[j], q)
        q = q - alpha[j] * glist[j]

    info(" @MINIMIZE: First L-BFGS loop recursion completed", verbosity.debug)

    # Two possiblities for scaling: using first or most recent
    # members of the gradient and position lists
    hk = np.dot(glist[bound1], qlist[bound1]) / np.dot(glist[bound1], glist[bound1])
    # hk = np.dot(glist[0], qlist[0]) / np.dot(glist[0], glist[0])
    xi = hk * q

    # Second loop
    for j in range(0, bound2, 1):
        beta[j] = rho[j] * np.dot(glist[j], xi)
        xi = xi + qlist[j] * (alpha[j] - beta[j])

    # Update direction xi
    xi = -xi.reshape(d0.shape)

    info(" @MINIMIZE: Second L-BFGS loop recursion completed", verbosity.debug)
    info(" @MINIMIZE: Updated search direction", verbosity.debug)

    return (x, fx, xi, qlist, glist)


def nichols(f0, d, dynmax, m3, big_step, mode=1):
    """Find new movement direction. JCP 92,340 (1990)
    IN    f       = physical + spring forces        (n,)
          d       = dynmax eigenvalues
          dynmax  = dynmax       (n x n-m) with m = # external modes
          m3      = mass vector
    OUT   DX      = displacement in cartesian basis

    INTERNAL
          ndim = dimension
          f    = forces
          d    = hessian eigenvalues
          w    = hessian eigenvector (in columns)
          g    = gradient in cartesian basis  #Note in trm g is the force vector
          gE   = gradient in eigenvector basis
          DX   = displacement in cartesian basis
          DXE  = displacement in eigenvector basis
    """

    # Resize
    ndim = f0.size
    shape = f0.shape
    f = (
        f0.reshape((1, ndim)) / m3.reshape((1, ndim)) ** 0.5
    )  # From cartesian base to mass-weighted base

    # Change of basis to eigenvector space
    d = d[:, np.newaxis]  # dimension nx1
    gEt = -np.dot(f, dynmax)  # Change of basis  #
    gE = gEt.T  # dimension (n-m)x1
    # The step has the general form:
    # d_x[j] =  alpha *( gE[j] )  / ( lambda-d[j] )

    if mode == 0:
        # Minimization
        alpha = 1.0
        lamb = 0.0

        d_x = alpha * (gE) / (lamb - d)

        if d[0] < 0 or np.vdot(d_x, d_x) > big_step**2:
            lamb = d[0] - np.absolute(gE[0] / big_step)
            d_x = alpha * (gE) / (lamb - d)

    elif mode == 1:
        if d[0] > 0:
            if d[1] / 2 > d[0]:
                alpha = 1
                lamb = (2 * d[0] + d[1]) / 4
            else:
                alpha = (d[1] - d[0]) / d[1]
                lamb = (
                    3 * d[0] + d[1]
                ) / 4  # midpoint between b[0] and b[1]*(1-alpha/2)

        elif d[1] < 0:  # Jeremy Richardson
            if d[1] >= d[0] / 2:
                alpha = 1
                lamb = (d[0] + 2 * d[1]) / 4
            else:
                alpha = (d[0] - d[1]) / d[1]
                lamb = (d[0] + 3 * d[1]) / 4

        # elif d[1] < 0:  #Litman for Second Order Saddle point
        #    alpha = 1
        #    lamb = (d[1] + d[2]) / 4
        #    print 'WARNING: We are not using the standar Nichols'
        #    print 'd_x', d_x[0],d_x[1]
        else:  # Only d[0] <0
            alpha = 1
            lamb = (d[0] + d[1]) / 4

        d_x = alpha * (gE) / (lamb - d)
    # Some check or any type of reject? ALBERTO

    DX = np.dot(dynmax, d_x)  # From ev base to mass-weighted base
    DX = DX.reshape(shape)
    DX = np.multiply(DX, m3 ** (-0.5))  # From mass-weighted base to cartesion base

    return DX


def Powell(d, Dg, H):
    """Update Cartesian Hessian using gradient.
    Input: d  = Change in position
           Dg = change in gradient
            H = Cartesian Hessian

    Output: H = Cartesian Hessian"""

    ddi = 1 / np.dot(d, d)
    y = Dg - np.dot(H, d)
    H += ddi * (np.outer(y, d) + np.outer(d, y) - np.dot(y, d) * np.outer(d, d) * ddi)
    return H
