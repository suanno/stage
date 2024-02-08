# Implicit Euler scheme for solving TDGL
$$\partial_t u = C u - u^3 + \Delta u$$
In Fourier space is
$$\partial_t U = CU - F[u^3] - q^2 U$$
If you use **Backward euler** formula for solving the equation:
$$\partial_t U(t) \simeq (U(t)-U(t-h))/h$$
or equivalently
$$U(t) \simeq U'(t)h + U(t-h)$$
Then, using the TDGL equation to explicitly write $U'(t) = \partial_t U(t)$
$$U(t) = -q^2U(t) + CU(t) - F[u^3](t) + U(t-h)$$
So bringing all terms $U(t)$ to the left side
$$U(t)(1+q^2-C) = U(t-h) - hF[u^3]$$
And we get the **Implicit Euler equation**
$$U(t)=\frac{U(t-h) - hF[u^3](t)}{1+q^2-C}$$

[WHY!!!]In the TUNG algorithm, $F[u^3](t)$ is changed with $F[u^3](t-h)$

## Explicit / Implicit Euler
The explicit algorithm is simply $\partial_t u(t) \simeq (u(t+h)-u(t))/h$ from which you get
$$u(t+h) = u(t) + u'(t)h$$

While the implicit is
$$u(t) = u'(t)h+u(t-h)$$

It is called implicit, because this expression cannot be exploited as it is, as it requires the knowledge of $u'(t)$ to calculate $u(t)$ and the expression that gives $u'(t)$, that is the equation, contains the value $u(t)$ that you want to find.

But from this formula you can still get an exploitable formula.

### Implicit Euler for a LINEAR differental equation
We say that a linear equation is an equation with this shape
$$\partial_t u(x,t) = L[\partial_x]u(x,t)$$
where $L$ is a linear operator, like $L = \partial_x^2 + C$.

Following the same steps done for the TDGL (that is NOT linear, due to the non-linear term $u^3$) you find out
$$U(t) = \frac{U(t-h)}{1-L[iq]h}$$
Where $U(t,q) = F[u(t,x)]$ is the Fourier transform.

## Advantage of Implicit methods: Stability!
If you write the linear differental equation in fourier space, you get
$$\dot{U}_q(t) = L[iq]U_q(t)$$
where $L[iq]$ is a matrix [in 1D a number] containing only numbers (not derivative symbols) and depends parametrically on $q$ (and not on the variable t).

So this is the equation of a **dynamical system with _continuos time_**.
If all the modes of the dynamical system are stable (or equivalently if any mode $q$ that solves the differential equation does not explode) then **all the eigenvalues** of the matrix $L[iq]$ are **negative**.

Now consider the _implicit euler formula_
$$U(t) = \frac{U(t-h)}{1-L[iq]h}$$

This is a **dynamical system with discrete time**
$$U_{n+1} = \frac{U_{n-1}}{1-L[iq]h} = MU_{n-1}$$
where the matrix of the dynamical system is $M = (1-L[iq]h)^{-1}$.

All the modes of the system are stable (so they do not explode) only if all the eigenvalues of $M$ have magnitude less than 1.
And this is so, as the eigenvalues $\lambda_n$ of matrix $L[iq]$ are all $\lambda_n < 0$ and so $|(1-\lambda_n h)|>0$ so its reciprocal is less than 1, **for any value of h**.

#### Interpretation of the last result
Using the implicit euler formula, we are approximating the dynamics of the continuos time system, with the one of the discrete time system. As the one of the discrete time system **has stable modes for any choice of h**, then 

_**"The numerical solution we get from the computer is a sum of modes that do not explode (as time goes on), INDEPENDENTLY on the choice of $h$ (even if $h$ is very LARGE)."**_

This is what we mean by **stable/non-stiff algorithm**.

The Implicit Euler method is stable (at least for **linear** differential equations), while the explicit method isn't stable, because we loose the stability if $h$ is large.

#### Example of the instability of Explicit Euler
The explicit Euler algorithm, in Fourier space is
$$U(q,t+h) = hL[iq]U(q,t)+U(q,t) = (hL[iq]+1)U(q,t)$$
so we have the dynamical system with discrete time
$$U_{n+1} = MU_n$$
where $M = (1+hL[iq])$.

Now you know that the eigenvalues of the matrix $L[iq]$ are negative, but this is not sufficient to tell that the eigenvalues of $M$ have magnitude less than 1. It is so only if you choose $h>-\frac{2}{\lambda_{min}}$ where $\lambda_{min}$ is the smallest eigenvalue of $L[iq]$ and remember that its eigenvalues are negative (so the $-$ sign gets the fraction positive).

This means that, for sufficiently large values of $h$, there are some modes of the discrete time dynamical system that explode!
**"This means that, if you use Explicit Euler with a large step $h$, the numerical solution you get will contain some modes that explode (as time goes on)."**

For this reason, we say that Explicit Euler is an **unstable/stiff algorithm**.

## Comments on stability
**Intuitively**, the reason why Explicit Euler is unstable is that the errors (truncation + roundoff) that are present at each step, accumulate in a somewhat "constructive" way. While it is not the case for an implicit algorithm.

A naive way of extimating /**the order of** an error associated to the value of the solution at **long time** $t$ is to assume that the error is the same at each step and so the error at a long time is
$$NO(h^n) = \frac{L}{h}O(h^n) = LO(h^{n-1})$$
so the error is 1 order of magnitude smaller.

Note that it makes sense to use this criteria only if the algorithm is stable (at least for the chosen value of $h$).