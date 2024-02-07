# 1D Allen–Cahn map w/out $\Delta\Phi$ term
$$\partial_t\phi = C\phi - \phi^3$$

## Stationary (or even called _equilibrium_) states
The stationary solutions of the equation ($\partial_t \Phi = 0$) are **homogeneous** and are:
- $\Phi = \pm\sqrt{C}$ if $C>0$. Then there is even $\Phi = 0$ but that state in unstable.
- $\Phi = 0$ if $C<0$ that is the only stationary solutiona and its stable.

**Stable** means that the system will (eventually) move towards one of this stationary states (two types of motivations for why this happens are given below).

We can say that thinking in two ways:
- The A.L. equation has the form $\partial_t \Phi = -\frac{dU}{dx}$ where $U(x) = -C\Phi^2/2 + \Phi^4/4$ and this potential has two minima if $C>0$ and one minima if $C<0$. And this kind of equation is typical of dissipative systems, so we expect the system to move closer to one of the minima.
- We ask $\partial_t\Phi = 0$ to find the stationary solutions and then we look at the sign of the derivative, finding out that it is positive if $\Phi < \sqrt{C}$ and negative if $\Phi > \sqrt{C}$.
So the state $\Phi(x) = \sqrt{C}\quad\forall x$ is stable.

## The role of the non-linear term $\Phi^3$
If we include the $\Delta \Phi$ term and then consider to have an **initial state** that is close to zero (zero everywhere), then the non-linear term can be neglected from the equation.

BUT, if you do so, you find out that the long-$\lambda$ modes are unstable: their amplitude explodes in a characteristic time $\sim \frac{1}{C}$.
So the term $\Phi^3$ is no more negligible.
If you consider this term when $t\gtrsim \frac{1}{C}$, you get a contribute to $\partial_t \Phi$ that is $-\Phi^3$ and so the _peaks emerging_ become _plateaus_.
The hight of plateaus is $\pm\sqrt{C}$ for fixed in time $C$. This is not strange as this is the value of $\Phi$ in the two _homogenoeus_ stationary states.

Up and down plateaus are formed, those are respectively called anti-kinks and _kinks_.

## At very long times
As the peaks become plateaus, $\Phi(x)$ becames **slowly varying in Space** and so the A.L. equation becomes well approximated by the one without the spatial derivative. In that case we know that (if $C>0$) the system has only two homogenous stable stationary states **and it moves towards one of them**.

So, **wheter if the system is extended (large) or small** ([?][for small you need to check experimentally what happens as you expect all modes to die and so to stay in $\Phi = 0$)] **it will _eventually_ reach one of the two homogeneous states**.
But this can took a very long time [in 1D the kinks should increase size as $\log t$].

## Initial, intermediate and equilibrium behaviour [?]
- Initial ($t < 1/C$): The system is in a random state but the amplitude of the modes contituting it are close to zero, so you can describe the system neglecting the non-linear term. **You have long-$\lambda$ modes rising exponentially**.
- Intermediate ($t\sim1/C$): The peaks grow a little bit and the $\Phi^3$ term is no more neligible. So you have plateaus formation, but you have still a profile varying in time.
- Equilibrium ($t>10/C$): Here the system is approaching one of the two homogeneous states [?] and the profile varies very slowly in time.

## Time scales
By controlling $C$ we can control **what are the unstable modes**. As the unstable modes explode with a characteristic time that is $\tau = 1/C$, then this is a characteristic time for the change of $C$ to be effective.
So we have to take into account this when changing $C(t)$ in time.