## Universality
There are a lot of physical systems that behave the same near a critical point, this concept is known as _universality_ and the physical systems behaving the same are called to be in the same _class_.
All systems that have/are
- up-down simmetry
- dissipative (strongly interacting w/ thermostat, so _overdamped_)
- scalar order parameter $\phi(x,t)$
- long wavelenght instability (meaning that the unstable modes, the ones which evolution _eplodes in time_, are the long $\lambda$ ones)

behave, _near a (the ?) critical point_, as the Ising model; so they behave the same.

Those systems, near the critical point, are described by the "Allen–Cahn" or "Time-dependent Gitzburg-Landau equation":

$$\partial_t\phi = -\frac{\delta F}{\delta \phi} = \epsilon\phi - \phi^3 + \Delta \phi$$
where $\epsilon$ is a parameter **that we assume to be able to control**, and for the ising model it depends on temperature $\epsilon \propto (T-Tc)$.

We are concerned more about what the equation tells us rather than the physics, so we are not concerned about free energy and we'll think to the RHS of the above eq.
[In some systems a Free energy can't be defined, but they follow this equation near the critical point]

### Why the evolution of the system is described by an equation of the kind $\dot{x} = -\frac{dU}{dx}$?
We are used to the equation $m\ddot{x} = -\frac{dU}{dx}-\nu \dot{x}$, but if you consider a very dissipative (so _overdamped_) system, you can neglect the inertial term $m\ddot{x}$ and you'll get $\dot{x} = -\frac{1}{\nu} \frac{dU}{dx}$.

## Equilibrium states and dynamics

The first two terms in the RHS determine the equilibrium states, that are **two** _homogeneous_ states (all spins up or all spins down) because of 
the absence of space derivatives.

But in order to describe the dynamics close to the critical point, you need to include the laplacian term. This enables space-dependent solutions where you have magnetic domains (Ising) and so smooth time dependent-solutions.

### Extendend systems [?]
You expect that the system, if you wait long enough, will reach one of the two equilibrium states. So we can say equation describes an **intermediate behaviour** of the system, where _intermediate_ means "long before reaching equilibrium".

[?] But the system eventually reaches an equilibrium states **only if it is small**. That's because **only the long wavelenght modes are unstable**.

### Why we talk about _long-wavelenght unstability_? [?]

In order to be an order parameter, so one that describes the phase transition, the two equilibrium solutions must be both $\phi^* = 0$ **at** the critical point (otherwise it makes no sense expanding the free energy in terms of the order parameter).

As a consequence you can tell that, when you work _close to_ the critical point, not only the deviation from the equilibrium state $\delta \phi = \phi(x,t) - \phi^*$ is small, but even $\phi(x,t)$ is small as $\phi^*$ is small.

Then, if you differentiate the equation A.C. then you can neglect the term $\sim \phi^3$ in favour of the one $\sim \phi$ [Why you do not do this in the A.C. equation???]

So you get $$\partial_t\delta \phi(x,t) = \epsilon\delta \phi - \Delta\delta \phi$$

that in Fourier domain is
$$\partial_t\delta \phi_q(q,t) = (\epsilon-q^2)\delta \phi$$
That has exponential solutions $\delta\phi_q(t) = e^{(\epsilon-q^2)t}$.
So the q-th mode of the scalar field $\phi(x,t)$
- dies zero: if $\epsilon < q^2$
- explodes: if $\epsilon > q^2$

So only small q (and so high $\lambda = \frac{1}{q}$) modes explode.
### Time scale and space scale of unstability
- The space scale of the unstability, meaning the $\lambda$-span of the unstable modes is given by $q^2 < \epsilon$, so $\Delta\lambda \sim \epsilon^{-\frac12}$ ($\Delta q\sim\epsilon^\frac12$).
- The time scale of the unstability, meaning the time it takes for the long-$\lambda$ modes to explode is $\tau = \frac{1}{\epsilon-q^2}\sim\epsilon^{-1}$ using $q^2$ is less relevant than $\epsilon$.


### Rescaling A.L. equation
...

### Kinks, attractors, phenomenology when t is not so far from characteristic $\tau$
...

### Control of 1D A.L. eq