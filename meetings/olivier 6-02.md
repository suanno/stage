# Behaviour of the system at long time $t$. There are onlt 3 stationary states!
At very long time $t$ the system **will eventually reach**
- **one of the two homogeneous state** $\Phi = \pm\sqrt{C}$ (if $C>0$) 
- or the state $\Phi(x) = \pm\sqrt{C}\tanh{(x/\sqrt{2})}$

wheter the system is small or extended (respect to the characteristic lenght, that is the smallest unstable wavelengh $\lambda_{min}=2\pi C$).

## Two domains stationary state
Those **3** states (actually 4 if you consider $\pm$ the 3rd state), for $C>0$, are **the only stationary solutions** of A.C. equation. And the "two domains" stationary state is no more a stationary state if you adopt PBC boundaries.

## More domains
There cannot be solutions with more than 2 infinite domains, because the **kinks and anti-kinks will eventually all disappear (except one if they are not all paired)** but it takes a lot of time to see this, as:

If you consider to have only one anti-kink and next one kink, they will come closer and eventually annichilate. But this process is slow as the speed they get closer is
        $$\dot{d(t)} \sim e^{-2\sqrt{2}d(t)}$$
        where $d$ is the kink and anti-kink pair distance. So small domains tend to _slowly_ disappear and they disappear **in order** from the smaller to the bigger, as the disappearing speed **decreases exponentially with the size**.

The calculation of the speed is done in the limit $d>>1$, but if you use the formula for any distance $d$, you get a time for the annichilation of the pair that is:

$$T\sim e^{2\sqrt{2}d(t=0)}$$
If you consider the solution $d(t)$ of the equation describing the velocity with which the pair kink-antikink become closer, you find that the distance, and so the dimensions of domains, changes size (it decreases or increases, it depends if you consider a small domain in a sea or a big domain near close opposite domains) **logaritmically**
$$d\sim\log (e^{2\sqrt{2}d(0)}-24t)$$
So the size of the domains changes very slowly, until the domain is very small, at that stage it changes fast (you can see it better from the expression of the derivative $\dot{d}$).