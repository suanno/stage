## Numerical evidence of $u(nT)\sim n^{-\frac12}$ asymptotic behavior

### Pre-asymptotic behavior
Analytically, we anticipate that 
$$\frac{1}{u^2(nT)} = \frac{1}{u_0^2} + 2nI_0(T)$$
where $I_0(T)=\int_0^T dt e^{2\sin{\frac{2\pi t}{T}}}$, and numerically evaluating the integral suggests linearity in $T$ 
$$I_0(T)\simeq 2.28T\text{,\quad if }C(t) = A\sin(2\pi t/T)\text{, with } A=1$$

As $n$ becomes large, the term with $n$ dominates the right-hand side, yielding the **asymptotic behavior**
$$u(nT)\simeq \frac{1}{2I_0(T)}n^{-\frac12}$$
Quantitatively, this behavior holds true when $2nI_0(T) >> \frac{1}{u_0^2}$, which occurs when $n>>\tilde{n}$ with $\tilde{n} = \frac{1}{2*2.28*T*u_0^2}$.

The pre-asymptotic behavior can be eliminated **by increasing $T$ or $u_0$**. However, significantly increasing $u_0$ could render the term $\frac{1}{u_0^2}$ negligible even at low times $t$, leading the analytical solution to suggest $u(nT)$ as the reciprocal of $2nI_0(T)$, which, for low $n$ (and thus $t$), approaches zero. Consequently, computational overflow errors are expected due to division by numbers close to zero.

Here we can see the pre-asymptotic behavior with $u_0 = 0.2$ ![here](../codes_tdgl/codes_tdgl/1D/Plots/at%20long%20times%20becames%20constant%202%20u0=0.2.png)

### Asymptotic behavior
The asymptotic behavior is observable in the numerical solutions ![here](../codes_tdgl/codes_tdgl/1D/Plots/at%20long%20times%20becames%20constant%202%20u0=10%20with%20analytical.png), where we adopted a large $u_0=10$ in order to get rid of the pre-asymptotic behaviour

**However:**
- The asymptotic behavior persists **until** a certain $n_{max}$, which decreases with the time step.
- In the region where the asymptotic behavior is visible, it is **up-shifted** in _log-log_ scale compared to the analytically expected curve.

### Post-asymptotic behavior
Following the region where the asymptotic behavior holds, $u(nT)$ converges to a constant value that increases with $dt$.
While $u(nT)$ stabilizes, the curve $u(t)$ evaluated at times different from $nT$ exhibits oscillations ![as shown here](../codes_tdgl/codes_tdgl/1D/Plots/long%20time%20oscillation%20when%20u(nT)%20is%20constant.png)

Increasing the period $T$ with fixed $dt$ (and so increasing the accuracy with which the control function $C(t)$ is encoded numerically) lowers the $n_{max}$ when the _saturation behaviour_ begins. ![here](../codes_tdgl/codes_tdgl/1D/Plots/varying%20T%20fixed%20dt.png) (where $u_0=10$).

It's worth noting that increasing the period $T$ makes $u(t)$ go faster to zero, as $u(t)\sim\frac{1}{I_0{T}}\sim\frac{1}{T}$ is the asymptotic behaviour expected.
So it's not strange that increasing $T$ anticipates the pre-asymptotic behaviour.

The _lenght_ of the _saturation tail_ seems to be infinite as shows this $tspan = 10^4$ evolution ($u_0 = 10$)
![here](../codes_tdgl/codes_tdgl/1D/Plots/long%20tail.png)

### Resolution of $C(t)$ encoding
Keeping $T/dt = 1000$ **fixed** you can see that the deviation of $C(nT)$ from the expected value (zero) is _the same_ independently on the couple ($T, dt$) you choose.
You see that, because the _curves are overlapping_.
![see here](../codes_tdgl/codes_tdgl/1D/Plots/C(t)%20resolution%20is%20T%20over%20dt.png)