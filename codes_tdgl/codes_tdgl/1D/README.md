# TUTORIAL: solve TDGL
- Execute "flatinit" to have a flat initial state or "datainit" to have a random state very close to $u=0$ .
        
        ./flatinit <lattice points> <flat value>
        ./datainit <lattice points>
    Both file will create a "fileinit.dat" which has the _seed_ and the _number of lattice points_ in the first line.
    In the following lines it contains the $u(x,0)$ values.

- Execute "tdglfd", it evolves the state with **Explicit euler** (_fd_ stands for _finite difference_) for a time $tspan$ with varying in time $C(t) = A\sin(2\pi t/T)$ (where we adopt the convention $T/2 = -1$ means that C is costant and $C=A$)
        
        ./tdglfd <tspan> A T/2
        

    - If the initial state is not smooth (it was generated by "datainit") you _have to_ execute this code _at least for a short period_, because "tdgl" needs a smooth initial state [cit. TUNG]
    - While if the initial state is flat, then you don't need to evolve the system with this code. Anyway you still need to execute it for a timespan equal to zero, in order to create the file "tdgl_result.dat" that _is necessary_ for running "tdgl".

    The code creates a file "tdgl_result.dat" that contains in the first row
    - the number of lattice points N
    - the time $t$ of the state (the following lines contain the state $(x, u(x,t))$)
    - dx and dt integration steps
    - the seed (this information can be used when reading the data to label a plot or to perform again the same evolution for _troubleshooting_)
    - the amplitude $A$ and the half-period $T/2$ of the $C(t)$ evolution
    
- Execute "tdgl", this code evolves **progressively** the state in Fourier space. The syntax is the same of "tdglfd"

        ./tdgl <tspan> A T/2
    It needs the existance of the file "tdgl_result.dat" to be run and the initial state must be smooth (because a sinusoidal decomposition must make sense).

    It **updates** the state recorded in "tdgl_result.dat", updating even the time $t$ of the state.
    The new values of $C(t)$ in time and of the space average of $u(x,t)$ are **appended** in the files "fileCout.dat" and "fileAveout.dat", so those files contains information about **the whole** simulation (from $t=0$).

    The files "fileCout.dat" and "fileAveout.dat" are emptied when "flatinit" or "datainit" are executed and the file "fileinit.dat" is overwritten.

## Plotting
- To plot the actual state of the system, you can use "plot_state.ipynb" that even saves the plot in a directory with the _seed_ name, while the name of the image is the _time t_ of the state.

- To plot the space average of $u(x,t)$ at times $t=nT$, where $T$ is the period of $C(t)$, you can use "flat_plots.ipynb".