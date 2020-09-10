# Modeling-Simulation

This is a simulation for fixed bed reactor of Oxidative dehydration ethane to ethylene by Matlab.

## Chemical process description

   Ethylene is usually produced by steam cracking of saturated hydrocarbons and is produced in recent developments by the cracking of catalytic fluidized bed and catalytic dehydrogenation of ethane,
although these processes require high reaction temperatures as well as problems with coke formation and uncontrollable side reactions. Now, the oxidative dehydrogenation technology has two major advantages:
1. Low Energy comsumption
2. Less COx sides produced and coke formation

But two of the leading challenges of this technology can be mentioned below:
1. Development of suitable catalyst capable of lower selectively ethane to produce ethylene at low temperature
2. Desgin of industrial reactor

### Oxidative dehydration reactions
![ODH Reactions](/img/Rxns.JPG)

Solve systems of time variante PDEs:
 Time domain | Space domain
------------ | -------------
  t variable | r ,z variables
ODE Runge-Kutta | Orthogonal Collocation 

### Gas Mass Balance
![Gas](/img/GasMass.JPG)

### Gas Energy Balance
![Gas](/img/GasEnergy.JPG)

### Solid Mass Balance
![Solid](/img/SolidMass.JPG)

### Solid Energy Balance
![Solid](/img/SolidEnergy.JPG)

### References
[Kinetic of ODH & Exp data](https://doi.org/10.1016/j.cej.2014.04.042)

[Reactor modeling](https://doi.org/10.1080/00986445.2017.1396538)

[Just help to code as an example](https://doi.org/10.1016/j.cherd.2019.04.028)
