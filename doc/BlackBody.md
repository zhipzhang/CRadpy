In the `Utils.h`, we introduced the method of `AddThermalTargetPhotons`. It will auto generate the photon number density from the Blackbody Radiation.

# BlackBody Radiation

In the textbook, we can find the following equation:
$$
u(\nu, T) = \frac{8\pi h\nu^3}{c^3}\frac{1}{exp(h\nu/kT) -1}
$$

- it represent the energy density per frequency per volume.
- The units of `u` is `erg/cm^3/Hz`

In order to convert the above equation to `1/cm^3/erg`, first we need understand:
$$
u(\nu)d\nu = u(E)dE\\
E = h\nu
$$
So we can get the $u(E)$:
$$
u(E) = u(\nu)/h
$$

- The units of `u(E)` is `erg/cm^3/erg`

So the last things we need is to divided the energy E :
$$
n(E)dE = u(E)d(E)/E = \frac{u(\nu)}{h^2 \nu}dE 
$$
So we can get the number of density :
$$
n(E) = u(\nu)/(h^2 \nu) = \frac{8\pi \nu^2}{hc^3} \frac{1}{exp(h\nu/kT) -1}
$$

## GreyBody Radiation

Considering the `graybody`, which have a emissivity less than 1. We can get the emissivity from the energy density:
$$
\epsilon = \frac{\rho_e \cdot 15 h^3 c^3}{8\pi^5 k^4 T^4}
$$
 for blackbody : $\epsilon = 1$

Using the $\epsilon$ and $n(E)$, we can get the distribution of GreyBody