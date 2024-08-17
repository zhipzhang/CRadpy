import CRadpy
import astropy.units as u
import numpy as np



energy = np.logspace(-3, 3, 400) * u.TeV
def function(energy):
    amplitude = 1e36 * u.Unit("1/eV")
    E0        = 1 * u.TeV
    index     = 2.1
    cutoff    = 13 * u.TeV
    return amplitude * np.power(energy/E0, -1 * index) * np.exp( -1 * np.power(energy/cutoff, 1))

energy = energy.to("erg")
density = function(energy).to("erg^-1")

rad = CRadpy.Radation()
rad.SetElectronDistribution(energy, density)
rad.AddBlackBodyPhotons(2.725)
