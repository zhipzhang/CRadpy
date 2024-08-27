import numpy as np
import astropy.units as u
import CRadpy
import matplotlib.pyplot as plt

energy = np.logspace(-3, 3, 400) * u.TeV
def function(energy):
    amplitude = 1e36 * u.Unit("1/eV")
    E0        = 1 * u.TeV
    index     = 2.1
    cutoff    = 13 * u.TeV
    return amplitude * np.power(energy/E0, -1 * index) * np.exp( -1 * np.power(energy/cutoff, 1))

density = function(energy).to("erg^-1")
energy = energy.to('erg')
distance = 1000 * CRadpy.constants.pc_to_cm;
rad = CRadpy.Radiation()
rad.SetElectronDistribution(energy, density)
rad.AddBlackBodyPhotons(2.725)
rad.SetB_Field(1e-4);
spectrum_energy = np.logspace(-1,14,500) * CRadpy.constants.eV_to_erg
rad.SetDistance(distance)
rad.CalculateDifferentialSpectrum(spectrum_energy)
inverse_comton = rad.GetICSpectrum()
synchrotron    = rad.GetSynSpectrum()
plt.ylim(1e-40, 1e-7)
#plt.loglog(spectrum_energy, synchrotron * spectrum_energy**2, label = "synchrotron")
plt.loglog(spectrum_energy, inverse_comton * spectrum_energy**2, label = "inverse_compton")
plt.legend()
plt.show()


