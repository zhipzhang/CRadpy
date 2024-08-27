#include "RadData.h"
#include "nanobind/nanobind.h"
#include "nanobind/eigen/dense.h"
#include "nanobind/nb_defs.h"
#include "Radiation.h"
#include "Eigen/Core"
#include "nanobind/stl/vector.h"
namespace  nb = nanobind;
using namespace nb::literals;
using CRad::Radiation;
using CRad::RadData;
NB_MODULE(CRadpy, m)  
{
    nb::class_<Radiation>(m, "Radiation")
        .def(nb::init())
        .def("SetProtonDistribution", &Radiation::SetProtonDistribution,"energy"_a, "density"_a, "set the spectra of proton ")
        .def("SetElectronDistribution", &Radiation::SetElectronDistribution,"energy"_a, "density"_a, "set the spectra of electron ")
        .def("AddBlackBodyPhotons", &Radiation::AddBlackBodyPhotons, "temperature"_a)
        .def("CalculateDifferentialSpectrum", &Radiation::CalculateDifferentialSpectrum)
        .def("GetTotalSpectrum", &Radiation::GetTotalSpectrum)
        .def("GetICSpectrum", &Radiation::GetICSpectrum)
        .def("GetSynSpectrum", &Radiation::GetSynSpectrum)
        .def("GetProtonSynSpectrum", &Radiation::GetProtonSynSpectrum)
        .def("SetDistance", &Radiation::SetDistance)
        .def("SetB_Field", &Radiation::SetB_Field)
        .def("SetSynApprox", &Radiation::SetSynApprox)
        .def("GetRadiationData", &Radiation::GetRadiationData, nb::rv_policy::reference);

        auto constants_module = nb::module_::import_("CRadpy_constants");
        m.attr("constants") = constants_module;

    nb::class_<RadData>(m, "RadData")
        .def_ro("electron_distribution", &RadData::electron_distribution)
        .def_ro("photon_distribution", &RadData::target_photons);
    
}


