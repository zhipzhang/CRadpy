#include "nanobind/nanobind.h"
#include "nanobind/eigen/dense.h"
#include "nanobind/nb_defs.h"
#include "Radiation.h"
#include "Eigen/Core"
namespace  nb = nanobind;
using namespace nb::literals;
using CRad::Radiation;
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
        .def("SetDistance", &Radiation::SetDistance);

        auto constants_module = nb::module_::import_("CRadpy_constants");
        m.attr("constants") = constants_module;
    
}


