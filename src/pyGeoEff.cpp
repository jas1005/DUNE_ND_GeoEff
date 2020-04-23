#include "geoEff.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(pyGeoEff, m) {
    m.doc() = "DUNE ND Geometric efficiency python module";
    
    py::class_<geoEff>(m, "geoEff")
      .def(py::init<int, bool>(), py::arg("seed"), py::arg("verbose")=false)
      .def("setNthrows", &geoEff::setNthrows)
      .def("setVertex", &geoEff::setVertex)
      .def("setHitSegEdeps", &geoEff::setHitSegEdeps)
      .def("setHitSegPoss", &geoEff::setHitSegPoss)
      .def("setRangeX", &geoEff::setRangeX)
      .def("setRangeY", &geoEff::setRangeY)
      .def("setRangeZ", &geoEff::setRangeZ)
      .def("setActiveX", &geoEff::setActiveX)
      .def("setActiveY", &geoEff::setActiveY)
      .def("setActiveZ", &geoEff::setActiveZ)
      .def("setOffsetX", &geoEff::setOffsetX)
      .def("setOffsetY", &geoEff::setOffsetY)
      .def("setOffsetZ", &geoEff::setOffsetZ)
      .def("setBeamDir", &geoEff::setBeamDir)
      .def("setVetoSizes", &geoEff::setVetoSizes)
      .def("setVetoEnergyThresholds", &geoEff::setVetoEnergyThresholds)
      .def("throwTransforms", &geoEff::throwTransforms)
      .def("getCurrentThrowTranslationsX", &geoEff::getCurrentThrowTranslationsX)
      .def("getCurrentThrowTranslationsY", &geoEff::getCurrentThrowTranslationsY)
      .def("getCurrentThrowTranslationsZ", &geoEff::getCurrentThrowTranslationsZ)
      .def("getCurrentThrowRotations", &geoEff::getCurrentThrowRotations)
      .def("getHadronContainment", &geoEff::getHadronContainment)
      .def("setSeed", &geoEff::setSeed);
}
