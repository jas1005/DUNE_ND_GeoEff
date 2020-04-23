#include "geoEff.h"

#include <pybind11/pybind11.h>
namespace py = pybind11;

PYBIND11_MODULE(pyGeoEff, m) {
    m.doc() = "DUNE ND Geometric efficiency python module";
    
    py::class_<geoEff>(m, "geoEff")
      .def(py::init<int, bool>(), py::arg("seed"), py::arg("verbose")=false);
}
