/*
###############################################################################
#
#   Test: interfacing from Python to C++
#   AUTHOR: Maciej_Bak
#   CONTACT: wsciekly.maciek@gmail.com
#
###############################################################################
*/

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/pytypes.h>

#include <Eigen/LU>
#include <mlpack.hpp>
#include <armadillo>
#include <boost/math/differentiation/finite_difference.hpp>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>

// simple function 1 (scalars)
double add2numbers(double i, double j) {
    return i + j;
}

// simple function 2 (vectorized)
double square_a_number(double x) {
    return x * x;
}

// Boost example:
// https://www.boost.org/
// https://www.boost.org/doc/libs/master/libs/math/doc/html/math_toolkit/diff.html
auto f = [](double x) { return 2 * std::exp(x); };
double boost_deriv(double x){
    return boost::math::differentiation::finite_difference_derivative(f, x);
}

// Eigen matrix from/to numpy array
// https://eigen.tuxfamily.org/dox/
// https://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html
Eigen::MatrixXd inverse(Eigen::MatrixXd m) {
    return m.inverse();
}

// Armadillo example:
// https://arma.sourceforge.net/
Eigen::MatrixXd lm(Eigen::MatrixXd XX, Eigen::MatrixXd yy) {
    arma::mat X(XX.data(), XX.rows(), XX.cols());
    arma::vec y(yy.data(), yy.size(), false);
    arma::colvec coeffs = arma::solve(X, y);
    arma::mat coeffmat = arma::mat(coeffs.memptr(), X.n_cols, 1);
    Eigen::MatrixXd eigen_coeffmat(coeffmat.n_rows, coeffmat.n_cols);
    for (int i = 0; i < coeffmat.n_rows; ++i) {
        for (int j = 0; j < coeffmat.n_cols; ++j) {
            eigen_coeffmat(i, j) = coeffmat(i, j);
        }
    }
    return eigen_coeffmat;
}

// GSL example:
// https://www.gnu.org/software/gsl/
// https://www.gnu.org/software/gsl/doc/html/usage.html
// https://www.gnu.org/software/gsl/doc/html/diff.html
double g (double x, void* params) {
  (void)(params); /* avoid unused parameter warning */
  return 2 * std::exp(x);
}
pybind11::tuple gsl_deriv() {
  gsl_function F;
  double result, abserr;
  F.function = &g;
  F.params = 0;
  gsl_deriv_central (&F, 0.0, 1e-8, &result, &abserr);
  return pybind11::make_tuple(result, abserr);
}

// mlpack example:
// https://www.mlpack.org/
double mlpack_example() {
  arma::mat data("0.339406815,0.843176636,0.472701471; \
                  0.212587646,0.351174901,0.81056695;  \
                  0.160147626,0.255047893,0.04072469;  \
                  0.564535197,0.943435462,0.597070812"); 
  data = data.t(); 
  mlpack::NeighborSearch<mlpack::NearestNeighborSort, mlpack::ManhattanDistance> nn(data);
  arma::Mat<size_t> neighbors;
  arma::mat distances;
  nn.Search(1, neighbors, distances);
  double avgdist = arma::accu(distances.row(0)) / 4;
  return avgdist;
}

// binding C++ code to Python
PYBIND11_MODULE(functions, m) {
    m.doc() = "Module: example C++ functions";
    m.def("add2numbers", &add2numbers, pybind11::arg("x"), pybind11::arg("y"));
    m.def("square_a_number", pybind11::vectorize(square_a_number));
    m.def("boost_deriv", &boost_deriv);
    m.def("inverse", &inverse);
    m.def("lm", &lm);
    m.def("gsl_deriv", &gsl_deriv);
    m.def("mlpack_example", &mlpack_example);
}
