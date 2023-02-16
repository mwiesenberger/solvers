#pragma once
// diag.h
#pragma once
#include "toefl.h"
#include "parameters.h"

namespace toefl
{

struct Variables
{
    Explicit<dg::x::CartesianGrid2d, dg::x::DMatrix, dg::x::DVec>& rhs;
    const dg::x::CartesianGrid2d& grid;
    const Parameters& p;
};

struct Record
{
    std::string name; // variable name in the output file
    std::string long_name; // longer description as an attribute
    std::function<void(dg::x::DVec&, Variables&)> function;
    // function that generates the data points for the variable
};

// time - independent output (only called once)
std::vector<Record> diagnostics2d_static_list = {
    { "xc", "x-coordinate in Cartesian coordinate system",
        []( dg::x::DVec& result, Variables& v ) {
            result = dg::evaluate( dg::cooX2d, v.grid);
        }
    },
    { "yc", "y-coordinate in Cartesian coordinate system",
        []( dg::x::DVec& result, Variables& v ) {
            result = dg::evaluate( dg::cooY2d, v.grid);
        }
    },
    { "weights", "Gaussian Integration weights",
        []( dg::x::DVec& result, Variables& v ) {
            result = dg::create::weights( v.grid);
        }
    }
};

// time - dependent output (called periodically)
std::vector<Record> diagnostics2d_list = {
    {"ne", "Electron density in 2d",
        []( dg::x::DVec& result, Variables& v) {
            dg::blas1::copy( v.rhs.var(0), result);
        }
    },
    {"gy", "Ion-gyro-center density in 2d",
        []( dg::x::DVec& result, Variables& v) {
            dg::blas1::copy( v.rhs.var(1), result);
        }
    },
    {"phi", "Electric potential",
        []( dg::x::DVec& result, Variables& v) {
            dg::blas1::copy(v.rhs.phi(0), result);
        }
    },
    {"psi", "Gyro-center potential",
        []( dg::x::DVec& result, Variables& v) {
            dg::blas1::copy(v.rhs.phi(1), result);
        }
    },
    {"lapNe", "+Delta ne",
        []( dg::x::DVec& result, Variables& v) {
            dg::blas2::gemv( -1., v.rhs.laplacianM(), v.rhs.var(0), 0., result);
        }
    },
    {"lapNi", "+Delta ni",
        []( dg::x::DVec& result, Variables& v) {
            dg::blas2::gemv( -1., v.rhs.laplacianM(), v.rhs.var(1), 0., result);
        }
    },
    {"lapPhi", "+Delta Phi",
        []( dg::x::DVec& result, Variables& v) {
            dg::blas2::gemv( -1., v.rhs.laplacianM(), v.rhs.phi(0), 0., result);
        }
    },
    {"Se", " ne ln ne",
        []( dg::x::DVec& result, Variables& v) {
            dg::blas1::transform( v.rhs.var(0), result, dg::LN<double>());
            dg::blas1::pointwiseDot(  1., v.rhs.var(0), result, 0., result);
        }
    },
    {"Si", " tau_i ni ln ni",
        []( dg::x::DVec& result, Variables& v) {
            dg::blas1::transform( v.rhs.var(1), result, dg::LN<double>());
            dg::blas1::pointwiseDot(  v.p.tau, v.rhs.var(1), result, 0., result);
        }
    },
    {"U", " 0.5 ni u_E^2",
        []( dg::x::DVec& result, Variables& v) {
            dg::blas1::pointwiseDot( 0.5, v.rhs.var(1),
                     v.rhs.uE2(), 0., result);
        }
    }
};


} //namespace toefl

