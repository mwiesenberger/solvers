#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <functional>

#include "dg/algorithm.h"
#include "dg/geometries/geometries.h"
#include "dg/file/json_utilities.h"

//global relative error in L2 norm is O(h^P)
//as a rule of thumb with n=4 the true error is err = 1e-3 * eps as long as eps > 1e3*err

const double lx = M_PI;
const double ly = 2.*M_PI;
dg::bc bcx = dg::DIR;
dg::bc bcy = dg::PER;

double initial( double x) {return 0.;}
double initial( double x, double y) {return 0.;}
//double amp = 0.9999; // LGMRES has problems here
double amp = 0.9;
double pol( double x) {return 1. + amp*sin(x); } //must be strictly positive
double pol( double x, double y) {return 1. + amp*sin(x)*sin(y); } //must be strictly positive
//double pol( double x, double y) {return 1.; }
//double pol( double x, double y) {return 1. + sin(x)*sin(y) + x; } //must be strictly positive

double rhs( double x) { return sin(x) + amp*sin(x)*sin(x) - amp*cos(x)*cos(x);}
double rhs( double x, double y) { return 2.*sin(x)*sin(y)*(amp*sin(x)*sin(y)+1)-amp*sin(x)*sin(x)*cos(y)*cos(y)-amp*cos(x)*cos(x)*sin(y)*sin(y);}
//double rhs( double x, double y) { return 2.*sin( x)*sin(y);}
//double rhs( double x, double y) { return 2.*sin(x)*sin(y)*(sin(x)*sin(y)+1)-sin(x)*sin(x)*cos(y)*cos(y)-cos(x)*cos(x)*sin(y)*sin(y)+(x*sin(x)-cos(x))*sin(y) + x*sin(x)*sin(y);}
double sol(double x)  { return sin( x);}
double sol(double x, double y)  { return sin( x)*sin(y);}
double derX(double x, double y)  { return cos( x)*sin(y);}
double derY(double x, double y)  { return sin( x)*cos(y);}
double vari(double x, double y)  { return pol(x,y)*pol(x,y)*(derX(x,y)*derX(x,y) + derY(x,y)*derY(x,y));}


std::function<void(const dg::DVec&, dg::DVec&)> create_solver(
    dg::file::WrappedJsonValue js, // the solver
    dg::CartesianGrid2d grid,
    unsigned& number,
    dg::Elliptic<dg::CartesianGrid2d,dg::DMatrix,dg::DVec>& elliptic
    )
{
    dg::DVec w2d = dg::create::weights( grid);
    unsigned max_iter = js["max_iter"].asUInt();
    std::string preconditioner = js["preconditioner"]["type"].asString();
    std::string solver = js["type"].asString();
    double eps = js["eps"].asDouble();
    dg::DVec precond = elliptic.precond();
    if ( preconditioner == "none")
        precond = dg::evaluate( dg::one, grid);
    else if( preconditioner == "diagonal")
        ;
    else
        throw dg::Error( dg::Message( _ping_) << "Precond "<<preconditioner<<" not recognized!\n");

    if( solver == "CG")
    {
        unsigned check_every = js["check-every"].asUInt();
        dg::PCG<dg::DVec > pcg( w2d, max_iter);
        return [&, pcg, precond, w2d, eps, check_every]( const auto& y, auto& x) mutable
        {
            number = pcg.solve( elliptic, x, y, precond, w2d, eps, 1., check_every);
        };
    }
    else if( solver == "LGMRES")
    {
        unsigned inner_m = js["inner_m"].asUInt(), outer_k = js["outer_k"].asUInt();
        unsigned max_iter = js["max_iter"].asUInt();
        dg::LGMRES<dg::DVec> lgmres( w2d, inner_m, outer_k, max_iter/inner_m);
        return [&, lgmres, precond, w2d, eps] ( const auto& y, auto& x) mutable
        {
            number = lgmres.solve( elliptic, x, y, precond, w2d, eps);
        };
    }
    else if( solver == "BICGSTABl")
    {
        unsigned l_input = js["l_input"].asUInt();
        unsigned max_iter = js["max_iter"].asUInt();
        dg::BICGSTABl<dg::DVec> bicg( w2d, max_iter, l_input);
        return [&, bicg, precond, w2d, eps] ( const auto& y, auto& x) mutable
        {
            number = bicg.solve( elliptic, x, y, precond, w2d, eps);
        };
    }
    else
        throw dg::Error( dg::Message( _ping_) << "Solver "<<solver<<" not recognized!\n");
}


int main( int argc, char* argv[])
{
    dg::file::WrappedJsonValue js( dg::file::error::is_throw);
    std::string input = argc==1 ? "input.json" : argv[1];
    dg::file::file2Json( input, js.asJson(), dg::file::comments::are_discarded);

    unsigned  n = js["grid"] ["n"].asUInt(), Nx = js["grid"]["Nx"].asUInt(),
             Ny = js["grid"]["Ny"].asUInt();

    amp = js["equations"]["amp"].asDouble();

    std::cout << "# Computation on: "<< n <<" x "<< Nx <<" x "<< Ny << std::endl;
    //std::cout << "# of 2d cells                 "<< Nx*Ny <<std::endl;
    std::stringstream sout;

	dg::CartesianGrid2d grid( 0, lx, 0, ly, n, Nx, Ny, bcx, bcy);
    dg::DVec w2d = dg::create::weights( grid);
    //create functions A(chi) x = b
    dg::DVec x =    dg::evaluate( initial, grid);
    dg::DVec b =    dg::evaluate( rhs, grid);
    dg::DVec chi =  dg::evaluate( pol, grid);
    dg::DVec chi_inv(chi);
    dg::blas1::transform( chi, chi_inv, dg::INVERT<double>());
    dg::DVec temp = x;
    //compute error
    const dg::DVec solution = dg::evaluate( sol, grid);
    const dg::DVec derivati = dg::evaluate( derX, grid);
    const dg::DVec variatio = dg::evaluate( vari, grid);
    const double norm = dg::blas2::dot( w2d, solution);
    dg::DVec error( solution);

    double jfactor = js["elliptic"]["jfactor"].asDouble();;
    std::string direction = js["elliptic"]["direction"].asString();
    auto dir = dg::str2direction(direction);

    dg::Elliptic<dg::CartesianGrid2d, dg::DMatrix, dg::DVec> elliptic( grid, dir, jfactor);
    elliptic.set_chi(chi);

    std::string solver = js["solver"]["type"].asString();
    dg::Timer t;
    x = temp;
    if( solver != "Multigrid")
    {
        unsigned number = 0;
        auto inverse_elliptic = create_solver( js["solver"], grid, number, elliptic);
        t.tic();
        dg::apply( inverse_elliptic, b, x);
        t.toc();
        sout << "time: "<<t.diff()<<"\n";
        sout << "iter: "<<number<<"\n";
    }

    else if( solver == "Multigrid-FAS")
    {
        unsigned num_stages = js["solver"]["num_stages"].asUInt();
        dg::NestedGrids<dg::CartesianGrid2d, dg::DMatrix, dg::DVec> nested( grid, num_stages);
        const std::vector<dg::DVec> multi_chi = nested.project( chi);
        std::vector<dg::DVec> multi_x = nested.project( x);
        std::vector<dg::DVec> multi_b = nested.project( b);
        std::vector<dg::Elliptic<dg::CartesianGrid2d, dg::DMatrix, dg::DVec> >
            multi_pol( num_stages);
        std::vector<std::function<void( const dg::DVec&, dg::DVec&)> >
            multi_inv_pol(num_stages);
        std::vector<unsigned> numbers(num_stages, 0);
        for( unsigned u=0; u<num_stages; u++)
        {
            multi_pol[u].construct( nested.grid(u), dir, jfactor);
            multi_pol[u].set_chi( multi_chi[u]);
            multi_inv_pol[u] = [=, &sout, inverse = create_solver( js["solver"]["solvers"][u], nested.grid(u),
                numbers[u], multi_pol[u])] (const auto& y, auto& x) mutable
            {
                dg::Timer t;
                t.tic();
                dg::apply( inverse, y, x);
                t.toc();
                sout << "stage: "<<u<<"\n";
                sout << "    time: "<<t.diff()<<"\n";
                sout << "    iter: "<<numbers[u]<<"\n";
            };
        }
        t.tic();
        dg::nested_iterations( multi_pol, x, b, multi_inv_pol, nested);
        t.toc();
        sout <<"time "<<t.diff()<<"\n";
    }
    //compute the error (solution contains analytic solution
    dg::blas1::axpby( 1.,x,-1., solution, error);

    //compute the L2 norm of the error
    double err = dg::blas2::dot( w2d, error);
    sout << "error: "<<sqrt(err/norm)<<"\n";
    sout << "error_abs: "<<sqrt(err)<<"\n";

    if( argc == 3)
    {
        std::fstream outf( argv[2]);
        outf << sout.str();
    }
    else
        std::cout << sout.str();

    return 0;
}

