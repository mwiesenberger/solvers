#include <iostream>

#include "dg/algorithm.h"
#include "dg/file/file.h"
#include "dg/matrix/matrix.h"


int main( int argc, char* argv[])
{
    dg::file::WrappedJsonValue js( dg::file::error::is_throw);
    try{
        std::string inputfile = "input/default.json";
        if( argc != 1) inputfile = argv[1];
        dg::file::file2Json( inputfile.c_str(), js.asJson(),
                dg::file::comments::are_discarded, dg::file::error::is_throw);
    } catch( std::exception& e) {
        std::cerr << "ERROR in input file "<<argv[1]<<std::endl;
        std::cerr << e.what()<<std::endl;
        dg::abort_program();
    }
    std::cout << js.asJson() << std::endl;

    //Construct grid
    unsigned n  = js["grid"]["n"].asUInt();
    unsigned Nx = js["grid"]["Nx"].asUInt();
    unsigned Ny = js["grid"]["Ny"].asUInt();
    dg::bc bcx = dg::str2bc(js["bc"][0].asString());
    dg::bc bcy = dg::str2bc(js["bc"][1].asString());
    dg::CartesianGrid2d grid( -1., 1., -1., 1., n, Nx, Ny, bcx, bcy);
    dg::direction elliptic_dir =  dg::str2direction( js["elliptic"]["direction"].asString());

    // assemble diffusion tensor
    dg::Elliptic<dg::CartesianGrid2d, dg::DMatrix, dg::DVec> elliptic( grid, bcx, bcy,
        elliptic_dir, 1);
    dg::SparseTensor<dg::DVec> bb;
    bb.idx(0,0) = 0, bb.idx(0,1) = 1;
    bb.idx(1,0) = 1, bb.idx(1,1) = 2;
    std::vector<dg::DVec> values( 3);
    values[0] = dg::evaluate( [](double x, double y){ return y*y;}, grid);
    values[1] = dg::evaluate( [](double x, double y){ return -x*y;}, grid);
    values[2] = dg::evaluate( [](double x, double y){ return x*x;}, grid);
    dg::blas1::scal( values, -1);//scale with -1 to counteract elliptic sign
    bb.values() = values;
    elliptic.set_chi( bb);
    // now make a matrix functional out of it
    double eps_rel = js["matrix-function"]["eps-rel"].asDouble();
    dg::mat::MatrixFunction<dg::DVec> rhs( elliptic, elliptic.weights(), eps_rel, 1, 500);

    // init vector

    dg::DVec y0 = dg::evaluate( [](double x, double y) {
        return 0.+10.*exp( -((x+0.6)*(x+0.6) + y*y)/0.04);}, grid);

    dg::mat::ExponentialStep<dg::DVec> rk(y0);
    double Tend = js["output"].get("tend", 1.0).asDouble();
    unsigned maxout = js["output"].get("maxout", 10).asUInt();
    double deltaT = Tend/(double)maxout;
    dg::SinglestepTimeloop<dg::DVec> timeloop ( rk, rhs, deltaT);

    dg::Timer t;
    t.tic();
    // open netcdf file
    std::string outputfile;
    if( argc == 1 || argc == 2)
        outputfile = "diffusion.nc";
    else
        outputfile = argv[2];
    // Create netcdf file
    dg::file::NC_Error_Handle err;
    int ncid=-1;
    try{
        err = nc_create(outputfile.c_str(), NC_NETCDF4|NC_CLOBBER, &ncid);
    }catch( std::exception& e)
    {
        std::cerr << "ERROR creating file "<<argv[1]<<std::endl;
        std::cerr << e.what() << std::endl;
        dg::abort_program();
    }
    std::map<std::string, std::string> att;
    att["title"] = "Output file of solvers/diffusion/diffusion.cpp";
    att["Conventions"] = "CF-1.8";
    ///Get local time and begin file history
    auto ttt = std::time(nullptr);

    std::ostringstream oss;
    ///time string  + program-name + args
    oss << std::put_time(std::localtime(&ttt), "%F %T %Z");
    for( int i=0; i<argc; i++) oss << " "<<argv[i];
    att["history"] = oss.str();
    att["comment"] = "Find more info in solvers/diffusion.ipynb";
    att["source"] = "FELTOR";
    att["git-hash"] = GIT_HASH;
    att["git-branch"] = GIT_BRANCH;
    att["compile-time"] = COMPILE_TIME;
    att["references"] = "https://github.com/feltor-dev/feltor";
    // Here we put the inputfile as a string without comments so that it can be read later by another parser
    att["inputfile"] = js.asJson().toStyledString();
    for( auto pair : att)
        DG_RANK0 err = nc_put_att_text( ncid, NC_GLOBAL,
            pair.first.data(), pair.second.size(), pair.second.data());

    int dim_ids[3], tvarID;
    int id3d;
    // the dimensions are the ones of grid_out!
    err = dg::file::define_dimensions( ncid, dim_ids, &tvarID, grid,
                    {"time", "y", "x"});
    err = nc_def_var( ncid, "u", NC_DOUBLE, 3, dim_ids, &id3d);
    size_t start = {0};
    size_t count = {1};
    dg::HVec transferH(y0);
    dg::file::put_vara_double( ncid, id3d, start,
                grid, transferH);
    double time = 0;
    err = nc_put_vara_double( ncid, tvarID, &start, &count, &time);
    err = nc_close( ncid);

    bool abort = false;
    for( unsigned u=1; u<=maxout; u++)
    {
        dg::Timer ti;
        ti.tic();
        try{
            timeloop.integrate( time, y0, u*deltaT, y0,
                              u < maxout ? dg::to::at_least : dg::to::exact);
        }catch ( std::exception& fail)
        {
            std::cerr << "ERROR in Timestepper\n";
            std::cerr << fail.what() << std::endl;
            std::cerr << "Writing last output and exit ..."<<std::endl;
            abort = true;
        }
        ti.toc();
        std::cout << "\n\t Time "<<time <<" of "<<Tend <<" with current timestep "<<timeloop.get_dt();
        std::cout << "\n\t Average time for one step: "<<ti.diff()<<"s\n\n"<<std::flush;
        start = u;
        err = nc_open(outputfile.c_str(), NC_WRITE, &ncid);
        // First write the time variable
        err = nc_put_vara_double( ncid, tvarID, &start, &count, &time);
        dg::assign( y0, transferH);
        dg::file::put_vara_double( ncid, id3d, start, grid, transferH);
        err = nc_close( ncid);
        if( abort) break;
    }
    ////////////////////////////////////////////////////////////////////
    t.toc();
    unsigned hour = (unsigned)floor(t.diff()/3600);
    unsigned minute = (unsigned)floor( (t.diff() - hour*3600)/60);
    double second = t.diff() - hour*3600 - minute*60;
    DG_RANK0 std::cout << std::fixed << std::setprecision(2) <<std::setfill('0');
    DG_RANK0 std::cout <<"Computation Time \t"<<hour<<":"<<std::setw(2)<<minute<<":"<<second<<"\n";
    DG_RANK0 std::cout <<"which is         \t"<<t.diff()<<"s\n";

#ifdef WITH_MPI
    MPI_Finalize();
#endif //WITH_MPI

    return 0;
}



