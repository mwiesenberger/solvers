#pragma once
#include <exception>

#include "dg/algorithm.h"
#include "parameters.h"

namespace toefl
{

template< class Geometry,  class Matrix, class Container >
struct Explicit
{
    Explicit( const Geometry& g, const Parameters& p );

    const Container& phi(unsigned i ) const { return m_phi[i];}
    const Container& var(unsigned i ) const { return m_ype[i];}
    const Container& uE2() const { return m_uE2;}

    dg::Elliptic<Geometry, Matrix, Container>& laplacianM( ) {
        return m_laplaceM;
    }

    dg::Helmholtz<Geometry, Matrix, Container >&  gamma_inv() {
        return m_multi_gamma1[0];
    }
    unsigned ncalls() const{ return m_ncalls;}

    void operator()( double t, const std::array<Container,2>& y,
            std::array<Container,2>& yp);

    void compute_psi( double t);
    void polarisation( double t, const std::array<Container,2>& y);
  private:
    //use chi and omega as helpers to compute square velocity in uE2
    Container m_chi, m_omega, m_uE2;
    Container m_temp0, m_temp1;

    const Container m_binv; //magnetic field

    std::array<Container,2> m_phi, m_dxphi, m_dyphi, m_ype;
    std::array<Container,2> m_lapy, m_v;
    Container m_gamma_n;
    Container m_fine_n, m_fine_phi, m_fine_dxphi, m_fine_dyphi, m_fine_v, m_fine_yp;
    std::array<Container,3> m_fine_temp;
    Container m_fine_binv;

    //matrices and solvers
    dg::Elliptic<Geometry, Matrix, Container> m_laplaceM;
    std::vector<dg::Elliptic<Geometry, Matrix, Container> > m_multi_pol;
    std::vector<dg::Helmholtz<Geometry,  Matrix, Container> > m_multi_gamma1;
    std::array<Matrix,2> m_dphi;
    dg::MultiMatrix<Matrix,Container> m_inter, m_project;
    std::array<Matrix,2> m_fine_forward, m_fine_backward, m_fine_dphi;
    dg::Advection< Geometry, Matrix, Container> m_adv;
    dg::ArakawaX< Geometry, Matrix, Container> m_arakawa;

    dg::MultigridCG2d<Geometry, Matrix, Container> m_multigrid;
    dg::Extrapolation<Container> m_old_phi, m_old_psi, m_old_gammaN;
    std::vector<Container> m_multi_chi;

    Parameters m_p;

    unsigned m_ncalls = 0;

};

template< class Geometry, class M, class Container>
Explicit< Geometry, M, Container>::Explicit( const Geometry& grid, const Parameters& p ):
    m_chi( evaluate( dg::zero, grid)), m_omega(m_chi), m_uE2(m_chi),
    m_temp0(m_chi), m_temp1(m_chi),
    m_binv( evaluate( dg::LinearX( p.kappa, 1.-p.kappa*p.posX*p.lx), grid)),
    m_phi( {m_chi, m_chi}), m_dxphi(m_phi), m_dyphi( m_phi), m_ype(m_phi),
    m_lapy(m_phi), m_v(m_phi),
    m_gamma_n(m_chi),
    m_laplaceM( grid,  p.diff_dir),
    m_adv( grid), m_arakawa(grid),
    m_multigrid( grid, p.num_stages),
    m_old_phi( 2, m_chi), m_old_psi( 2, m_chi), m_old_gammaN( 2, m_chi),
    m_p(p)
{
    m_multi_chi= m_multigrid.project( m_chi);
    for( unsigned u=0; u<p.num_stages; u++)
    {
        m_multi_pol.push_back({ m_multigrid.grid(u),  p.pol_dir, 1.});
        m_multi_gamma1.push_back({-0.5*p.tau, { m_multigrid.grid(u), p.pol_dir}});
    }
    m_dphi = {dg::create::dx( grid, m_p.bcx, m_p.dphi_dir),
                  dg::create::dy( grid, m_p.bcy, m_p.dphi_dir )};
    if( "div-upwind" == m_p.advection)
    {
        Geometry fine_grid = grid;
        fine_grid.set( 2*grid.n(), grid.Nx(), grid.Ny());
        //theoretically we only need 2n-1 but it isn't wrong to take more
        m_inter = dg::create::fast_interpolation( grid, 2, 1, 1);
        m_project = dg::create::fast_projection( fine_grid, 2, 1, 1);

        m_fine_dphi[0] = dg::create::dx( fine_grid, m_p.bcx, m_p.dphi_dir);
        m_fine_dphi[1] = dg::create::dy( fine_grid, m_p.bcy, m_p.dphi_dir);
        m_fine_forward[0] = dg::create::dx( fine_grid, dg::inverse( m_p.bcx), dg::forward);
        m_fine_forward[1] = dg::create::dy( fine_grid, dg::inverse( m_p.bcy), dg::forward);
        m_fine_backward[0] = dg::create::dx( fine_grid, dg::inverse( m_p.bcx), dg::backward);
        m_fine_backward[1] = dg::create::dy( fine_grid, dg::inverse( m_p.bcy), dg::backward);
        m_fine_phi = dg::evaluate( dg::zero, fine_grid);
        m_fine_dxphi = m_fine_dyphi = dg::evaluate( dg::zero, fine_grid);
        m_fine_n = m_fine_yp = m_fine_v = dg::evaluate( dg::zero, fine_grid);
        m_fine_temp[0] = dg::evaluate( dg::zero, fine_grid);
        m_fine_temp[1] = dg::evaluate( dg::zero, fine_grid);
        m_fine_temp[2] = dg::evaluate( dg::zero, fine_grid);
        m_fine_binv = dg::evaluate( dg::LinearX( p.kappa, 1.-p.kappa*p.posX*p.lx), fine_grid);
    }
}

template< class G, class M, class Container>
void Explicit<G, M, Container>::compute_psi( double t)
{
    if (m_p.tau == 0.) {
        dg::blas1::axpby( 1.,m_phi[0], 0.,m_phi[1]); //chi = N_i - 1
    }
    else {
        m_old_psi.extrapolate( t, m_phi[1]);
        m_multigrid.set_benchmark( true, "Gamma Phi   ");
        m_multigrid.solve( m_multi_gamma1, m_phi[1], m_phi[0], m_p.eps_gamma);
        m_old_psi.update( t, m_phi[1]);
    }
    //compute (nabla phi)^2
    m_multi_pol[0].variation(m_phi[0], m_uE2);
    //compute psi
    dg::blas1::pointwiseDot( 1., m_binv, m_binv, m_uE2, 0., m_uE2);
    dg::blas1::axpby( -0.5, m_uE2, 1., m_phi[1]);
}


//computes and modifies expy!!
template<class G, class M, class Container>
void Explicit<G, M, Container>::polarisation( double t,
        const std::array<Container,2>& y)
{
    //compute chi
    dg::blas1::evaluate( m_chi, dg::equals(), []DG_DEVICE
            ( double nt, double binv){
                return (nt+1.)*binv*binv;
            }, y[1], m_binv);
    m_multigrid.project( m_chi, m_multi_chi);
    for( unsigned u=0; u<3; u++)
        m_multi_pol[u].set_chi( m_multi_chi[u]);
    //compute polarisation
    if (m_p.tau == 0.) {
        dg::blas1::axpby( 1., y[1], 0., m_gamma_n); //chi = N_i - 1
    }
    else {
        m_old_gammaN.extrapolate(t, m_gamma_n);
        m_multigrid.set_benchmark( true, "Gamma N     ");
        m_multigrid.solve( m_multi_gamma1, m_gamma_n, y[1], m_p.eps_gamma);
        m_old_gammaN.update(t, m_gamma_n);
    }
    dg::blas1::axpby( -1., y[0], 1., m_gamma_n, m_omega); //omega = a_i\Gamma n_i - n_e
    //invert

    m_old_phi.extrapolate(t, m_phi[0]);
    m_multigrid.set_benchmark( true, "Polarisation");
    m_multigrid.solve( m_multi_pol, m_phi[0], m_omega, m_p.eps_pol);
    m_old_phi.update( t, m_phi[0]);
}

template< class G, class M, class Container>
void Explicit<G, M, Container>::operator()( double t,
        const std::array<Container,2>& y, std::array<Container,2>& yp)
{
    m_ncalls ++ ;
    //y[0] = N_e - 1
    //y[1] = N_i - 1 || y[1] = Omega

    polarisation( t, y);
    compute_psi( t);

    ///////////////////////////////////////////////////////////////////////
    dg::blas1::transform( y, m_ype, dg::PLUS<double>(1.));
    std::array<double, 2> tau = {-1., m_p.tau};
    for( unsigned u=0; u<2; u++)
    {
        // ExB + Curv advection with updwind scheme
        if ( "upwind" == m_p.advection)
        {
            dg::blas2::symv( m_dphi[0], m_phi[u], m_dxphi[u]);
            dg::blas2::symv( m_dphi[1], m_phi[u], m_dyphi[u]);
            dg::blas1::pointwiseDot( -1., m_binv, m_dyphi[u], 0., m_v[0]);
            dg::blas1::pointwiseDot( +1., m_binv, m_dxphi[u], 0., m_v[1]);
            dg::blas1::plus( m_v[1], -tau[u]*m_p.kappa);
            m_adv.upwind( -1., m_v[0], m_v[1], y[u], 0., yp[u]);
            // Div ExB velocity
            dg::blas1::pointwiseDot( m_p.kappa, m_ype[u], m_dyphi[u], 1., yp[u]);
        }
        else if ( "div-upwind" == m_p.advection)
        {
            dg::blas2::symv( m_inter, m_ype[u], m_fine_n);
            dg::blas2::symv( m_inter, m_phi[u], m_fine_phi);
            dg::blas1::copy( 0., m_fine_yp);
            dg::blas2::symv( m_fine_dphi[0], m_fine_phi, m_fine_dxphi);
            dg::blas2::symv( m_fine_dphi[1], m_fine_phi, m_fine_dyphi);
            //dx ( nv_x)
            dg::blas1::pointwiseDot( -1., m_fine_binv, m_fine_dyphi, 0., m_fine_v);
            //v_x
            dg::blas1::pointwiseDot( m_fine_n, m_fine_v, m_fine_temp[0]); //f_x
            dg::blas2::symv( m_fine_forward[0], m_fine_temp[0], m_fine_temp[1]);
            dg::blas2::symv( m_fine_backward[0], m_fine_temp[0], m_fine_temp[2]);
            dg::blas1::evaluate( m_fine_yp, dg::minus_equals(), dg::Upwind(), m_fine_v, m_fine_temp[2], m_fine_temp[1]);
            //dy ( nv_y)
            dg::blas1::pointwiseDot( +1., m_fine_binv, m_fine_dxphi, 0., m_fine_v);
            dg::blas1::plus( m_fine_v, -tau[u]*m_p.kappa);
            //v_y
            dg::blas1::pointwiseDot( m_fine_n, m_fine_v, m_fine_temp[0]); //f_y
            dg::blas2::symv( m_fine_forward[1], m_fine_temp[0], m_fine_temp[1]);
            dg::blas2::symv( m_fine_backward[1], m_fine_temp[0], m_fine_temp[2]);
            dg::blas1::evaluate( m_fine_yp, dg::minus_equals(), dg::Upwind(), m_fine_v, m_fine_temp[2], m_fine_temp[1]);
            dg::blas2::symv( m_project, m_fine_yp, yp[u]);
        }
    }

    for( unsigned u=0; u<2; u++)
    {
        if( m_p.nu > 0)
        {
            dg::blas1::copy( y[u], m_temp0);
            for( unsigned s=0; s<m_p.diff_order; s++)
            {
                using std::swap;
                swap( m_temp0, m_temp1);
                dg::blas2::symv( 1., m_laplaceM, m_temp1, 0., m_temp0);
            }
            dg::blas1::axpby( -m_p.nu, m_temp0, 1, yp[u]);
        }
    }
}

}//namespace dg
