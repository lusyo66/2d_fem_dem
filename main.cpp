#include <iostream>
#include <fstream>
#include "routines.h"
#include "const.h"
#include "realtypes.h"
#include <Eigen/Dense>

using namespace dem;

int main()
{

  std::cout<<"Start 2d-FEM"<<std::endl;
  // Soil Material Parameters
  REAL rho = 198750;
//  REAL shear_vel = 250;
//  REAL mu = rho * shear_vel * shear_vel;
  REAL nu = 0.25;
  REAL Emod = 4.5e+10;
  REAL lambda = Emod * nu/ (1 + nu)/ (1 - 2*nu);
  REAL mu = Emod / 2/ (1 + nu);
  Eigen::Matrix3d D;
  D << lambda + 2* mu, lambda, 0,
       lambda, lambda + 2* mu, 0,
       0,      0,             mu;
  std::cout<<"D = "<<std::endl<<D<<std::endl;

  // Bedrock Material Parameters
  REAL shear_vel_rock = 760;
  REAL rho_rock = 2400;

  // Soil Rayleigh damping coefficients
//  REAL zeta = 0.0;
//  REAL omega1 = 2* Pi *0.2;
//  REAL omega2 = 2* Pi *20;
  REAL a_R = 0.0;
  REAL b_R = 0.0;

  // Wavelength Parameters
//  REAL fMax = 100;
//  REAL wavelength = shear_vel / fMax;
//  int nel_wave = 10;
  REAL el_height = 0.1;

  // gravitation
  REAL grav = 0.0;
  REAL grav_t_ramp = 0.0;

  // Mesh parameters
  REAL soil_height = 0.1;
  int nel = round(soil_height / el_height);
  int neldof = 8;
  int ncoords = 8;
  int numips = 4;
  int nstress = 3;
  Eigen::MatrixXd coords(nel, ncoords);
  std::vector<Eigen::MatrixXd> Kk(nel);
  std::vector<Eigen::MatrixXd> Mm(nel);
  Eigen::MatrixXd Fff(neldof, nel);
  REAL el_width = 0.1;

  std::cout<<"Parameters done"<<std::endl;
  std::cout<<"nel = "<<nel<<std::endl;

  // Coordinates, stiffness, mass matrices and body force vectors
  for (int el = 0; el < nel; el++){
//    std::cout<<"element "<<el<<std::endl;
    coords.row(el) << 0.0, el_width, el_width, 0.0, el*el_height, el*el_height, (el+1)*el_height, (el+1)*el_height;
//    std::cout<<"coords = "<<std::endl<<coords.row(el)<<std::endl;
    Kk[el] = el_kd_g2int(coords.row(el), D, neldof, numips);
//    std::cout<<"Kk ="<<std::endl<<Kk[el]<<std::endl;
    Mm[el] = el_md_g1int(coords.row(el), rho, neldof);
//    std::cout<<"Mm ="<<std::endl<<Mm[el]<<std::endl;
    Fff.col(el) = el_f_g4int(coords.row(el), rho, grav, neldof);
//    std::cout<<"Fff ="<<std::endl<<Fff.col(el)<<std::endl<<std::endl;
  }
  
  std::cout<<"Coordinates, stiffness, mass, body force done"<<std::endl;
  
  // Dashpot coefficient
  REAL c_dashpot = 0;

  int ndof = 2 * nel +1;
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(ndof, ndof);
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(ndof, ndof);
  Eigen::VectorXd Ff = Eigen::VectorXd::Zero(ndof);
//  Eigen::MatrixXd g(neldof,nel);

  Eigen::MatrixXi LM(neldof, nel);
  Eigen::RowVectorXi LM_row(neldof);
  LM_row << 1, 0, 1, 0, 2, 3, 2, 3;
  LM.col(0) = LM_row.transpose();

  for (int el = 1; el < nel; el++) {
    LM_row << 2*el, 2*el+1, 2*el, 2*el+1, 2*el+2, 2*el+3, 2*el+2, 2*el+3;
    LM.col(el) = LM_row.transpose();
  }

  std::cout<<"LM = "<<std::endl<<LM<<std::endl;

  Eigen::VectorXi LM_col;
  int I, J;
  for (int el = 0; el < nel; el++) {
    for (int i = 0; i < neldof; i++) {
      I = LM(i,el);
      if (I > 0) {
        Ff(I-1) = Ff(I-1) + Fff(i,el);
        for (int j = 0; j < neldof; j++) {
          J = LM(j,el);
          if (J > 0) {
            K(I-1, J-1) = K(I-1, J-1) + Kk[el](i,j);
            M(I-1, J-1) = M(I-1, J-1) + Mm[el](i,j);
          }
        }
      }
    }
  }
  std::cout<<"Ff = "<<Ff<<std::endl;
  std::cout<<"K ="<<std::endl<<K<<std::endl;
  std::cout<<"M ="<<std::endl<<M<<std::endl;
  Eigen::MatrixXd C = a_R * M + b_R * K;
  C(0,0) = C(0,0) + c_dashpot;

  // time integration parameters
  REAL dt = 5e-3;
  REAL gamma = 0.5;
  REAL beta = 0.25;
  Eigen::MatrixXd G = M + gamma * dt * C + beta * dt * dt * K;

  Eigen::VectorXd v = Eigen::VectorXd::Zero(ndof);
  Eigen::VectorXd d = Eigen::VectorXd::Zero(ndof);
  Eigen::VectorXd a = Eigen::VectorXd::Zero(ndof);

  REAL t = 0;
  REAL totaltime = grav_t_ramp + 1.0;
  int numsteps = round(totaltime / dt);
//  int numsteps = 1;

  Eigen::VectorXd asolve = Eigen::VectorXd::Zero(numsteps+1);
  Eigen::VectorXd dsolve = Eigen::VectorXd::Zero(numsteps+1);
  Eigen::VectorXd tsolve = Eigen::VectorXd::Zero(numsteps+1);

  // base force
  REAL cFactor = el_width * rho_rock * shear_vel_rock;
  REAL omega = 2 * Pi * 0.5;
  REAL vel0 = 1.0;

  std::ofstream d_file, strain_file, stress_file, stress_dem_file;
  d_file.open("d.txt");
  strain_file.open("strain.txt");
  stress_file.open("stress.txt");
  stress_dem_file.open("stress_dem.txt");
  d_file<<d.transpose()<<std::endl;
  strain_file<<Eigen::Matrix3d::Zero()<<std::endl;
  stress_file<<Eigen::Matrix3d::Zero()<<std::endl;

  Eigen::VectorXd dtilde, vtilde;
  REAL grav_factor;
  Eigen::MatrixXd FF = Eigen::MatrixXd::Zero(ndof, numsteps + 1);
  Eigen::VectorXd Ftotal, RHS;
  Eigen::MatrixXd d_el = Eigen::MatrixXd::Zero(nel, neldof);
  std::vector<Eigen::MatrixXd> stress_el(nel);
  std::vector<Eigen::MatrixXd> strain_el(nel);
  std::vector<Eigen::MatrixXd> stress_el_dem(nel);
  std::vector<Eigen::MatrixXd> strain_el_dem(nel);
  Eigen::MatrixXd zero_strain = Eigen::MatrixXd::Zero(numips, nstress);
  std::vector<Eigen::MatrixXd> strain_el_last(nel, zero_strain);
  Eigen::MatrixXd stress(numips, nstress);
  
  std::cout<<"Start solving"<<std::endl;
  //Solve
  for (int n = 0; n < numsteps; n++) {
    std::cout<<"t = "<<t<<std::endl;
    t = t + dt;
    dtilde = d + dt * v + (dt * dt / 2) * (1 - 2 * beta) * a;
    vtilde = v + dt * (1 - gamma) *a;
    if (t < grav_t_ramp)
      grav_factor = t / grav_t_ramp;
    else
      grav_factor = 1;
    
    FF(0,n+1) = cFactor * vel0 * sin(omega * t);
    Ftotal = grav_factor * Ff + FF.col(n+1);
    RHS = Ftotal - C * vtilde - K * dtilde;
//    std::cout<<"RHS done"<<std::endl;
    a = G.ldlt().solve(RHS);
//    std::cout<<"a done"<<std::endl;
    d = dtilde + beta * dt * dt * a;
    v = vtilde + gamma * dt * a;

    for (int el = 0; el < nel; el++) {
      for (int i = 0; i < neldof; i++) {
        I = LM(i,el);
        if (I > 0)
          d_el(el,i) = d(I-1);
      }
      el_stress_isv(coords.row(el), d_el.row(el), numips, nstress, D, strain_el[el], stress_el[el]);
      el_stress_dem(coords.row(el), d_el.row(el), numips, nstress, D, strain_el_dem[el], stress_el_dem[el], strain_el_last[el],n);
    }
//    std::cout<<"stress done"<<std::endl;
    d_file<<d.transpose()<<std::endl;
    strain_file<<strain_el[0]<<std::endl;
    stress_file<<stress_el[0]<<std::endl;
    stress_dem_file << stress_el_dem[0] <<std::endl;
    strain_el_last = strain_el;
  }  

  d_file.close();
  strain_file.close();
  stress_file.close();
  stress_dem_file.close();

  return 0;            
}
