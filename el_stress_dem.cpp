#include "routines.h"
#include <fstream>
#include <string>
#include <iostream>
#include <stdio.h>
#include <sstream>

void el_stress_dem(Eigen::RowVectorXd coordsx,
                   Eigen::RowVectorXd d,
                   int numips,
                   int nstress,
                   Eigen::MatrixXd D,
                   Eigen::MatrixXd & strain_el,
                   Eigen::MatrixXd & stress_el,
		   Eigen::MatrixXd strain_el_last,
                   int n) {
  REAL x1, x2, x3, x4, y1, y2, y3, y4;
  REAL xi, eta;
  REAL dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi;
  REAL dN1_deta, dN2_deta, dN3_deta, dN4_deta;
  REAL dx_dxi, dx_deta, dy_dxi, dy_deta;
  REAL j;
  REAL dStrainRate, s11, s12, s13, s21, s22, s23, s31, s32, s33;

  int dem_steps = 1000;

  std::string line, filename;
  std::stringstream dump1, dump2, dump3;

  Eigen::Matrix2d Je, Jeinv;
  Eigen::MatrixXd B(3,8);
  Eigen::RowVector2d dN1, dN2, dN3, dN4;
  Eigen::RowVector2d dN1_dx_vect, dN2_dx_vect, dN3_dx_vect, dN4_dx_vect;

  x1 = coordsx(0);
  x2 = coordsx(1);
  x3 = coordsx(2);
  x4 = coordsx(3);
  y1 = coordsx(4);
  y2 = coordsx(5);
  y3 = coordsx(6);
  y4 = coordsx(7);

  strain_el = Eigen::MatrixXd::Zero(numips, nstress);
  stress_el = Eigen::MatrixXd::Zero(numips, nstress);

  REAL const0 = 1/sqrt(3);
  Eigen::MatrixXd xi_vect(4,2);
  xi_vect << -const0, -const0,
              const0, -const0,
              const0,  const0,
             -const0,  const0;
  Eigen::Vector4d weight(1.0, 1.0, 1.0, 1.0);
  std::fstream copy_old, stress_output;
  std::ofstream dem_input, copy_new;


  for (int i=0; i<4; i++){

    xi = xi_vect(i,0);
    eta = xi_vect(i,1);

    dN1_dxi = -0.25 * (1 - eta);
    dN2_dxi = -dN1_dxi;
    dN3_dxi = 0.25 * (1 + eta);
    dN4_dxi = -dN3_dxi;

    dN1_deta = -0.25 * (1 - xi);
    dN2_deta = -0.25 * (1 + xi); 
    dN3_deta = -dN2_deta;
    dN4_deta = -dN1_deta;

    dx_dxi = dN1_dxi * x1 + dN2_dxi * x2 + dN3_dxi * x3 + dN4_dxi * x4;
    dx_deta = dN1_deta * x1 + dN2_deta * x2 + dN3_deta * x3 + dN4_deta * x4;
    dy_dxi = dN1_dxi * y1 + dN2_dxi * y2 + dN3_dxi * y3 + dN4_dxi * y4;
    dy_deta = dN1_deta * y1 + dN2_deta * y2 + dN3_deta * y3 + dN4_deta * y4;
    Je << dx_dxi, dx_deta,
          dy_dxi, dy_deta;
    j = Je.determinant();
    Jeinv = Je.inverse();

    dN1 << dN1_dxi, dN1_deta;
    dN1_dx_vect = dN1 * Jeinv;
    B.block(0,0,3,2) << dN1_dx_vect(0), 0,
                         0, dN1_dx_vect(1),
                         dN1_dx_vect(1), dN1_dx_vect(0);
    dN2 << dN2_dxi, dN2_deta;
    dN2_dx_vect = dN2 * Jeinv;
    B.block(0,2,3,2) << dN2_dx_vect(0), 0,
                         0, dN2_dx_vect(1),
                         dN2_dx_vect(1), dN2_dx_vect(0);

    dN3 << dN3_dxi, dN3_deta;
    dN3_dx_vect = dN3 * Jeinv;
    B.block(0,4,3,2) << dN3_dx_vect(0), 0,
                         0, dN3_dx_vect(1),
                         dN3_dx_vect(1), dN3_dx_vect(0);

    dN4 << dN4_dxi, dN4_deta;
    dN4_dx_vect = dN4 * Jeinv;
    B.block(0,6,3,2) << dN4_dx_vect(0), 0,
                         0, dN4_dx_vect(1),
                         dN4_dx_vect(1), dN4_dx_vect(0);

    strain_el.row(i) = B * d.transpose();
    dStrainRate = (strain_el(i,2) - strain_el_last(i,2)) / dem_steps;

    copy_old.open("input_org.txt");
    dem_input.open("input.txt");
    while(getline(copy_old, line)){
      dem_input << line <<std::endl;
    }
    dem_input << std::endl << "shearRate "<< dStrainRate;
    dem_input.close();
    copy_old.close();


    filename = "ini_boundary_" + std::to_string(i);
    rename(filename.c_str(), "ini_boundary");
    filename = "ini_particle_" + std::to_string(i);
    rename(filename.c_str(), "ini_particle");


//    std::cout<<"Call paraEllip3d"<<std::endl;

    std::system("./paraEllip3d input.txt");

    stress_output.open("simpleShear_tensor_001");
    while( getline (stress_output, line) ){
//      std::cout<<line<<std::endl;
      if (*line.begin() == 's'){
        getline(stress_output, line);
        dump1 << line;
	dump1 >> s11 >> s12 >> s13;
        getline(stress_output, line);
        dump2 << line;
        dump2 >> s21 >> s22 >> s23;
        getline(stress_output, line);
        dump3 << line;
        dump3>> s31 >> s32 >> s33;
        break;
      }
    }
    stress_el(i,2) = s32;
//    std::cout<<"stress = "<<std::endl<<s11<<' '<<s12<<' '<<s13<<std::endl<<s21<<' '<<s22<<' '<<s23<<std::endl<<s31<<' '<<s32<<' '<<s33<<std::endl;
    

    filename = "output/simpleShear_particle_" + std::to_string(i) + "_" + std::to_string(n);
    copy_old.open("simpleShear_particle_end");
    copy_new.open(filename);
    while(getline(copy_old, line)){
      copy_new << line <<std::endl;
    }
    copy_old.close();
    copy_new.close();

    filename = "output/simpleShear_boundary_" + std::to_string(i) + "_" + std::to_string(n);
    copy_old.open("simpleShear_boundary_end");
    copy_new.open(filename);
    while(getline(copy_old, line)){
      copy_new << line <<std::endl;
    }
    copy_old.close();
    copy_new.close();

    filename = "ini_boundary_" + std::to_string(i);
    rename("simpleShear_boundary_end", filename.c_str());

    filename = "ini_particle_" + std::to_string(i);
    rename("simpleShear_particle_end", filename.c_str());

    dem_input.close();
    stress_output.close();
  }
}
