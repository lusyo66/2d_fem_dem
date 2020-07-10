#include "routines.h"


void el_stress_isv(Eigen::RowVectorXd coordsx,
                   Eigen::RowVectorXd d,
                   int numips,
                   int nstress,
                   Eigen::MatrixXd D,
                   Eigen::MatrixXd & strain_el,
                   Eigen::MatrixXd & stress_el) {
  REAL x1, x2, x3, x4, y1, y2, y3, y4;
  REAL xi, eta;
  REAL dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi;
  REAL dN1_deta, dN2_deta, dN3_deta, dN4_deta;
  REAL dx_dxi, dx_deta, dy_dxi, dy_deta;
  REAL j;

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
    stress_el.row(i) = D * B * d.transpose();
  }
}
