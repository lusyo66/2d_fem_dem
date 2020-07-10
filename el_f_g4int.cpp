#include "routines.h"

Eigen::VectorXd el_f_g4int( Eigen::RowVectorXd coordsx, REAL rho, REAL grav, int neldof) {

  REAL x1, x2, x3, x4, y1, y2, y3, y4;
  REAL xi, eta;
  REAL N1, N2, N3, N4;
  REAL dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi;
  REAL dN1_deta, dN2_deta, dN3_deta, dN4_deta;
  REAL dx_dxi, dx_deta, dy_dxi, dy_deta;
  REAL j;

  Eigen::Matrix2d Je;
  Eigen::MatrixXd N(2,8);
  Eigen::Vector2d g(0.0, -grav);

  x1 = coordsx(0);
  x2 = coordsx(1);
  x3 = coordsx(2);
  x4 = coordsx(3);
  y1 = coordsx(4);
  y2 = coordsx(5);
  y3 = coordsx(6);
  y4 = coordsx(7);

  Eigen::VectorXd f_e = Eigen::VectorXd::Zero(neldof);

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

    N1 = 0.25 * (1 - xi) * (1 - eta);
    N2 = 0.25 * (1 + xi) * (1 - eta);
    N3 = 0.25 * (1 + xi) * (1 + eta);
    N4 = 0.25 * (1 - xi) * (1 + eta);

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

    N << N1, 0, N2, 0, N3, 0, N4, 0,
         0, N1, 0, N2, 0, N3, 0, N4;

    f_e = f_e + rho * N.transpose() * g * j * weight(i);
  }

  return f_e;
}
