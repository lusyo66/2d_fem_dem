#include <Eigen/Dense>
#include "realtypes.h"

Eigen::MatrixXd el_kd_g2int(Eigen::RowVectorXd coordsx, Eigen::MatrixXd D, int neldof, int numips);
Eigen::MatrixXd el_md_g1int(Eigen::RowVectorXd coordsx, REAL rho, int neldof);
Eigen::VectorXd el_f_g4int(Eigen::RowVectorXd coordsx, REAL rho, REAL grav, int neldof);
void el_stress_isv(Eigen::RowVectorXd coordsx, Eigen::RowVectorXd d, int numips, int nstress, Eigen::MatrixXd D, Eigen::MatrixXd & strain_el, Eigen::MatrixXd & stress_el );
void el_stress_dem(Eigen::RowVectorXd coordsx, Eigen::RowVectorXd d, int numips, int nstress, Eigen::MatrixXd D, Eigen::MatrixXd & strain_el, Eigen::MatrixXd & stress_el, Eigen::MatrixXd strain_el_last, int n);
