#include <iostream>
#include <string>
#include <math.h>

#include "Eigen/Dense"
using namespace Eigen;

#include "fileio.h"
using namespace std;

#include<vector>
typedef vector <Matrix2d> M3d;

Vector2d gravity(0.0, -1.0);

class FEM{

  double E, nu, mu, lam, ball_radius;
  Vector2d ball_pos;

  MatrixXd pos, f2v;
  int NF, NV;

  MatrixXd vel, force;
  M3d B, F;
  VectorXd V, phi;
  double U, node_mass;

  void init_mesh(string meshdir); //mesh initializer
  Matrix2d calcP(int t); //calc P function

public:

  FEM(); //constructor equivalent to __init__
  void update_U(); //potential energy updater(?)
  void advance(double dt); //advance the simul
  void save_state(int itr, string statedir); //write state to file

};

//constructor definition

FEM::FEM() {
  E  = 1e3; nu = 0.2; //Young's modulus & Poisson's ratio
  mu = E /  2 / (1 + nu);
  lam= E * nu / (1 + nu) / (1 - 2 * nu);

  ball_pos(0) = 0.4; ball_pos(1) = 0.0;
  ball_radius = 0.5;
  init_mesh("./circle/"); //change this to change the mesh being used

  vel = MatrixXd::Zero(NV, 2); force = MatrixXd::Zero(NV, 2);
  B.resize(NF, Matrix2d(2,2)); F.resize(NF, Matrix2d(2,2));
  //V is area/volume and phi is potential energy of each face
  V = VectorXd::Zero(NF); phi = VectorXd::Zero(NF);
  U = 0.0; //total potential energy

  for (int i=0; i<NF; i++) {
    int ia = f2v(i,0); int ib = f2v(i,1); int ic = f2v(i,2);

    Matrix2d t1; t1.row(0) = pos.row(ia); t1.row(1) = pos.row(ib);
    Matrix2d t2; t2.row(0) = pos.row(ic); t2.row(1) = pos.row(ic);

    Matrix2d B_i_inv = t1 - t2;
    B[i] = B_i_inv.inverse().transpose();
  }
}

void FEM::init_mesh(string meshdir) {

  string meshopdir = meshdir + "mesh_output/";

  Vector2i pars = get_params(meshopdir + "mesh_params.dat");
  pos = read_data(meshopdir + "mesh_node.dat", pars(0), 2);
  f2v = read_data(meshopdir + "mesh_tnode.dat", pars(1), 3);

  NV = pars(0); NF = pars(1);

  for (int i=0; i<NF; i++) {
    double temp = f2v(i,1);
    //swapping to make clockwise
    f2v(i,1) = f2v(i,2); f2v(i,2) = temp;
    //converting 1-indexing to 0-indexing
    for (int j=0; j<3; j++) f2v(i,j) = int(f2v(i,j) - 1);
  }
}

void FEM::update_U() {

  for (int i=0; i<NF; i++) {
    int ia = f2v(i,0); int ib = f2v(i,1); int ic = f2v(i,2);

    Vector3d a(pos(ia,0), pos(ia,1), 0);
    Vector3d b(pos(ib,0), pos(ib,1), 0);
    Vector3d c(pos(ic,0), pos(ic,1), 0);

    Vector3d crossprod = (a-c).cross(b-c);
    V(i) = abs(crossprod(2)); //k component

    Matrix2d t1; t1.row(0) = pos.row(ia); t1.row(1) = pos.row(ib);
    Matrix2d t2; t2.row(0) = pos.row(ic); t2.row(1) = pos.row(ic);
    Matrix2d D_i = t1 - t2;

    F[i] = D_i.transpose()*B[i];
  }

  for (int i=0; i<NF; i++) {

    Matrix2d F_i = F[i];
    double J_i = max(F_i.determinant(), 0.01);
    double log_J_i = log(J_i);

    double phi_i = mu / 2 * ((F_i.transpose()*F_i).trace() - 2);
    phi_i -= mu * log_J_i;
    phi_i += lam / 2 * log_J_i * log_J_i;

    phi(i) = phi_i; U += V(i)*phi_i;
  }
}

Matrix2d FEM::calcP(int i) {

  Matrix2d F_i = F[i];
  Matrix2d F_T = F_i.inverse().transpose();
  double J = max(F_i.determinant(), 0.01);

  Matrix2d P = (mu*(F_i-F_T)) + (lam*log(J)*F_T);
  return P;
}

void FEM::advance(double dt = 1e-4) {

  node_mass = 1.0;
  for (int i=0; i<NV; i++) force.row(i) = node_mass*gravity;

  for (int i=0; i<NF; i++) {
    Matrix2d P = calcP(i);
    Matrix2d H = (-1*V(i)) * (P*B[i].transpose());

    Vector2d h1(H(0,0), H(1,0)); Vector2d h2(H(0,1), H(1,1));
    int ia = f2v(i,0); int ib = f2v(i,1); int ic = f2v(i,2);
    force.row(ia) += h1; force.row(ib) += h2;
    force.row(ic) += -1*(h1+h2);
  }

  for (int i=0; i<NV; i++) {
    Vector2d acc = force.row(i)/node_mass;
    vel.row(i) += dt*acc;
    vel.row(i) *= exp(-1.0*dt);
  }

  for (int i=0; i<NV; i++) {
    Vector2d disp = pos.row(i); disp -= ball_pos;
    double disp2 = disp.dot(disp);

    if (disp2 <= pow(ball_radius,2)) {
      double NoV = vel.row(i).dot(disp);
      if (NoV < 0) vel.row(i) -= NoV * disp / disp2;
    }

    bool cond0 = ( (pos(i,0)<0) && (vel(i,0)<0) )
              || ( (pos(i,0)>1) && (vel(i,0)>0) );
    bool cond1 = ( (pos(i,1)<0) && (vel(i,1)<0) )
              || ( (pos(i,1)>1) && (vel(i,1)>0) );

    if(cond0) vel(i,0) = 0; if(cond1) vel(i,1) = 0;
  }

  for (int i=0; i<NV; i++) pos.row(i) += dt*vel.row(i);
}

void FEM::save_state(int itr, string statedir) {
  string fname = statedir + "state" + to_string(itr) + ".txt";
  write_state(fname, pos, f2v, ball_pos, ball_radius, NV, NF);
}
