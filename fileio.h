#include <iostream>
#include <fstream>
#include <string>

#include "Eigen/Dense"
using namespace Eigen;
using namespace std;

Vector2i get_params(string fname) {

  Vector2i pars(2,0); //index 0 for pos, index 1 for f2v

  fstream f; f.open(fname, ios::in);
  if (f.is_open()) {

    string tp;
    while(getline(f, tp)) {

      stringstream iss(tp); string temp_str;
      int temp_int, store_int;

      while(!iss.eof()) {
        iss >> temp_str; //take words into temp_str one by one
        if(stringstream(temp_str) >> temp_int) store_int = temp_int;
        else {
          if (temp_str == "NumPhysNodes") pars(0) = store_int;
          if (temp_str == "NumPhysElems") pars(1) = store_int;
        }
        temp_str = ""; //clear temp string
      }
    }
  }

  else {cout << "ERROR: could not open file " << fname << endl;}
  return pars;
}

MatrixXd read_data(string fname, int numrows, int numcols) {

  MatrixXd result(numrows, numcols); int rowind = 0;
  fstream f; f.open(fname, ios::in);

  if (f.is_open()){
     string tp;

     while(getline(f, tp)){

        std::istringstream iss(tp);
        long double value; int colind = 0;

        while ( iss >> value ) result(rowind, colind++) = value;
        rowind++;
        if(rowind==numrows) break;
     }
     f.close(); //close the file object.
  }
  else {cout << "ERROR: could not open file " << fname << endl;}
  return result;
}

void write_state(string fname, MatrixXd pos, MatrixXd f2v,
                 Vector2d ball_pos, double ball_radius, int NV, int NF) {

   // open a file in write mode.
   ofstream outfile; outfile.open(fname, ios::out);

   outfile << "pos" << endl;

   for (int i=0; i<NV; i++) outfile << "(" << pos(i,0) << " " << pos(i,1) << " " << 0 << ")" << endl;
   outfile << "f2v" << endl;
   //outfile << "NV is " << NV << endl;

   for (int i=0; i<NF; i++) outfile << "3(" << f2v(i,0) << " " << f2v(i,1) << " " << f2v(i,2) << ")" << endl;
   //outfile << "NF is " << NF << endl;
   outfile << "ball_pos" << endl; outfile << ball_pos << endl;
   outfile << "ball_radius" << endl; outfile << ball_radius << endl;

   outfile.close();
}
