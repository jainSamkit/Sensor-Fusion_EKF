#include <iostream>
#include "tools.h"
#include <cmath>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	Eigen::VectorXd rmse(4);
	rmse << 0,0,0,0;

	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	for(int i=0;i<estimations.size();++i)
	{
		Eigen::VectorXd residual(4);
		residual = ground_truth[i]-estimations[i];
		residual = residual.array() * residual.array();
		rmse=rmse+residual;
	}
	// cout<<rmse<<endl;
	rmse=rmse/estimations.size();
	// cout<<rmse<<endl;
	rmse=rmse.array().sqrt();
	// cout<<rmse<<endl;
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
	Eigen::MatrixXd Hj_(3,4);
	float px = x_state[0];
    float py = x_state[1];
    float vx = x_state[2];
    float vy = x_state[3];

    float c1 = px*px+py*py;
    float c2 = sqrt(c1);
    float c3 = (c1*c2);

    // cout<<"poiu"<<endl;

    Hj_ << (px/c2), (py/c2), 0, 0,
      -py/c1, px/c1, 0, 0,
      py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

    return Hj_;
}

float Tools::normalize_phi(float phi)
{
	float pi = M_PI;
	float negative_pi=-pi;
	if(phi<negative_pi || phi > pi)
	{
		if(phi<negative_pi)
		{
			while(phi<negative_pi)
			{
				phi+=2*pi;
			}
		}
		if(phi>pi)
		{
			while(phi>pi)
			{
				phi-=2*pi;
			}
		}
	}
	// cout<<phi<<endl;
	return phi;
}
