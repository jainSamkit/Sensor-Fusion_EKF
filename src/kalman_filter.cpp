#include "kalman_filter.h"
#include "tools.h"
#include<iostream>
#include "Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in,MatrixXd &H_, MatrixXd &R_in,MatrixXd &R_, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_laser_ = H_in;
  Hj_=H_;
  R_radar_= R_;
  R_laser_ = R_in;
  Q_ = Q_in;

}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_*x_;
  Eigen::MatrixXd F_t(4,4);
  F_t=F_.transpose();
  P_ = F_ * P_ * F_t + Q_;
  // cout<<"predicted x"<<endl;
  // cout<<x_<<endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_laser_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_laser_.transpose();
  MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  x_=x_+(K*y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_=(I-K*H_laser_) * P_;
  // cout<<"samitnvjnv"<<endl;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  // cout<<"Radar Update step"<<endl;
  // cout<<"Printing present state"<<endl;
  // cout<<x_<<endl;
  float px= x_[0];
  float py= x_[1];
  float vx= x_[2];
  float vy= x_[3];

  float rho,phi,rho_;
  VectorXd z_pred = VectorXd(3);

  if(fabs(px)<0.0001 || fabs(py)<0.001)
  {
    if(fabs(px)<0.0001)
    {
      px=0.0001;
      // cout<<"Too small px"<<endl;
    }

    if(fabs(py)<0.0001)
    {
      py=0.0001;
      // cout<<"Too small py"<<endl;
    }
    // cout<<"small"<<endl;
    rho=sqrt(px*px+py*py);
    phi=0;
    rho_=0;
  }
  else
  {
    // cout<<"big enough"<<endl;
    rho=sqrt(px*px+py*py);
    phi=atan2(py,px);
    rho_=(px*vx+py*vy)/rho;
  }
  z_pred<<rho,phi,rho_;
  VectorXd y = z - z_pred;
  
  y[1] = tools1.normalize_phi(y[1]);
  
  MatrixXd Ht = Hj_.transpose();
  MatrixXd S = Hj_ * P_ * Ht + R_radar_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  // cout<<"X before update"<<endl;
  // cout<<x_<<endl;
  x_=x_+(K*y);
  // cout<<"Printing updated x"<<endl;
  // cout<<x_<<endl;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_=(I-K*Hj_) * P_;
  // cout<<"samitnvjnv"<<endl;
}

