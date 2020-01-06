#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"
#include<math.h>



using Eigen::MatrixXd;
using Eigen::VectorXd;

/**f
 * Initializes Unscented Kalman filter
 */
//namespace plt = matplotlibcpp;

UKF::UKF()
{
   is_initialized_ = false;


  use_laser_ = true;  // if this is false, laser measurements will be ignored (except during init)


  use_radar_ = true;  // if this is false, radar measurements will be ignored (except during init)


  x_ = VectorXd(5);   // initial state vector
  x_.fill(0);


  P_ = MatrixXd(5, 5); // initial covariance matrix

  Q_ = MatrixXd(2,2);  // Process Noise Covariance


  std_a_ = 2;          //ms-1  // Process noise standard deviation of  longitudinal acceleration in m/s^2

                       // Process noise standard deviation of yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;      // rad/s^2   vary Normalised Innovation Squared
  
  std_laspx_ = 0.15;   // Laser measurement noise standard deviation position1 in m

  std_laspy_ = 0.15;   // Laser measurement noise standard deviation position2 in m

  std_radr_ = 0.3;     // Radar measurement noise standard deviation radius in m

  std_radphi_ = 0.03;  // Radar measurement noise standard deviation angle in rad
  
  std_radrd_ = 0.3;    // Radar measurement noise standard deviation radius change in m/s

  n_aug_ = 7;

  n_x_ = 5 ;

  lambda_ = 3 - n_x_ ;

  n_z_ = 3 ;

  time_us_ = 0;

  previous_timestamp_ = 0.0;

  X_sig = MatrixXd(5, 15);
  X_sig.fill(0.0);

  dt = 0;


  S_R = MatrixXd(3 , 3);
  S_R.fill(0);                               // RADAR Measurement Covariance

  S_L = MatrixXd(2,2);
  S_L.fill(0);                              // LIDAR Measurement Covariance


  weights = VectorXd(2* n_aug_ +1) ;

  cycle = 0;


}

UKF::~UKF() {}


void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{

	if(!is_initialized_)
	{
		if(use_laser_ == true && meas_package.sensor_type_ == MeasurementPackage::LASER)
		{
	      // Update with Laser Measurements
		 VectorXd inter = meas_package.raw_measurements_ ;
		 x_(0) = inter(0);
		 x_(1) = inter(1);
		}
		else if(use_radar_ == true && meas_package.sensor_type_ == MeasurementPackage::RADAR)
		{
	       // Update with Radar measurements
			double rho =  meas_package.raw_measurements_(0);
			double phi = meas_package.raw_measurements_(1);
			double rhodot = meas_package.raw_measurements_(2);
			x_(0) = rho * cos(phi) ;
			x_(1) = rho * sin(phi) ;
			double vx = rhodot * cos(phi);
			double vy = rhodot * sin(phi);
			x_(2) = sqrt(vx*vx + vy*vy);
			x_(3) = 0;
			x_(4) = 0;
		}

       P_ = MatrixXd::Identity(n_x_, n_x_);  // Covariance matrix set to Identity


       is_initialized_ = true;

       cycle = 0;


       previous_timestamp_ = meas_package.timestamp_ ;

       Q_ << std_a_ , 0,
      			     0 , std_yawdd_ ;

       std::cout << std::endl;
       std::cout << "1     $$$$$$$$$$$$$$$$$$$$$$   Initialization Done!   $$$$$$$$$$$$$$$$$$$$$$" << std::endl;
       std::cout << std::endl << std::endl;

       return ;

	}

	dt  = (meas_package.timestamp_ - previous_timestamp_)/1000000.0 ;
	previous_timestamp_ = meas_package.timestamp_ ;

	std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$  CYCLE : " << cycle << "    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl;

	cycle++;

	Prediction(dt);  //  Sigma points generated , augmented , mapped and x_ and P_ is triangulated .

	// If next meas is LIDAR or RADAR meas Update with the state and Covariance using Kalman Gain
    if(cycle >= 497){

    	std::cout << "**********  PLOT HERE ? ******************" << std::endl;
        //plt::plot(NIS_R_vec);
        //plt::savefig("NIS_R_vec.pdf");
    	std::cout << "NIS_R_vec "<<std::endl ;
    	for(int i =0 ; i < NIS_R_vec.size() ; ++i)
    	{
    		std::cout << NIS_R_vec[i] <<std::endl;
    	}
    }

	if(use_laser_ == true && meas_package.sensor_type_ == MeasurementPackage::LASER)
	{
        UpdateLidar(meas_package);
	}
	else if(use_radar_ == true && meas_package.sensor_type_ == MeasurementPackage::RADAR)
	{
       UpdateRadar(meas_package);
	}

}


void UKF::Prediction(double delta_t)
{

     std::cout << "$$$$$$$$$$$$$$$$$$$$$$$    PREDICTION    $$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl;
     std::cout << std::endl;



     VectorXd x_aug_ = VectorXd(n_aug_);
     x_aug_.head(n_x_) = x_ ;
     x_aug_(5) = 0 ;
     x_aug_(6) = 0 ;


     MatrixXd P_aug = MatrixXd(n_aug_,  n_aug_);

     P_aug.topLeftCorner(n_x_ , n_x_) = P_;

     P_aug(5,5) = std_a_ * std_a_;

     P_aug(6,6) = std_yawdd_ * std_yawdd_ ;

     MatrixXd A = MatrixXd(n_aug_ , n_aug_);

     A = P_aug.llt().matrixL();  // CHOLESKY DECOMPOSITION returning a PSD matrix

     MatrixXd Xsig_aug_ = MatrixXd(7 , 15);

     lambda_ = 3 - n_x_;

     Xsig_aug_.col(0) = x_aug_ ;

     for(int i =0 ; i < n_aug_ ; ++i)
     {

   	  Xsig_aug_.col(i + 1) = x_aug_ + sqrt(lambda_ + n_aug_)* A.col(i) ;

   	  Xsig_aug_.col(i + n_aug_ + 1) = x_aug_ - sqrt(lambda_ + n_aug_)* A.col(i);

     }


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     	 	 	 	 	 	 	 	 	 	        //   PREDICTION   //
                                      //   MAPPING THROUGH A NON LINEAR PROCESS   //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




    for(int i =0 ; i < 2*n_aug_ +1 ; ++i)
    {
  	  double px = Xsig_aug_.col(i)[0];
  	  double py = Xsig_aug_.col(i)[1];
  	  double v = Xsig_aug_.col(i)[2];
  	  double yaw = Xsig_aug_.col(i)[3];
  	  double yaw_rate = Xsig_aug_.col(i)[4];
  	  double acc_noise = Xsig_aug_.col(i)[5];
  	  double yawrate_noise = Xsig_aug_.col(i)[6];
      double px_m;
      double py_m;
      double v_m;
      double yaw_m;
      double yaw_rate_m;

  	  if(fabs(yaw_rate) > 0.001 )
  	  {
  		  px_m = px + ((v/yaw_rate)*( sin(yaw + yaw_rate*dt) - sin(yaw) )  ) + ((1/2)*(dt*dt)* cos(yaw)*acc_noise);
  		  py_m = py + ((v/yaw_rate) * (-cos(yaw + yaw_rate*dt) + cos(yaw) )) + ((1/2)*(dt*dt)* cos(yaw)*acc_noise);
  		  v_m = v + 0 + dt*acc_noise;
  		  yaw_m = yaw + yaw_rate*dt+ (1/2)*(dt*dt)*yawrate_noise;
  		  yaw_rate_m = yaw_rate + 0 + dt*yawrate_noise;

  	  }
  	  else
  	  {
   		 px_m = px + (v*cos(yaw)*dt) + ((1/2)*(dt*dt)* cos(yaw)*acc_noise);
   		 py_m = py + (v*cos(yaw)*dt) + ((1/2)*(dt*dt)* cos(yaw)*acc_noise);
   		 v_m = v + 0 + dt*acc_noise;
   		 yaw_m = yaw + 0 + (1/2)*(dt*dt)*yawrate_noise;
   		 yaw_rate_m = yaw_rate + 0 + dt*yawrate_noise;
  	  }

  	    X_sig(0 , i) = px_m;
  	    X_sig(1 , i) = py_m;
  	    X_sig(2 , i) = v_m;
  	    X_sig(3 , i) = yaw_m;
  	    X_sig(4 , i) = yaw_rate_m;
    }


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  ///////////////////////////   TRIANGULATE MEAN AND COVARIANCE    ////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




      for(int i =0 ; i < 2*n_aug_+1 ; ++i)
      {
     	   if(i == 0)
     	   {
     		   weights(i) = lambda_/(lambda_ + n_aug_);
     	   }
     	   else
     	   {
               weights(i) = 1/(2*(lambda_ + n_aug_));
     	   }
      }


	  x_.fill(0.0);

	  for(int j= 0 ; j < 2*n_aug_ +1 ; ++j)
	  {
	   x_ = x_ + weights(j)*X_sig.col(j);
	  }

	  P_.fill(0.0);

	  for(int k = 0 ; k < 2*n_aug_+1 ; ++k)
	  {

	     VectorXd x_diff = X_sig.col(k) - x_ ;

		 while(x_diff(3) > M_PI)
	     {
			x_diff(3) -= 2.*M_PI ;

		 }

		 while(x_diff(3) < -M_PI){
			 x_diff(3) += 2.*M_PI ;

		 }

		 P_ = P_ + weights(k)*(x_diff)* x_diff.transpose() ;

	  }

}



void UKF::UpdateLidar(MeasurementPackage meas_package)
{


	std::cout << "################     UPDATING WITH  LIDAR MEASUREMENTS    ##################"<< std::endl  << std::endl;

	VectorXd inter_ = meas_package.raw_measurements_ ;

	VectorXd z_ = VectorXd(2);
    z_(0) = inter_(0);
    z_(1) = inter_(1);

    MatrixXd Z_sig = MatrixXd(2 , 2*n_aug_+1);

    for(int i = 0 ; i < 2*n_aug_+1 ; ++i)
    {
    	Z_sig.col(i)(0) = X_sig.col(i)(0);
    	Z_sig.col(i)(1) = X_sig.col(i)(1);
    }

    ///////////  MEAN PREDICTED MEASUREMENT /////////////////
    VectorXd z_pred = VectorXd(2);

    lambda_ = 3 - n_x_;

    z_pred.fill(0.0);

    for(int i =0 ; i < 2*n_aug_+1 ; ++i)
    {
   	   if(i == 0)
   	   {
   		   weights(i) = lambda_/(lambda_ + n_aug_);
   	   }
   	   else
   	   {
             weights(i) = 1/(2*(lambda_ + n_aug_));
   	   }
    }


    for(int i = 0 ; i< 2 * n_aug_ +1 ; ++i){

    	z_pred = z_pred + weights(i) * Z_sig.col(i);

    }


    ///////////  MEASUREMENT COVARIANCE MATRIX ///////////////
	
    VectorXd z_diff = VectorXd(2);

    for(int i =0 ; i < 2*n_aug_+1 ; ++i)
    {
   	   if(i == 0)
   	   {
   		   weights(i) = lambda_/(lambda_ + n_aug_);
   	   }
   	   else
   	   {
             weights(i) = 1/(2*(lambda_ + n_aug_));
   	   }
    }


    for(int i =0 ; i < 2*n_aug_+ 1 ; ++i)
   	{
    	z_diff =  Z_sig.col(i) - z_pred ;
    	S_L = S_L + weights(i)* z_diff * z_diff.transpose();
   	}


    MatrixXd R = MatrixXd(2 , 2);

    R << std_laspx_*std_laspx_ , 0 ,
    		                     0,  std_laspy_*std_laspy_ ;

    S_L = S_L + R ;


	MatrixXd T = MatrixXd(5,2);

	T.fill(0.0);

    VectorXd x_diff = VectorXd(n_x_);

	for(int i = 0 ; i < 2* n_aug_ +1 ; ++i)
	{

		x_diff = X_sig.col(i) - x_ ;

		while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
		while(x_diff(3) < M_PI) x_diff(3) += 2.*M_PI;

		z_diff = Z_sig.col(i) - z_pred ;

		T = T + weights(i) * x_diff * z_diff.transpose();

	}

	MatrixXd K = MatrixXd(n_x_ , n_x_);

	K.fill(0.0);

	K = T * S_L.inverse();

	double NIS_val_RADAR = z_diff.transpose() * S_L.inverse() * z_diff ;

	NIS_L_vec.push_back(NIS_val_RADAR);

	x_ = x_ + K * (z_ - z_pred);  // z_ -> SENSOR MEASUREMENT    z_pred -> PREDICTED MEASUREMETN

	P_ = P_ - K * S_L * K.transpose();

}


void UKF::UpdateRadar(MeasurementPackage meas_package) {


	 std::cout << "$$$$$$$$$$$$  UPDATING WITH RADAR MEASUREMENTS $$$$$$$$$$$$$$$" << std::endl;

     VectorXd z = meas_package.raw_measurements_;


     int n_z = 3;


     MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	  /// MAP points into measurement space

	  for (int i = 0; i < 2 * n_aug_ + 1; i++)
	  {

		double p_x = X_sig(0, i);
		double p_y = X_sig(1, i);
		double v   = X_sig(2, i);
		double yaw = X_sig(3, i);

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;


		Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);
		Zsig(1, i) = atan2(p_y, p_x);
		Zsig(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);
	  }

	  //mean predicted measurement
	  VectorXd z_pred = VectorXd(n_z);

	  for(int i =0 ; i < 2*n_aug_+1 ; ++i)
	  {
		   if(i == 0)
		   {
			   weights(i) = lambda_/(lambda_ + n_aug_);
		   }
		   else
		   {
			   weights(i) = 1/(2*(lambda_ + n_aug_));
		   }
	  }

	  z_pred.fill(0.0);
	  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + weights(i) * Zsig.col(i);
	  }



	  // MEASUREMENT COVARINCE matrix
	  S_R.fill(0.0);
	  for (int i = 0; i < 2 * n_aug_ + 1; i++)
	  {

		VectorXd z_diff = Zsig.col(i) - z_pred;


		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

		S_R = S_R + weights(i) * z_diff * z_diff.transpose();
	  }


	  MatrixXd R = MatrixXd(n_z, n_z);
	  R << std_radr_*std_radr_,                       0,                     0,
							 0, std_radphi_*std_radphi_,                     0,
							 0,                       0, std_radrd_*std_radrd_;
	  S_R = S_R + R;


	  MatrixXd Tc = MatrixXd(n_x_, n_z);


	  Tc.fill(0.0);

	  for (int i = 0; i < 2 * n_aug_ + 1; i++)
	  {

		VectorXd z_diff = Zsig.col(i) - z_pred;

		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;


		VectorXd x_diff = X_sig.col(i) - x_;

		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		Tc = Tc + weights(i) * x_diff * z_diff.transpose();

	  }


	  MatrixXd K = Tc * S_R.inverse();

	  VectorXd z_diff = z - z_pred;

	  //angle normalization
	  while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
	  while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

	  //calculate NIS

	  double NIS_value_RADAR = z_diff.transpose() * S_R.inverse() * z_diff;

	  NIS_R_vec.push_back(NIS_value_RADAR); // k length at kth  cycle



	  //update state mean and covariance matrix
	  x_ = x_ + K * z_diff;
	  P_ = P_ - K*S_R*K.transpose();


}




