#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"
#include "measurement_package.h"
#include "tools.h"
#include<vector>
#include<string>


class UKF
{
 public:


		UKF();       // CONSTRUCTOR

		virtual ~UKF(); // DESTRUCTOR


		void ProcessMeasurement(MeasurementPackage meas_package); // the latest measurement data

		/**
		 * Prediction Predicts sigma points, the state, and the state covariance
		 * matrix
		 * @param delta_t Time between k and k+1 in s
		 */


		void Prediction(double delta_t);

		  /**
		   * Updates the state and the state covariance matrix using a laser measurement
		   * @param meas_package The measurement at k+1
		   */
		void UpdateLidar(MeasurementPackage meas_package);

		  /**
		   * Updates the state and the state covariance matrix using a radar measurement
		   * @param meas_package The measurement at k+1
		   */
	    void UpdateRadar(MeasurementPackage meas_package);



	    bool is_initialized_;   // initially set to false, set to true in first call of ProcessMeasurement


	    bool use_laser_;        // if this is false, laser measurements will be ignored (except for init)


	    bool use_radar_;        // if this is false, radar measurements will be ignored (except for init)


	    Eigen::VectorXd x_;     // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad


	    Eigen::MatrixXd P_;  		  // state covariance matrix

	    Eigen::MatrixXd Q_;          // Process Noise Covarinace MAtrix

		long long time_us_;    // time when the state is true, in us

		float dt ;

		long long previous_timestamp_ ;


		double std_a_;         // Process noise standard deviation longitudinal acceleration in m/s^2


		double std_yawdd_;    // Process noise standard deviation yaw acceleration in rad/s^2


		double std_laspx_;    // Laser measurement noise standard deviation position1 in m


	    double std_laspy_;    // Laser measurement noise standard deviation position2 in m


		double std_radr_;     // Radar measurement noise standard deviation radius in m


	    double std_radphi_;    // Radar measurement noise standard deviation angle in rad


	    double std_radrd_ ;    // Radar measurement noise standard deviation radius change in m/s


	    Eigen::VectorXd weights_;   // Weights of sigma points


	    int n_x_;   // State dimension

	    int n_aug_; // Augmented state dimension

	    int n_z_;  	    //RADAR state dimension

	    double lambda_;  // Sigma point spreading parameter


	    Eigen::MatrixXd X_sig;  // predicted sigma points matrix

	    Eigen::MatrixXd S_R ;      // Measurement Covariance RADAR


	    Eigen::MatrixXd S_L ;      // Measurement Covaraince LIDAR

	    Eigen::VectorXd weights;

	    int cycle;


};

#endif  // UKF_H
