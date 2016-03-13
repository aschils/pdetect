/*
 * Charge.hpp
 *
 *  Created on: 14 févr. 2016
 *      Author: aschils
 */

#pragma once

#include "Constants.hpp"

class Charge {

public:

	double get_mobility_saturation(Tensor<1, 2> electric_field) {

		double electric_field_y = fabs(electric_field[1]);
		//double beta_e = beta_0_e * pow((Temp / 300), beta_exp_e);
		//double beta_h = beta_0_h * pow((Temp / 300), beta_exp_h);
		double beta = beta_0;
		double inv_beta = 1.0 / beta;

		double mob_over_velocity = mobility_300K / saturation_velocity_300K;
		double mobility = mobility_300K
				/ pow((1 + pow((mob_over_velocity * electric_field_y), beta)),
						inv_beta);

		/*
		 *
		 * mu = mu_300K / (1+ (mu * |E_y|/saturation_velocity)^beta)^(1/beta)
		 *
		 */


		return mobility;
	}

	double get_mobility_300K(){
		return mobility_300K;
	}

	virtual int get_charge_sign() = 0;

protected:

	double mobility_300K;
	double SILICON_MOBILITY, HELIUM_MOBILITY;
	double SILICON_SATURATION_VELOCITY_300K, HELIUM_SATURATION_VELOCITY_300K;
	unsigned type_material;
	double beta_0;
	double saturation_velocity_300K;

	Charge(){}

	Charge(unsigned type_material) {
		this->type_material = type_material;
	}

	void set_mobility() {
		switch (type_material) {
		case TYPE_SILICON:
			mobility_300K = SILICON_MOBILITY;
			saturation_velocity_300K = SILICON_SATURATION_VELOCITY_300K;
			break;
		case TYPE_HELIUM:
			mobility_300K = HELIUM_MOBILITY;
			saturation_velocity_300K = HELIUM_SATURATION_VELOCITY_300K;
			break;
		default:
			mobility_300K = SILICON_MOBILITY;
			saturation_velocity_300K = SILICON_SATURATION_VELOCITY_300K;
		}
	}
};

class Hole: public Charge {

	//static const double MASS = 1;
public:

	Hole(){}

	Hole(unsigned type) :
			Charge(type) {
		this->SILICON_MOBILITY = 4.5e10; //(microm)^2/(Vs)
		this->HELIUM_MOBILITY = 1e7;
		this->beta_0 = 1.213;
		this->SILICON_SATURATION_VELOCITY_300K = 8.37e10;
		this->HELIUM_SATURATION_VELOCITY_300K = 5e7; //microm/s
		//Vitesse de sat de 5cm / micros
		//E/N 10 V/m^2
		//mob: 10^-4 m2 atm /(Vs)
		set_mobility();
	}

	int get_charge_sign() {
		return 1;
	}
};

class Electron: public Charge {

public:

	Electron(){}

	Electron(unsigned type) :
			Charge(type) {
		this->SILICON_MOBILITY = 1.35e11; //(microm)^2/(Vs)
		this->HELIUM_MOBILITY = 1e10;
		this->beta_0 = 1.109;
		this->SILICON_SATURATION_VELOCITY_300K = 1.07e11; //Vitesse 10^5 m/s pr Helium
		this->HELIUM_SATURATION_VELOCITY_300K = 1e11;
		set_mobility();
	}

	int get_charge_sign() {
		return -1;
	}
};
