/*
 * Charge.hpp
 *
 *  Created on: 14 févr. 2016
 *      Author: aschils
 */

#pragma once

#define TYPE_SILICIUM 1

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
		return mobility;
	}

	double get_mobility_300K(){
		return mobility_300K;
	}

	virtual int get_charge_sign() = 0;

protected:

	double mobility_300K;
	double SILICON_MOBILITY;
	unsigned type_material;
	double beta_0;
	double saturation_velocity_300K;


	Charge(unsigned type_material) {
		this->type_material = type_material;
	}

	void set_mobility() {
		switch (type_material) {
		case TYPE_SILICIUM:
			mobility_300K = SILICON_MOBILITY;
			break;
		default:
			mobility_300K = SILICON_MOBILITY;
		}
	}
};

class Hole: public Charge {

	//static const double MASS = 1;
public:
	Hole(unsigned type) :
			Charge(type) {
		this->SILICON_MOBILITY = 4.5e10; //(microm)^2/(Vs)
		this->beta_0 = 1.213;
		this->saturation_velocity_300K = 8.37e10;
		set_mobility();
	}

	int get_charge_sign() {
		return 1;
	}
};

class Electron: public Charge {

public:
	Electron(unsigned type) :
			Charge(type) {
		this->SILICON_MOBILITY = 1.35e11; //(microm)^2/(Vs)
		this->beta_0 = 1.109;
		this->saturation_velocity_300K = 1.07e11;
		set_mobility();
	}

	int get_charge_sign() {
		return -1;
	}
};
