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
	//static const double MASS;

	double get_mobility() {
		return mobility;
	}

	virtual int get_charge_sign() = 0;

protected:

	double mobility;
	double SILICIUM_MOBILITY;
	unsigned type;

	Charge(unsigned type) {
		this->type = type;
	}

	void set_mobility() {
		switch (type) {
		case TYPE_SILICIUM:
			mobility = SILICIUM_MOBILITY;
			break;
		default:
			mobility = SILICIUM_MOBILITY;
		}
	}
};

class Hole: public Charge {

	//static const double MASS = 1;
public:
	Hole(unsigned type) :
			Charge(type) {
		this->SILICIUM_MOBILITY = 4.5e10; //(microm)^2/(Vs)
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
		this->SILICIUM_MOBILITY = 1.35e11; //(microm)^2/(Vs)
		set_mobility();
	}

	int get_charge_sign() {
		return -1;
	}
};
