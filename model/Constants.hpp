/*
 * Constants.hpp
 *
 *  Created on: 15 févr. 2016
 *      Author: aschils
 */

#pragma once

//Physical constants
const double ELECTRON_CHARGE = 1.602176e-13; // microCoulomb
const double GAS_CONSTANT = 8.3144598;// J/(K mol)
const double ATMOSPHERIC_PRESSURE = 1.01325e5; // Pascal
const double MOLAR_MASS_HELIUM = 0.0040026; // kg/mol
const double LIGHT_SPEED = 299792458e6; //microm/s

#define TYPE_SILICON 1
#define TYPE_HELIUM 2

//Energy per length deposited by particle passing through matter
const double HELIUM_DE_DX = 0.032; //eV/microm

//Average energy to produce an ion pair
const double HELIUM_ION_PAIR_E = 24.5; //eV
