#include <iostream>
#include "src/constants.hpp"
#include "src/model.hpp"
#include "src/integrator.hpp"

int sundial_test() {// в пределах дня

	sundial_model model(rad(55), rad(37), get_JDN(2024, 3, 15, 0, 0, 0));
	DormandPrinceIntegrator integrator(1e-16l);
	integrator.run(model);
	model.load_res2file("res1.txt");

	return 0;
}

int blag_test() { // год
	blag_time_model model;
	DormandPrinceIntegrator integrator(1e-16l);
	integrator.run(model);
	model.load_res2file("res2.txt");

	return 0;
}

int main() {

	// Запуск теста для модели солнечных часов
	std::cout << "Running sundial test" << std::endl;
	sundial_test();
	std::cout << "Sundial test completed" << std::endl;

	// Запуск теста для анализа продолжительности светового дня
	std::cout << "Running daylight duration analysis" << std::endl;
	blag_test();
	std::cout << "Daylight duration analysis completed" << std::endl;

	return 0;
}