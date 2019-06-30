#include "bilayer.h"

int main(int argc, char *argv[]){
	Model model;
	Configuration config;
	Observables ob;
	int looptimes;
	
	initialization(model, config, argc, argv);
	warmup(model, config);
	looptimes = computelooptimes(model, config);
	display(model, config);
	simulation(model, config,ob, looptimes);
	writetofile(model, ob);
	return 0;
}