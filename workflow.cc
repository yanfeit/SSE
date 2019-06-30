#include "bilayer.h"

void initialization(Model & model, Configuration & config, int argc, char **argv){
	
	readfile(model, argc, argv);
	makelattice(model);
	calweight(model);
	caltable(model);
	askspace(model, config);
	
	cout << "Successful loading files!\nNow begin Monte Carlo Warm up process!\n";
	
}

void warmup(const Model & model, Configuration & config){
	
	for (int i = 0; i < model.warmsteps; ++i){
		int looptimes = 0;
		diagupdate(model, config);
		linkvertices(model, config);
		while (looptimes < config.msl)
			looptimes += loopupdate(model, config);
		updatestate(model, config);
		adjustcutoff(model, config);
		if (i%(model.warmsteps/10) == 0)
			cout << i << " warm up steps finished!\n";
	}
	cout << "Now begin to compute loop times per Monte Carlo steps\n";	
}

int computelooptimes(const Model & model, Configuration & config){
	
	int times = 0, turns = 100;
	for (int i = 0; i < turns; ++i){
		int looptimes = 0;
		diagupdate(model, config);
		linkvertices(model, config);
		while (looptimes < config.msl){
			looptimes += loopupdate(model, config);
			times++;
		}
		updatestate(model, config);
	}
	if (times/turns >= 1)
		return times/turns;
	else
		return 1;		
}

void simulation(const Model &model, Configuration &config, Observables & ob, int looptimes){
	
	cout << looptimes << " times needed per Monte Carlo steps\n";
	for (int i = 0; i < model.mcsteps; ++i){
		
		diagupdate(model, config);
		linkvertices(model, config);
		for (int j = 0; j < looptimes; ++j)
			loopupdate(model, config);
		updatestate(model, config);
		vector<double> winding = computewind(model, config);
		ob.density.push_back(computedensity(config));
		ob.energy.push_back(computeenergy(model, config));
		ob.rhosTT.push_back(winding[0]);
		ob.rhosTB.push_back(winding[1]);
		ob.rhos.push_back(winding[2]);
		if (i%(model.mcsteps/10) == 0){
			cout << i << " simulation finished!\n";
		}
		
	}
	cout << "Now ends.";
}

void display(const Model & model, const Configuration & config){
	
	cout << "The cut off length is: " << config.msl << endl;
	
}

// need to be modified here...
void writetofile(const Model & model, const Observables &ob){
	
	// string file1;
	string file2;
	ostringstream os[5];
	os[0] << model.lx;
	os[1] << model.t_parallel;
	os[2] << model.t_perp;
	os[3] << model.chempot;
	os[4] << model.beta;
	
	/* file1 = "lx" + os[0].str() + "t0" + os[1].str() + "t1" + os[2].str() + "mu" + os[3].str()
		+ "beta" + os[4].str() + "raw.dat";
	*/
	file2 = "lx" + os[0].str() + "t0" + os[1].str() + "t1" + os[2].str() + "mu" + os[3].str()
		+ "beta" + os[4].str() + ".dat";
	/*
	ofstream ofile(&file1[0]);
	if (ofile.is_open()){
		for (int i = 0; i < int(ob.density.size()); ++i){
			ofile << ob.density[i] << '\t' << ob.energy[i] << '\t' << ob.rhosTT[i] << '\t' <<
				ob.rhosTB[i] << '\t' << ob.rhos[i] << endl;
		}		
	}
	*/
	
	vector<double> quantity[10];
	quantity[0] = binning(ob.density);
	quantity[1] = binning(ob.energy);
	quantity[2] = binning(ob.rhosTT);
	quantity[3] = binning(ob.rhosTB);
	quantity[4] = binning(ob.rhos);
	quantity[5] = computek(model, ob.density);
	
	ofstream result(&file2[0]);
	if (result.is_open()){
		
		
		result << "lx\t" << model.lx << endl;
		result << "t_parallel\t" << model.t_parallel << endl;
		result << "t_perp\t" << model.t_perp << endl;
		result << "interaction\t" << model.interaction << endl;
		result << "mu\t" << model.chempot << endl;
		result << "mpps\t" << model.mpps << endl;
		result << "beta\t" << model.beta << endl;
		result << "thermalization\t" << model.warmsteps << endl;
		result << "simulation\t" << model.mcsteps << endl;
		result << "seed\t" << model.seed << endl;
		
		result << "density\t" << quantity[0][0] << '\t' << quantity[0][1] << endl;
		result << "energy\t" << quantity[1][0] << '\t' << quantity[1][1] << endl;
		result << "rhosTT\t" << quantity[2][0] << '\t' << quantity[2][1] << endl;
		result << "rhosTB\t" << quantity[3][0] << '\t' << quantity[3][1] << endl;
		result << "rhos\t" << quantity[4][0] << '\t' << quantity[4][1] << endl;
		result << "compressibility\t" << quantity[5][0] << '\t' << quantity[5][1] << endl;		
	}
	
	
}














