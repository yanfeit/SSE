#include "bilayer.h"

double computedensity(const Configuration & config){
	
	return accumulate(config.state.begin(), config.state.end(), 0)
		/double(config.state.size());	
}

double computeenergy(const Model & model, const Configuration & config){
	
	return (-config.nh/model.beta+model.nb*model.c)/model.ns;
	
}

vector<double> computewind(const Model & model, const Configuration & config){
	
	int op, b, site, sum = 0;
	// NL[0] number of boson hops to the left in the bottom layer
	// NL[1] number of boson hops to the left in the top layer
	int NL[2] = {0, 0}, NR[2]= {0, 0}, NU[2]= {0, 0}, ND[2]= {0, 0};
	
	// rhos[0] is the superfluid density for top layer top layer coupling,
	// rhos[1] is the superfluid density for top layer bottom layer coupling.
	// rhos[2] is the conventional way to compute superfluid denstiy.
	vector<double> rhos(3, 0.0);
	
	int period[4]; // bond [0, lxly) in bottom layer, bond [lxly, 2lxly) in top layer,
	               // bond [2lxly, 3lxly) in bottom layer, bond [3lxly, 4lxly) in top layer,
	
	vector<int> coupling(model.lx*model.ly, 0);
	
	period[0] = model.lx*model.ly;
	period[1] = 2*period[0];
	period[2] = 3*period[0];
	period[3] = 4*period[0];
	
	
	
	for (int i = 0; i < config.msl; ++i){		
		op = config.opstring[i];		
		if ( op%2 == 1){
			b = op/2 - 1;
			if (b < period[0] || (b >= period[1] && b < period[2] )){
	            if ( ((b+1)%model.lx == 0) && (b < period[0]) ) {
	                if (config.leg[i][0] == config.leg[i][2] - 1) {
	                    ++NL[0];
	                }
	                else if ( config.leg[i][0] == config.leg[i][2] + 1 ){
	                    ++NR[0];
	                }
	            }
	            else if ( b >= period[2] - model.lx){
	                if (config.leg[i][0] == config.leg[i][2] - 1) {
	                    ++ND[0];
	                }
	                else if ( config.leg[i][0] == config.leg[i][2] + 1 ){
	                    ++NU[0];
	                }
	            }
			}
			
			else if ( (b >= period[0] && b< period[1]) || (b>= period[2] && b < period[3]) ){
	            if ( ((b+1)%model.lx == 0) && (b < period[1]) ) {
	                if (config.leg[i][0] == config.leg[i][2] - 1) {
	                    ++NL[1];
	                }
	                else if ( config.leg[i][0] == config.leg[i][2] + 1 ){
	                    ++NR[1];
	                }
	            }
	            else if ( b >= period[3] - model.lx){
	                if (config.leg[i][0]  == config.leg[i][2] - 1) {
	                    ++ND[1];
	                }
	                else if ( config.leg[i][0] == config.leg[i][2] + 1 ){
	                    ++NU[1];
	                }
	            }
			}
			
			else{
				
				site = model.bsites[b][0];
                if (config.leg[i][0] == config.leg[i][2] - 1) {
                    --coupling[site];
                }
                else if ( config.leg[i][0] == config.leg[i][2] + 1 ){
                    ++coupling[site];
                }
			}			
		}		
	}
	for (int i =0; i < period[0]; ++i)
		sum += i*coupling[i];
	/*
	rhos[0] = double(model.lx*model.lx*(NL[1] - NR[1])*(NL[1] - NR[1]) + model.ly*model.ly*(NU[1] - ND[1])*(NU[1] - ND[1]) -
		(model.lx*(NL[1] - NR[1]) + model.ly*(NU[1] - ND[1]))*sum);
	rhos[1] = double(model.lx*model.lx*(NL[1] - NR[1])*(NL[0] - NR[0]) + model.ly*model.ly*(NU[1] - ND[1])*(NU[0] - ND[0]) +
		(model.lx*(NL[0] - NR[0]) + model.ly*(NU[0] - ND[0]))*sum);
	*/
	rhos[0] = double(model.lx*model.lx*(NL[1] - NR[1])*(NL[1] - NR[1]) + model.ly*model.ly*(NU[1] - ND[1])*(NU[1] - ND[1]));
		
	rhos[1] = double(model.lx*model.lx*(NL[1] - NR[1])*(NL[0] - NR[0]) + model.ly*model.ly*(NU[1] - ND[1])*(NU[0] - ND[0]));
		
	rhos[2] = double( model.lx*model.lx*(NL[0]+NL[1]-NR[0]-NR[1])*(NL[0]+NL[1]-NR[0]-NR[1]) + 
		model.ly*model.ly*(NU[0]+NU[1]-ND[0]-ND[1])*(NU[0]+NU[1]-ND[0]-ND[1]) );
	
	rhos[0] /= (4.0*model.t_parallel*model.beta*period[0]);
	rhos[1] /= (4.0*model.t_parallel*model.beta*period[0]);
	rhos[2] /= (4.0*model.t_parallel*model.beta*period[0]);
	
	
	
	return rhos;
	
}

vector<double> jackknife(const vector<double> &pool){
	
	vector<double> re(2, 0.0);
	vector<double> sample;
	vector<double> diff;
	int block = 10;
	double sum;
	
	sum = accumulate(pool.begin(), pool.end(), 0.0);
	re[0] = sum/double(pool.size());
	
	for (int i = 0; i < block; ++i){
		vector<double>::const_iterator begin, end;
		double partialsum = 0;
		begin = pool.begin() + (pool.end() - pool.begin())/block*i;
		end = pool.begin() + (pool.end() - pool.begin())/block*(i+1);
		for(vector<double>::const_iterator it = begin; it != end; ++it)
			partialsum += *it;
		sample.push_back((sum-partialsum)/(block-1)/(pool.size()/block));
		diff.push_back(sample[i] - re[0]);
		
	}
	re[1] = sqrt(inner_product(diff.begin(), diff.end(), diff.begin(), 0.0)*(block-1)/block);
	return re;
	
}


vector<double> binning(const vector<double> &pool){
	
	vector<double> re(2, 0.0);
	vector<double> sample;
	vector<double> diff;
	int block = 10;
	double sum;
	
	sum = accumulate(pool.begin(), pool.end(), 0.0);
	re[0] = sum/double(pool.size());
	
	for (int i = 0; i < block; ++i){
		vector<double>::const_iterator begin, end;
		double partialsum = 0;
		begin = pool.begin() + (pool.end() - pool.begin())/block*i;
		end = pool.begin() + (pool.end() - pool.begin())/block*(i+1);
		for(vector<double>::const_iterator it = begin; it != end; ++it)
			partialsum += *it;
		sample.push_back(partialsum/(pool.size()/block));
		diff.push_back(sample[i] - re[0]);
		
	}
	re[1] = sqrt(inner_product(diff.begin(), diff.end(), diff.begin(), 0.0)/block);
	return re;
	
}

vector<double> computek(const Model &model, const vector<double> &pool){
	
	vector<double> re(2, 0.0);	
	vector<double> kappa;
	int block = 10;
	double ave;	
	ave = accumulate(pool.begin(), pool.end(), 0.0)/pool.size();
	
	for (int i = 0; i < block; ++i){
		vector<double>::const_iterator begin, end;
		vector<double> sample(pool.size()/block);
		vector<double> diff;
		begin = pool.begin() + (pool.end() - pool.begin())/block*i;
		end = pool.begin() + (pool.end() - pool.begin())/block*(i+1);
		copy(begin, end, sample.begin());
		
		for (int j = 0; j < int(sample.size()); ++j){
			diff.push_back(sample[j]-ave);
		}
		kappa.push_back(model.beta*model.ns*inner_product(diff.begin(), diff.end(), diff.begin(), 0.0)/(pool.size()/block));
				
	}	
	re[0] = accumulate(kappa.begin(), kappa.end(), 0.0)/kappa.size();
	for (int i = 0; i < int(kappa.size()); ++i){
		kappa[i] = kappa[i] -re[0];
	}
	re[1] = sqrt(inner_product(kappa.begin(), kappa.end(), kappa.begin(), 0.0)/kappa.size());	
	return re;
	
}












