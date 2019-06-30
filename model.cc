#include "bilayer.h"

void readfile(Model &model,int argc, char **argv){
	
	if (argc != 2){
		cout << "No specifying input file!\n";
		cout << "Usage: " << argv[0] << " <filename>\n";
		// modify here should exit program 
		// insert code...
		abort();
	}
	else{
		ifstream ifile(argv[1]);
		if (ifile.is_open()){
			ifile >> model.lx >> model.ly;
			ifile >> model.beta >> model.mpps >> model.warmsteps >> model.mcsteps;
			ifile >> model.t_parallel >> model.t_perp >> model.interaction >> model.chempot;
			ifile >> model.seed;
			// insert code here display the model's parameter, 
			// Do I need to display the modle's parameter ??
		}
	}
	model.chempot += 1e-7;
	
}

void makelattice(Model &model){
	
	int lx = model.lx;
	int ly = model.ly;
	int ns = lx*ly*2;
	int nb = lx*ly*5;
	int s, x, y;
	model.bsites.resize(5*lx*ly, vector<int> (2, 0));	
	
	model.ns = ns;
	model.nb = nb;
	// for horizontal bonds
	for (int k = 0; k < 2; ++k){
		for (int j = 0; j < ly; ++j){
			for (int i = 0; i < lx; ++i){
				s = i + j*lx + k*lx*ly;
				model.bsites[s][0] = s; // index of the site on the left side of the bond
				x = (i+1)%lx;
				model.bsites[s][1] = x + j*lx + k*lx*ly; // index of the site on the right side of the bond
				model.bsites[s+ns][0] = s; // index of the site on the front side of the bond
				y = (j+1)%ly;
				model.bsites[s+ns][1] = i + y*lx + k*lx*ly; // index of the site on the back side of the bond				
			}
		}
	}
	// for vertical bonds
	for (int j = 0; j < lx; ++j){
		for (int i = 0; i < ly; ++i){
			s = i + j*lx;
			model.bsites[s+2*ns][0] = s; // index of the site on the down side of the bond
			model.bsites[s+2*ns][1] = i + j*lx + lx*ly; // index of the site on the upper side of the bond
		}
	}	
}

void calweight(Model & model){
	// import parameter from Model
	int mpps = model.mpps;
	double v = model.interaction;
	double t0 = model.t_parallel;
	double t1 = model.t_perp;
	double mu = model.chempot;
	
	double small = numeric_limits<double>::max();
	double z = 4.0;
	double epislon = 0.000001;
	double np, nq, t, c;
	int sub; // index of the weight
	
	model.weight.resize(2*3*(mpps+1)*(mpps+1), 0.0);
		
	for (int p = 0; p < mpps+1; ++p){
		np = p;
		for (int q = 0; q < mpps+1; ++q){
			nq = q;
			sub = p*(mpps+1) + q;
			model.weight[sub] = -v*np*nq + mu*(np+nq)/z - 0.5/z*(np*(np-1) + nq*(nq-1));
			small = (((small > model.weight[sub]) ? model.weight[sub] : small));
		} 
	}
	c = -small + epislon;
	model.c = c;
	for (int kind = 0; kind < 2; ++kind){
		for (int type = 0; type < 3; ++type){
			for (int p = 0; p < mpps+1; ++p){
				np = p;
				for (int q = 0; q < mpps+1; ++q){
					nq = q;
					sub = kind*3*(mpps+1)*(mpps+1) + type*(mpps+1)*(mpps+1) + p*(mpps+1) + q;
					t = ((kind == 0) ? t0 : t1);
					switch (type){
						case 0:
						model.weight[sub] += c;
						break;
						case 1:
						if (p != mpps && q != 0)
							model.weight[sub] = t*sqrt(nq*(np+1.0));
						else
							model.weight[sub] = 0.0;
						break;
						case 2:
						if (p != 0 && q != mpps)
							model.weight[sub] = t*sqrt(np*(nq+1.0));
						else
							model.weight[sub] = 0.0;
						break;
					}
					//cout << model.weight[sub] << endl;
				}
			}
		}
	}
}

void askspace(const Model & model, Configuration & config){
	
	boost::random::mt19937 gen(model.seed);
	boost::random::uniform_int_distribution<> dist(0, model.mpps);
	config.msl = 20;
	config.nh = 0;
	
	for (int i = 0; i < model.ns; ++i)
		config.state.push_back(dist(gen));
	
	config.opstring.assign(config.msl, 0);
	config.leg.resize(config.msl, vector<int> (4, -1));
	config.vertexlist.resize(4*config.msl);
	config.frststateop.resize(model.ns);
	config.laststateop.resize(model.ns);
}

void diagupdate(const Model & model, Configuration & config){
	
	static boost::random::mt19937 gen(model.seed);
	boost::random::uniform_int_distribution<> dist(0, model.nb - 1);
	boost::random::uniform_real_distribution<double> rnd(0, 1);
	
	int op, b, np, nq, kind, sub;
	
	for (int i = 0; i < config.msl; ++i){
		
		op = config.opstring[i];
		
		if (op == 0){
			// add an diagonal operator on a random bond
			b = dist(gen);
			np = config.state[model.bsites[b][0]];
			nq = config.state[model.bsites[b][1]];
			kind = (( b < 2*model.ns) ? 0 : 1);
			sub = kind*3*(model.mpps+1)*(model.mpps+1) + np*(model.mpps+1) + nq;
			if ( rnd(gen)*(config.msl - config.nh) <= model.nb*model.beta*model.weight[sub] ){
				config.opstring[i] = 2*b + 2;
				config.nh++;
				config.leg[i][0] = np; config.leg[i][1] = nq;
				config.leg[i][2] = np; config.leg[i][3] = nq;
			}
		}
		else if (op%2 == 0){
			// decrease an operator
			b = op/2 - 1;
			kind = (( b < 2*model.ns) ? 0 : 1);
			np = config.state[model.bsites[b][0]];
			nq = config.state[model.bsites[b][1]];
			sub = kind*3*(model.mpps+1)*(model.mpps+1) + np*(model.mpps+1) + nq;
			if ( rnd(gen)*model.weight[sub]*model.nb*model.beta <= double(config.msl - config.nh +1.0)){
				config.opstring[i] = 0;
				config.nh--;
				config.leg[i][0] = -1; config.leg[i][1] = -1;
				config.leg[i][2] = -1; config.leg[i][3] = -1;
			}
		}
		else{
			// an off-diagonal operator to propagate the state
			b = (op - 3)/2;
			config.state[model.bsites[b][0]] = config.leg[i][2];
			config.state[model.bsites[b][1]] = config.leg[i][3];
		}		
	}	
}

void adjustcutoff(const Model & model, Configuration &config){
	
	// length of opstirng, leg, vertexlist need to be modified
	int new_msl = config.nh + config.nh/3;
	if (new_msl > config.msl){
		config.opstring.resize(new_msl);
		config.vertexlist.resize(4*new_msl);
		config.leg.resize(new_msl, vector<int> (4, -1));
		config.msl = new_msl;
	}	
}

void linkvertices(const Model & model, Configuration & config){
	
	int op, b, s1, s2, v1, v2;
	
	fill(config.frststateop.begin(), config.frststateop.end(), -1);
	fill(config.laststateop.begin(), config.laststateop.end(), -1);
	
    for (int v0 = 0; v0 < 4*config.msl; v0 = v0 + 4){
        op = config.opstring[v0/4];
        if (op != 0){
            b = op/2 - 1;
            s1 = model.bsites[b][0];
            s2 = model.bsites[b][1];
            v1 = config.laststateop[s1];
            v2 = config.laststateop[s2];
            if (v1 != -1){
                config.vertexlist[v1] = v0;
                config.vertexlist[v0] = v1;
            }
            else{
                config.frststateop[s1] = v0;
            }
            if (v2 != -1){
                config.vertexlist[v2] = v0 + 1;
                config.vertexlist[v0 + 1] = v2;
            }
            else{
                config.frststateop[s2] = v0 + 1;
            }
            config.laststateop[s1] = v0 + 2;
            config.laststateop[s2] = v0 + 3;
        }
        else{
            config.vertexlist[v0] = -1;
            config.vertexlist[v0 + 1] = -1;
            config.vertexlist[v0 + 2] = -1;
            config.vertexlist[v0 + 3] = -1;
        }
    }
    
    for (int s1 = 0; s1 < model.ns; ++s1){
        v1 = config.frststateop[s1];
        if (v1 != -1){
            v2 = config.laststateop[s1];
            config.vertexlist[v2] = v1;
            config.vertexlist[v1] = v2;
        }
    }	
}

// crement == 0 increase a particle, crement == 1 decrese a particle 
double bounce(const Model & model, int kind, int type, int np, int nq, int crement){
	int sub = kind*3*(model.mpps+1)*(model.mpps+1) + type*(model.mpps+1)*(model.mpps+1) + np*(model.mpps+1) + nq;
	return model.weight[sub];
}

// modify code here, throw exception probably, overflow of array in this case
double reverse(const Model &model, int kind, int type, int np, int nq, int crement){
	if (crement == 0){
		++np; --nq;
		if (type == 0)
			type = 2;
		else if (type == 1)
			type = 0;
		else
			type = 3;
	}
	else{
		--np, ++nq;
		if (type == 0)
			type = 1;
		else if (type == 1)
			type = 3;
		else
			type = 0;
	}
	int sub = kind*3*(model.mpps+1)*(model.mpps+1) + type*(model.mpps+1)*(model.mpps+1) + np*(model.mpps+1) + nq;
	if (np < 0 || nq < 0 || np > model.mpps || nq > model.mpps || type >2)
		return 0.0;
	else
		return model.weight[sub];	
}

double straight(const Model &model, int kind, int type, int np, int nq, int crement){
	(crement == 0) ? ++np : -- np;
	int sub = kind*3*(model.mpps+1)*(model.mpps+1) + type*(model.mpps+1)*(model.mpps+1) + np*(model.mpps+1) + nq;
	if (np < 0 || nq < 0 || np > model.mpps || nq > model.mpps || type >2)
		return 0.0;
	else
		return model.weight[sub];
} 

double jump(const Model &model, int kind, int type, int np, int nq, int crement){
	if (crement == 0){
		++np;
		if (type == 0)
			type = 2;
		else if (type == 1)
			type = 0;
		else
			type = 3;
	}
	else{
		--np;
		if (type == 0)
			type = 1;
		else if (type == 1)
			type = 3;
		else
			type = 0;
	}
	int sub = kind*3*(model.mpps+1)*(model.mpps+1) + type*(model.mpps+1)*(model.mpps+1) + np*(model.mpps+1) + nq;
	if (np < 0 || nq < 0 || np > model.mpps || nq > model.mpps || type >2)
		return 0.0;
	else
		return model.weight[sub];	
}

void caltable(Model & model){
	
	const int Kind = 2;
	const int Type = 3;
	const int Mpps = model.mpps + 1;
	const int Crement = 2;
	const int Exitleg = 4;
	
	//ofstream ofile("table.dat");
	
	model.table.resize(Kind*Crement*Type*Mpps*Mpps*Exitleg, 0.0);
	
	for (int kind = 0; kind < 2; ++kind){
		for (int crement = 0; crement < Crement; ++crement){
			for (int type = 0; type < Type; ++type){
				for (int np = 0; np < Mpps; ++np){
					for (int nq = 0; nq < Mpps; ++nq){
						int sub = kind*Crement*Type*Mpps*Mpps*Exitleg + crement*Type*Mpps*Mpps*Exitleg +
							type*Mpps*Mpps*Exitleg + np*Mpps*Exitleg + nq*Exitleg;
						vector<double> w;
						dMatrix a;
						w.push_back(bounce(model, kind, type, np, nq, crement));
						w.push_back(reverse(model, kind, type, np, nq, crement));
						w.push_back(straight(model, kind, type, np, nq, crement));
						w.push_back(jump(model, kind, type, np, nq, crement));
						
						if (accumulate(w.begin(), w.end(), 0.0) == 0){
							model.table[sub] = 0;
							model.table[sub+1] = 0;
							model.table[sub+2] = 0;
							model.table[sub+3] = 0;
						}
						else{
							a = solution(w);
							model.table[sub] = a[0][0];
							model.table[sub+1] = a[0][1];
							model.table[sub+2] = a[0][2];
							model.table[sub+3] = a[0][3];
							
						}
						//ofile << model.table[sub] << '\t' << model.table[sub+1] << '\t'  
						//	<< model.table[sub+2] << '\t' << model.table[sub+3] << '\n';
					}
				}
			}
		}
	}
	
}



















