#include "bilayer.h"

//      type 0           type 1             type 2
//    np      nq      np+1    nq-1       np-1    nq+1  
//    |--------|       |--------|         |--------|
//    |--------|       |--------|         |--------|
//    np      nq      np       nq         np      nq
//    

int determinetype(const int x, const int y){
	
	if (x == y)
		return 0;
	else if (x == y - 1)
		return 1;
	else if (x == y + 1)
		return 2;
	else{
		cout << "Error :: unexpected operator type!" << endl;
		abort();
	}
	
} 

// (0,0)-> 0; (0,1)-> 1; (0,2)-> 2; (0,3)-> 3;
// (1,0)-> 1; (1,1)-> 0; (1,2)-> 3; (1,3)-> 2;
// (2,0)-> 2; (2,1)-> 3; (2,2)-> 0; (2,3)-> 1;
// (3,0)-> 3; (3,1)-> 2; (3,2)-> 1; (3,3)-> 0;
int rotation(int x, int y){
	
	return (x+int(pow(-1, x%2))*y+4)%4;
	
}

int loopupdate(const Model &model, Configuration &config){
	
	const int Type = 3;
	const int Mpps = model.mpps + 1;
	const int Crement = 2;
	const int Exitleg = 4;
	
	int v0, v1, v2, crement, p, entraleg, exitleg, op, b, a, np, nq, sub, kind;
	int type, indexcre, y, z, looplength = 0;
	double w[4], rndnum;
	
	static boost::random::mt19937 gen(model.seed);
	boost::random::uniform_real_distribution<double> rnd(0, 1);
	
	v0 = rnd(gen)*4*config.msl;
	v1 = v0;
	
	if (config.vertexlist[v0] != -1){
		
		crement = ( rnd(gen) < 0.5 ? -1 : 1);
		
		do{
			
			p = v1/4;
			entraleg = v1%4;
			op = config.opstring[p];
			b = op/2 - 1;
			a = op%2 + 1;
			
			// determine the type and give value to np and nq
			y = (entraleg+2)%4;
			z = ( (entraleg%2) ? entraleg-1: entraleg+1);
			np = config.leg[p][entraleg];
			nq = config.leg[p][z];
			type = determinetype(np, config.leg[p][y]);
			
			indexcre = ( (crement == 1) ? 0 : 1 );
			kind = ( (b < model.ns*2) ? 0 : 1 );
			
			sub = kind*Crement*Type*Mpps*Mpps*Exitleg + indexcre*Type*Mpps*Mpps*Exitleg +
							type*Mpps*Mpps*Exitleg + np*Mpps*Exitleg + nq*Exitleg;
			
			for (int i = 0; i < 4; ++i)
				w[i] = model.table[sub+i];
			
			rndnum = rnd(gen);
			if (rndnum < w[0])
				exitleg = 0;
			else if (rndnum < w[0]+w[1])
				exitleg = 1;
			else if (rndnum < w[0]+w[1]+w[2])
				exitleg = 2;
			else if (rndnum < w[0]+w[1]+w[2]+w[3])
				exitleg = 3;
			else{
				cout << "Warning:: Can't find out the exitleg" << endl;
				abort();
			}
			
			exitleg = rotation(entraleg, exitleg);
			
			if (exitleg != entraleg)
				looplength++;
			// can be rewritten and make it short
			// definitely needed to be revised 
			switch (entraleg){
				case 0:
				switch (exitleg){
					case 0:
					crement = -crement;
					break;
					case 1:
					config.leg[p][0] += crement;
					config.leg[p][1] -= crement;
					crement = -crement;
					a = !(a-1) + 1;
					break;
					case 2:
					config.leg[p][0] += crement;
					config.leg[p][2] += crement;
					break;
					case 3:
					config.leg[p][0] += crement;
					config.leg[p][3] += crement;
					a = !(a-1) + 1;
					break;
				}
				break;
				
				case 1:
				switch (exitleg){
					case 0:
					config.leg[p][1] += crement;
					config.leg[p][0] -= crement;
					crement = -crement;
					a = !(a-1) + 1;
					break;
					case 1:
					crement = -crement;
					break;
					case 2:
					config.leg[p][1] += crement;
					config.leg[p][2] += crement;
					a = !(a-1) + 1;
					break;
					case 3:
					config.leg[p][1] += crement;
					config.leg[p][3] += crement;
					break;
				}
				break;
				
				case 2:
				switch (exitleg){
					case 0:
					config.leg[p][2] += crement;
					config.leg[p][0] += crement;
					break;
					case 1:
					config.leg[p][2] += crement;
					config.leg[p][1] += crement;
					a = !(a-1) + 1;
					break;
					case 2:
					crement = -crement;
					break;
					case 3:
					config.leg[p][2] += crement;
					config.leg[p][3] -= crement;
					crement = -crement;
					a = !(a-1) + 1;
					break;
				}
				break;
				
				case 3:
				switch (exitleg){
					case 0:
					config.leg[p][3] += crement;
					config.leg[p][0] += crement;
					a = !(a-1) + 1;
					break;
					case 1:
					config.leg[p][3] += crement;
					config.leg[p][1] += crement;
					break;
					case 2:
					config.leg[p][3] += crement;
					config.leg[p][2] -= crement;
					crement = -crement;
					a = !(a-1) + 1;
					break;
					case 3:
					crement = -crement;
					break;
				}
				break;
			}
			
			config.opstring[p] = 2*b + a + 1;
			v2 = 4*p + exitleg;
			if (v2 == v0)
				break;
			v1 = config.vertexlist[v2];
			
		}while(v0 != v1);
		
	}
	return looplength;
}

void updatestate(const Model & model, Configuration & config){
	
	static boost::random::mt19937 gen(model.seed);
	boost::random::uniform_int_distribution<> dist(0, model.mpps);
	
	int op, b, s1, s2, v1, v2;
	fill(config.frststateop.begin(), config.frststateop.end(), -1);
	fill(config.laststateop.begin(), config.laststateop.end(), -1);
	fill(config.state.begin(), config.state.end(), model.mpps + 1);
	
	for (int i = 0; i < config.msl; ++i){
		op = config.opstring[i];
		if (op != 0){
			b = op/2 - 1;
			s1 = model.bsites[b][0];
			s2 = model.bsites[b][1];
			v1 = config.frststateop[s1];
			v2 = config.laststateop[s2];
			if (v1 == -1){
				config.state[s1] = config.leg[i][0];
				config.frststateop[s1] = 1;
			}
			if (v2 == -1){
				config.state[s2] = config.leg[i][1];
				config.laststateop[s2] = 1;
			}
		}
	}
	for (int i = 0; i < model.ns; ++i){
		if (config.state[i] == model.mpps + 1)
			config.state[i] = dist(gen);
	}
}
















