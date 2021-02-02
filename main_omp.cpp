//
//  map_align
//

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <iterator>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <omp.h>

using namespace std;
typedef vector<char> char_1D;
typedef vector<bool> bool_1D;
typedef vector<bool_1D> bool_2D;
typedef vector<int> int_1D;
typedef vector<int_1D> int_2D;
typedef vector<double> double_1D;
typedef vector<double_1D> double_2D;
typedef vector<double_2D> double_3D;
typedef vector<string> string_1D;
typedef vector<string_1D> string_2D;

class Opt{
public:
	string file_a;
	string file_b;
	string file_list_a;
	string file_list_b;
	string file_out;
	bool_1D range_a;
	bool_1D range_b;
	bool use_gap_ss = false;
	double gap_ss_w = 2;
	bool use_prf = false;
	double prf_w = 1;
	double gap_open = -1;
	double gap_ext = -0.01;
	double gap_ext_w;
	int sep_cutoff = 3;
	double sco_cutoff = 0.0;
	int iter = 20;
	bool silent = false;
	bool get_tot = true;
	int cpu = omp_get_max_threads();
	void get(string_1D &arg);
};

class Map{
public:
	string name;
	int_1D m2n;
	int_1D n2m;
	double_2D mtx;
	int_1D vec;
	int_1D vec_div;
	int_2D vec_i;
	double_2D prf;
	char_1D aa;
	char_1D ss;
	double_1D gap;
	int size;
	void load(string file, int sep_cutoff, double sco_cutoff, bool_1D range);
	void mod_gap(double gap_ss_w);
};
typedef vector<Map> Map_1D;

inline bool exists (const std::string& name) {ifstream f(name);return f.good();}

double gaussian(double mean, double stdev, double x){return exp(-pow((x - mean),2)/(2*(pow(stdev,2))));}
int_1D align(double_1D &gap_a, double_1D &gap_b, double &gap_e, double_2D &sco_mtx, double_2D &p_sco_mtx);
double Falign(double *sco_mtx, int rows, int cols);


void do_a_vs_b(Map_1D &maps_a, Map_1D &maps_b, Opt opt);
string do_it(Map &map_a, Map &map_b, Opt opt);

double sepw(double sep){if(sep <= 4){return 0.50;}else if(sep == 5){return 0.75;}else{return 1.00;}}
void ini_SCO(double sep_x, double sep_y, double_2D &SCO, Map &map_a, Map &map_b);
void ini_prf_SCO(double_2D &P_SCO, double &prf_w, Map &map_a,Map &map_b);
int_1D mod_SCO(int iter, double gap_e, double_2D &SCO, double_2D &P_SCO, Map &map_a, Map &map_b);
void chk(double &gap_e_w,double &con_sco,double &gap_sco,double &prf_sco,int_1D &a2b,double_2D &P_SCO,Map &map_a,Map &map_b);
Map load_map(string file, bool_1D range, Opt opt);
void load_maps(Map_1D &maps, string file_list, Opt opt);
void get_tot(double_1D &tot, double_1D &tot_aln, int &aln_len, Map &map_a, Map &map_b, int_1D &a2b_max);


int main(int argc, const char * argv[]){
	// parse input arguments
	Opt opt; string_1D arg(argv+1,argv+argc); opt.get(arg);
	opt.gap_ext_w = fabs(opt.gap_ext)/fabs(opt.gap_open);

	// set number of threads
	omp_set_num_threads(opt.cpu);

	// load maps
	Map_1D maps_a;
	Map_1D maps_b;
	if(!opt.file_a.empty()){maps_a.push_back(load_map(opt.file_a,opt.range_a,opt));}
	if(!opt.file_b.empty()){maps_b.push_back(load_map(opt.file_b,opt.range_b,opt));}
	if(!opt.file_list_a.empty()){load_maps(maps_a,opt.file_list_a,opt);}
	if(!opt.file_list_b.empty()){load_maps(maps_b,opt.file_list_b,opt);}

	if(maps_a.size() == 0){cout << "ERROR: no 'a' map(s) loaded" << endl; exit(1);}
	if(maps_b.size() == 0){cout << "ERROR: no 'b' map(s) loaded" << endl; exit(1);}

	if(maps_a.size() > maps_b.size()){do_a_vs_b(maps_b,maps_a,opt);}
	else{do_a_vs_b(maps_a,maps_b,opt);}

	return 0;
}

void do_a_vs_b(Map_1D &maps_a, Map_1D &maps_b, Opt opt){
	// get results
	string_2D logs(maps_a.size(),string_1D(maps_b.size()));
	for(int a = 0; a < maps_a.size(); a++){
		#pragma omp parallel for
		for(int b = 0; b < maps_b.size(); b++){
			logs[a][b] = do_it(maps_a[a],maps_b[b],opt);
		}
	}
	// save results
	ofstream out(opt.file_out);
	if(out.is_open()){
		for(int a = 0; a < maps_a.size(); a++){
			for(int b = 0; b < maps_b.size(); b++){out << logs[a][b];}
		}
		out.close();
	}
}

string do_it(Map &map_a, Map &map_b, Opt opt){
	stringstream out;

	// if use_prf on, initialize profile SCO matrix
	double_2D P_SCO; if(opt.use_prf){ini_prf_SCO(P_SCO,opt.prf_w,map_a,map_b);}

	// STARTING ALIGNMENT!!!
	// keeping track of the BEST alignment
	int max_sep_x = 0;
	int max_sep_y = 0;
	int max_g_e = 0;
	double con_max = -1;
	double gap_max = 0;
	double prf_max = 0;
	int_1D a2b_max;

	// try different sep (sequence seperation difference) penalities
	double_1D sep_x_steps {0,1,2}; // (constant, linear, quadratic)
	for(int sx = 0; sx < sep_x_steps.size(); sx++){
		double sep_x = sep_x_steps[sx];

		//try different scaling factors for sep penalities
		double_1D sep_y_steps {1,2,4,8,16,32};
		for(int sy = 0; sy < sep_y_steps.size(); sy++){
			double sep_y = sep_y_steps[sy];

			// Get initial score matrix
			double_2D C_SCO(map_a.size,double_1D(map_b.size,0));
			ini_SCO(sep_x,sep_y,C_SCO,map_a,map_b);
			// try different gap_ext penalities!
			double_1D gap_e_steps {5,10,100,1000};
			for(int g_e = 0; g_e < gap_e_steps.size(); g_e++){
				double gap_e = 1.0/gap_e_steps[g_e];

				// restart SCO matrix
				double_2D SCO = C_SCO;

				// get alignment (a2b mapping) after X iterations
				int_1D a2b = mod_SCO(opt.iter,gap_e,SCO,P_SCO,map_a,map_b);

				// compute number of contacts/gaps made
				double con_sco = 0;
				double gap_sco = 0;
				double prf_sco = 0;
				chk(opt.gap_ext_w,con_sco,gap_sco,prf_sco,a2b,P_SCO,map_a,map_b);

				// print info
				/*
				 if(!opt.silent){
					out
					<< "TMP " << sep_x
					<< "_" << sep_y
					<< "_" << g_e
					<< " a: " << map_a.name
					<< " b: " << map_b.name
					<< " con: " << con_sco
					<< " gap: " << gap_sco;
					if(opt.use_prf){out << " prf: " << prf_sco;}
					out << endl;
				 }
				 */
				// save if BEST!
				if(con_sco+gap_sco+prf_sco > con_max+gap_max+prf_max){
					max_sep_x = sep_x;
					max_sep_y = sep_y;
					max_g_e = g_e;
					con_max = con_sco;
					gap_max = gap_sco;
					prf_max = prf_sco;
					a2b_max = a2b;
				}
			}
		}
	}
	// Report the BEST score
	out
	<< "best_params: " << max_sep_x << "_" << max_sep_y << "_" << max_g_e
	<< " id_a: " << map_a.name
	<< " id_b: " << map_b.name
	<< " len_a: " << map_a.vec.size()
	<< " len_b: " << map_b.vec.size()
	<< " con: " << con_max
	<< " gap: " << gap_max;
	if(opt.use_prf){out << " prf: " << prf_max;}
	if(opt.get_tot){
		double_1D tot(2,0.0);
		double_1D tot_aln(2,0.0);
		int aln_len = 0;
		get_tot(tot,tot_aln,aln_len,map_a,map_b,a2b_max);
		out
		<< " aln_len: " << aln_len
		<< " tot_a: " << tot[0]
		<< " tot_b: " << tot[1]
		<< " tot_aln_a: " << tot_aln[0]
		<< " tot_aln_b: " << tot_aln[1];
	}
	for(int a = 0; a < map_a.size; a++){
		int b = a2b_max[a];
		if(b != -1){
			out << " " << map_a.m2n[a] << ":" << map_b.m2n[b];
		}
	}
	out << endl;
	return(out.str());
}
int_1D align(double_1D &gap_a, double_1D &gap_b, double &gap_e, double_2D &sco_mtx, double_2D &p_sco_mtx){
	// LOCAL_ALIGN
	// Start	0
	// [A]lign	1
	// [D]own	2
	// [R]ight	3

	double max_sco = 0;
	int rows = sco_mtx.size();
	int cols = sco_mtx[0].size();

	bool add_prf = false;if(p_sco_mtx.size() == rows){add_prf = true;}

	int_1D a2b(rows,-1);

	double_2D sco(rows+1,double_1D(cols+1,0));
	int_2D label(rows+1,int_1D(cols+1,0));

	int max_i = 0;int max_j = 0;
	for (int i = 1; i <= rows; i++){
		for (int j = 1; j <= cols; j++){
			double A = sco[i-1][j-1] + sco_mtx[i-1][j-1]; if(add_prf == true){A += p_sco_mtx[i-1][j-1];}
			double D = sco[i-1][j];
			double R = sco[i][j-1];

			if(label[i-1][j] == 1){D += gap_b[j-1];}else{D += gap_b[j-1] * gap_e;}
			if(label[i][j-1] == 1){R += gap_a[i-1];}else{R += gap_a[i-1] * gap_e;}

			if(A <= 0 and D <= 0 and R <= 0){label[i][j] = 0;sco[i][j] = 0;}
			else{
				if(A >= R){if(A >= D){label[i][j] = 1;sco[i][j] = A;}else{label[i][j] = 2;sco[i][j] = D;}}
				else{if(R >= D){label[i][j] = 3;sco[i][j] = R;}else{label[i][j] = 2;sco[i][j] = D;}}
				if(sco[i][j] > max_sco){max_i = i;max_j = j;max_sco = sco[i][j];}
			}
		}
	}
	int i = max_i;int j = max_j;
	while(1){
		if(label[i][j] == 0){break;}
		else if(label[i][j] == 1){a2b[i-1] = j-1;i--;j--;}
		else if(label[i][j] == 2){i--;}
		else if(label[i][j] == 3){j--;}
	}
	return(a2b);
}
double Falign(double *sco_mtx, int rows, int cols){
	double max_sco = 0;
	double sco[rows+1][cols+1]; memset(sco, 0, sizeof(sco));
	for (int i = 1; i <= rows; i++){
		for (int j = 1; j <= cols; j++){
			double A = sco[i-1][j-1] + sco_mtx[(i-1)*cols+(j-1)];
			double D = sco[i-1][j];
			double R = sco[i][j-1];

			if(A >= R){if(A >= D){sco[i][j] = A;}else{sco[i][j] = D;}}
			else{if(R >= D){sco[i][j] = R;}else{sco[i][j] = D;}}

			if(sco[i][j] > max_sco){max_sco = sco[i][j];}
		}
	}
	return(max_sco);
}
void Map::mod_gap(double gap_ss_w){
	for(int i = 0; i < ss.size()-1; i++){
		if((ss[i] == 'H' && ss[i+1] == 'H') || (ss[i] == 'E' && ss[i+1] == 'E')){gap[i] *= gap_ss_w;}
	}
}
// INITIATE SCORE MATRIX: function for populating the initial similarity matrix
void ini_SCO(double sep_x, double sep_y, double_2D &SCO, Map &map_a, Map &map_b){
	// Get initial score matrix
	for(int i=0; i < map_a.vec.size(); i++){ // go through columns (map_a.vec) in map_a that has contacts
		int ai = map_a.vec[i];
		for(int j=0; j < map_b.vec.size(); j++){ // go through columns (map_b.vec) in map_b that has contacts
			int bi = map_b.vec[j];
			int A[2] = {(int)map_a.vec_div[ai],(int)(map_a.vec_i[ai].size()-map_a.vec_div[ai])};
			int B[2] = {(int)map_b.vec_div[bi],(int)(map_b.vec_i[bi].size()-map_b.vec_div[bi])};
			for(int k=0; k <= 1; k++){ // left and right of diagonal
				if(A[k] > 0 and B[k] > 0){
					double M[A[k]*B[k]];
					for(int n=0; n < A[k]; n++){
						int nn = n; if(k == 1){nn += map_a.vec_div[ai];}
						int aj = map_a.vec_i[ai][nn];
						int sep_a = abs(ai-aj);
						for(int m=0; m < B[k]; m++){
							int mm = m; if(k == 1){mm += map_b.vec_div[bi];}
							int bj = map_b.vec_i[bi][mm];
							int sep_b = abs(bi-bj);
							int sep_D = abs(sep_a-sep_b);
							double sep_M = min(sep_a,sep_b);
							double sep_std = sep_y*(1+pow(sep_M-2,sep_x));
							if(sep_D/sep_std < 6){
								M[n*B[k]+m] = map_a.mtx[ai][aj] * map_b.mtx[bi][bj] * sepw(sep_M) * gaussian(0,sep_std,sep_D);
							}else{M[n*B[k]+m] = 0;}
						}
					}
					SCO[ai][bi] += Falign(M,A[k],B[k]);
				}
			}
		}
	}
}

// MODIFY SCORE MATRIX: function for modifying the initial similarity matrix
int_1D mod_SCO(int iter, double gap_e, double_2D &SCO, double_2D &P_SCO, Map &map_a, Map &map_b){
	// iterate
	int_1D a2b_tmp;
	for(int it=0; it < iter; it++){
		// align
		a2b_tmp = align(map_a.gap,map_b.gap,gap_e,SCO,P_SCO);

		// update similarity matrix
		double IT = (double)it + 1;
		double s1 = (IT/(IT+1)); double s2 = (1/(IT+1));
		for(int a=0; a < map_a.vec.size(); a++){ // go through columns (map_a.vec) in map_a that has contacts
			int ai = map_a.vec[a];
			for(int b=0; b < map_b.vec.size(); b++){ // go through columns (map_b.vec) in map_b that has contacts
				int bi = map_b.vec[b];
				double sco_contact = 0;
				for(int n=0; n < map_a.vec_i[ai].size(); n++){ // go through contacts in map_a.vec
					int aj = map_a.vec_i[ai][n];
					int bj = a2b_tmp[aj]; // get mapping
					if(bj != -1){ // if mapping exists
						if((ai > aj and bi > bj) or (ai < aj and bi < bj)){ // if ai-aj in same direction as bi-bj
							double sep_M = min(abs(ai-aj),abs(bi-bj));
							sco_contact += map_a.mtx[ai][aj] * map_b.mtx[bi][bj] * sepw(sep_M);
						}
					}
				}
				SCO[ai][bi] = s1*SCO[ai][bi] + s2*sco_contact;
			}
		}
	}
	return(a2b_tmp);
}
// CHK: compute number of contacts/gaps made
void chk(
		 double &gap_e_w,
		 double &con_sco,
		 double &gap_sco,
		 double &prf_sco,
		 int_1D &a2b,
		 double_2D &P_SCO,
		 Map &map_a,
		 Map &map_b){

	bool use_prf = false; if(P_SCO.size() == map_a.size){use_prf = true;}

	int a = 0;int b = 0;
	for(int ai = 0; ai < map_a.size; ai++){
		int bi = a2b[ai];
		if(bi != -1){
			if(use_prf){prf_sco += P_SCO[ai][bi];}
			if(a > 0){ // compute number of gaps
				double num_gap_a = ((ai-a)-1); if(num_gap_a > 0){gap_sco += map_a.gap[ai] + map_a.gap[ai] * gap_e_w * (num_gap_a-1);}
				double num_gap_b = ((bi-b)-1); if(num_gap_b > 0){gap_sco += map_b.gap[bi] + map_b.gap[bi] * gap_e_w * (num_gap_b-1);}
			}
			for(int m=0; m < map_a.vec_div[ai]; m++){ // compute number of contacts
				int aj = map_a.vec_i[ai][m];
				int bj = a2b[aj];
				if(bj != -1){
					double sep_M = min(abs(ai-aj),abs(bi-bj));
					con_sco += map_a.mtx[ai][aj] * map_b.mtx[bi][bj] * sepw(sep_M);
				}
			}
			a = ai;b = bi;
		}
	}
	gap_sco /= 2;
}

// compute profile similarity matrix
void ini_prf_SCO(double_2D &P_SCO, double &prf_w, Map &map_a,Map &map_b){
	int size_a = map_a.prf.size();
	int size_b = map_b.prf.size();

	P_SCO.resize(size_a,double_1D(size_b,0));

	// compute background frequencies
	double_1D pb(20,0);
	int prf_size = map_a.prf[0].size();
	double pb_size = 0;
	for(int ai = 0; ai < size_a; ai++){
		if(map_a.aa[ai] != 'X'){ // ignore positions that have no identity
			for(int p=0; p < prf_size; p++){pb[p] += map_a.prf[ai][p];}
			pb_size += 1;
		}
	}
	for(int bi = 0; bi < size_b; bi++){
		if(map_b.aa[bi] != 'X'){ // ignore positions that have no identity
			for(int p=0; p < prf_size; p++){pb[p] += map_b.prf[bi][p];}
			pb_size += 1;
		}
	}
	for (int i=0; i < size_a; i++){
		for (int j=0; j < size_b; j++){
			if(map_a.aa[i] == 'X' or map_b.aa[j] == 'X'){P_SCO[i][j] = 0;} // if no identity, return score of 0
			else{
				// profile comparison calculation, similar to HHsuite from Soeding.
				double tmp_sco = 0;
				for(int p=0; p < prf_size; p++){tmp_sco += (map_a.prf[i][p]*map_b.prf[j][p])/(pb[p]/pb_size);}
				P_SCO[i][j] = log2(tmp_sco)/5 * prf_w;
			}
		}
	}
}
// compute the expected score over the full and aligned regions
void get_tot(double_1D &tot, double_1D &tot_aln, int &aln_len, Map &map_a, Map &map_b, int_1D &a2b_max){
	int_1D b2a_max(map_b.size,-1);
	for(int ai = 0; ai < map_a.size; ai++){
		int bi = a2b_max[ai];
		if(bi != -1){
			b2a_max[bi] = ai;
			aln_len++;
		}
		for(int m=0; m < map_a.vec_div[ai]; m++){
			int aj = map_a.vec_i[ai][m];
			int bj = a2b_max[aj];
			if(bi != -1 and bj != -1){
				int sep = min(abs(ai-aj),abs(bi-bj));
				double con = pow(map_a.mtx[ai][aj],2) * sepw(sep);
				tot[0] += con;
				tot_aln[0] += con;
			}
			else{
				int sep = abs(ai-aj);
				double con = pow(map_a.mtx[ai][aj],2) * sepw(sep);
				tot[0] += con;
			}
		}
	}
	for(int bi = 0; bi < map_b.size; bi++){
		int ai = b2a_max[bi];
		for(int m=0; m < map_b.vec_div[bi]; m++){
			int bj = map_b.vec_i[bi][m];
			int aj = b2a_max[bj];
			if(ai != -1 and aj != -1){
				int sep = min(abs(bi-bj),abs(ai-aj));
				double con = pow(map_b.mtx[bi][bj],2) * sepw(sep);
				tot[1] += con;
				tot_aln[1] += con;
			}
			else{
				int sep = abs(bi-bj);
				double con = pow(map_b.mtx[bi][bj],2) * sepw(sep);
				tot[1] += con;
			}
		}
	}
}
Map load_map(string file, bool_1D range, Opt opt){
	Map map; map.load(file,opt.sep_cutoff,opt.sco_cutoff,range);
	map.gap.resize(map.size,opt.gap_open);if(opt.use_gap_ss){map.mod_gap(opt.gap_ss_w);}
	return(map);
}

void load_maps(Map_1D &maps, string file_list, Opt opt){
	string line;
	ifstream in(file_list);
	while(getline(in,line)){
		istringstream is(line);
		string file;
		is >> file;
		if(file[0] != '#'){
			if(exists(file) == 1){
				if(!opt.silent){cout << "Adding: " << file << endl;}
				bool_1D null;
				maps.push_back(load_map(file,null,opt));
			}else if(!opt.silent){
				cout << "Error: '" << file << "' not found!" << endl;
			}
		}
	}
	in.close();
}
void Map::load(string file, int sep_cutoff, double sco_cutoff, bool_1D range){
	string line;
	ifstream in(file);
	while(getline(in,line)){
		istringstream is(line);
		string label;
		is >> label;
		if(label == "LEN" or label == "SIZE"){
			int size; is >> size;

			if(range.size() > 0){range.resize(size,0);} // if range previously defined fill rest with "0"
			else{range.resize(size,1);} // else set all to "1"

			int m = 0;n2m.resize(size,-1);
			for(int n = 0; n < size; n++){
				if(range[n] == 1){
					n2m[n] = m;
					m2n.push_back(n);
					m++;
				}
			}

			mtx.resize(m,vector<double>(m,0));
			prf.resize(m,vector<double>(20,0));
			aa.resize(m,'X');
			ss.resize(m,'X');
		}
		else if(label == "CON"){
			int i, j; is >> i >> j;
			if(range[i] == 1 and range[j] == 1){
				double sco;
				if(abs(j-i) >= sep_cutoff){
					if(is >> sco){}else{sco = 1;}
					if(sco >= sco_cutoff){
						mtx[n2m[i]][n2m[j]] = sco;
						mtx[n2m[j]][n2m[i]] = sco;
					}
				}
			}
		}
		else if(label == "PRF"){
			int i; is >> i;
			if(range[i] == 1){
				char tmp;
				is >> tmp; aa[n2m[i]] = tmp;
				is >> tmp; ss[n2m[i]] = tmp;
				double val;
				int j = 0;
				while(is >> val){prf[n2m[i]][j] = val;j++;}
			}
		}
	}
	in.close();
	for(int i=0; i < mtx.size(); i++){
		vec_i.push_back(vector<int>());
		for(int j=0; j < mtx.size(); j++){
			if(i == j){
				if(vec_i[i].empty()){vec_div.push_back(0);}
				else{vec_div.push_back(vec_i[i].size());}
			}
			if(mtx[i][j] > 0){
				vec_i[i].push_back(j);
			}
		}
		if(vec_i[i].size() > 0){vec.push_back(i);}
	}
	size = mtx.size();
	name = file;
}
void Opt::get(string_1D &opt){
	for (int a = 0; a < opt.size(); a++){
		string arg = opt[a];
		if (arg[0] == '-'){
			if(arg == "-a"){file_a = opt[a+1]; a++;}
			else if(arg == "-list_a" or arg == "-a_list"){file_list_a = opt[a+1]; a++;}
			else if(arg == "-b"){file_b = opt[a+1]; a++;}
			else if(arg == "-list_b" or arg == "-b_list"){file_list_b = opt[a+1]; a++;}
			else if(arg == "-out" or arg == "-o"){file_out = opt[a+1]; a++;}
			else if(arg == "-range_a"){
				while(a+1 < opt.size() && opt[a+1].substr(0,1) != "-"){
					string r = opt[a+1];
					int i = stoi(r.substr(0,r.find('-')));
					int j = stoi(r.substr(r.find('-')+1));

					if(j+1 > range_a.size()){range_a.resize(j+1,0);}
					for(int n = i; n <= j; n++){range_a[n] = 1;}
					a++;
				}
			}
			else if(arg == "-range_b"){
				while(a+1 < opt.size() && opt[a+1].substr(0,1) != "-"){
					string r = opt[a+1];
					int i = stoi(r.substr(0,r.find('-')));
					int j = stoi(r.substr(r.find('-')+1));
					if(j+1 > range_b.size()){range_b.resize(j+1,0);}
					for(int n = i; n <= j; n++){range_b[n] = 1;}
					a++;
				}
			}
			else if(arg == "-use_prf"){
				if(a+1 < opt.size() and opt[a+1][0] != '-'){use_prf = stoi(opt[a+1]); a++;}
				else{use_prf = true;}
			}
			else if(arg == "-prf_w"){prf_w = stod(opt[a+1]); a++;}
			else if(arg == "-use_gap_ss"){
				if(a+1 < opt.size() and opt[a+1][0] != '-'){use_gap_ss = stoi(opt[a+1]); a++;}
				else{use_gap_ss = true;}
			}
			else if(arg == "-gap_ss_w"){gap_ss_w = stod(opt[a+1]); a++;}
			else if(arg == "-gap_o"){gap_open = stod(opt[a+1]); a++;}
			else if(arg == "-gap_e"){gap_ext  = stod(opt[a+1]); a++;}
			else if(arg == "-sep_cut"){sep_cutoff = stoi(opt[a+1]); a++;}
			else if(arg == "-sco_cut"){sco_cutoff = stof(opt[a+1]); a++;}
			else if(arg == "-iter"){iter = stoi(opt[a+1]); a++;}
			else if(arg == "-silent"){
				if(a+1 < opt.size() and opt[a+1][0] != '-'){silent = stoi(opt[a+1]); a++;}
				else{silent = true;}
			}
			else if(arg == "-get_tot"){
				if(a+1 < opt.size() and opt[a+1][0] != '-'){get_tot = stoi(opt[a+1]); a++;}
				else{get_tot = true;}
			}
			else if(arg == "-cpu"){cpu = stoi(opt[a+1]); a++;}
		}
	}
	if(
	   (file_a.empty() && file_list_a.empty())||
	   (file_b.empty() && file_list_b.empty())||
	   (exists(file_a) == 0 && exists(file_list_a) == 0)||
	   (exists(file_b) == 0 && exists(file_list_b) == 0)||
	   file_out.empty()
	   ){
		cout << "---------------------------------------------------------------------\n";
		cout << "                         MAP_ALIGN - Parallel                        \n";
		cout << "---------------------------------------------------------------------\n";
		cout << "  -a             contact map A                [REQUIRED -a OR -list_a]   \n";
		cout << "  -b             contact map B                [REQUIRED -b OR -list_b]   \n";
		cout << "  -o             save results to              [REQUIRED]                 \n";
		cout << "  -cpu           number of threads to use     [Default=" << cpu        << "]\n";
		cout << "  -gap_o         gap opening penalty          [Default=" << gap_open   << "]\n";
		cout << "  -gap_e         gap extension penalty        [Default=" << gap_ext    << "]\n";
		cout << "  -sep_cut       seq seperation cutoff        [Default=" << sep_cutoff << "]\n";
		cout << "  -sco_cut       prob cutoff                  [Default=" << sco_cutoff << "]\n";
		cout << "  -iter          number of iterations         [Default=" << iter       << "]\n";
		cout << "  -silent                                     [Default=" << silent     << "]\n";
		cout << "---------------------------------------------------------------------\n";
		cout << " Advanced options                                                    \n";
		cout << "---------------------------------------------------------------------\n";
		cout << "  -list_a        list of contact maps                                \n";
		cout << "  -list_b        list of contact maps                                \n";
		cout << "  -range_a       trim map A to specified range(s) (eg. 0-20 50-100)  \n";
		cout << "  -range_b       trim map B to specified range(s)                    \n";
		cout << "---------------------------------------------------------------------\n";
		cout << " Experimental features\n";
		cout << "---------------------------------------------------------------------\n";
		cout << "  -use_gap_ss    penalize gaps at secondary   [Default=" << use_gap_ss << "]\n";
		cout << "                 structure elements(SSE)\n";
		cout << "  -gap_ss_w      gap penality weight at SSE   [Default=" << gap_ss_w   << "]\n";
		cout << "  -use_prf       use sequence profile         [Default=" << use_prf    << "]\n";
		cout << "  -prf_w         profile weight               [Default=" << prf_w      << "]\n";
		cout << "  -get_tot       compute expected total score [Default=" << get_tot    << "]\n";
		cout << "---------------------------------------------------------------------\n";
		exit(1);
	}
	else if(!silent){
		cout << "OPT -------------------------------------------------------------------\n";
		cout << "OPT                       MAP_ALIGN - Parallel                         \n";
		cout << "OPT -------------------------------------------------------------------\n";
		if(!file_a.empty()){     cout << "OPT   -a          " << file_a         << endl;}
		if(!file_list_a.empty()){cout << "OPT   -list_a     " << file_list_a    << endl;}
		if(!file_b.empty()){     cout << "OPT   -b          " << file_b         << endl;}
		if(!file_list_b.empty()){cout << "OPT   -list_b     " << file_list_b    << endl;}
		cout << "OPT   -o          " << file_out     << endl;
		cout << "OPT   -cpu        " << cpu          << endl;
		cout << "OPT   -gap_o      " << gap_open     << endl;
		cout << "OPT   -gap_e      " << gap_ext      << endl;
		cout << "OPT   -sep_cut    " << sep_cutoff   << endl;
		cout << "OPT   -sco_cut    " << sco_cutoff   << endl;
		cout << "OPT   -iter       " << iter         << endl;
		cout << "OPT   -silent     " << silent       << endl;
		cout << "OPT -------------------------------------------------------------------\n";
		cout << "OPT   -use_gap_ss  " << use_gap_ss  << endl; if(use_gap_ss == true){cout << "OPT   -gap_ss_w    " << gap_ss_w << endl;}
		cout << "OPT   -use_prf     " << use_prf     << endl; if(use_prf    == true){cout << "OPT   -prf_w       " << prf_w    << endl;}
		cout << "OPT   -get_tot     " << get_tot     << endl;
		cout << "OPT -------------------------------------------------------------------\n";
	}
	if(gap_open > 0 || gap_ext > 0){
		cout << "ERROR: gap penality should be < 0\n";
		exit(1);
	}
}
