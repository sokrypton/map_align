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
#include <tuple>
#include <string>

using namespace std;

typedef vector<vector<int> > mtx_int;
typedef vector<vector<double> > mtx_double;
typedef vector<int> vec_int;
typedef vector<double> vec_double;
typedef vector<char> vec_char;

double exp_fast(double x){
    // WARNING fails if |x| > 1024
    //https://codingforspeed.com/using-faster-exponential-approximation/
    x = 1 + x/1024;
    x *= x; x *= x; x *= x; x *= x;
    x *= x; x *= x; x *= x; x *= x;
    x *= x; x *= x;
    return x;
}

inline bool exists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

double gaussian(double mean, double stdev, double x){return exp_fast(-pow((x - mean),2)/(2*(pow(stdev,2))));}

vec_int align(double gap_o, double gap_e, mtx_double &sco_mtx);
double Falign(double *sco_mtx, int rows, int cols);

double sepw(double sep){if(sep <= 4){return 0.50;}else if(sep == 5){return 0.75;}else{return 1.00;}}

void ini_SCO(double sep_x, double sep_y, mtx_double &SCO, vec_int &vec_a_div,vec_int &vec_b_div,vec_int &vec_a,vec_int &vec_b,mtx_int &vec_a_i,mtx_int &vec_b_i,mtx_double &mtx_a,mtx_double &mtx_b);
void ini_prf_SCO(mtx_double &P_SCO, mtx_double &prf_a, vec_char &aa_a, mtx_double &prf_b, vec_char &aa_b);
vec_int mod_SCO(double do_it, double gap_o, double gap_e, mtx_double &SCO, vec_int &vec_a_div,vec_int &vec_b_div,vec_int &vec_a, vec_int &vec_b,mtx_int &vec_a_i, mtx_int &vec_b_i,mtx_double &mtx_a, mtx_double &mtx_b);
void chk (double gap_o, double gap_e, double& con_sco,double& gap_sco,vec_int& vec_a_div,mtx_int& vec_a_i,mtx_double& mtx_a,mtx_double& mtx_b,vec_int& a2b);

void add_mtx (mtx_double &A, mtx_double &B){
    for(int i = 0; i < A.size(); i++){
        for(int j = 0; j < A[i].size(); j++){
            A[i][j] += B[i][j];
        }
    }
}
void load_data (string file, int ss_cutoff, mtx_double &mtx, vec_int &vec_div, vec_int &vec, mtx_int &vec_i, mtx_double &prf, vec_char &aa, vec_char &ss);
int main(int argc, const char * argv[])
{
    
    string file_a;
    string file_b;
    bool use_prf = false;
    double gap_open = -1.00;
    double gap_ext = -0.01;
    int ss_cutoff = 3;
    ////////////////////////////////////////////////////////////////////////////////
    for (int a = 1; a < argc; a++)
    {
        string arg = argv[a];
        if (arg.substr(0,1) == "-")
        {
                 if(arg == "-a"){file_a = argv[a+1]; a++;}
            else if(arg == "-b"){file_b = argv[a+1]; a++;}
            else if(arg == "-prf"){use_prf = true;}
            else if(arg == "-gap_o"){gap_open = stod(argv[a+1]); a++;}
            else if(arg == "-gap_e"){gap_ext  = stod(argv[a+1]); a++;}
            else if(arg == "-ss_cut"){ss_cutoff  = stoi(argv[a+1]); a++;}
        }
    }
    if(file_a.empty() || file_b.empty() || exists(file_a) == 0 || exists(file_b) == 0)
    {
        cout << "-----------------------------------------------------\n";
        cout << "                      MAP_ALIGN                      \n";
        cout << "-----------------------------------------------------\n";
        cout << "  -a      contact map A             [REQUIRED]\n";
        cout << "  -b      contact map B             [REQUIRED]\n";
        cout << "  -prf    use sequence profile      [Default=" << use_prf   << "]\n";
        cout << "  -gap_o  gap opening penalty       [Default=" << gap_open  << "]\n";
        cout << "  -gap_e  gap extension penalty     [Default=" << gap_ext   << "]\n";
        cout << "  -ss_cut seq seperation cutoff     [Default=" << ss_cutoff << "]\n";
        cout << "-----------------------------------------------------\n";
        exit(1);
    }
    else
    {
        cout << "----------------------------------------------------\n";
        cout << "                      MAP_ALIGN                     \n";
        cout << "----------------------------------------------------\n";
        cout << "  -a      " << file_a    << endl;
        cout << "  -b      " << file_b    << endl;
        cout << "  -prf    " << use_prf   << endl;
        cout << "  -gap_o  " << gap_open  << endl;
        cout << "  -gap_e  " << gap_ext   << endl;
        cout << "  -ss_cut " << ss_cutoff << endl;
        cout << "----------------------------------------------------\n";
    }
    ////////////////////////////////////////////////////////////////////////////////
    
    // load data from contact map A
    mtx_double mtx_a; vec_int vec_a; vec_int vec_a_div; mtx_int vec_a_i; mtx_double prf_a; vec_char aa_a; vec_char ss_a;
    load_data(file_a,ss_cutoff,mtx_a,vec_a_div,vec_a,vec_a_i,prf_a,aa_a,ss_a);
    int size_a = mtx_a.size();
    
    // load data from contact map B
    mtx_double mtx_b; vec_int vec_b; vec_int vec_b_div; mtx_int vec_b_i; mtx_double prf_b; vec_char aa_b; vec_char ss_b;
    load_data(file_b,ss_cutoff,mtx_b,vec_b_div,vec_b,vec_b_i,prf_b,aa_b,ss_b);
    int size_b = mtx_b.size();
    
    mtx_double P_SCO;if(use_prf == true){ini_prf_SCO(P_SCO,prf_a,aa_a,prf_b,aa_b);}
    
    // keeping track of the BEST alignment
    
    int max_sep_x = 0;
    int max_sep_y = 0;
    int max_mode = 0;
    int max_g_e = 0;
    double con_max = -1;
    double gap_max = 0;
    vec_int a2b_max;
    
    vec_double gap_e_steps {5,10,100,1000};
    
    // try different sep (sequence seperation difference) penalities
    vec_double sep_x_steps {0,1,2}; // (constant, linear, quadratic)
    for(int sx = 0; sx < sep_x_steps.size(); sx++){double sep_x = sep_x_steps[sx];
        
        //try different scaling factors for sep penalities
        vec_double sep_y_steps {1,2,4,8,16,32};
        for(int sy = 0; sy < sep_y_steps.size(); sy++){double sep_y = sep_y_steps[sy];
            
            // Get initial score matrix
            mtx_double C_SCO(size_a,vector<double>(size_b,0));
            ini_SCO(sep_x,sep_y,C_SCO,vec_a_div,vec_b_div,vec_a,vec_b,vec_a_i,vec_b_i,mtx_a,mtx_b);
            
            // try adding profile info (mode 0 = with; mode 1 = without)
            int start_mode = 1; if(use_prf == true){start_mode = 0;}
            for(int mode = start_mode; mode <= 1; mode++){
                
                // try different gap_ext penalities!
                for(int g_e = 0; g_e < gap_e_steps.size(); g_e++){double gap_e = gap_open/gap_e_steps[g_e];
                    
                    // restart SCO matrix
                    mtx_double SCO = C_SCO;if(mode == 0){add_mtx(SCO,P_SCO);}
                    
                    // get alignment (a2b mapping)
                    vec_int a2b = mod_SCO(20,gap_open,gap_e,SCO,vec_a_div,vec_b_div,vec_a,vec_b,vec_a_i,vec_b_i,mtx_a,mtx_b);
                    
                    // compute number of contacts/gaps made
                    double con_sco = 0;double gap_sco = 0; chk(gap_open,gap_ext,con_sco,gap_sco,vec_a_div,vec_a_i,mtx_a,mtx_b,a2b);
                    
                    // print info
                    cout << "TMP " << " " << sep_x << "_" << sep_y << "_" << mode << "_" << g_e << " " << con_sco << " " << gap_sco << " " << con_sco+gap_sco << endl;
                    
                    // save if BEST!
                    if(con_sco+gap_sco > con_max+gap_max){
                        max_sep_x = sep_x;
                        max_sep_y = sep_y;
                        max_g_e = g_e;
                        max_mode = mode;
                        con_max = con_sco;
                        gap_max = gap_sco;
                        a2b_max = a2b;
                    }
                }
            }
        }
    }
    int aln_len = 0;for(int ai = 0; ai < size_a; ai++){int bi = a2b_max[ai];if(bi != -1){aln_len++;}}
    // Report the BEST score
    cout << "MAX " << max_sep_x << "_" << max_sep_y << "_" << max_mode << "_" << max_g_e << " " << file_a << " " << file_b << " " << con_max+gap_max << " " << aln_len;
    for(int a = 0; a < size_a; a++){int b = a2b_max[a];if(b != -1){cout << " " << a << ":" << b;}}
    cout << endl;
    return 0;
}
vec_int align(double gap_o, double gap_e, mtx_double &sco_mtx)
{
    // LOCAL_ALIGN
    // Start	0
    // [A]lign	1
    // [D]own	2
    // [R]ight	3
    
    double max_sco = 0;
    int rows = sco_mtx.size();
    int cols = sco_mtx[0].size();
    
    vec_int a2b(rows,-1);
    
    mtx_double sco(rows+1,vector<double>(cols+1,0));
    mtx_int label(rows+1,vector<int>(cols+1,0));
    
    int max_i = 0;int max_j = 0;
    for (int i = 1; i <= rows; i++){
        for (int j = 1; j <= cols; j++){
            
            double A = sco[i-1][j-1] + sco_mtx[i-1][j-1];
            double D = sco[i-1][j];
            double R = sco[i][j-1];
            if(label[i-1][j] == 1){D += gap_o;}else{D += gap_e;}
            if(label[i][j-1] == 1){R += gap_o;}else{R += gap_e;}
            if(A <= 0 and D <= 0 and R <= 0){label[i][j] = 0;sco[i][j] = 0;}
            else
            {
                if(A >= R){if(A >= D){label[i][j] = 1;sco[i][j] = A;}else{label[i][j] = 2;sco[i][j] = D;}}
                else{if(R >= D){label[i][j] = 3;sco[i][j] = R;}else{label[i][j] = 2;sco[i][j] = D;}}
                if(sco[i][j] > max_sco){max_i = i;max_j = j;max_sco = sco[i][j];}
            }
        }
    }
    int i = max_i;int j = max_j;
    while(1){
        if(label[i][j] == 0){break;}
        if(label[i][j] == 1)
        {
            a2b[i-1] = j-1;
            i--;j--;
        }
        else if(label[i][j] == 2){i--;}
        else if(label[i][j] == 3){j--;}
    }
    return(a2b);
}
double Falign(double *sco_mtx, int rows, int cols)
{
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
void load_data (string file, int ss_cutoff, mtx_double &mtx, vec_int &vec_div, vec_int &vec, mtx_int &vec_i, mtx_double &prf, vec_char &aa, vec_char &ss)
{
    string line;
    ifstream in(file.c_str());
    while(getline(in,line)){
        istringstream is(line);
        string label;
        is >> label;
        if(label == "LEN" or label == "SIZE"){
            int size; is >> size;
            for (int i=0; i < size; i++){
                mtx.push_back(vector<double>());
                prf.push_back(vector<double>());
                for (int j=0; j < size; j++){mtx[i].push_back(0);}
            }
        }
        else if(label == "CON"){
            int i, j; double sco; is >> i >> j;
            if(abs(j-i) >= ss_cutoff){
                if(is >> sco){}else{sco = 1;}
                mtx[i][j] = sco;
                mtx[j][i] = sco;
            }
        }
        else if(label == "PRF"){
            int i; is >> i;
            char tmp;
            is >> tmp; aa.push_back(tmp);
            is >> tmp; ss.push_back(tmp);
            double val;
            while(is >> val){prf[i].push_back(val);}
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
}
// INITIATE SCORE MATRIX: function for populating the initial similarity matrix
void ini_SCO(double sep_x, double sep_y, mtx_double &SCO,
             vec_int &vec_a_div,vec_int &vec_b_div,vec_int &vec_a,vec_int &vec_b,mtx_int &vec_a_i,mtx_int &vec_b_i,mtx_double &mtx_a,mtx_double &mtx_b)
{
    // Get initial score matrix
    for(int i=0; i < vec_a.size(); i++){ // go through columns (vec_a) in map_a that has contacts
        int ai = vec_a[i];
        for(int j=0; j < vec_b.size(); j++){ // go through columns (vec_b) in map_b that has contacts
            int bi = vec_b[j];
            int A[2] = {(int)vec_a_div[ai],(int)(vec_a_i[ai].size()-vec_a_div[ai])};
            int B[2] = {(int)vec_b_div[bi],(int)(vec_b_i[bi].size()-vec_b_div[bi])};
            for(int k=0; k <= 1; k++){ // left and right of diagonal
                if(A[k] > 0 and B[k] > 0){
                    double M[A[k]*B[k]];
                    for(int n=0; n < A[k]; n++){
                        int nn = n; if(k == 1){nn += vec_a_div[ai];}
                        int aj = vec_a_i[ai][nn];
                        int sep_a = abs(ai-aj);
                        for(int m=0; m < B[k]; m++){
                            int mm = m; if(k == 1){mm += vec_b_div[bi];}
                            int bj = vec_b_i[bi][mm];
                            int sep_b = abs(bi-bj);
                            int sep_D = abs(sep_a-sep_b);
                            double sep_M = min(sep_a,sep_b);
                            double sep_std = sep_y*(1+pow(sep_M-2,sep_x));
                            if(sep_D/sep_std < 6){
                                M[n*B[k]+m] = mtx_a[ai][aj] * mtx_b[bi][bj] * sepw(sep_M) * gaussian(0,sep_std,sep_D);
                            }else{M[n*B[k]+m] = 0;}
                        }
                    }
                    SCO[ai][bi] += Falign(M,A[k],B[k]);
                }
            }
        }
    }
}
void ini_prf_SCO(mtx_double &P_SCO, mtx_double &prf_a, vec_char &aa_a, mtx_double &prf_b, vec_char &aa_b)
{
    int size_a = prf_a.size();
    int size_b = prf_b.size();
    
    // compute profile similarity
    P_SCO.resize(size_a,vector<double>(size_b,0));
    vec_double pb(20,0); int prf_size = prf_a[0].size();
    double pb_size = 0;
    for(int ai = 0; ai < size_a; ai++)
    {
        if(aa_a[ai] != 'X')
        {
            for(int p=0; p < prf_size; p++){pb[p] += prf_a[ai][p];}
            pb_size += 1;
        }
    }
    for(int bi = 0; bi < size_b; bi++)
    {
        if(aa_b[bi] != 'X')
        {
            for(int p=0; p < prf_size; p++){pb[p] += prf_b[bi][p];}
            pb_size += 1;
        }
    }
    for (int i=0; i < size_a; i++){
        for (int j=0; j < size_b; j++){
            if(aa_a[i] == 'X' or aa_b[j] == 'X'){P_SCO[i][j] = 0;}
            else{
                double tmp_sco = 0;
                for(int p=0; p < prf_size; p++){tmp_sco += (prf_a[i][p]*prf_b[j][p])/(pb[p]/pb_size);}
                P_SCO[i][j] = log2(tmp_sco)/5;
            }
        }
    }
    
}
// MODIFY SCORE MATRIX: function for modifying the initial similarity matrix
vec_int mod_SCO(double do_it, double gap_o, double gap_e, mtx_double &SCO,
                vec_int &vec_a_div,vec_int &vec_b_div,vec_int &vec_a, vec_int &vec_b,mtx_int &vec_a_i, mtx_int &vec_b_i,mtx_double &mtx_a, mtx_double &mtx_b)
{
    // iterate
    vec_int a2b_tmp;
    for(int it=1; it < do_it; it++)
    {
        // align
        a2b_tmp = align(gap_o,gap_e,SCO);
        double IT = (double)it;
        double s1 = (IT/(IT+1)); double s2 = (1/(IT+1));
        for(int a=0; a < vec_a.size(); a++){ // go through columns (vec_a) in map_a that has contacts
            int ai = vec_a[a];
            for(int b=0; b < vec_b.size(); b++){ // go through columns (vec_b) in map_b that has contacts
                int bi = vec_b[b];
                double sco_contact = 0;
                for(int n=0; n < vec_a_i[ai].size(); n++){ // go through contacts in vec_a
                    int aj = vec_a_i[ai][n];
                    int bj = a2b_tmp[aj]; // get mapping
                    if(bj != -1){ // if mapping exists
                        if((ai > aj and bi > bj) or (ai < aj and bi < bj)){ // if ai-aj in same direction as bi-bj
                            double sep_M = min(abs(ai-aj),abs(bi-bj));
                            sco_contact += mtx_a[ai][aj] * mtx_b[bi][bj] * sepw(sep_M);
                        }
                    }
                }
                SCO[ai][bi] = s1*SCO[ai][bi] + s2*sco_contact;
            }
        }
    }
    return(a2b_tmp);
}
// CHK: check number of contact and gaps made and compute alignment score
void chk (double gap_o, double gap_e, double& con_sco,double& gap_sco,vec_int& vec_a_div,mtx_int& vec_a_i,mtx_double& mtx_a,mtx_double& mtx_b,vec_int& a2b)
{
    
    int size_a = mtx_a.size();
    // compute number of contacts/gaps made
    int a = 0;int b = 0;
    for(int ai = 0; ai < size_a; ai++){
        int bi = a2b[ai];
        if(bi != -1){
            if(a > 0){ // compute number of gaps
                double gap = ((ai-a)-1) + ((bi-b)-1);
                if(gap > 0){gap_sco += (gap_o + gap_e*(gap-1))/2;}
            }
            for(int m=0; m < vec_a_div[ai]; m++){ // compute number of contacts
                int aj = vec_a_i[ai][m];
                int bj = a2b[aj];
                if(bj != -1){
                    double sep_M = min(abs(ai-aj),abs(bi-bj));
                    con_sco += mtx_a[ai][aj] * mtx_b[bi][bj] * sepw(sep_M);
                    
                    
                }
            }
            a = ai;b = bi;
        }
    }
}