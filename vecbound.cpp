//SEE APPENDIX c FOR A DETAILED DISCUSSION OF THIS FILE
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")
#include <bits/stdc++.h>
using namespace std;
#define ll long long

//macros to output a variable/vector to terminal/file
#define out(x) cout<<#x " "<<x<<endl
#define vout(x) cout<<#x " "; for(const auto& vectorelem: x){cout <<vectorelem << ' ';}cout<<endl
#define ffout(x) fout<<#x " "<<x<<endl
#define fvout(x) fout<<#x " "; for(const auto& vectorelem: x){fout <<vectorelem << ' ';}fout<<endl

#define pw(x) (((ll)1)<<x)

const int t = 9; //an upper bound for v(F')
//If v(F') < 9, we effectively add isolated vertices to the weighted graph until it reaches 9 vertices. This does not affect the distributions of V_F, E_F, C_F.

const int q = 3; //how many possible colors per vertex

const int E = (t*t-t)/2; //how many possible edges for the F'
//We represent unweighted graphs as long long integers, which are bitmasks of length E.

const int totalcolors = pow(q,t-1); //how many possible colorings for all vertices. Since the individual colors are identical, we may fix the coloring of a singular vertex, hence the total is pow(q,t-1) rather than pow(q,t).

const int vtxtotal = 14; //how many vertices in F
//change this to vtxtotal = 15 to verify 15-minimal hypergraphs.

const vector<ll> vmask = {255, 32513, 2064642, 65044996, 1008796680, 7587629072, 26986418208, 45382901824, 56406065280};
// vmask[i] represents a bitmask for the edges that contain vertex i.
// this allows us to do graph&vmask[i]!=0 to see if vertex i is isolated in graph or not.

//upper bound for V, E, C respectively (+1 because these determine the size of the arrays used to represent p.m.f.s, which are 0-indexed)
const int vmax = 7+1;
const int emax = 6+1;
const int cmax = 4+1;

// For an edge (a,b), tobit(a,b) is the index of edge (a,b) in the bitmask
// So graph[tobit(a,b)] is 1/0 depending on whether or not (a,b) is in graph or not respectively.
inline int to_bit(int a, int b){
    if (b < a) swap(a,b);
    return (a*(2*t-1-a))/2+(b-a)-1;
}



// Joint probability distributions of V, E, C can be represented via a 3-D vector of size [vmax][emax][cmax]. 
// To improve efficiency, we unroll this into a 1-D vector of size datalen=vmax*emax*cmax.
const int emaxtimescmax = emax*cmax;
const int datalen = vmax*emax*cmax;

//function to get the 1-D vector index from the values of V, E, and C
inline int tind(int v,int e, int c){
    return c+e*cmax+v*emaxtimescmax;
}
//functions to get the value of V, E, or C respectively from the 1-D vector index
inline int tv(int a){
    return a/emaxtimescmax;
}
inline int te(int a){
    return (a%emaxtimescmax)/cmax;
}
inline int tc(int a){
    return a%cmax;
}


//adds two polynomial distributions a and b
//for example, if G and H are disjoint and we have the probability distributions for both,
//then we can use this function to get the distirbution for G \cup H.
vector<double> polymult(vector<double> const &a, vector<double> const &b){
    vector<double> res(datalen,0);
    for (int i=0; i<datalen; i++){
        if (a[i] == 0){
            continue;
        }
        for (int j=0; j<datalen; j++){
            if (b[j] == 0){
                continue;
            }
            int newc = max(tc(i),tc(j));
            int newe = te(i)+te(j);
            int newv = tv(i)+tv(j);
            if (newc < cmax && newe < emax && newv < vmax){
                res[tind(newv,newe,newc)] += a[i]*b[j];
            }
        }
    }
    return res;
}

//adds two polynomial distributions, in the cases where the maximum codegree is additive rather than maximized.
//we'll use this, for example, if a represented the distribution for a single edge, and b represented
//a distribution for a single edge that shared 2 vertices with the edge from a, so in this scenario the codegree would be additive rather than maximized.
vector<double> cpolymult(vector<double> const &a, vector<double> const &b){
    vector<double> res(datalen,0);
    for (int i=0; i<datalen; i++){
        if (a[i] == 0){
            continue;
        }
        for (int j=0; j<datalen; j++){
            if (b[j] == 0){
                continue;
            }
            int nc = tc(i)+tc(j);
            int ne = te(i)+te(j);
            int nv = tv(i)+tv(j);
            if (nc < cmax && ne < emax && nv < vmax){
                res[tind(nv,ne,nc)] += a[i]*b[j];
            }
        }
    }
    return res;
}

int main(){
    time_t start_time = time(0);
    ios_base::sync_with_stdio(false);cin.tie(0);cout.precision(18);
    ofstream fout; //output file
    fout.open ("bounding\\vecboundout"+to_string(vtxtotal)+".txt");
    fout.precision(18);
    ifstream coefreader; //file containing our c_H
    coefreader.open("results\\coefs.txt");
    ll inputgraphbitmask; double inputgraphcoefficient;
    
    //given our c_H, these vectors store the values of the outer max term in the RHS of lemma 5.8 as pvals[i][j][k] and identically for nvals
    vector<vector<vector<double>>> pvals(vmax,vector<vector<double>>(emax,vector<double>(cmax,0)));
    vector<vector<vector<double>>> nvals(vmax,vector<vector<double>>(emax,vector<double>(cmax,0)));
    inputgraphbitmask = 2147483665; inputgraphcoefficient = 4.0; int inputvertexcount = 7; //handle the one c_H on 7 vtx first
    do{ //read our c_H and populate pvals & nvals
        int inputgraphpopcount = __popcount(inputgraphbitmask); //number of edges
        int inputgraphmaxedges = inputvertexcount*(inputvertexcount-1)*(inputvertexcount-2)/6;
        vector<int> bmmap(inputgraphmaxedges,-1);
        int cz = 0;
        for (int i=0; i < inputgraphmaxedges; i++){
            if (inputgraphbitmask&pw(i)){
                bmmap[i] = cz;
                cz++;
            }
        }
        for (int bm = 0; bm < (1<<inputgraphpopcount); bm++){ //iterate over all subgraphs (since the inner maximum is taken H' \subseteq H)
            int ind = 0;
            vector<int> deg(7,0); //degrees
            vector<int> cdeg(49,0); //codegrees
            for (int i=0; i<inputvertexcount; i++){ //computes degrees & codegrees of the subgraph
                for (int j=i+1; j<inputvertexcount; j++){
                    for (int k=j+1; k<inputvertexcount; k++){
                        if (inputgraphbitmask&(((ll)1)<<ind)){
                            if (bm&pw(bmmap[ind])){
                                deg[i]++;deg[j]++;deg[k]++;
                                cdeg[7*i+j]++;cdeg[7*i+k]++;cdeg[7*j+k]++;
                            }
                        }
                        ind++;
                    }
                }
            }
            int v = 0; //vertices in subgraph
            for (int i:deg){
                if (i > 0){
                    v++;
                }
            }
            int e = __popcount(bm); //edges in subgraph
            int c = 0; //maximum codegree of subgraph
            for (int i:cdeg){
                c = max(c,i);
            }
            if (inputgraphcoefficient > 0){ //populate pvals with this subgraph if applicable
                for (int i=0; i<=min(vmax-1,v); i++){
                    pvals[i][e][c] = max(pvals[i][e][c],inputgraphcoefficient);
                }
            } else{ //populate nvals with this subgraph if applicable
                for (int i=0; i<=min(vmax-1,v); i++){
                    nvals[i][e][c] = max(nvals[i][e][c],-inputgraphcoefficient);
                }
            }
        }
        inputvertexcount = 6; //after handling the only 7 vertex coefficient, switch to handling 6 vertex coefficients
    } while (coefreader>>inputgraphbitmask>>inputgraphcoefficient); //read from file
    
     
    vector<double> pcoef(datalen,0); //unroll pvals & nvals from 3D vectors to 1D vectors for greater efficiency
    vector<double> ncoef(datalen,0);
    for (int i=0; i<vmax; i++){
        for (int j=0; j<emax; j++){
            for (int k=0; k<cmax; k++){
                pcoef[k+j*cmax+i*cmax*emax] = pvals[i][j][k];
                ncoef[k+j*cmax+i*cmax*emax] = nvals[i][j][k];
            }
        }
    }

    map<pair<int,int>,vector<double>> extras; //extras[{i,j}] represents the joint p.m.f. of V, E, C for i 3-leaf and 2-leaf edges
    //So a pmf for (F,e_2,e_3) can be obtained by adding extras[{e_3,e_2}] to the pmf for F (calculated later).
    vector<double> v3(datalen,0); //distribution for a singular 3-leaf edge
    vector<double> v2(datalen,0); //distribution for a singular 2-leaf edge
    v3[0] = 4.0/9; v2[0] = 4.0/9;
    v3[tind(3,1,1)] = 5.0/9; v2[tind(2,1,1)] = 5.0/9;
    for (int i=0; i<100; i++){
        for (int j=0; j<100; j++){
            if (3*i+2*j <= vtxtotal){
                vector<double> res(datalen,0);
                res[0] = 1;
                for (int k=0; k<i; k++){
                    res = polymult(res,v3);
                }
                for (int k=0; k<j; k++){
                    res = polymult(res,v2);
                }
                extras[{i,j}] = res;
            }
        }
    }

    vector<ll> qpow = {1}; //powers of q
    for (int i=0; i<t+1; i++){
        qpow.push_back(q*qpow.back());
    }

    unordered_map<ll, int> colordict; //possible colorings of V(F'), in the format {bitmask of non-monochromatic edges, # of occurences (if the same bitmask happens for multiple colorings)}.
    for (int cl=0; cl<totalcolors; cl++){
        ll res = 0;
        int ind = 0;
        for (int i=0; i<t; i++){
            for (int j=i+1; j<t; j++){
                if ((cl%qpow[i+1])/qpow[i] != (cl%qpow[j+1])/qpow[j]){ //if color of i != color of j
                    res += ((ll)1)<<(ind);
                }
                ind++;
            }
        }
        if (colordict.find(res) == colordict.end()){
            colordict[res] = 0;
        }
        colordict[res]++;
    }


    //transfer unordered_map to vectors for more efficiency
    vector<ll> colors;
    vector<int> colorcts;
    for (const auto& elem: colordict){
        colors.push_back(elem.first);
        colorcts.push_back(elem.second);
    }

    cout << "Generated " << colordict.size() << " colorings"<< endl;

    vector<pair<int,int>> edges;
    for (int i=0; i<t; i++){
        for (int j=i+1; j<t; j++){
            edges.push_back({i,j});
        }
    }

    //joint pmf for a pair of monochromatic vertices in F_2 with codegree x is represented with onethirdexp[x]
    vector<vector<double>> onethirdexp = {vector<double>(datalen,0)};
    onethirdexp[0][0] = 1;
    vector<double> onethird(datalen,0);
    onethird[0] = 2.0/3;
    onethird[tind(1,1,1)] = 1.0/3;
    for (int i=0; i<vtxtotal; i++){
        onethirdexp.push_back(cpolymult(onethirdexp.back(),onethird));
    }

    //joint pmf for a pair of non-monochromatic vertices in F_2 with codegree x is represented with twothirdexp[x]
    vector<vector<double>> twothirdexp = {vector<double>(datalen,0)};
    twothirdexp[0][0] = 1;
    vector<double> twothird(datalen,0);
    twothird[0] = 1.0/3;
    twothird[tind(1,1,1)] = 2.0/3;
    for (int i=0; i<vtxtotal; i++){
        twothirdexp.push_back(cpolymult(twothirdexp.back(),twothird));
    }

    ifstream greader;
    //read all equivalence classes from file
    greader.open("bounding\\equivclasses"+to_string(vtxtotal)+".txt"); 
    char instr;
    while (greader>>instr){
        char junk;
        ll fprimebitmask;
        greader>>fprimebitmask; //read F' bitmask
        vector<int> edgeweights;
        greader>>junk;
        int edgecount = __popcount(fprimebitmask); //e(F')
        for (int i=0; i<edgecount; i++){
            int a; greader>>a; edgeweights.push_back(a); //read edge weights
        }
        int esub3; int esub2;
        greader>>junk>>esub3>>esub2; //read e_3 and e_2
        int vc = 0; //v(F')
        for (auto vm : vmask){
            if (fprimebitmask&vm){
                vc++;
            }
        }
        vector<int> edgetoweightlocationmap(E,-1); //maps an edge to the location in edgeweights that stores its weight
        int c = 0;
        for (int i=0; i<E; i++){
            if (fprimebitmask&(((ll)1)<<i)){
                edgetoweightlocationmap[i] = c;
                c++;
            }
        }
        vector<int> inverseedgetoweightlocationmap(E,-1);
        for (int i=0; i<E; i++){
            if (edgetoweightlocationmap[i] != -1){
                inverseedgetoweightlocationmap[edgetoweightlocationmap[i]] = i;
            }
        }
        vector<double> fprimepmf(datalen,0); //pmf for F'
        if (fprimebitmask == 0){
            fprimepmf[0] = 1;
        } else{
            for (int iii=0; iii<colors.size(); iii++){
                ll nonmonochromaticbitmask = colors[iii]&fprimebitmask; //bitmask for nonmonochromatic edges under this coloring
                vector<vector<double>> edgedistribs(edgecount,vector<double>()); //list of pmfs for each individual pair of vertices in F_2
                int index = 0;
                for (int i=0; i<t; i++){
                    for (int j=i+1; j<t; j++){
                        if (fprimebitmask&(((ll)1)<<index)){ //if codegree is not 0 in F
                            if (nonmonochromaticbitmask&(((ll)1)<<index)){ // diff color
                                edgedistribs[edgetoweightlocationmap[index]] = twothirdexp[edgeweights[edgetoweightlocationmap[index]]];
                            } else{ //same color
                                edgedistribs[edgetoweightlocationmap[index]] = onethirdexp[edgeweights[edgetoweightlocationmap[index]]];
                            }
                        }
                        index++;
                    }
                }
                vector<double> colorpmf(datalen,0); //final pmf GIVEN the coloring of F_2
                for (ll Q=0; Q<(((ll)1)<<edgecount); Q++){ //iterate over all cases Q
                    vector<double> colorandcasepmf(datalen,0);
                    colorandcasepmf[0] = 1;
                    for (int i=0; i<edgedistribs.size(); i++){
                        vector<double> tomult;
                        if (Q&(((ll)1)<<i)){ //codegree is non-zero in Q
                            tomult = edgedistribs[i];
                            tomult[0] = 0;
                        } else{ //codegree is 0 in Q
                            tomult = vector<double>(datalen,0);
                            tomult[0] = edgedistribs[i][0];
                        }
                        colorandcasepmf = polymult(colorandcasepmf,tomult);
                    }
                    vector<bool> vexist(t,false);
                    for (int i=0; i<edgedistribs.size(); i++){
                        if (Q&(((ll)1)<<i)){ //mark vertices that are contained some pair in Q
                            vexist[edges[inverseedgetoweightlocationmap[i]].first] = true;
                            vexist[edges[inverseedgetoweightlocationmap[i]].second] = true;
                        }
                    }
                    int alpha_Q = 0;
                    for (bool i:vexist){
                        if (i){
                            alpha_Q++;
                        }
                    }
                    for (int i=0; i<datalen; i++){
                        if (tv(i)+alpha_Q < vmax){ //shift distrib by alpha_Q
                            colorpmf[tind(tv(i)+alpha_Q,te(i),tc(i))] += colorandcasepmf[i];
                        }
                    }
                }
                for (int i=0; i<colorpmf.size(); i++){
                    fprimepmf[i] += colorpmf[i]*colorcts[iii];
                }
            }
            for (int i=0; i<fprimepmf.size(); i++){
                fprimepmf[i] /= totalcolors;
            }
        }
        vector<double> finaldistrib = polymult(fprimepmf,extras[{esub3,esub2}]); //Algorithm 5.13
        double posmubound = 0;
        double negmubound = 0;
        for (int i=0; i<datalen; i++){ //Lemma 5.8
            posmubound += pcoef[i]*finaldistrib[i];
            negmubound += ncoef[i]*finaldistrib[i];
        }
        if (max(posmubound,negmubound) > 1){ //equivalence class fails to be bounded by Lemma 5.8, output to file
            out(fprimebitmask);
            ffout(fprimebitmask);
            vout(edgeweights);
            fvout(edgeweights);
            cout << "extra " << esub3 << " " << esub2 << endl;
            fout << "extra " << esub3 << " " << esub2 << endl;
            out(posmubound);
            out(negmubound);
            ffout(posmubound);
            ffout(negmubound);
            cout << endl;
            fout << endl;
        }
    }
    
    fout << "Check finished in " << difftime(time(0), start_time) << " seconds" <<  endl;
    cout << "Check finished in " << difftime(time(0), start_time) << " seconds" <<  endl;
    fout.close();
    return 0;
}