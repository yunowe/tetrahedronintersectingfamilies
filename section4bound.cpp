#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")
#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define out(x) cout<<#x " "<<x<<endl
#define vout(x) cout<<#x " "; for(const auto& vectorelem: x){cout <<vectorelem << ' ';}cout<<endl

const int t = 10; //number of vertices in the hypergraphs. must be run seperately for different values of V(K)
const int edgecheck = 9; //number of edges in the hypergraphs we wish to check
const int q = 3; //number of colors
const int E = (t*(t-1)*(t-2))/6; //number of possible edges on a 3-uniform hypergraph w/ t vertices.
const int totalcolors = pow(q,t-1); //how many possible colorings for all vertices. Since the individual colors are identical, we may fix the coloring of a singular vertex, hence the total is pow(q,t-1) rather than pow(q,t).
const string b64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
const int fterm = 2 + (-1+(2*t-3))/2 - ((t-3)*(t-2)*(t-1))/6; //used in to_bit


int to_bit(int x, int y, int z){ //given an edge e={x,y,z}, maps it to the index of that edge in the bitmask
    if (z < y) swap(y,z);
    if (y < x) swap(x,y);
    if (z < y) swap(y,z);
    return z + (-y*y+(2*t-3)*y)/2 - ((t-3-x)*(t-2-x)*(t-1-x))/6 - fterm;
}

bool last7alltrue(vector<bool> v){
    for (int i=v.size()-1; i>=v.size()-7; i++){
        if (!v[i]){
            return false;
        }
    }
    return true;
}


int main(){
    ios_base::sync_with_stdio(false);cin.tie(0);cout.precision(18);
    ofstream fout;
    fout.open ("sec4results\\sec4out"+to_string(t)+".txt"); //output file
    fout.precision(18);
    cout << "Program started" << endl;
    fout << "Program started" << endl;
    unordered_map<bitset<E>,double> pcoef; //isomorphism lookup table, values of the max term used in the computation of M+(K) and M-(K)
    unordered_map<bitset<E>,double> ncoef;
    int maxpc = 0;
    int cnum = 0;
    vector<bitset<E>> vmask;
    // vmask[i] represents a bitmask for the edges that contain vertex i.
    // this allows us to do hypergraph&vmask[i]!=0 to see if vertex i is isolated in hypergraph or not.

    for (int vtx=0; vtx<t; vtx++){
        int ind = 0;
        bitset<E> mask;
        for (int i=0; i<t; i++){
            for (int j=i+1; j<t; j++){
                for (int k=j+1; k<t; k++){
                    if (i==vtx || j==vtx || k == vtx){
                        mask.set(ind);
                    }
                    ind++;
                }
            }
        }
        vmask.push_back(mask);
    }
    vector<bool> combi(t, false); //for all coefficients c_H except for one, H is on at most 6 vertices.
    //its most efficient to choose which 6 vertices will be used by H then permute within those 6 vertices rather than permuting directly.
    //this takes (t choose 6)*6! iterations rather than t!
    vector<vector<int>> edgeperms; //the resulting edge mappings of all possible permutations of vertices
    fill(combi.begin(),combi.begin()+6,true);
    do {
        vector<int> perm;
        for (int i=0; i<combi.size(); i++){
            if (combi[i]){
                perm.push_back(i);
            }
        }
        do{
            edgeperms.push_back(vector<int>(E,0));
            for (int i=0; i<6; i++){
                for (int j=i+1; j<6; j++){
                    for (int k=j+1; k<6; k++){
                        edgeperms.back()[to_bit(i,j,k)] = to_bit(perm[i],perm[j],perm[k]);
                    }
                }
            }
        }while (next_permutation(perm.begin(),perm.end()));
    } while (prev_permutation(combi.begin(), combi.end()));
    cout << "Permutations generated" << endl;
    fout << "Permutations generated" << endl;
    ifstream coefreader;
    coefreader.open("results\\monotonecoefs.txt"); //contains the coefficients used when computing M+(K) and M-(K)
    ll inp1; double inp2; double inp3;
    while (coefreader>>inp1>>inp2>>inp3){ 
        int ind = 0;
        bitset<E> coefvec;
        for (int i=0; i<6; i++){
            for (int j=i+1; j<6; j++){
                for (int k=j+1; k<6; k++){
                    if (inp1&(((ll)1)<<ind)){
                        coefvec.set(to_bit(i,j,k));
                    }
                    ind++;
                }
            }
        }
        cnum++;
        maxpc = max(maxpc,static_cast<int>(coefvec.count())); //maximum number of edges in H with nonzero c_H
        for (int j=0; j<edgeperms.size(); j++){ //iterate through permutations
            bitset<E> res = 0; //resulting graph after vertex permutation
            for (int h=0; h<E; h++){
                if (coefvec.test(h)){
                    res.set(edgeperms[j][h]);
                }
            }
            pcoef[res] = abs(inp2); //populate isomorphism lookup table
            ncoef[res] = abs(inp3);
        }
        cout << cnum << '\r';
        cout.flush();
    }
    coefreader.close();
    ifstream extrareader; //we have 1 H on 7 vertices so we handle it seperately, the process is identical to those with 6 vertices
    extrareader.open("results\\monotoneextracoefs.txt");
    int ptct = 0;
    while (extrareader>>inp1>>inp2>>inp3){
        int ind = 0;
        bitset<E> coefvec;
        for (int i=0; i<7; i++){
            for (int j=i+1; j<7; j++){
                for (int k=j+1; k<7; k++){
                    if (inp1&(((ll)1)<<ind)){
                        coefvec.set(to_bit(i,j,k));
                    }
                    ind++;
                }
            }
        }
        fill(combi.begin(),combi.end(),false);
        fill(combi.begin(),combi.begin()+7,true);
        do {
            vector<int> perm;
            for (int i=0; i<combi.size(); i++){
                if (combi[i]){
                    perm.push_back(i);
                }
            }
            do{
                vector<int> ep(E,0);
                for (int i=0; i<7; i++){
                    for (int j=i+1; j<7; j++){
                        for (int k=j+1; k<7; k++){
                            ep[to_bit(i,j,k)] = to_bit(perm[i],perm[j],perm[k]);
                        }
                    }
                }
                bitset<E> res = 0;
                for (int h=0; h<E; h++){
                    if (coefvec.test(h)){
                        res.set(ep[h]);
                    }
                }
                pcoef[res] = abs(inp2);
                ncoef[res] = abs(inp3);
                ptct++;
            }while (next_permutation(perm.begin(),perm.end()));
            cout << ptct << '\r';
            cout.flush();
        } while (prev_permutation(combi.begin(), combi.end()));
    }
    extrareader.close();
    //end processing of 7 vertex H

    cout << "Generated "<< pcoef.size() << " permutations of H with non-zero c_H" << endl;
    fout << "Generated "<< pcoef.size() << " permutations of H with non-zero c_H" << endl;

    vector<ll> qpow = {1}; //powers of q
    for (int i=0; i<t+1; i++){
        qpow.push_back(q*qpow.back());
    }

    unordered_map<bitset<E>, int> colordict;
    //generate bitsets of which edges are deleted under each coloring.
    //since some colorings produce the same bitset, colordict is structured as {bitmask,frequency} for greater efficiency
    for (int cl=0; cl<totalcolors; cl++){
        bitset<E> res = 0;
        int ind = 0;
        for (int i=0; i<t; i++){
            for (int j=i+1; j<t; j++){
                for (int k=j+1; k<t; k++){
                    int ci = (cl%qpow[i+1])/qpow[i];
                    int cj = (cl%qpow[j+1])/qpow[j];
                    int ck = (cl%qpow[k+1])/qpow[k];
                    if ((ci != cj) && (cj != ck) && (ck != ci)){ //check the two ways for a edge to remain (pairwise different colors OR sum to 2 mod 3)
                        res.set(ind);
                    } else if ((ci+cj+ck)%3 == 2){
                        res.set(ind);
                    }
                    ind++;
                }
            }
        }
        if (colordict.find(res) == colordict.end()){
            colordict[res] = 0;
        }
        colordict[res]++;
    }

    cout << "Generated " << colordict.size() << " colorings"<< endl;
    fout << "Generated " << colordict.size() << " colorings"<< endl;

    //convert colordict to pair of vectors for more efficiency
    vector<bitset<E>> colors;
    vector<int> colorcts;
    for (const auto& elem: colordict){
        colors.push_back(elem.first);
        colorcts.push_back(elem.second);
    }

    vector<bitset<E>> hypergraphlist;
    ifstream isocreader;
    isocreader.open("sec4graphlists\\newb64graphs"+to_string(t)+".txt");
    string isocinp;

    //read hypergraphs K from file
    
    while (isocreader>>isocinp){
        bitset<E> toadd;
        for (int i=0; i<isocinp.size(); i++){
            int val = b64.find(isocinp[i]);
            for (int j=0; j<6; j++){
                if (val&(1<<j)){
                    toadd.set(6*i+j);
                }
            }
        }
        hypergraphlist.push_back(toadd);
    }

    for (int ii=0; ii<hypergraphlist.size(); ii++){ //check each K individually
        bitset<E> K = hypergraphlist[ii];
        if (K.count() != edgecheck){ //only check hypergraphs with correct number of edges
            continue;
        }
        double Mpos = 0; 
        double Mneg = 0;
        int vertexcount = 0;
        for (auto i:vmask){
            if ((i&K).any()){
                vertexcount++;
            }
        }
        if (vertexcount != t){ //check if correct number of vertices
            continue;
        }

        for (int i=0; i<colors.size(); i++){ //for each color,
            bitset<E> Kstar = K&colors[i]; //varphi(K)
            if (Kstar.count() <= maxpc){
                if (pcoef.find(Kstar) != pcoef.end()){ //compute term of summation
                    Mpos += pcoef[Kstar]*colorcts[i];
                    Mneg += ncoef[Kstar]*colorcts[i];
                }
            }
        }
        Mpos /= totalcolors;
        Mneg /= totalcolors;

        if (Mpos > 1 || (-Mpos) > 1 || Mneg > 1 || (-Mneg) > 1){ //if fail to bound, output to log
                cout << "BAD " << K << ' ' << Mpos << ' ' <<Mneg << endl;
                fout << "BAD " << K << ' ' << Mpos << ' ' <<Mneg << endl;
        }

        if (ii%1000 == 0){ //status update
            cout << ii << '\r';
            cout.flush();
        }

    }

    cout << "Check finished" << endl;
    fout << "Check finished" << endl;

    fout.close();
    return 0;
}