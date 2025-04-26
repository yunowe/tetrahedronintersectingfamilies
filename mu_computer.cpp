//Please refer to the comments in section4bound.cpp for an explanation of how this file works
//Recall that M(K) is computed essentially identically to mu(G), and the file that computes M, section4bound.cpp, has explanations and comments in the code.
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")
#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define out(x) cout<<#x " "<<x<<endl
#define vout(x) cout<<#x " "; for(const auto& vectorelem: x){cout <<vectorelem << ' ';}cout<<endl

const int t = 9;
const int q = 3;
const int E = (t*(t-1)*(t-2))/6;
const int totalcolors = pow(q,t-1);
const string b64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";


const int fterm = 2 + (-1+(2*t-3))/2 - ((t-3)*(t-2)*(t-1))/6;
int to_bit(int x, int y, int z){ 
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
    cout << "Program started" << endl;
    unordered_map<bitset<E>,double> coef;
    int maxpc = 0;
    vector<bool> combi(t, false);
    vector<vector<int>> edgeperms;
    int cnum = 0;
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
    ifstream coefreader;
    coefreader.open("results\\coefs.txt");
    ll inp1; double inp2;
    while (coefreader>>inp1>>inp2){
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
        maxpc = max(maxpc,static_cast<int>(coefvec.count()));
        for (int j=0; j<edgeperms.size(); j++){
            bitset<E> res = 0;
            for (int h=0; h<E; h++){
                if (coefvec.test(h)){
                    res.set(edgeperms[j][h]);
                }
            }
            coef[res] = inp2;
        }
        cout << cnum << '\r';
        cout.flush();
    }
    coefreader.close();
    ifstream extrareader;
    extrareader.open("results\\extracoefs.txt");
    int ptct = 0;
    while (extrareader>>inp1>>inp2){
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
                coef[res] = inp2;
                ptct++;
            }while (next_permutation(perm.begin(),perm.end()));
            cout << ptct << '\r';
            cout.flush();
        } while (prev_permutation(combi.begin(), combi.end()));
    }
    extrareader.close();

    cout << "Generated "<< coef.size() << "permutations of H with non-zero c_H" << endl;

    vector<ll> qpow = {1};
    for (int i=0; i<t+1; i++){
        qpow.push_back(q*qpow.back());
    }

    unordered_map<bitset<E>, int> colordict;
    for (int cl=0; cl<totalcolors; cl++){
        bitset<E> res = 0;
        int ind = 0;
        for (int i=0; i<t; i++){
            for (int j=i+1; j<t; j++){
                for (int k=j+1; k<t; k++){
                    int ci = (cl%qpow[i+1])/qpow[i];
                    int cj = (cl%qpow[j+1])/qpow[j];
                    int ck = (cl%qpow[k+1])/qpow[k];
                    if ((ci != cj) && (cj != ck) && (ck != ci)){
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

    vector<bitset<E>> colors;
    vector<int> colorcts;
    for (const auto& elem: colordict){
        colors.push_back(elem.first);
        colorcts.push_back(elem.second);
    }

    vector<bitset<E>> isoclasses;
    ifstream isocreader;
    isocreader.open("data\\graphs"+to_string(t)+".txt");
    string isocinp;
    
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
        isoclasses.push_back(toadd);
    }

    cout << "Checking " << isoclasses.size() << " graphs" << endl;
    for (int ii=0; ii<isoclasses.size(); ii++){
        bitset<E> iso = isoclasses[ii];
        double mu = 0;

        for (int i=0; i<colors.size(); i++){
            bitset<E> G = iso&colors[i];
            if (G.count() <= maxpc){
                if (coef.find(G) != coef.end()){
                    mu += coef[G]*colorcts[i];
                }
            }
        }
        mu /= totalcolors;

        if (mu > 1 || (-mu) > 1){
                cout << "BAD " << iso << ' ' << mu << endl;
        }

        if (ii%1000 == 0){
            cout << ii << '\r';
            cout.flush();
        }

    }

    cout << "Check finished" << endl;

    return 0;
}