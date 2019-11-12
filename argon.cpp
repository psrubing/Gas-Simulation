#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include "argon.h"
#include <cstdlib>


using namespace std;

void savefile(State &state,char *oFile,int N,string rp,bool b,double t){


    //  Zapis wartości XYZ
    if (rp == "r") {
        ofstream ofile;
        if (b) {
            ofile.open(oFile,std::ios_base::app);
        }
        else {
            ofile.open(oFile);
        }

        if(ofile.is_open()){
            ofile<<"125"<<"\n"<<endl;
            for (int i=0;i<N;i++)ofile<<"Ar "<<state.r[i][0]<<' '<<state.r[i][1]<<' '<<state.r[i][2]<<endl;
            ofile.close();
        }
        else {
            cout << "Unable to save file!"<<endl;
            exit(1);
        }
    }

    //  Zapis wartości Px Py Pz
    if (rp == "p") {
        ofstream ofile;
        if (b) {
            ofile.open(oFile,std::ios_base::app);
        }
        else {
            ofile.open(oFile);
        }
        if(ofile.is_open()){
            for (int i=0;i<N;i++)ofile<<state.p[i][0]<<'\t'<<state.p[i][1]<<'\t'<<state.p[i][2]<<endl;
            ofile<<"\n"<<endl;
            ofile.close();
        }
        else {
            cout << "Unable to save file!"<<endl;
            exit(1);
        }
    }

    //  Zapis wartości Fx Fy Fz
    if (rp == "f") {
        ofstream ofile;
        if (b) {
            ofile.open(oFile,std::ios_base::app);
        }
        else {
            ofile.open(oFile);
        }

        if(ofile.is_open()){
            for (int i=0;i<N;i++)ofile<<state.F[i][0]<<'\t'<<state.F[i][1]<<'\t'<<state.F[i][2]<<endl;
            ofile<<"\n"<<endl;
            ofile.close();
        }
        else {
            cout << "Unable to save file!"<<endl;
            exit(1);
        }
    }

    //  Zapis wartości stanu
    if (rp == "s") {
        ofstream ofile;
        if (b) {
            ofile.open(oFile,std::ios_base::app);
        }
        else {
            ofile.open(oFile);
        }

        if(ofile.is_open()){
            if(b==0){
                ofile << "TIME" << '\t' << '\t' << "H" << '\t' << '\t' << "V" << '\t' << '\t' << "T" << '\t' << '\t' << "P" << endl;
                ofile <<t<<'\t'<< '\t' << state.H << '\t'<< '\t' << state.V << '\t'<< '\t' << state.T << '\t'<< '\t' << state.P << endl;
            }
            else{
                ofile <<t<<'\t'<< '\t' << state.H << '\t'<< '\t' << state.V << '\t'<< '\t' << state.T << '\t'<< '\t' << state.P << endl;
            }
            ofile.close();
        }
        else {
            cout << "Unable to save file!"<<endl;
            exit(1);
        }
    }

    // Zapis wartości XYZ wraz z Ekin
    if (rp == "e") {
        ofstream ofile;
        if (b) {
            ofile.open(oFile,std::ios_base::app);
        }
        else {
            ofile.open(oFile);
        }

        if(ofile.is_open()){
            for (int i=0;i<N;i++)ofile<<state.r[i][0]<<'\t'<<state.r[i][1]<<'\t'<<state.r[i][2]<<'\t'<<state.Ekin<<endl;
            ofile<<"\n"<<endl;
            ofile.close();
        }
        else {
            cout << "Unable to save file!"<<endl;
            exit(1);
        }
    }


}
void readfile(Parameters &parameters,char *oFile){

    //  Odczyt wartosci parametrow z pliku
    ifstream ifile;
    ifile.open(oFile);
    if (ifile.is_open()) {
        ifile>>parameters.n>>parameters.m>>parameters.a>>parameters.f>>parameters.R>>parameters.e>>parameters.L>>parameters.T0>>parameters.tao>>parameters.So>>parameters.Sd>>parameters.Sout>>parameters.Sxyz;
    }
    else {
        cout << "Unable to read file!" << endl;
        exit(1);
    }
}
void count_position(State &state, Parameters &parameters,double *b0, double *b1, double *b2){

    //  Liczenie początkowych położen w krysztale
    int i = 0;
    for(int i0=0;i0<parameters.n;i0++){
        for(int i1=0;i1<parameters.n;i1++){
            for(int i2=0;i2<parameters.n;i2++){

                state.r[i][0] = (i0 - (parameters.n - 1.0) / 2.0) * b0[0] + (i1 - (parameters.n - 1.0) / 2.0) * b1[0] +
                                (i2 - (parameters.n - 1.0) / 2.0) * b2[0];
                state.r[i][1] = (i0 - (parameters.n - 1.0) / 2.0) * b0[1] + (i1 - (parameters.n - 1.0) / 2.0) * b1[1] +
                                (i2 - (parameters.n - 1.0) / 2.0) * b2[1];
                state.r[i][2] = (i0 - (parameters.n - 1.0) / 2.0) * b0[2] + (i1 - (parameters.n - 1.0) / 2.0) * b1[2] +
                                (i2 - (parameters.n - 1.0) / 2.0) * b2[2];
                i++;

            }
        }
    }
}
void count_pendulum(State &state, Parameters &parameters,int N){

    //  Liczenie początkowych pędów wraz z ich wyśrodkowaniem
    double ekinx = 0;
    double ekiny = 0;
    double ekinz = 0;
    double lambdax = 0;
    double lambday = 0;
    double lambdaz = 0;
    double Px = 0;
    double Py = 0;
    double Pz = 0;


    for(int j=0;j<N;j++) {

    //  Obliczanie pedow poczatkowych
    lambdax = drand48();
    lambday = drand48();
    lambdaz = drand48();

    ekinx = (-1.0) * 0.5 * parameters.kb * parameters.T0 * log(lambdax);
    ekiny = (-1.0) * 0.5 * parameters.kb * parameters.T0 * log(lambday);
    ekinz = (-1.0) * 0.5 * parameters.kb * parameters.T0 * log(lambdaz);

    lambdax = drand48();
    lambday = drand48();
    lambdaz = drand48();

    (lambdax > 0.5) ? lambdax = -1 : lambdax = 1;
    (lambday > 0.5) ? lambday = -1 : lambday = 1;
    (lambdaz > 0.5) ? lambdaz = -1 : lambdaz = 1;

    state.p[j][0] = lambdax * sqrt(2.0 * parameters.m * ekinx);
    state.p[j][1] = lambday * sqrt(2.0 * parameters.m * ekiny);
    state.p[j][2] = lambdaz * sqrt(2.0 * parameters.m * ekinz);

    //  Wysrodkowanie pedow

    Px += state.p[j][0];
    Py += state.p[j][1];
    Pz += state.p[j][2];

    state.p[j][0] -= Px * 1.0 / N;
    state.p[j][1] -= Py * 1.0 / N;
    state.p[j][2] -= Pz * 1.0 / N;
    }
}
void simulation(State &state, Parameters &parameters,int N){

    double ri = 0;
    double rij = 0;
    double riL = 0;
    double lri = 0;
    double R_rij = 0;
    double R_12 = 0;
    double R_6 = 0;
    double v_chwil = 0;

    //  Zerowanie wszystkich sił
    for(int i=0;i<N;i++){
        state.F[i][0] = 0.0;
        state.F[i][1] = 0.0;
        state.F[i][2] = 0.0;
    }
    //  Zerowanie wartości stanu kryształu
    state.V = 0.0;
    state.P = 0.0;
    state.H = 0.0;

    for(int i=0;i<N;i++){

        ri = pow((pow(state.r[i][0],2) + pow(state.r[i][1],2) + pow(state.r[i][2],2)),0.5);

        //  Obliczanie potencjału oraz ciśnienia wraz z ich akumulacją do stanu kryształu
        if (ri < parameters.L){
            state.V += 0;
            state.F[i][0] += 0;
            state.F[i][1] += 0;
            state.F[i][2] += 0;

        }
        else {
            riL = ri-parameters.L;
            lri = parameters.L-ri;

            state.V += 0.5 * parameters.f * pow(riL,2);
            state.F[i][0] += parameters.f * state.r[i][0] * lri / ri;
            state.F[i][1] += parameters.f * state.r[i][1] * lri / ri;
            state.F[i][2] += parameters.f * state.r[i][2] * lri / ri;
        }

        //  Obliczanie ciśnienia chwilowego
        state.P += sqrt(state.F[i][0] * state.F[i][0] + state.F[i][1] * state.F[i][1] + state.F[i][2] * state.F[i][2]) / (4.0 * parameters.Pi * parameters.L * parameters.L);

        for(int j=0;j<i;j++){

            rij = sqrt((pow((state.r[i][0] - state.r[j][0]),2) + pow((state.r[i][1] - state.r[j][1]),2) + pow((state.r[i][2] - state.r[j][2]),2)));
            R_rij = parameters.R / rij;
            R_12 = pow(R_rij, 12);
            R_6 = pow(R_rij, 6);

            v_chwil = parameters.e * (R_12 - 2*R_6);
            state.V += v_chwil;

            state.F[i][0] += 12 * (R_12-R_6) * (state.r[i][0] - state.r[j][0]) / pow(rij,2);
            state.F[i][1] += 12 * (R_12-R_6) * (state.r[i][1] - state.r[j][1]) / pow(rij,2);
            state.F[i][2] += 12 * (R_12-R_6) * (state.r[i][2] - state.r[j][2]) / pow(rij,2);

            state.F[j][0] -= 12 * (R_12-R_6) * (state.r[i][0] - state.r[j][0]) / pow(rij,2);
            state.F[j][1] -= 12 * (R_12-R_6) * (state.r[i][1] - state.r[j][1]) / pow(rij,2);
            state.F[j][2] -= 12 * (R_12-R_6) * (state.r[i][2] - state.r[j][2]) / pow(rij,2);

        }
        //  Obliczanie członu związanego z energia kinecztyczna i akumalcja do energi całkowitej
        state.H += (pow(state.p[i][0],2) + pow(state.p[i][1],2) + pow(state.p[i][2],2)) / (2 * parameters.m);
    }
    //  Akumulacja członu zwiazanego z potencjałem do energi całkowitej
    state.H += state.V;
}

int main(int argc, char *arg[]){

//WARUNKI POCZATKOWE
    Parameters parameters = Parameters();
    State state = State();
//  Wczytywanie wartości parametrow z pliku
    readfile(parameters,arg[1]);
    int N = pow(parameters.n,3);
    double krok_t = 0;

//  WEKTORY WLASNE
    double b0 [3]={parameters.a,0.0,0.0};
    double b1 [3]={parameters.a / (2.0),parameters.a * sqrt(3.0) / (2.0),0.0};
    double b2 [3]={parameters.a / (2.0),parameters.a * sqrt(3.0) / (6.0),parameters.a * sqrt((2.0) / (3.0))};

//  OBLICZANIE POLOZEN POCZATKOWYCH
    count_position(state,parameters,b0,b1,b2);

//  OBLICZANIE PEDOW
    count_pendulum(state,parameters,N);

//  ZAPIS DO PLIKOW TEKSTOWYCH
    savefile(state,arg[3],N,"r",0,krok_t);
    savefile(state,arg[2],N,"s",0,krok_t);
    //savefile(state,arg[4],N,"e",0,krok_t);


// ALGORYTM 2

    simulation(state, parameters,N);
    cout<<"V poczatkowe: "<<state.V<<endl;
    cout<<"P poczatkowe: "<<state.P<<endl;

//ALGORYTM 1

//  Zerowanie wartości uśrednionych
    state.H_srd=0.0;
    state.P_srd=0.0;
    state.T_srd=0.0;

    for (int s = 0; s < (parameters.So+parameters.Sd);s++) {

        for (int i = 0; i < N; i++) {

            //  Modyfikacja pędów
            state.p[i][0] += 0.5 * state.F[i][0] * parameters.tao;
            state.p[i][1] += 0.5 * state.F[i][1] * parameters.tao;
            state.p[i][2] += 0.5 * state.F[i][2] * parameters.tao;

            //  Modyfikacja położeń
            state.r[i][0] += (1.0 / parameters.m) * state.p[i][0] * parameters.tao;
            state.r[i][1] += (1.0 / parameters.m) * state.p[i][1] * parameters.tao;
            state.r[i][2] += (1.0 / parameters.m) * state.p[i][2] * parameters.tao;

        }

        //Obliczanie nowych sil
        simulation(state, parameters, N);

        state.Ekin = 0.0;

        for (int i = 0; i < N; i++) {
            //  Modyfikacja pędów
            state.p[i][0] += 0.5 * state.F[i][0] * parameters.tao;
            state.p[i][1] += 0.5 * state.F[i][1] * parameters.tao;
            state.p[i][2] += 0.5 * state.F[i][2] * parameters.tao;
            //  Obliczanie energi kinetycznej chwilowej
            state.Ekin += (pow(state.p[i][0],2) + pow(state.p[i][1],2) + pow(state.p[i][2],2)) / (2 * parameters.m);

        }
        //  Obliczanie temperatury chwilowej
        state.T = 2.0 / (3.0 * N * parameters.kb) * state.Ekin;
        //  Obliczanie całkowitej energii chwilowej
        state.H = state.Ekin + state.V;

        cout<<"P chwilowe: "<<state.P<<endl;
        cout<<"T chwilowe: "<<state.T<<endl;
        cout<<"Ekin chwilowa : "<<state.Ekin<<endl;
        cout<<"V chwilowy: "<<state.V<<endl;
        cout<<"H chwilowy: "<<state.H<<endl;

        savefile(state,arg[3],N,"r",1,krok_t);
        // Jezeli s jest wielkrotnoscia Sout to zapis stamu krysztalu do pliku
        if(s%(int)parameters.Sout==0){

            savefile(state,arg[2],N,"s",1,krok_t);
        }
        // Jezeli s jest wielkrotnoscia Sxyz to zapis wartosci xyz wraz z Ekin do pliku
        if(s%(int)parameters.Sxyz==0){

            //savefile(state,arg[4],N,"e",1,krok_t);
        }
        // Jezeli s jest wieksze niz So to akumulacja wartości uśrednianych
        if(s>=parameters.So) {
            state.T_srd += state.T;
            state.P_srd += state.P;
            state.H_srd += state.H;
        }
        krok_t += parameters.tao;

    }
    state.T_srd = state.T_srd / parameters.Sd;
    state.P_srd = state.P_srd / parameters.Sd;
    state.H_srd = state.H_srd / parameters.Sd;

    cout<<"Temperatura srednia:"<<state.T_srd<<endl;
    cout<<"Cisnienie srednie:"<<state.P_srd<<endl;
    cout<<"Energia srednia:"<<state.H_srd<<endl;

    return 0;
}