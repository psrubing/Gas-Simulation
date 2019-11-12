#ifndef ARGON_ARGON_H
#define ARGON_ARGON_H
#define MAX_ATOMS 5000

using namespace std;

struct Parameters{

    //parametery
    int n;
    double m;
    double a;
    double f;
    double R;
    double e;
    double L;
    double T0;

    //parametry dynamiki
    double tao;
    int So;
    int Sd;

    //parametry wyjscia
    int Sout;
    int Sxyz;

    const double kb = 0.00861733;
    const double Pi = 3.14;
};

struct State{

    //  parametry stanu układu
    double V = 0.0;
    double P = 0.0;
    double H = 0.0;
    double T = 0.0;
    double Ekin = 0.0;

    double T_srd = 0.0;
    double P_srd = 0.0;
    double H_srd = 0.0;

    //tablice
    double r[MAX_ATOMS][3];
    double p[MAX_ATOMS][3];
    double F[MAX_ATOMS][3];

};

void savefile(State&,char*,int,string,bool,double);
void readfile(Parameters&,char*);
void simulation(State&, Parameters&,int);
void count_position(State&, Parameters&, double*, double*, double*);
void count_pendulum(State&, Parameters&, int);



#endif //ARGON_ARGON_H
