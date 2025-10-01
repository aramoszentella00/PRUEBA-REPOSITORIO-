#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

const double Rg = 0.518;
const double Pc = 4600.0;
const double Tc = 191.0;
const double T = 233.0;
const double P = 65000.0;

double calcA(){
    return 0.42748 * (Rg*Rg*pow(Tc,2.5)) / Pc;
}
double calcB(){
    return 0.08664 * (Rg*Tc) / Pc;
}

double fun(double v,double A,double B){
    return P - (Rg*T)/(v-B) + A/(sqrt(T)*v*(v+B));
}

double dFun(double v,double A,double B){
    return (Rg*T)/((v-B)*(v-B)) - A*(2*v+B)/(sqrt(T)*v*v*(v+B)*(v+B));
}

double newton(double A,double B,double x0){
    double x=x0;
    for(int k=0;k<150;k++){
        double y=fun(x,A,B);
        if(fabs(y)<1e-8) return x;
        double dy=dFun(x,A,B);
        if(fabs(dy)<1e-14) break;
        x=x-y/dy;
        if(x<=B) x=B+1e-6;
    }
    return x;
}

int main(){
    cout<<setprecision(6)<<fixed;
    cout<<"--- Calculo usando Redlich-Kwong ---\n";

    double A=calcA();
    double B=calcB();

    cout<<"Valores de constantes:\n";
    cout<<"A="<<A<<"  B="<<B<<"\n";

    double vGasIdeal=(Rg*T)/P;
    double v0=vGasIdeal;
    if(v0<=B) v0=B+0.00001;

    cout<<"\nAproximacion gas ideal: "<<vGasIdeal<<" m3/kg\n";
    cout<<"Usando v0="<<v0<<" como arranque\n";

    double vSol=newton(A,B,v0);

    cout<<"\nResultado aproximado:\n";
    cout<<"Volumen especifico="<<setprecision(9)<<vSol<<" m3/kg\n";

    double tanque=3.0;
    double m=tanque/vSol;
    cout<<"Masa contenida en "<<tanque<<" m3: "<<setprecision(6)<<m<<" kg\n";

    cout<<"(Chequeo f(v)="<<fun(vSol,A,B)<<")\n";

    return 0;
}
