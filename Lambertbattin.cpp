//
// Created by josorian on 23/03/2023.
//

#include <cstring>
#include "Lambertbattin.h"
#include "seebattk.h"
#include "seebatt.h"
void
lambartbattin(vector<double> ro, vector<double> r, const char *dm, double Dtsec, vector<double> *vo, vector<double> *v) {
    double small=0.000001;
    double mu=3.986004418e14;
    double y1=0;
    double pi=3.14;
    double magr= norm(r,3);
    double magro= norm(ro,3);
    double  CosDeltaNu= dot(ro,r,3,3);
    vector<double> rcrossr= cross(ro, r);
    double magrcrossr = norm(rcrossr,3);
    double SinDeltaNu;
    double DNu;
    double Ror;
    double eps;
    double tan2w;
    double rp;
    double L;
    double m;
    double x;
    double xn;
    double chord;
    double s;
    double  lim1;
    int Loops;
    double tempx;
    double Denom;
    double h1;
    double h2;
    double b;
    double u;
    double k2;
    double y;
    double a;
    if (strcmp(dm, reinterpret_cast<const char *>('pro')) == 0){
        SinDeltaNu=magrcrossr/(magro*magr);
    }
    else{
        SinDeltaNu=-magrcrossr/(magro*magr);
    }
    DNu= atan2(CosDeltaNu,SinDeltaNu);
    if (DNu<0.0){
        DNu=2.0*M_PI+DNu;
    }
    Ror=magr/magro;
    eps=Ror-1.0;
    tan2w=0.25*eps*eps/(sqrt(Ror)+Ror*(2.0+ sqrt(Ror)));
    rp= sqrt(magro*magr)*(pow(cos(DNu*0.25),2)+tan2w);
    if(DNu<M_PI){
        L=(pow(sin(DNu*0.25),2)+tan2w)/((pow(sin(DNu*0.25),2))+tan2w+ cos(DNu*0.5));
    }
    else{
        L=(pow(cos(DNu*0.25),2)+tan2w- cos(DNu*0.5))/((pow(cos(DNu*0.25),2))+tan2w);
    }
    m=mu*Dtsec*Dtsec /(8.0*rp*rp*rp);
    x=10.0;
    xn=L;
    chord=sqrt( magro*magro + magr*magr - 2.0*magro*magr*cos( DNu ) );
    s    = ( magro + magr + chord )*0.5;
    lim1 = sqrt(m/L);
    Loops=1;
    while(1){
        x=xn;
        tempx= seebatt(x);
        Denom= 1.0 / ( (1.0+2.0*x+L) * (4.0*x + tempx*(3.0+x) ) );
        h1   = pow( L+x ,2) * ( 1.0+ 3.0*x + tempx )*Denom;
        h2   = m*( x - L + tempx )*Denom;
        b = 0.25*27.0*h2 / (pow(1.0+h1,3) );
        if(b<-1.0){
            xn=1.0-2.0*L;
        }
        else{
            if(y1>lim1){
                xn=xn*(lim1/y1);
            }
            else{
                u = 0.5*b / ( 1.0 + sqrt( 1.0 + b ) );
                k2 = seebattk(u);
                y = ( ( 1.0+h1 ) / 3.0 )*( 2.0 + sqrt( 1.0+b )/( 1.0+2.0*u*k2*k2 ) );
                xn= sqrt( pow( (1.0-L)*0.5 ,2) + m/(y*y) ) - ( 1.0+L )*0.5;
            }

        }
        Loops = Loops + 1;
        y1=sqrt(m/((L+x)*(1.0+x)) );
        if((fabs(xn-x) < small)&&(Loops>30)){
            break;
        }
    }
    a=  mu*Dtsec*Dtsec / (16.0*rp*rp*xn*y*y );

}
