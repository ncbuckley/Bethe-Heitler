#include "file2.h"
#include "EnergyLookup.h"
#include <math.h>
#include <TVector3.h>



double dSigmaBorn(double x, double phi1, double phi2, double theta1, double theta2,double E0){
  double alpha = 0.0072973525664;
  int Z = 1;
  int xi3 = 1;
  double me = 5.109989461e-4;

  
  TVector3 u1(sin(theta1)*cos(phi1),sin(theta1)*sin(phi1),cos(theta1));
  TVector3 u2(sin(theta2)*cos(phi2),sin(theta2)*sin(phi2),cos(theta2));

  TVector3 k1(u1.X()*(sqrt(pow(x,2)*pow(E0,2)-pow(me,2))),u1.Y()*(sqrt(pow(x,2)*pow(E0,2)-pow(me,2))),u1.Z()*(sqrt(pow(x,2)*pow(E0,2)-pow(me,2))));
  TVector3 k2(u2.X()*(sqrt(pow((1-x),2)*pow(E0,2)-pow(me,2))), u2.Y()*(sqrt(pow((1-x),2)*pow(E0,2)-pow(me,2))), u2.Z()*(sqrt(pow((1-x),2)*pow(E0,2)-pow(me,2))));

  TVector3 p1(k1.X(),k1.Y(),0);
  TVector3 p2(k2.X(),k1.Y(),0);

  TVector3 q(k1.X()+ k2.X(), k1.Y()+ k2.Y(),0);


  double c1 = p1.X()*p1.X()+p1.Y()*p1.Y() + pow(me,2);
  double c2 = p2.X()*p2.X()+p2.Y()*p2.Y() + pow(me,2);

  double S = 1/c1 - 1/c2;
  TVector3 T(p1.X()/c1+p2.X()/c2,p1.Y()/c1 + p2.Y()/c2,0);
  
  double p12(p1.X()*p1.X()+p1.Y()*p1.Y());
  double p22(p2.X()*p2.X()+p2.Y()*p2.Y());
  double T2(T.X()*T.X()+T.Y()*T.Y());
  double q2(q.X()*q.X()+q.Y()*q.Y());

  double dsig1 = ((2 * pow(alpha,3) * pow(Z,2) * pow(E0,4) * pow(x,2) * pow((1-x),2))/(pow(M_PI,2) * pow(q2,2)));
  
  double dsig2 = pow(me,2) * pow(S,2) + (pow(x,2) + pow((1-x),2)) * (T2) - 2*x*(1-x)*xi3 *( (p12/(pow(c1,2)))*cos(2*phi1) + (p22/(pow(c2,2)))*cos(2*phi1) + 2*(p1.Mag()/c1)*(p2.Mag()/c2)*cos(phi1+phi2));
  double dsig = dsig1 * dsig2;


  const double alpha1 = 0.1;
  const double alpha2 = 0.55;
  const double alpha3 = 0.35;
  const double b1 = 6.0;
  const double b2 = 1.2;
  const double b3 = 0.3;
  const double mu1 = (me * pow(Z,1/3)) * b1;
  const double mu2 = (me * pow(Z,1/3)) * b2;
  const double mu3 = (me * pow(Z,1/3)) * b3;
  double AtomicFF = 1 - q.Mag2() * ( alpha1/( pow(mu1,2) + q.Mag2() )+ alpha2/( pow(mu2,2) + q.Mag2() )+ alpha3/( pow(mu3,2) + q.Mag2() ));
  double FullOnFF = abs( pow( 1 + (q.Mag2() /0.71) , -2) - AtomicFF );

  double hbarc = 3.893793e5; //nbn GeV^2
  double adj_dsig =dsig * pow(FullOnFF,2) * abs(pow(theta1,2) * sin(theta1)) * abs(pow(theta2,2) * sin(theta2)) * hbarc * EnergyLookup(E0);
  
  return(adj_dsig);
		  

}
