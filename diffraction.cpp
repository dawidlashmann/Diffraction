#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
//////////////////////////////////////////////////////////////////////////////
// ZMIENNE DO MODYFIKACJI PRZEZ UZYTKOWNIKA
#define Y_MIN_1 -1.5
#define Y_MAX_1 1.5  //szeroko�� 1 szczeliny
#define Y_LEFT -100000.0
#define Y_RIGHT 100000.0 //szeroko�� rysowania wynikow na obrazie (lepiej dawac wieksze niz szcelina)
#define N 1000
#define N_INT 200
#define LAMBDA 1.0 //d�ugo�c fali
#define AMPLITUDE 10.0 //amplituda
#define D 200000.0 //odleg�o�� od ekranu
#define L 10.0 //odleglosc pomiedzy szczelinami
#define Y_MIN_2 -2.5
#define Y_MAX_2 2.5  // szerokosc 2 szczeliny
//////////////////////////////////////////////////////////////////////////////
// ZMIENNE PROGRAMOWE
double k = 2.0 * 3.1415 / LAMBDA;
double h = (Y_RIGHT - Y_LEFT)/(N-1);

double y[N];
double complex_amplitude[N];

fstream plik;
//////////////////////////////////////////////////////////////////////////////
// FUNKCJE
double Real(double x, double y){
  double r = sqrt( D*D+(y-x)*(y-x) );
  double amplitude = AMPLITUDE * cos(k*r)/sqrt(r);

  return amplitude;
}

double Imag(double x, double y){
  double r = sqrt( D*D+(y-x)*(y-x) );
  double amplitude = AMPLITUDE * sin(k*r)/sqrt(r);

  return amplitude;
}

double Simpson(double x_min, double x_max, double y, double n,  double (*Function)(double, double)){
  double h_ = (x_max - x_min) / (2.0 * n );
  double simpson = 0.0;
  double x = x_min+h_;

  for(int i=0; i<n; i++){
    simpson = simpson + h_*( Function(x+h_, y) + 4.0*Function(x, y) + Function(x-h_, y) ) / 3.0;
    x = x+2.0*h_;
  }

  return simpson;
}

double Simpson(double len1, double len2, double l, double y, double n,  double (*Function)(double, double), int number_of_slits){
  double h_1 = len1 / (2.0 * n );
  double h_2 = len2 / (2.0 * n );
  double simpson = 0.0; 
  
  bool even = 1;
  if(number_of_slits % 2 == 1){
    double x = -len1 / 2 + h_1;
    while(x <= len1 / 2){
      simpson += h_1*( Function(x+h_1, y) + 4.0*Function(x, y) + Function(x-h_1, y) ) / 3.0;
      x += 2.0*h_1;
    }
    even = 0;
  }

  for(int i = -number_of_slits / 2; i < 0; i++){
    double anchor;
    if(even){
      anchor = ((i + 1) * l) + (i * len1) - (l / 2);
    }else{
      anchor = (i * l) + (i * len1) - (len1 / 2);
    }

    double x = anchor + h_1;
    while(x <= anchor + len1){
      simpson += h_1*( Function(x+h_1, y) + 4.0*Function(x, y) + Function(x-h_1, y) ) / 3.0;
      x += 2.0*h_1;
    }
  }

  for(int i = 0; i < number_of_slits / 2; i++){
    double anchor;
    if(even){
      anchor = (i * l) + (i * len1) + (l / 2);
    }else{
      anchor = ((i + 1) * l) + (i * len1) + (len1 / 2);
    }

    double x = anchor + h_1;
    while(x <= anchor + len2){
      simpson += h_1*( Function(x+h_1, y) + 4.0*Function(x, y) + Function(x-h_1, y) ) / 3.0;
      x += 2.0*h_1;
    }
  }

  return simpson;
}

/////////////////////////////////////////////////////////////////////////////
// FUNKCJA GLOWNA

int main(){
  double real, imag;
  double maximum = 0.0;
  for(int i=0; i<N; i++){ y[i] = i*h + Y_LEFT; }

//single slit
  for(int i=0; i<N; i++){
    real = Simpson(Y_MIN_1, Y_MAX_1, y[i], N_INT, Real);
    imag = Simpson(Y_MIN_1, Y_MAX_1, y[i], N_INT, Imag);

    complex_amplitude[i] = real*real + imag*imag;

    maximum = (complex_amplitude[i] > maximum) ? complex_amplitude[i] : maximum;
  }

  plik.open("wave.data", ios::out);

  for(int i=0; i<N; i++){
    plik << y[i] << "   "
	 << complex_amplitude[i] / maximum
	 << endl;
  }
////////////////////////////////////////////////////////////////

fstream diffraction_grating("diff_grating.data", ios::out);
// diffraction grate
double diff_grat_amplitude[N];
const int number_of_slits = 10;

for(int i = 0; i < N; i++){
    real = Simpson(Y_MAX_1 - Y_MIN_1, Y_MAX_1 - Y_MIN_1, L, y[i], N_INT, Real, number_of_slits);
    imag = Simpson(Y_MAX_1 - Y_MIN_1, Y_MAX_1 - Y_MIN_1, L, y[i], N_INT, Imag, number_of_slits);

    diff_grat_amplitude[i] = real*real + imag*imag;
    
    maximum = (diff_grat_amplitude[i] > maximum) ? diff_grat_amplitude[i] : maximum;
  }


  for(int i=0; i<N; i++){
      diffraction_grating << y[i] << "   "
    << diff_grat_amplitude[i] / maximum
    << endl;
  }


//////////////////////////////////////////////////////////////
  plik.close();
  diffraction_grating.close();

  return 0;
}
