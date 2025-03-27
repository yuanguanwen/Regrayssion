// preliminar version of the code, z* = 1090
// da capire segmentation falut con z = 0 (int)

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>

#define N 200
// 200*Z --> the integral is insensitive to the actual value of z* provided z* >> 10 (500 is ok)
#define N1 100000

// LCDM: neglect radiation and assume spatially flat universe
 
#define H0 72
#define H0_s 2.33e-18		// hubble consant in seconds
#define Ocdm 0.212		// DM density parameter
#define Om 0.2552		// Matter density parameter
#define z_rec 1090
#define Z 500			// maximum redshift
#define M 			// BH mass

using namespace std;

int import_data(int);
void flux_rescaling(double z);
double H(double z);
void integral(int);

string mass;
double de;
double flux[N], energy[N], temp_energy[N];
double flux1[N1], energy1[N1], temp_flux[N1];
int names[3] = {3, 15, 45};

int main(){
	
	mass = "1e15";
	
	//cycle over the reg. parameter
	for(int n = 0; n < 3; n++){
		
		import_data(n);
		de = energy[N-1]/Z/N;
		
		for(int i = 0; i < N1; i++){
			energy1[i] = de*i;
			flux1[i] = 0;
		}
			
		integral(n);
		
	}
	return 0;
}

int import_data(int n){
	
	int j = 0;
	string filename;
	filename = "a0" + to_string(names[n]) + "rh_" + mass + ".txt";
	
	ifstream file(filename);
	if (!file.is_open()) {
		cerr << "Errore nell'apertura del file." << endl;
	        return 1;	
	}
	while ( j < N && file >> energy[j] >> flux[j]) {
        	j++;
    	}
    	file.close();
       	return 0;
}

void flux_rescaling(double z){
	double m, q;
	int j = 1, i = 0;
	
	for(int ii = 0; ii < N1; ii++){
		flux1[ii] = 0;
	}
	// rescale to redshift z
	for(int ii = 0; ii < N; ii++){
		temp_energy[ii] = energy[ii]/(1+z);
		//cout << temp_energy[ii] << endl;
	}
	
	
	while(energy1[i] < temp_energy[0]){
		flux1[i] = 0;
		i++;
		
	}
	 
	while(energy1[i] < temp_energy[N-1] && i < N1){
		if(energy1[i] < temp_energy[j]){
			m = (flux[j] - flux[j-1]) / (temp_energy[j] - temp_energy[j-1]);
			q = flux[j-1];
			flux1[i] = m*(energy1[i] - temp_energy[j-1]) + flux[j-1];
			i++;
		}
		else if (energy1[i] == temp_energy[j]){
			flux1[i] = flux[j];
			j++;
			i++;
		}
		else {
			m = (flux[j+1] - flux[j]) / (temp_energy[j+1] - temp_energy[j]);
			flux1[i] = m*(energy1[i]-temp_energy[j-1]) + flux[j-1];
			i++;
			j++;			
		} 
	}
	
	//debugging
	/*
	ofstream file_resc("rescaling.txt");
	for(int ii = 0; ii < N1; ii++){
		file_resc << energy1[ii] << "    " << flux1[ii] << endl;
	}*/
}

double H(double z){
	return H0_s*sqrt(Om*pow((1+z),3.0) + (1.0 - Om));
}

void integral(int n){

	double z;
	for(int i = 0; i < N1; i++){
		temp_flux[i] = 0;
	}
	for(int i = 0; i < 3270; i++){
		
		z = (double)i/3;
		//cout << z << "    " << H(z) << endl;
		
		flux_rescaling(z);
		for(int j = 0; j < N1; j++){
			temp_flux[j] += 1/H(z) * flux1[j] ; 
		}

	}
	// "./yourdirectory/filename"
	string filename2 = "./" + mass + "/Bsum_"+ mass + "_a0" + to_string(names[n]) + "rh.txt";
	ofstream file2(filename2);
	for(int jj = 0; jj < N1; jj+= 300){
		file2 << energy1[jj] << "    " << energy1[jj]*temp_flux[jj]/3 << endl;
	}
	file2.close();
}

