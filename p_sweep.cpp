#include "pqd.h"

int main(){

	// Make sure your inputs are correct.
	// If you want a single p, bind both p_start and p_end to it and divisions to 1. 

	int grid_size;
	float p_start;
	float p_end;
	int divisions;
	float feedback_strength;
	float death_probability;
	int initializations;
	int realizations;
	int lag;
	float how_many_occupied_initially;
	int updates_per_site;
	int collect_frames;

	cout << "Enter grid size: ";
    cin >> grid_size;

    cout << "Enter starting p: ";
    cin >> p_start;

    cout << "Enter final p: ";
    cin >> p_end;

    cout << "Enter number of divisions: ";
    cin >> divisions;

    cout << "Enter feedback strength q: ";
    cin >> feedback_strength;

    cout << "Enter death probability d: ";
    cin >> death_probability;

    cout << "How many initializations: ";
    cin >> initializations;

    cout << "Enter number of realizations for each initialization: ";
    cin >> realizations;

    cout << "Enter lag between realizations: ";
    cin >> lag;

    cout << "How many occupied initially: ";
    cin >> how_many_occupied_initially;


    cout << "Enter time to equilibrium: ";
    cin >> updates_per_site;

    cout << "Do you want to collect frames? Answer with 1 (yes) or 0 (no): ";
    cin >> collect_frames;

    auto start = high_resolution_clock::now(); 

    float densities[divisions][initializations*realizations];

    vector<double> birth_probabilities = linspace(p_start, p_end, divisions);

    for (int i=0; i < divisions; i++){

		//cout << i << endl; 

		float these_densities[initializations*realizations];

		equilibrium_density_pqd(grid_size, birth_probabilities[i], feedback_strength, death_probability, realizations, lag, initializations, how_many_occupied_initially, these_densities, updates_per_site, collect_frames);

		for (int j=0; j < initializations*realizations; j++){
			densities[i][j] = these_densities[j];
		}
	}

	ofstream outdata;

	outdata.open("dump/"+std::to_string(grid_size)+"_"+std::to_string(p_start)+"_"+std::to_string(p_end)+".txt");
	if( !outdata ) {
		cerr << "File error, try again." << endl;
		exit(1);
	}
	for (int i=0; i< divisions; i++){
		outdata << birth_probabilities[i];
		for (int j=0; j<initializations*realizations; j++){
			outdata << " " << densities[i][j];
		}
		outdata << endl;
	}
	outdata.close();

	auto stop = high_resolution_clock::now(); 
	auto duration = duration_cast<seconds>(stop - start); 
	cout << "CPU Time: " << duration.count() << " seconds" << endl;
}
