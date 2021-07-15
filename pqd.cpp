#include "pqd.h"

//---------------------Some Primitives -------------------------------------------//

void seed_twister(int seed){
	rng.seed(seed);
}

double random_real(double initial, double last) {

    std::uniform_real_distribution<double> distribution(initial, last);
    return distribution(rng);  // Use rng as a generator
}

int random_int(int initial, int last) {

    std::uniform_int_distribution<int> distribution(initial, last);
    return distribution(rng);  // Use rng as a generator
}

float mean_of_array(float array[],int size){
	
	float sum = 0.0;
	
	for (int i =0; i<size; ++i){
		sum += array[i];
	}
	return sum/(float)size;
}

float standard_deviation_of_array(float array[],int size){
	
	float mean = mean_of_array(array,size);
	float sum = 0.0;
	
	for (int i = 0; i < size; ++i){
		sum += (array[i]-mean)*(array[i]-mean);
	}
	
	float variance = sum/(float)size;

	return sqrt(variance);
}

float mean_of_vector(vector<float> array,int size){
	
	float sum = 0.0;
	
	for (int i =0; i<size; ++i){
		sum += array[i];
	}
	return sum/(float)size;
}

void random_frame(int frame[], int grid_size) {      

	for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {

            frame[i*grid_size + j] = random_int(0, 1);
        		
        }
    }
}

void random_frame_of_density(float density, int frame[], int grid_size) {      

	for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
        	
        	if (random_real(0,1)<=density){
        		frame[i*grid_size + j] = 1;
        	}
        	else{
        		frame[i*grid_size + j] = 0;
        	}
        }
    }
}

void randomly_occupied_frame(int how_many_occupied, int frame[], int grid_size) {      

	zeros(frame, grid_size);

	for (int i = 0; i < how_many_occupied; ++i) {

		int site = random_int(0, grid_size*grid_size);
		frame[site] = 1;
    }
}

void print_frame(int frame[], int grid_size){
	for (int i = 0; i < grid_size; i++){
		for (int j =0; j < grid_size; j++){
			cout << frame[i*grid_size+j];
		}
		cout << endl;
	}
}

void print_vector(std::vector<double> vec){
  
  std::cout << "size: " << vec.size() << std::endl;
  for (double d : vec)
    std::cout << d << " ";
  std::cout << std::endl;
}

void save_frame(int frame[], int grid_size, string name){
	
	//Credits: Sumithra Sankaran

	ofstream output;  
	string filename = "dump/" + name + ".txt";
	output.open("dump/" + name + ".txt");

	for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            output << frame[i*grid_size + j];
        }
        output << '\n';
    }
	output.close();
}

void zeros(int frame[], int grid_size) {      

	for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            frame[i*grid_size + j] = 0;
        }
    }
}

float calculate_density(int frame[], int grid_size){

	float occupancy = 0;

	for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            occupancy += frame[i*grid_size + j];
        }
    }

    float density = occupancy/(grid_size*grid_size);
	return density;
}

// -------------------------------pqd model---------------------------------------//

coordinates select_neighbor_of_site(coordinates site, int grid_size){

	// Selects at random one of the four neighbours in von-neumann radius of 1 with periodic boundary conditions 

	int x = site.x;
	int y = site.y;

	coordinates neighbours[4];

	neighbours[0].x = (x-1)%grid_size, neighbours[0].y = y; //top
	neighbours[1].x = (x+1)%grid_size, neighbours[1].y = y; //bottom
	neighbours[2].x = x, neighbours[2].y = (y-1)%grid_size; //left
	neighbours[3].x = x, neighbours[3].y = (y+1)%grid_size; //right
 
	int who = random_int(0,3);

	coordinates neighbour = neighbours[who];

	if (neighbour.x < 0){
		neighbour.x = grid_size-1; //correction if selected site is top of the first row
	}

	if (neighbour.y < 0){
		neighbour.y = grid_size-1; //correction if selected site is left of the first column
	}

	return neighbour;
}

coordinates select_neighbor_of_pair(coordinates site, coordinates neighbour, int grid_size){

	// Selects at random one of the six neighbours of a pair of sites. Here too we have periodic boundary conditions. 
	
	int diff_x = site.x - neighbour.x;
	int diff_y = site.y - neighbour.y; 

	coordinates neighbours_of_pair[6];

	if ((diff_x ==1||diff_x==-(grid_size-1)) && diff_y ==0){ // neighbour is to the left of focal site 

		neighbours_of_pair[0].x = (site.x+1)%grid_size, neighbours_of_pair[0].y = site.y;
		neighbours_of_pair[1].x = (neighbour.x-1)%grid_size, neighbours_of_pair[1].y = site.y;
		neighbours_of_pair[2].x = site.x, neighbours_of_pair[2].y = (site.y+1)%grid_size;
		neighbours_of_pair[3].x = site.x, neighbours_of_pair[3].y = (site.y-1)%grid_size;
		neighbours_of_pair[4].x = neighbour.x, neighbours_of_pair[4].y = (neighbour.y+1)%grid_size;
		neighbours_of_pair[5].x = neighbour.x, neighbours_of_pair[5].y = (neighbour.y-1)%grid_size;
	}
	if ((diff_x ==-1 ||diff_x==(grid_size-1)) && diff_y ==0){ // neighbour is to the right of focal site 

		neighbours_of_pair[0].x = (site.x-1)%grid_size, neighbours_of_pair[0].y = site.y;
		neighbours_of_pair[1].x = (neighbour.x+1)%grid_size, neighbours_of_pair[1].y = site.y;
		neighbours_of_pair[2].x = site.x, neighbours_of_pair[2].y = (site.y+1)%grid_size;
		neighbours_of_pair[3].x = site.x, neighbours_of_pair[3].y = (site.y-1)%grid_size;
		neighbours_of_pair[4].x = neighbour.x, neighbours_of_pair[4].y = (neighbour.y+1)%grid_size;
		neighbours_of_pair[5].x = neighbour.x, neighbours_of_pair[5].y = (neighbour.y-1)%grid_size;
	}
	if (diff_x ==0 && (diff_y ==1 ||diff_y==-(grid_size-1))){ // neighbour is below the focal site 

		neighbours_of_pair[0].x = site.x, neighbours_of_pair[0].y = (site.y+1)%grid_size;
		neighbours_of_pair[1].x = neighbour.x, neighbours_of_pair[1].y = (neighbour.y-1)%grid_size;
		neighbours_of_pair[2].x = (site.x+1)%grid_size, neighbours_of_pair[2].y = site.y;
		neighbours_of_pair[3].x = (site.x-1)%grid_size, neighbours_of_pair[3].y = site.y;
		neighbours_of_pair[4].x = (neighbour.x+1)%grid_size, neighbours_of_pair[4].y = neighbour.y;
		neighbours_of_pair[5].x = (neighbour.x-1)%grid_size, neighbours_of_pair[5].y = neighbour.y;
	}
	if (diff_x ==0 && (diff_y ==-1 ||diff_y==(grid_size-1))){ // neighbour is above the focal site 

		neighbours_of_pair[0].x = site.x, neighbours_of_pair[0].y = (site.y-1)%grid_size;
		neighbours_of_pair[1].x = neighbour.x, neighbours_of_pair[1].y = (neighbour.y+1)%grid_size;
		neighbours_of_pair[2].x = (site.x+1)%grid_size, neighbours_of_pair[2].y = site.y;
		neighbours_of_pair[3].x = (site.x-1)%grid_size, neighbours_of_pair[3].y = site.y;
		neighbours_of_pair[4].x = (neighbour.x+1)%grid_size, neighbours_of_pair[4].y = neighbour.y;
		neighbours_of_pair[5].x = (neighbour.x-1)%grid_size, neighbours_of_pair[5].y = neighbour.y;
	}

	int who = random_int(0,5);

	coordinates neighbour_of_pair = neighbours_of_pair[who];

	if (neighbour_of_pair.x < 0){
		neighbour_of_pair.x = grid_size-1; //correction if selected site is the left of first column 
	}

	if (neighbour_of_pair.y < 0){
		neighbour_of_pair.y = grid_size-1; //correction if selected site is the top of first row
	}

	return neighbour_of_pair;
}

void pqd_update(int frame[], int grid_size, float birth_probability,  float feedback_strength, float death_probability){

	coordinates site; 
	site.x = random_int(0, grid_size-1);
	site.y = random_int(0, grid_size-1);
	double chance;
	double another_chance;

	if (frame[site.x*grid_size+site.y]==1){ // selected site is occupied 

		coordinates neighbour = select_neighbor_of_site(site,grid_size);

		if (frame[neighbour.x*grid_size+neighbour.y]==0){ // randomly selected neighbour unoccupied

			chance = random_real(0, 1);

			if (chance < birth_probability){ 
				frame[neighbour.x*grid_size+neighbour.y] = 1; // birth at empty neighbour with probability p
			}

			chance = random_real(0, 1);

			if (chance < death_probability){
				frame[site.x*grid_size+site.y]=0; // death of selected site with probability d
			}

		}
		else{ // randomly selected neighbour occupied

			coordinates neighbour_of_pair = select_neighbor_of_pair(site, neighbour, grid_size);

			chance = random_real(0,1);
			another_chance = random_real(0,1);

			if (chance < feedback_strength){

				frame[neighbour_of_pair.x*grid_size+neighbour_of_pair.y] = 1;
			}
			else if (another_chance < death_probability){
				frame[site.x*grid_size+site.y]=0; // death of selected site with probability d(1-q)
			}
		}
	}
}

void simulate_pqd(int frame[], int grid_size, float birth_probability, float feedback_strength, float death_probability, int updates_per_site){

	// Takes a frame and simulates pqd for the specified number of updates per site with the specified parameters

	while (updates_per_site > 0) {

		for (int i = 0; i < grid_size*grid_size; i++){ // loop performs on average a single update per site

			pqd_update(frame, grid_size, birth_probability, feedback_strength, death_probability); 
		}
		updates_per_site -= 1;
	}
}

void equilibrium_density_pqd(int grid_size, float birth_probability, float feedback_strength, float death_probability, int realizations, int lag, int initializations, float how_many_occupied_initially, float densities[], int updates_per_site, int collect_frames){
	
	// Simulates pqd with specified parameters and populates the array densities. Capture the frames from which vegetation cover is calculated by binding collect_frames to 1.

    #pragma omp parallel for
    for (int i=0; i < initializations; i++){
    	
    	int frame[grid_size*grid_size];
    	randomly_occupied_frame(how_many_occupied_initially,frame,grid_size);
    	simulate_pqd(frame, grid_size, birth_probability,  feedback_strength, death_probability, updates_per_site);

    	for (int j=0; j<realizations; j++) {
    	
    		simulate_pqd(frame, grid_size, birth_probability, feedback_strength, death_probability, lag);
    		densities[i*realizations+j] = calculate_density(frame,grid_size);

	    	if (collect_frames ==1){ // conditional for collecting frames 

	    		ofstream outdata;  

				stringstream p, q, d, census, time_to_equilibrium, init;
				p << birth_probability;
				q << feedback_strength;
				d << death_probability;
				time_to_equilibrium << updates_per_site;
				init << i;

				string name = "pqd_frame_"+std::to_string(grid_size)+"_"+p.str()+"_"+q.str()+"_"+d.str()+"_"+time_to_equilibrium.str()+"_"+init.str();
				save_frame(frame, grid_size, name);
    		}
    	}
	}
}
