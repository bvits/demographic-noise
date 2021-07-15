#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <cstdlib>
#include <random>
#include <vector> 
#include <algorithm>
static std::mt19937_64 rng(std::random_device{}());
#include <chrono> 
using namespace std::chrono; 
#include<bits/stdc++.h> 
using namespace std; 
#include <sys/resource.h>
#include <cmath>
#include <iostream> 
#include<string>
using namespace std; 

//---------------------Some Primitives -------------------------------------------------//

struct coordinates { 
   int x;
   int y; 
}; 

void seed_twister(int seed);
double random_real(double initial, double last);
int random_int(int initial, int last);
float mean_of_array(float array[],int size);
float standard_deviation_of_array(float array[],int size);
float mean_of_vector(vector<float> array,int size);
void random_frame(int frame[], int grid_size);
void random_frame_of_density(float density, int frame[], int grid_size);
void randomly_occupied_frame(int how_many_occupied, int frame[], int grid_size);
void print_frame(int frame[], int grid_size);
void print_vector(std::vector<double> vec);
void save_frame(int frame[], int grid_size, string name);
void zeros(int frame[], int grid_size);
float calculate_density(int frame[], int grid_size);

template<typename T>std::vector<double> linspace(T start_in, T end_in, int num_in){

  //Source: https://stackoverflow.com/questions/27028226/python-linspace-in-c/27030598#27030598
  //Credits: https://stackoverflow.com/users/1078084/akavall

  //Equivalent of numpy's linspace method. 

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); 
  return linspaced;
};

// -------------------------------PQD Model-------------------------------------------//


coordinates select_neighbor_of_site(coordinates site, int grid_size);
coordinates select_neighbor_of_pair(coordinates site, coordinates neighbour, int grid_size);
void pqd_update(int frame[], int grid_size, float birth_probability, float feedback_strength, float death_probability);
void simulate_pqd(int frame[], int grid_size, float birth_probability, float feedback_strength, float death_probability, int updates_per_site);
void equilibrium_density_pqd(int grid_size, float birth_probability, float feedback_strength, float death_probability, int realizations, int lag, int initializations, float how_many_occupied_initially, float denities[], int updates_per_site=100000, int collect_frames=0);