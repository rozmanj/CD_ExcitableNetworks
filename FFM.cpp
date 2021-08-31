#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <random>
#include <iomanip>

using namespace std;

//**********************************************************************
//******************FUNCTIONS*******************************************
//**********************************************************************
int max_val (vector<int> clone_sizes)
{
    int max_v = 0;
    for (int i = 0; i < clone_sizes.size(); i++)
    {
        if (max_v < clone_sizes[i]) max_v = clone_sizes[i];
    }
    return max_v;
}

double total_val (vector<double> value_list)
{
    if (value_list.size() == 0) return 0.0;
    double avr_value = 0.0;
    for (int i = 0; i < value_list.size(); i ++)
    {
        avr_value += value_list[i];
    }
    return avr_value;
}

double total_val (vector<int> value_list)
{
    if (value_list.size() == 0) return 0.0;
    double avr_value = 0.0;
    for (int i = 0; i < value_list.size(); i ++)
    {
        avr_value += value_list[i];
    }
    return avr_value;
}

double avr (vector<double> value_list)
{
    if (value_list.size() == 0) return 0.0;
    double avr_value = 0.0;
    for (int i = 0; i < value_list.size(); i ++)
    {
        avr_value += value_list[i];
    }
    return avr_value / (double) value_list.size();
}

double avr (vector<int> value_list)
{
    if (value_list.size() == 0) return 0.0;
    double avr_value = 0.0;
    for (int i = 0; i < value_list.size(); i ++)
    {
        avr_value += value_list[i];
    }
    return avr_value / (double) value_list.size();
}

double standard_deviation (vector<double> value_list)
{
    if (value_list.size() == 0) return 0.0;
    double avr_value = avr(value_list);
    double std_value = 0.0;
    for (int i = 0; i < value_list.size(); i++)
    {
        std_value += pow(value_list[i] - avr_value, 2.0);
    }
    return sqrt(std_value / (double) value_list.size());
}

int consiscenty_check(vector<vector<int>> cells)
{
    int all_fine = 1;
    int cell_found;
    vector<int> current_cell;
    vector<int> other_cell;
    for (int i = 0; i < cells.size(); i++){
        current_cell = cells[i];
        for (int j = 2; j < current_cell.size(); j++)
        {
            other_cell = cells[current_cell[j]];
            cell_found = 0;
            for (int k = 2; k < other_cell.size(); k++)
            {
                if (other_cell[k] == i) {cell_found = 1;}
            }
            if (cell_found != 1) {all_fine = 0;}
        }
    }
    return all_fine;
}

int long_consistency_check(vector<vector<int>> cells)
{
    //checks a) if all links go in both direction, and b) if all linked cells in the same clone
    int linked_cell_index;
    if (consiscenty_check(cells) != 1 ) return 1;
    for (int cell_index = 0; cell_index < cells.size(); cell_index++)
    {
        for (int link_index = 2; link_index < cells[cell_index].size(); link_index++)
        {
            linked_cell_index = cells[cell_index][link_index];
            if (cells[cell_index][1] != cells[linked_cell_index][1]) return 1;
        }
    }
    return 0;
}

string double_to_string(double dbl)
{
    ostringstream strs;
    strs << dbl;
    return strs.str();;
}

string int_to_string(int integer)
{
    ostringstream strs;
    strs << integer;
    return strs.str();;
}


//**********************************************************************
//******************DIVERSITY MEASURES**********************************
//**********************************************************************
double gini (vector<int> clone_sizes)
{
    int G_n = 0;
    int G_d = 0;
    for (int i = 0; i < clone_sizes.size(); i++){
        G_d += clone_sizes[i];
        for (int j = 0; j < clone_sizes.size(); j++){
            G_n += abs(clone_sizes[i] - clone_sizes[j]);
        }
    }
    return (double) G_n / (2 * clone_sizes.size() * G_d);
}   

double hoover (vector<int> clone_sizes)
{
	double H_up = 0.;
	double H_down = 0.;
	double x0 = avr(clone_sizes);
	double x1;
	for (int i = 0; i < clone_sizes.size(); i ++ )
	{
		x1 = clone_sizes[i];
		H_down += x1;
		H_up   += fabs(x1 - x0);
	}
	return 0.5 * H_up / H_down;
}

double shannon(vector<int> clone_sizes)
{
	double NN = (double) total_val(clone_sizes);
	double Sh = 0.;
	double ni;
	for (int i = 0; i < clone_sizes.size(); i ++ )
	{
		ni = clone_sizes[i];
		Sh += -1 * (ni / NN) * log(ni / NN);
	}
	return Sh;
}

double simpson_r(vector<int> clone_sizes)
{
	double NN = (double) total_val(clone_sizes);
	double Sr = 0.;
	double ni;
	for (int i = 0; i < clone_sizes.size(); i ++ )
	{
		ni = clone_sizes[i];
		Sr += (ni / NN) * (ni / NN);
	}
	return Sr;
}

double simpson(vector<int> clone_sizes)
{
	double NN = (double) total_val(clone_sizes);
	double Sr = 0.;
	double ni;
	for (int i = 0; i < clone_sizes.size(); i ++ )
	{
		ni = clone_sizes[i];
		Sr += (ni / NN) * ((ni - 1) / (NN - 1));
	}
	return Sr;
}

double moment_mean(vector<int> clone_sizes)
{
	double sum_val = (double) total_val(clone_sizes);
	double count_val = (double) (clone_sizes.size());
	return sum_val/count_val;
}

double moment_variance(vector<int> clone_sizes)
{
	double NN = (double) (clone_sizes.size());
	double out_val = 0.;
    double mu = moment_mean(clone_sizes);
	for (int i = 0; i < clone_sizes.size(); i ++ )
	{
		out_val += pow(((double) clone_sizes[i] - mu),2.);
	}
	return out_val/NN;
}

double moment_skewness(vector<int> clone_sizes)
{
	double NN = (double) (clone_sizes.size());
	double out_val = 0.;
    double mu = moment_mean(clone_sizes);
    double sig1 = sqrt(moment_variance(clone_sizes));
	for (int i = 0; i < clone_sizes.size(); i ++ )
	{
		out_val += pow(((double) clone_sizes[i] - mu)/sig1,3.);
	}
	return out_val/NN;
}

double moment_kurtosis(vector<int> clone_sizes)
{
	double NN = (double) (clone_sizes.size());
	double out_val = 0.;
    double mu = moment_mean(clone_sizes);
    double sig1 = sqrt(moment_variance(clone_sizes));
	for (int i = 0; i < clone_sizes.size(); i ++ )
	{
		out_val += pow(((double) clone_sizes[i] - mu)/sig1,4.);
	}
	return out_val/NN;
}

//**********************************************************************
//******************MAIN************************************************
//**********************************************************************
int main(int argc, char *argv[])
{
	
	//***********SIMULATION PARAMETERS********************
	double p_i, p_t, p_r, p_d;
	p_i = std::stod(argv[1]);
	p_t = std::stod(argv[2]);
	p_r = std::stod(argv[3]);
	p_d = std::stod(argv[4]);
	std::cout << p_i << " " << p_t << " " << p_r << " " << p_d << endl;

   

    //**********MAIN VARIABLES****************************
    int N_cells;
    int time;
    vector<vector<int>> cells;    
    vector<int> cell_colors;
    vector<int> clone_sizes;
	int simulation_error = 0;

	//************FOR DIVERSITY************************
	vector<vector<vector<double>>> diversity_store = {};
    for (int i = 0; i < 14; i ++ ) 
    {
    	diversity_store.push_back({});
    	for (int j = 0; j < 1200; j ++ )
    	{
    		diversity_store[i].push_back({});
    	}
    }

    //***************RND************************************
    mt19937 generator;
    uniform_real_distribution<double> distribution(0.0,1.0);

    //***************FOR INICIALIZATION*********************
    int clone_blue;
    int size_of_clone; 
    int is_blue;
    vector<vector<int>> all_clone_sizes = {
        {4,4,4,3,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1},
        {5,5,3,3,3,3,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1},
        {8,4,2,2,2,2,2,4,2,1,1,1,1,1,1,1,1,1,1,1,1},
        {2,2,2,2,5,4,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
        {6,4,3,3,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1},
        {7,3,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1},
        {6,3,3,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1},
        {4,3,3,2,2,2,2,2,2,2,2,2,1,1,1,1,1},
        {8,5,5,4,3,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1},
        {4,4,3,3,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
        {7,4,3,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}
    };
    //***************FOR PROPAGATION**********************
    vector<vector<int>> new_cell_network;
    vector<int> updated_network;
    vector<int> cells_to_divide;
    vector<int> new_cell;
    vector<int> current_cell;
    vector<int> other_cell;

    //***************FOR DIVISION************************
    vector<int> dividing_cell;
    vector<int> old_daughter_cell;
    vector<int> new_daughter_cell;
    vector<int> original_links;
    int rnd_H_value;
    int found_other_link;

    //**************FOR ANALYSIS*************************
    vector<int> cell_visited;
    vector<int> component_sizes;
    int current_component_size;
    int observed_cell;
    vector<int> cell_queue;
    double Sh_val;
    double Sh_max;
    std::vector<double> max_fractions = {};
    std::vector<double> clone_count = {};
    std::vector<double> blue_fractions = {};
    double blue_count = 0;


    //**********************************************************************
    //******************RUN*************************************************
    //**********************************************************************
    for (int gini_repeats = 0; gini_repeats < 2200; gini_repeats ++ )
    {	
    	cells = {};
    	cell_colors = {};
        clone_sizes =  all_clone_sizes[(int) gini_repeats / 200];
        
        N_cells = 0;
        time = 0;

        //**********************************************************************
        //******************STARTING CELL NETWORK*******************************
        //**********************************************************************
        for (int i = 0; i < clone_sizes.size(); i++)
        {
            size_of_clone = clone_sizes[i];
            if (distribution(generator) < 0.5) is_blue = 0;
            else is_blue = 1;
            for (int j = 0; j < size_of_clone; j++)
            {
                N_cells += 1;
                new_cell = {0, i}; //creates a cell belonging to the i-th clone
                if (j == 0 && size_of_clone != 1)
                {
                    new_cell.push_back(N_cells);
                }
                else if (j == size_of_clone - 1 && size_of_clone != 1)
                {
                    new_cell.push_back(N_cells - 2);
                }
                else if (size_of_clone != 1)
                {
                    new_cell.push_back(N_cells);
                    new_cell.push_back(N_cells - 2);
                }
                cells.push_back(new_cell);
                cell_colors.push_back(is_blue);
            }

        }
        if (cells.size() != N_cells) {simulation_error += 1;}
        //**********************************************************************
        //******************SIMULATION******************************************
        //**********************************************************************
        while (N_cells <1000)
        {	
        	time += 1;
        	new_cell_network = {};
            cells_to_divide = {};
            if (N_cells != cells.size()) {simulation_error += 1;}
            
            //**********************************************************************
        	//******************PROPAGATE STATES************************************
        	//**********************************************************************
			for (int i = 0; i < N_cells; i++)
            {
                current_cell = cells[i];
                updated_network = {current_cell[0], current_cell[1]}; //makes an updated network of same cell type and clone, without links

                //***************TREE**********************
    			if (current_cell[0] == 0) 
                {
                	//***************SPONTANEOUS**********************
    				if (distribution(generator)< p_i) //checks if cell catches fire on its own
                	{
                		updated_network[0] = 1;
                	}
                	//**************FIRE TRANSFER*********************
                	else
                	{
                		for (int j = 2; j < current_cell.size(); j++)
                        {
                             other_cell = cells[current_cell[j]];
                             if (other_cell[0] == 1)
                             {
                                 if (distribution(generator) < p_t) 
                                 {
                                     updated_network[0] = 1;
                                 }
                             }
                         }
                	}
				}
      			//***************FIRE**********************
    			else if (current_cell[0] == 1)
                 {
                     updated_network[0] = 2;
                     cells_to_divide.push_back(i);
                 }
				//***************REFRACTORY****************
    			else if (current_cell[0] == 2)
                {
                    if (distribution(generator) < p_r) {updated_network[0] = 0;}
                }
                //**************ADD LINKS******************
                for (int j = 2; j < current_cell.size(); j++)
                {
                    updated_network.push_back(current_cell[j]);
                }
                new_cell_network.push_back(updated_network);
            }
            //**********************************************************************
        	//******************DIVIDE**********************************************
        	//**********************************************************************
			for (int i = 0; i < cells_to_divide.size(); i++)
            {
            	dividing_cell = new_cell_network[cells_to_divide[i]];
            	original_links = {};
                for (int j = 2; j < dividing_cell.size(); j++)
                {
                    original_links.push_back(dividing_cell[j]);
                }
               
                //*****************MAKE NEW CELLS*****************************
                new_daughter_cell = {2, dividing_cell[1]};
                new_cell_network[cells_to_divide[i]] = {2, dividing_cell[1]}; 

                //*****************ASSIGN ORIGINAL LINKS**********************
                for (int j = 0; j < original_links.size(); j++)
                {
                     
                    if (distribution(generator) > 0.5)
                    {
                        new_daughter_cell.push_back(original_links[j]);
                        found_other_link = 0;
                        for (int k = 2; k < new_cell_network[original_links[j]].size(); k++)
                        {
                            if (new_cell_network[original_links[j]][k] == cells_to_divide[i]) 
                            {
                            	new_cell_network[original_links[j]][k] = N_cells;
                            	found_other_link += 1;
                            }
                        }
                        if (found_other_link != 1) {simulation_error += 1;}
                    }
                    else
                    {
                        new_cell_network[cells_to_divide[i]].push_back(original_links[j]);
                    }
                    
                }
                
                //***************CONNECT DAUGHTER CELLS**********************
                if (distribution(generator) > p_d)
                {
                	new_cell_network[cells_to_divide[i]].push_back(N_cells);
                	new_daughter_cell.push_back(cells_to_divide[i]);
                }

                //***************COMPLETE DIVISION**************************
                new_cell_network.push_back(new_daughter_cell);
                N_cells += 1;
                clone_sizes[new_daughter_cell[1]] += 1;
                cell_colors.push_back(cell_colors[cells_to_divide[i]]);
            }

            //**********************************************************************
        	//******************UPDATE**********************************************
        	//**********************************************************************
			cells = new_cell_network;
            
            //**********************************************************************
        	//******************ANALYZE*********************************************
        	//**********************************************************************

			if (cells_to_divide.size() != 0)
            {
				//****************COMPONENTS********************************************
	        	component_sizes = {};

	        	cell_visited = {};
	        	for (int i = 0; i < N_cells; i ++)
	        	{
	        		cell_visited.push_back(0);
	        	}


	        	cell_queue = {};
	        	for (int i = 0; i < N_cells; i ++)
	        	{
	        		if (cell_visited[i] != 1)
	        		{
	        			cell_queue.push_back(i);
	        			current_component_size = 0;
	        			while(cell_queue.size() != 0)
	        			{
	        				observed_cell = cell_queue[0];
	        				if (cell_visited[observed_cell] == 1) {simulation_error += 1;}
	        				if (cells[observed_cell].size() < 2) {simulation_error += 1;}
	        				for (int j = 2; j < cells[observed_cell].size(); j ++ )
	        				{
	        					if (cell_visited[cells[observed_cell][j]] == 0)
	        					{
	        						cell_queue.push_back(cells[observed_cell][j]);
	        					}
	        				}
	        				cell_visited[observed_cell] = 1;
	        				current_component_size += 1;
	        				cell_queue.erase(cell_queue.begin());
	        			}
	        			component_sizes.push_back(current_component_size);
	        		}
	        	}

	        	//***************DIVERSITY MEASURES**************************************
	        	if (N_cells < 1150)
                {
                    Sh_val = shannon(component_sizes);
                    Sh_max = log(component_sizes.size());

                    diversity_store[0][N_cells].push_back(gini(component_sizes));
                    diversity_store[1][N_cells].push_back(Sh_val);
                    diversity_store[2][N_cells].push_back(exp(Sh_val));
                    diversity_store[3][N_cells].push_back(Sh_val/Sh_max);
                    diversity_store[4][N_cells].push_back(Sh_max - Sh_val);
                    diversity_store[5][N_cells].push_back(simpson_r(component_sizes));
                    diversity_store[6][N_cells].push_back(simpson(component_sizes));
                    diversity_store[7][N_cells].push_back(((double) max_val(component_sizes)) / (double) N_cells);
                    diversity_store[8][N_cells].push_back(hoover(component_sizes));
                    diversity_store[9][N_cells].push_back(moment_mean(component_sizes));
                    diversity_store[10][N_cells].push_back(moment_variance(component_sizes));
                    diversity_store[11][N_cells].push_back(moment_skewness(component_sizes));
                    diversity_store[12][N_cells].push_back(moment_kurtosis(component_sizes));
                    diversity_store[13][N_cells].push_back(component_sizes.size());
                    simulation_error += long_consistency_check(cells);
                }

        	}
        }

        //******************CHAMBER COMPLETE*************************************************
        if (total_val(component_sizes) != N_cells) {simulation_error += 1;}
        
        max_fractions.push_back(((double) max_val(component_sizes)) / (double) N_cells);
        clone_count.push_back(component_sizes.size());
        
        blue_fractions.push_back(((double) total_val(cell_colors)) / (double) N_cells);
        blue_count = 0.;
        for (int i = 0; i < cell_colors.size(); i ++ ) {if (cell_colors[i] == 1) blue_count += 1;}
        if (blue_count != total_val(cell_colors)) {simulation_error += 1;}
        if (cell_colors.size() != N_cells) {simulation_error += 1;}
        
        std::cout << gini_repeats << " " << gini(component_sizes) << std::endl;
 
    }
	
	//**********************************************************************
    //******************OUTPUT*************************************************
    //**********************************************************************
    		
    //******************MAX*************************************************
   	ofstream file_max("output//max//max_pi_" + double_to_string(p_i) +"_pt_" + double_to_string(p_t) +"_pr_" + double_to_string(p_r) +"_pd_" + double_to_string(p_d) + "_.dat");
    for (int i = 0; i < 2200; i ++)
   	{
   		file_max << max_fractions[i] << " " << clone_count[i] << endl;
   	}
   	file_max.close();

    //******************BLUE*************************************************
   	ofstream file_blue("output//bl//blue_pi_" + double_to_string(p_i) +"_pt_" + double_to_string(p_t) +"_pr_" + double_to_string(p_r) +"_pd_" + double_to_string(p_d) + "_.dat");
    for (int i = 0; i < 2200; i ++)
   	{
   		file_blue << blue_fractions[i] << endl;
   	}
   	file_blue.close();

    //******************DIVERSITY*************************************************
   	ofstream file_diversity("output//diversity//diversity_pi_" + double_to_string(p_i) +"_pt_" + double_to_string(p_t) +"_pr_" + double_to_string(p_r) +"_pd_" + double_to_string(p_d) + "_.dat");   
    for (int j = 0; j < diversity_store[0].size(); j ++ )
    {
    	file_diversity << j;
    	for (int i = 0; i < diversity_store.size(); i ++ )
    	{
    		if (diversity_store[i][j].size() > 1)
    		{
    			file_diversity << "\t" << avr(diversity_store[i][j]) << "\t" << standard_deviation(diversity_store[i][j]);
    		}
    		else
    		{
    			file_diversity << "\t-1\t-1";
    		}
    	}
    	file_diversity << "\t" << diversity_store[0][j].size() << "\n";
    }
    file_diversity.close();

   	//******************COMPLETE*************************************************
    ofstream file_complete("output//cl//complete_pi_" + double_to_string(p_i) +"_pt_" + double_to_string(p_t) +"_pr_" + double_to_string(p_r) +"_pd_" + double_to_string(p_d) + "_.dat");
    file_complete << simulation_error << endl;
    file_complete.close();
}
