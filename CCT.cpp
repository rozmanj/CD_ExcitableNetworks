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

//**********************************************************************
//******************MAIN************************************************
//**********************************************************************
int main(int argc, char *argv[])
{
	
	//***********SIMULATION PARAMETERS********************
	double t0_1, t0_2, t0_ratio, sig_1, sig_2, clone_prob, p_d;
	
    t0_1 = std::stod(argv[1]);
	t0_ratio = std::stod(argv[2]);
    t0_2 = t0_ratio * t0_1;
	
    sig_1 = std::stod(argv[3]);
	sig_2 = std::stod(argv[4]);
    
    clone_prob = std::stod(argv[5]);
    p_d = std::stod(argv[6]);
	std::cout << t0_1 << " " << sig_1 << " " << t0_2 << " " << sig_2 << " " << t0_ratio << " "<< clone_prob << " " << p_d << endl;

   

    //**********MAIN VARIABLES****************************
    int N_cells;
    int time;
    vector<vector<int>> cells;
    vector<double> cell_division_times;
    vector<int> cell_fast;
    vector<int> cell_colors;

    vector<int> clone_sizes;
	int simulation_error = 0;
    int time_error = 0;

	//************FOR DIVERSITY************************
	vector<vector<vector<double>>> diversity_store = {};
    for (int i = 0; i < 10; i ++ ) 
    {
    	diversity_store.push_back({});
    	for (int j = 0; j < 1200; j ++ )
    	{
    		diversity_store[i].push_back({});
    	}
    }

    //***************RND************************************
    mt19937 generator_gauss;
    normal_distribution<double> distribution_gauss(0.0,1.0);

    mt19937 generator_uniform (1);
    uniform_real_distribution<double> distribution_uniform(0.0,1.0);

    //***************FOR INICIALIZATION*********************
    int clone_blue;
    int size_of_clone; 
    int is_blue;
    int is_fast;
    int fast_clone_count;
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

    double t0_current;
    double sig_current;

    double current_time;
    double next_division_time;
    double current_cycle_time;
    double current_cycle_sig;

    int min_cell_index;
    double min_cell_time;

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
    std::vector<double> starting_blue_fractions = {};
    std::vector<double> starting_fast_fractions = {};
    std::vector<double> fast_fractions = {};
    std::vector<double> blue_fractions = {};
    std::vector<double> fin_times = {};
    double blue_count = 0;
    double fast_count = 0;



    //**********************************************************************
    //******************RUN*************************************************
    //**********************************************************************
    for (int gini_repeats = 0; gini_repeats < 2200; gini_repeats ++ )
    {	
    	cells = {};
    	cell_colors = {};
        cell_fast = {};
        cell_division_times = {};

        clone_sizes =  all_clone_sizes[(int) gini_repeats / 200];

        N_cells = 0;
        time = 0;

        //**********************************************************************
        //******************STARTING CELL NETWORK*******************************
        //**********************************************************************
        blue_count = 0;
        fast_count = 0;
        for (int i = 0; i < clone_sizes.size(); i++)
        {
            size_of_clone = clone_sizes[i];
            if (distribution_uniform(generator_uniform) < 0.5) is_blue = 1;
            else is_blue = 0;

            if (distribution_uniform(generator_uniform) < clone_prob) is_fast = 1;
            else is_fast = 0;

            
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
                cell_fast.push_back(is_fast);

                blue_count += is_blue;
                fast_count += is_fast;

                if (is_fast == 1) {t0_current = t0_2; sig_current = sig_2;}
                else {t0_current = t0_1; sig_current = sig_1;}
                next_division_time = -1.;
                while (next_division_time <= 0)
                {
                    next_division_time = t0_current + sig_current * distribution_gauss(generator_gauss);
                }
                cell_division_times.push_back(next_division_time);
                


            }

        }
        starting_blue_fractions.push_back(blue_count);
        starting_fast_fractions.push_back(blue_count);
        if (cells.size() != N_cells) {simulation_error += 1;}
        if (cell_division_times.size() != N_cells) {simulation_error += 1;}
        //**********************************************************************
        //******************SIMULATION******************************************
        //**********************************************************************
        current_time = 0.;
        while (N_cells <1000)
        {	
        	new_cell_network = {};
            cells_to_divide = {};
            if (N_cells != cells.size()) {/*cout << "ERROR: N_cells mismatch: " << N_cells << " " << cells.size() << std::endl; */simulation_error += 1;}
            
           //**********************************************************************
        	//******************MAKE UPDATED NETWORK********************************
        	//**********************************************************************
			for (int i = 0; i < N_cells; i++)
            {
                current_cell = cells[i];
                updated_network = {current_cell[0], current_cell[1]}; //makes an updated network of same cell type and clone, without links

                
                //**************ADD LINKS******************
                for (int j = 2; j < current_cell.size(); j++)
                {
                    updated_network.push_back(current_cell[j]);
                }
                new_cell_network.push_back(updated_network);
            }
            //**********************************************************************
            //******************FIND DIVIDING CELL**********************************
            //**********************************************************************
            min_cell_index = 0;
            min_cell_time = cell_division_times[0];
            for (int i = 0; i < N_cells; i++)
            {
                if (cell_division_times[i] < min_cell_time)
                {
                    min_cell_time = cell_division_times[i];
                    min_cell_index = i;
                }
            }

            if (current_time >= min_cell_time) time_error += 1;
            current_time = min_cell_time;
            cells_to_divide.push_back(min_cell_index);

            

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
                     
                    if (distribution_uniform(generator_uniform) > 0.5)
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
                if (distribution_uniform(generator_uniform) > p_d)
                {
                	new_cell_network[cells_to_divide[i]].push_back(N_cells);
                	new_daughter_cell.push_back(cells_to_divide[i]);
                }

                //***************COMPLETE DIVISION**************************
                new_cell_network.push_back(new_daughter_cell);
                N_cells += 1;
                
                clone_sizes[new_daughter_cell[1]] += 1;
                cell_colors.push_back(cell_colors[cells_to_divide[i]]);
                cell_fast.push_back(cell_fast[cells_to_divide[i]]);

                //***************SET DIVISION TIMES************************
                if (cell_fast[cells_to_divide[i]] == 1) {t0_current = t0_2; sig_current = sig_2;}
                else {t0_current = t0_1; sig_current = sig_1;}
                next_division_time = -1.;
                while (next_division_time <= 0)
                {
                    next_division_time = t0_current + sig_current * distribution_gauss(generator_gauss);
                }
                cell_division_times[cells_to_divide[i]] = current_time + next_division_time;

                next_division_time = -1.;
                while (next_division_time <= 0)
                {
                    next_division_time = t0_current + sig_current * distribution_gauss(generator_gauss);
                }
                cell_division_times.push_back(current_time + next_division_time);
                
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
                    diversity_store[9][N_cells].push_back(component_sizes.size());
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

        fast_fractions.push_back(((double) total_val(cell_fast)) / (double) N_cells);
        fin_times.push_back(current_time);
        
        std::cout << gini_repeats << " " << gini(component_sizes) <<  std::endl;

    }
	
	//**********************************************************************
    //******************OUTPUT**********************************************
    //**********************************************************************
    		
   	//******************MAX*************************************************
    ofstream file_max("output//max//max_t1_" + double_to_string(t0_1) +"_sig1_" + double_to_string(sig_1) +"_t0r_" + double_to_string(t0_ratio) + "_sig2_" + double_to_string(sig_2) + "_prob_" + double_to_string(clone_prob) +"_pd_" + double_to_string(p_d) + "_.dat");
    for (int i = 0; i < 2200; i ++)
   	{
   		file_max << max_fractions[i] << " " << clone_count[i] << endl;
   	}
   	file_max.close();

    //******************TIME*************************************************
    ofstream file_time("output//time//time_t1_" + double_to_string(t0_1) +"_sig1_" + double_to_string(sig_1) +"_t0r_" + double_to_string(t0_ratio) + "_sig2_" + double_to_string(sig_2) + "_prob_" + double_to_string(clone_prob) +"_pd_" + double_to_string(p_d) + "_.dat");
    for (int i = 0; i < 2200; i ++)
    {
        file_time << fin_times[i] << endl;
    }
    file_time.close();

    //******************BLUE*************************************************
    ofstream file_blue("output//bl//blue_t1_" + double_to_string(t0_1) +"_sig1_" + double_to_string(sig_1) +"_t0r_" + double_to_string(t0_ratio) + "_sig2_" + double_to_string(sig_2) + "_prob_" + double_to_string(clone_prob) +"_pd_" + double_to_string(p_d) + "_.dat");
    for (int i = 0; i < 2200; i ++)
   	{
   		file_blue << starting_blue_fractions[i] << " "  << blue_fractions[i] << " " << starting_fast_fractions[i] << " " << fast_fractions[i] << endl;
   	}
   	file_blue.close();
    
   	//******************DIVERSITY*************************************************
    ofstream file_diversity("output//diversity//diversity_t1_" + double_to_string(t0_1) +"_sig1_" + double_to_string(sig_1) +"_t0r_" + double_to_string(t0_ratio) + "_sig2_" + double_to_string(sig_2) + "_prob_" + double_to_string(clone_prob) +"_pd_" + double_to_string(p_d) + "_.dat");
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
    ofstream file_complete("output//cl//complete_t1_" + double_to_string(t0_1) +"_sig1_" + double_to_string(sig_1) +"_t0r_" + double_to_string(t0_ratio) + "_sig2_" + double_to_string(sig_2) + "_prob_" + double_to_string(clone_prob) +"_pd_" + double_to_string(p_d) + "_.dat");
    file_complete << simulation_error << " " << time_error << endl;
    file_complete.close();
}
