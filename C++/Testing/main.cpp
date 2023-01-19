#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <fstream>
using namespace std;


int main(){

    int THRESH = 0.00002; 

    string line; 
    ifstream file("pdf_2000.txt");
    int count = 0; vector<double> probs; vector<int> x_states; vector<int> y_states; vector<int> z_states; 
    while (getline(file, line, ' ')) {
        if (count % 4 == 0){
            probs.push_back(stod(line));
        }else if (count % 4 == 1){
            x_states.push_back(stoi(line));
        }else if (count % 4 == 2){
            y_states.push_back(stoi(line));
        }else{
            z_states.push_back(stoi(line));
        }
            
        count++;
    }
    
    vector<vector<int>> states; vector<int> state(3,0);

    for (int i = 0; i < x_states.size(); i++){
        state[0] = x_states[i]; 
        state[1] = y_states[i]; 
        state[2] = z_states[i]; 
        states.push_back(state);
    }

    int n = states.size(); int num; 

    for(int i = 0; i < n; i++){
        cout << "State: " << i+1 << endl; 
        vector<int> current_state = states[i]; vector<vector<int>> neighbors; 

        for (int i = current_state[1]-1; i <= current_state[1]+1; i++){
            for (int j = current_state[0]-1; j <= current_state[0]+1; j++){
                for (int k = current_state[2]-1; k <= current_state[2]+1; k++){
                    vector<int> new_state = {j, i, k}; neighbors.push_back(new_state);
                }
            }
        }

        int m = neighbors.size(); int no_neighbors = 1; 
        for (int j = 0; j < m; j++){
            vector<int> neighbor = neighbors[j]; 
            for (int k = 0; k < n; k++){
                if(neighbor == states[i]){
                    if(probs[i]>THRESH){
                        no_neighbors = 0; 
                        break; 
                    }
                }
            }
        }

        if(no_neighbors){
            num++; 
            cout << "GOTCHA" << endl; 
        }
    }

    cout << "Number of states without big neighbors: " << num << endl;

    return 0;
}

