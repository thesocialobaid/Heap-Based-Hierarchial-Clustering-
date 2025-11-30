// Muhammad Obaidullah [24L-0509] Assignment 03  - Visual Grouping Tool (Hierarchial Clustering)
// SIMPLIFIED EXPLANATION: 
// 1. We treat every point as a group initially 
// 2. We look at our "Distance Matrix " to find the two closest groups. 
// 3. We merge them into one group
// 4. We update the matrix to reflect the new distances 
// 5. We repeat until we hit the target number of groups (M)

#include<iostream>
#include <vector>
#include <list> 
#include <cmath>
#include <cfloat> 
#include <fstream>
#include <iomanip> 

using namespace std; 

// We are defining a point structure which holds the values for id, the x and the y-cordinate for a particular point. 
struct Point{
    int id; 
    double x; 
    double y; 
};

// Global Constant for "Infinity " --> this is used for distances to self 
const double INF = DBL_MAX;     // this defines the maximum value for the double. 

// --- HELPER FUNCTIONS ---- 

// 1. Calculating the simple distance between two raw points 
double getDistance(Point p1, Point p2){ 
    double dx = p1.x - p2.x; 
    double dy = p1.y - p2.y; 

    return sqrt(dx*dx + dy*dy); 
} 

struct HeapNode { 
    double dist; 
    int u; 
    int v; 
}; 

class MinHeap{ 
    private: 
    vector<HeapNode> tree; 
    int size; 
    vector<vector<int>>* posMatrix; 

    void swapNodes(int i, int j){ 
        HeapNode temp = tree[i]; 
        tree[i] = tree[j]; 
        tree[j] = temp; 

        (*posMatrix)[tree[i].u][tree[i].v] = i;     // For these lines we are trying to swap the indexes for the two values inside the Matrix. 
        (*posMatrix)[tree[i].v][tree[i].u] = i;

        (*posMatrix)[tree[j].u][tree[j].v] = j;
        (*posMatrix)[tree[j].v][tree[j].u] = j;
    }

    void siftUp(int i){     // reheaping the values up 
        while(i>0){ 
            int parent = (i-1) / 2; 
            if(tree[i].dist < tree[parent].dist){
                swapNodes(i,parent);
                i = parent; 
            } else{ 
                break; 
            }
        }
    }

    void siftDown(int i){    // reheaping the values down 
        int smallest = i; 
        int left = 2*i + 2; 
        int right = 2*i + 2; 

        if(left < size && tree[left].dist < tree[smallest].dist)
             smallest = left; 
        if (right < size && tree[right].dist < tree[smallest].dist){ 
            smallest = right; 
        }

        if(smallest != i){ 
            swapNodes(i, smallest); 
            siftDown(smallest); 
        }
    }

    public: 
    MinHeap(int capacity, vector<vector<int>>* matrixRef){ 
        tree.resize(capacity); 
        size = 0; 
        posMatrix = matrixRef; 

    }

    void insert(double d, int u, int v ){ 
        tree[size].dist = d; 
        tree[size].u = v; 
        tree[size].v = u; 
        (*posMatrix)[u][v] = size; 
        (*posMatrix)[v][u] = size; 
        siftUp(size);
        size++; 
    }

    HeapNode extractMin(){ 
        if(size == 0 ) return {INF, -1, -1}; 
        HeapNode minNode = tree[0]; 
        tree[0] = tree[size-1]; 

        // Updating the matrix for the moved node 
        (*posMatrix)[tree[0].u][tree[0].v] = 0; 
        (*posMatrix)[tree[0].v][tree[0].u] = 0; 

        //Invalidating the extracted Node 
        (*posMatrix)[minNode.u][minNode.v] = -1; 
        (*posMatrix)[minNode.v][minNode.u] = -1; 

        size--; 
        siftDown(0); 
        return minNode; 
    }

    void updateKey(int index, double newDist){ 
        if(index <0 || index >= size){      // This is the safety check to ensure that there is no loss of information 
            return; 
        }

        double oldDist = tree[index].dist; 
        tree[index].dist = newDist; 

        if(newDist < oldDist){ 
            siftUp(index);
        }
        else{ 
            siftDown(index); 
        }
    }

    // Helper functiont to allow us to peek at the distances for Single Linkage Logic 
    double getDistAt(int index){ 
        if(index <0 || index >= size) return INF; 
        return tree[index].dist; 
    }

    bool isEmpty() { return size == 0; }

}; 

// 2. Finding the two closest active groups in the matrix 
// APPROACH 
// First we scan the matrix and use a nested loop to look at every possible pair of groups.
// The outer loop (i) picks the first group and the inner loop (j) picks the second group. 
// The inner group starts at j = i+1 which is efficient because the distance from Group 1 to Group 2 is the same as Group 2 is the same as Group 2 to Group 1, so we only need to check one side ( the upper triangle of the matrix)

// Checking the active status allowing us to ensure that before checking the distance, we classify if the group still exists or not. 
// If a group was previously merged into another, it's active flag is false. The loop immediately skips it so we don't accidently merge a 'dead group' 

// We maintain a variable minVal ( initially set to infinity) as check each pair, we compare their distance to this minimum 
// By the time the loop finishes, we are having two variables u and v that will hold the IDs of the specific pair of groups that are closest to each other. 

void findClosestGroups(int N, const vector<vector<double>>& matrix, const vector<bool>& active, int&u, int& v, double& minVal){ 
    minVal = INF; 
    u = -1; 
    v = -1; 

    // Loop through every row (i) and column (j) 
    for(int i=0; i<N; i++){
        if(!active[i]) continue;   // Skipping the groups that we already merged 

        for(int j=i+1; j < N; j++){
            if(!active[i]) continue; // Skip groups we already merged 

            // If this pair is closer than anything we've seen so far, we must find the minimum value 
            if(matrix[i][j] < minVal){ 
                minVal = matrix[i][j]; 
                u = i; 
                v= j; 
            }
        }
    }
} 


// 3. Merging the group V into Group U 
// To merge the groups together we simply unhook the Group V chain and hook it to the end of the Group U chain 
// 
void mergeGroups(int u, int v, vector<list<Point>>& groups, vector<bool>& active){ 
   //splice moves elements from one list to another efficently. 
   groups[u].splice(groups[u].end(), groups[v]); // Groups is a vector of Linked Lists. Splice is a special command for Linked List. It transfers all elements from one list (groups[v]) to another (groups[u])
   // groups[u].end() tells the computer to put the new items at the very beggining of Group U. 
   // We use this because of the time complexity of 0(1) 

   // Deactivate v (it is now part of u) 
   active[v] = false; 
}

//4. Updating the Distance matrix after the merge. 
// - The new distance to the merged group is the MINIMUM of the old distances. 
// - logic: dist(MergedGroup, K) = min(dist(U,K), dist(V,K))
// The distance update is handled by the updatematrix function and it is based on the Single Linkage rule, which means the distance 
// between two groups is defined by the distance between their two closest members. 

void updateMatrix(int N, int u, int v, vector<vector<double>>& matrix, const vector<bool>& active){ 
    // We are merging Group V into Group U. U is the new, combined group (U+V)
    //Looping through every other group K in the entire matrix 
    for(int k=0; k < N; ++k){
        // 1. Guarding against updating distances to groups that don't exist 
        // (Groups already merged) or to Group U itself. 
        if(active[k] && k!= u){ 
            // 2. Getting the OLD Distances for the two components of the new  group 
            double distFromU = matrix[u][k]; // Distance from original U to K 
            double distFromV = matrix[v][k]; // Distance from original V to K 

            // Applying the Single Linkage Rule (Finding the minimum distance ) and the distance from the new merged group (U+V) to K is the 
            // minimum of the two old distances 
            double newDist = min(distFromU, distFromV); 

            // Updating the matrix, we update the row/column for U with the new minimum distance 
            // We don't touch row/columns V since it is now inactive 
            matrix[u][k] = newDist; 
            matrix[k][u] = newDist; // Since the matrix is symmetric we are going to be updating both the sides.
        }
    }
}


// -----------------------------------------------------------------------------------------------------------
// FROM HERE OWNWARDS WE ARE GOING TO BE WORKING ON THE IMPLEMENTATION OF QUESTION 2 
// In the first question, we scanned the entire table everytime, which is slow. A min Heap solves this by keeping the smallest number at always the top. We can grap the 
// smallest distance instantly. However, heaps are hard to search. If Group 4 and Group 7 merge, we need to update the distance between "Group 4" and "Group 9 ". But 
// where is "4 vs 9" inside the Heap Array. It could be at index 5, or index 100, we do not know. 

//  In order to resolve this issue we tend to create a map (Matrix) that tells us exactly where in the Heap Array a specific pair is located. 
// Matrix[i][j] = K means "The distance between Groups i and Group j is stored at Index K inside the heap array" 

// The Algorithmic Logic 
// 1. Initializing by calculating all the distances and putting them in the heap, we record their heap index in the Matrix 
// 2. Pop: we take the top of the Heap (closest pair, say U and V)
// 3. Merge: Moving the points from V to U 
// 4. Updating: We need to calculate the new distance from the merged U to every other group K. 
// New Distance = min(old_dist(U,K), old_dist(V,K))
// We look at Matrix[U][K] to find where that specific node is in the heap. We go there, change the value, and then "Fix" the heap (move the node up if the value got smaller)
// We effectively delete V's connections by setting them to infiity. 



int main(){ 
// // 1. OPEN FILE AND THE IMPLEMENTATION OF EXCEPTION HANDLING TO ENSURE THAT THERE ARE NO ERRORS WHILE WORKING WITH FILES 
// ifstream inFile("inputData.txt");
// if(!inFile){ 
//     cout << "Error: Could not open the input data file " << endl; 
//     return 1 ;  
// }


// // Reading the values for N= Total Points , M = Target Groups from the file 
// int N, M; 
// inFile >> N >> M; 

// cout << "Processing " << N << "points down to " << M << "groups .... " << endl; 

// // 2. SETUP STORAGE 
// vector<list<Point>> groups(N); 
// vector<bool> active(N,true); 
// vector<Point> rawPoints(N); 

// // Matrix initialization - The distance Matrix here is used to store the distance between the Point 1 and Point 2 and Point 3, etc. 
// // It then saves those distances in the distMatrix grid. 
// // The distance from a point to itself is set to INF in order to ensure that the code doesn't try to merge a point with itself. 

// vector<vector<double>> distMatrix(N, vector<double>(N)); 

// //Reading the data 
// for(int i=0; i< N; ++i){ 
//     double x,y; 
//     inFile >> x >> y; 
//     Point p = {i+1, x,y}; 

//     rawPoints[i] = p; 
//     groups[i].push_back(p); // Initially, every point is its own group 

// }
// inFile.close(); 

// // 3. INITIAL MATRIX CALCULATION 
// for(int i=0; i<N; i++){ 
//     for(int j=0; j<N; j++){ 
//         if(i == j ) distMatrix[i][j] = INF; 
//         else distMatrix[i][j] = getDistance(rawPoints[i], rawPoints[j]); 
//     }
// }


// // 4. THE MAIN LOOP ( Run until we have M groups left)
// // This keeps on running a cycle until the number of groups drop down to the target (M). 

// int currentGroups = N; 

// while (currentGroups > M)
// {
//     int u,v; 
//     double minDist; 

//     //A. Finding the closest pair 
//     findClosestGroups(N,distMatrix, active, u,v, minDist); // It looks for the smallest number in the whole grid (ignoring groups that have already been switched "OFF")
//                                                            // It finds that Group A and Group B are only 1.5 units apart. This is the smallest distance, so u becomes A and v becomes B 
                                                    

//     if(u == -1) break; 

//     cout << "Merging Group " << (v+1) << "into " << (u+1) << "(Distance: " << fixed << setprecision(2) << minDist<< ")" << endl; 

//     // Merging the groups together 
//     mergeGroups(u,v,groups,active);   // It takes all the points inside Group B(v) and moves them in Group A (u). Here the .splice() is used to ensure that instead of copying the data, we rewire the list pointers and snap 
//                                       // the list of B onto the end of A instantly. 
//                                       // Then it sets active[v] = false. Group B is now "dead". It exists only inside Group A now. 

//     // Updating the math (Distance Matrix) 
//     updateMatrix(N,u, v, distMatrix, active);  // We now have a new, bigger Group, we need to update the distances on our map. The code updates 
//                                                // the distance matrix with these new, shorter distances 
                                            

//     currentGroups--; 
    
// }

// // 5. Once the while loop is finished, the code stops merging and it loops through the active list one last time. 
// // If a group is still marked "ON" (true), it prints out "Group X " and lists all the points currently inside it . 

// cout << "\nFinal Groups: " << endl;
// cout << "*******************"<<endl; 
// int displayindex = 1; 
// for(int i=0; i<N; ++i){ 
//     if(active[i]) {
//         cout << "Group " << displayindex++ << ": "; 
//         for(const auto& p : groups[i]){
//             cout << "(" << p.x << ", " << p.y << ")"; 
//         }
//         cout << endl; 
        
//     }
// }


ifstream inFIle("InputData.txt"); 
if(!inFIle) { 
    cout << "Error: Create inputData.txt first! " <<endl; 
    return 1; 
}

int N,M; 
char comma; 
// Handling formatting " 3,9 " or "9,3 " 
if(!(inFIle >> N)) return 0; 

char nextchar; 
inFIle >> nextchar; 
if(nextchar == ','){ 
   inFIle >> M; 
} else { 
    inFIle.putback(nextchar); 
    inFIle>> M; 
}

vector<Point> rawPoints(N); 
vector<list<Point>> groups(N); 
vector<bool> active(N,true); 

// We are going to be initializing the position matrix first with -1 
vector<vector<int>> posMatrix(N, vector<int>(N-1)); 

MinHeap heap(N*N, &posMatrix); 

// 1. Reading the input 
for(int i=0; i < N; i++){ 
    for(int j= i+1 ; j < N; j++){ 
        double d = getDistance(rawPoints[i], rawPoints[j]); 
        heap.insert(d,i,j); 
    }
}

// 3. Main Loop 
int currentGroups = N; 
while (currentGroups> M)
{
    // Extracting the Minimum Element 
    HeapNode minNode = heap.extractMin(); 
    int u = minNode.u; 
    int v = minNode.v; 

    //Skipping if this edge involves already merged groups (lazy deletion check)
    if(!active[u] || !active[v]) continue; 

    // B. Merging V into U 
    // Keeping the smaller ID active (convention) or just keeping U . 
    // Letting explicitly keep min(u,v) active to be consisent 
    int keeper = (u<v) ? u:v; 
    int merger = (u<v) ? v:u; 

    u = keeper; 
    v = merger; 

    // Merging The Linked lists in O(1) time.
    groups[u].splice(groups[u].end(), groups[v]); 
    active[v] = false; 
    currentGroups--; 
    
    // Updating the Distances 
    for(int k=0; k< N; k++){
        if(active[k] && k!=u){ 
            int idxU = posMatrix[u][k]; 
            int idxV = posMatrix[v][k]; 

        if(idxU != -1 && idxV != -1 ){ 
            double d1 = heap.getDistAt(idxU); 
            double d2 = heap.getDistAt(idxV); 

            //Single Linkage: New Distance is the minimum of the two 
            double newDist = min(d1,d2); 

            //Updating U's distance to K with the new value 
            heap.updateKey(idxU, newDist); 

            // Effectively removing V's distance to K by setting to INF
            // It will sink to bottom and be ignored by active checks 
            heap.updateKey(idxV, INF); 
        }
        }
    }
}

int gNum = 1;
    for (int i = 0; i < N; i++) {
        if (active[i]) {
            cout << "Group " << gNum++ << ": ";
            for (const auto& p : groups[i]) {
                cout << "(" << p.x << "," << p.y << ") ";
            }
            cout << endl;
        }
    }



    return 0; 
} 
