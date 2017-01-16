#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
#include <algorithm>




using namespace std;

//Each city has a class with its name, its longitude and latitude
class city{
public:
    string name;
    double lati;
    double longi;
    int position;

    city(string n,double la,double lo,int pos){
        name = n;
        lati = la;
        longi = lo;
        position = pos;
    }
};


//this function is useful when we read the .csv file and we want to split each line to instantiate each city class
void Split(vector<string> &vectorr, string line, char separator)
{
	vectorr.clear();

	string::size_type stTemp = line.find(separator);

	while(stTemp != string::npos)
	{
		vectorr.push_back(line.substr(0, stTemp));
		line = line.substr(stTemp + 1);
		stTemp = line.find(separator);
	}

	vectorr.push_back(line);
}

//Function to read the .csv file and return a vector with all the cities and their GPS coordinates
vector<city> fileReader(string file){

    //Let's create the list of cities :

    vector<city> cities;

    //We create a string that is going to be used as the cursor of the file :
    string line;
    int index = 0;
    vector<string> vectorr;

    ifstream myfile (file.c_str());
    //Let's open the file :
    if (myfile.is_open())
    {
        while(getline(myfile,line))
        {
            //We don't want to get the first line :
            if(line != "Title,latitude,longitude"){
            //The split function separates the line and put each parts in the vector
            Split(vectorr,line,',');
            //vectorr is a vector of string. We need to convert vectorr[1] and vectorr[2] in double
            double lati = atof(vectorr[1].c_str());
            double longi = atof(vectorr[2].c_str());
            string cityName = vectorr[0];
            //Then we create a class for each line :
            city A(cityName,lati,longi,index);

            //Then we push this city in the vector that contains all the cities :
            cities.push_back(A);
            index++;
            }
        }
        myfile.close();
        //We close the file as we just finished reading it
    }
    return cities;

}

//Given two cities, this function return the distance between them
double Distance(city A,city B){
    const double earthRadius = 6372.795477598; //In kilometers
    double Alati = A.lati*(M_PI/180);//In radian
    double Alongi = A.longi*(M_PI/180);
    double Blati = B.lati*(M_PI/180);
    double Blongi = B.longi*(M_PI/180);

    double dist = pow(sin((Blati - Alati)/2),2);
    dist = dist + cos(Alati)*cos(Blati)*pow(sin((Blongi-Alongi)/2),2);
    dist = 2*earthRadius*asin(sqrt(dist));

    return dist;
}

//Function that fills the matrix with distances between each city and put it the edge in an adjacency list only if distance is greater than 100km
void AddDistanceMatrix(double** matrixDistance, vector<city> cities){

    double distance;
    for(int i=0; i<cities.size();i++){
        for(int j=0 ; j<cities.size() ; j++){
                distance = Distance(cities[i],cities[j]);
                if(distance > 100){
                    matrixDistance[i][j] = distance;
                }
                else if(distance < 100){
                    //For distance smaller than 100km, we put the larger double that we can
                    matrixDistance[i][j] = LONG_MAX;
                }
        }
    }
}

//Function to free the memory once we're done with the program
void DeleteMatrix(double** distance,int taille){
    for(int i = 0; i < taille; i++) {
        delete (distance[i]);
    }
    delete (distance);
}

//function to return the total distance of a path in km
double GetDistancePath(vector<int> path, double** distance){
    double pathLength = 0;
    for(int i = 1 ;i<path.size() ; i++){
        pathLength += distance[path[i-1]][path[i]];
    }
    return pathLength;
}

//This function generates a first path in order to have a kind of reference to compare with in the BruteForceAlgorithm
vector<int> PathGenerator(string cityDeparture, vector<city> cities){

    vector<int> Path;

    //We put the starting city into the solution, we just put into the vector the position of the city
    for( int i=0 ;i< cities.size(); i++){
        if(cities[i].name == cityDeparture){
            Path.push_back(i);
        }
    }
    //The path contains now the first city
    //We now add the others cities' index
    for (int i= 0 ; i<cities.size(); i++){
        if(cities[i].name != cityDeparture){
            Path.push_back(i);
        }
    }

    //We add the starting city at the end to complete the path
    Path.push_back(Path[0]);

    return Path;
}

//Algorithm that returns the longest path
vector<int> BruteForceAlgorithm(double **distance, string cityDeparture, vector<city> cities){

    //We need a first path and its distance to compare with :
    vector<int> pathReference = PathGenerator(cityDeparture,cities);
    double distancePath = GetDistancePath(pathReference,distance);

    //We create a vector that will contains the longest path :
    vector<int> longestPath;

    //Then we generate all the permutations of this path and we keep the longest
    //The library algorithm has a boolean function that makes all the permutations and returns false when it has done all of them
    //We need all the permutations from the second city to the penultimate
    while( next_permutation(pathReference.begin()+1 , pathReference.end()-1 )){
        //If the distance of the permuted path is greater than the
        if(GetDistancePath(pathReference,distance) > distancePath){
            longestPath = pathReference;
            distancePath = GetDistancePath(pathReference,distance);
        }
    }
    return longestPath;
}

//Void that displays path
void DisplayPathName(vector<int> path, vector<city> cities){
    cout << "Path : ";
    for (int i=0 ; i< path.size()-1 ; i++){
        cout<< cities[path[i]].name << "->" ;
    }
    cout << cities[path[path.size()-1]].name <<" \n";
}

//Function that returns true if the city named by its position is in the path
bool IsAlreadyInPath(vector<int> path,int n){
    bool IsInThePath = false;
    for(int i = 0; i<path.size(); i++){
        if(path[i] == n){
            IsInThePath = true;
        }
    }
    return IsInThePath;
}

//Local Heuristic algorithm that approximates the solution by finding the longest distance between each city
vector<int> FurthestAlgorithm(double** distance, string cityDeparture, vector<city> cities){

    //Vector that will contain the solution
    vector<int> solutionPath;

    //We insert the first city in the vector
    for(int i=0 ; i<cities.size(); i++){
        if(cities[i].name == cityDeparture){
            solutionPath.push_back(cities[i].position);
        }
    }

    double distanceMax;
    int cityDistanceMax;

    //We add n-1 cities. 1 has been inserted before and the end city will be inserted after that :
    for(int i =0 ; i<cities.size()-1; i++){
        distanceMax = 0;
        cityDistanceMax =0;
        //We search the furthest distance from the last city :
        for(int j = 0; j<cities.size() ;j++){
            if( distance[solutionPath[solutionPath.size()-1]][j]>distanceMax && IsAlreadyInPath(solutionPath,j) == false && distance[solutionPath[solutionPath.size()-1]][j]< LONG_MAX ){
                distanceMax = distance[solutionPath[solutionPath.size()-1]][j];
                cityDistanceMax = j;
            }
        }
        solutionPath.push_back(cityDistanceMax);
    }

    //We add the last city (which is the first) at the end of the solutionPath
    solutionPath.push_back(solutionPath[0]);

    return solutionPath;
}

//Function that gives back the position of the city given its name
int GetPositionCity(string cityName,vector<city> cities){
    int position = 0;
    for(int i=0 ;i<cities.size() ;i++){
        if(cities[i].name == cityName){
            position = i;
        }
    }
    return position;
}

//function to find the vertex with minimum key value from the set of vertices not yet included in the MST
int minKey(double key[],bool isInMst[],int nbVertices){
    double mini = LONG_MAX;
    int indexMini;
    for(int i=0 ; i < nbVertices ; i++){
        if(isInMst[i]==false && key[i]<mini){
            mini = key[i];
            indexMini = i;
        }
    }
    return indexMini;
}

//Function to print the MST
void printMST(int parent[], int nbVertices, double** distance,vector<city> cities)
{
    double distanceTotal = 0;
    cout << "Minimum Spanning Tree : \n";

    for (int i = 1; i < nbVertices; i++){
        cout << cities[parent[i]].name << " - "<< cities[i].name << " : "<< distance[i][parent[i]] << " \n";
        distanceTotal += distance[i][parent[i]];
    }

    cout << " \nTotal distance : " << distanceTotal << " kilometers. \n";
}


void PrimsAlgorithm(double** distance, string city1,string city2, vector<city> cities){
    int nbVertices = cities.size();

    int parent[nbVertices]; //Store the MST
    double key[nbVertices];    //Key value to pick the minimum distance
    bool mstSet[nbVertices];    //To know if a vertex is already in the MST

    //All keys are set INFINITE at the beginning
    for(int i=0 ; i<nbVertices ; i++){
        key[i] = LONG_MAX;
        mstSet[i] = false; //No vertex in the MST
    }

    //We put a first vertex in the MST :
    key[0] = 0;
    parent[0] = -1; //The first vextex has not parent

    int iterat =1;
    while(iterat<= nbVertices){
        int u = minKey(key,mstSet,nbVertices);  //u is the index of the min in key array
        mstSet[u] = true;

        //We now update key value and parent index of the adjacent vertices of the picked vertex
        //We only take vertices that are not yet in the MST :

        for(int i =0; i<nbVertices ; i++){
            if(mstSet[i]==false && distance[u][i]<key[i]){
                parent[i] = u;
                key[i] = distance[u][i];
            }
        }
        iterat++;
    }
    //At the end we change the parent of the city1 by city2
    parent[GetPositionCity(city2,cities)] = GetPositionCity(city1,cities);

    //Once we have the MST, we print it :
    printMST(parent,nbVertices,distance,cities);
}


int main()
{
    cout << "Welcome to the ADSA project ! " << " \n";
    cout << " \n";

    //Reading of the .csv file that contains only a few cities in order to test the brute force algorithm :
    vector<city> shortListCities = fileReader("shortcities.csv");
    int nShortList = shortListCities.size(); //Gives the number of cities

    //Initialization of the pointer of pointer which is the matrix that will contain all the distance
    double** distanceShortList = new double*[nShortList];
    for (int i=0 ;i<nShortList ;i++){
        distanceShortList[i] = new double[nShortList];
    }
    //We put the distances for each cities
    AddDistanceMatrix(distanceShortList,shortListCities);

//Question 1a) - Brute Force :

    vector<int> longestPathBruteForce = BruteForceAlgorithm(distanceShortList,"ANGERS",shortListCities);
    //this vector contains the longest path, cities are represented by their index (position)
    cout<< "Brute force algorithm distance : " <<GetDistancePath(longestPathBruteForce,distanceShortList) << " kilometers \n";
    DisplayPathName(longestPathBruteForce, shortListCities);

    //At the end we clean up the pointer of pointer :
    DeleteMatrix(distanceShortList,nShortList);


//Question 1c) - Local Algorithm :

    //We read the .csv file that contains all the cities :
    vector<city> cities = fileReader("Cites.csv");
    int n = cities.size();

    //this pointer of pointer will contain all the distance to connect cities with each others :
    double** distance = new double*[n];
    for(int i = 0; i < n; ++i)
        distance[i] = new double[n];

    AddDistanceMatrix(distance,cities);
    //Now the matrix distance is filled out with all the distance between each city

    vector<int> longestPathFurthestAlgo2 = FurthestAlgorithm(distance,"PARIS",cities);
    cout << "\n \nHeuristic local algorithm distance : "<<GetDistancePath(longestPathFurthestAlgo2,distance) << " kilometers \n";
    DisplayPathName(longestPathFurthestAlgo2,cities);

//Minimum Spanning Tree :

    cout << "\n \n ";
    PrimsAlgorithm(distance,"PARIS","SAINT GEORGES",cities);

    //At the end we clean up the pointer of pointer :
    DeleteMatrix(distance, n);


    return 0;
}
