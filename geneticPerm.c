/***************************************************************
  Paul Lewis
  File Name: geneticPerm.c
  Permutations and Genetic Algorithms

  A program for running 2 algoritms for solving the Traveling Salesman Problem.
  The first is a brute force permutation algorithm which check each possible permutation.
  The second is a genetic algorithm.
***************************************************************/

#include "geneticPerm.h"

/*
 *  Global array to contain values for cities to permute
 */
int s[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};

int main(void) {

    int numCities;
    int genSize;
    int numGen;
    int numMutants;
    int numElites;
    double opt;
    double graph[GRAPH_SIZE][GRAPH_SIZE];
    
    populateGraph(graph);
    getPresets(&numCities, &genSize, &numGen, &numMutants, &numElites);

    /* run brute force algorithm */
    opt = perm2(numCities,graph);
    /* run genetic algorithm */
    genetic(graph, numCities, genSize, numGen, numMutants, numElites, opt);

    return 0;
}
/* 
 * Function to populate graph with values from cities.txt
 *
 * @param double graph[][GRAPH_SIZE], the graph
 *
 * @local int i, j; counters
 * @local line[LINE_SIZE], buffer
 * @local FILE *fp, file pointer
 */
void populateGraph(double graph[][GRAPH_SIZE]) {

    int i, j;
    char line[LINE_SIZE];
    FILE *fp;    

    if((fp = fopen("cities.txt","r")) == NULL) {
        perror("Error opening cities.txt\n");
        exit(1);       
    }
    
    for(i=0;i<GRAPH_SIZE;i++) {
        for(j=0;j<GRAPH_SIZE;j++) {
            if(i == j)
                graph[i][j] = 0.0;
            else {
                fgets(line,sizeof(line),fp);
                graph[i][j] = atof(line);
            } 
        }
    }
    
    fclose(fp);
}

/*
 * Function to parse preset information from presets.txt
 *
 * @param int *numCities, reference to integer containing number of cities
 * @param int *genSize, reference to integer containing size of generation
 * @param int *numGen, reference to integer containing number of generations to run
 * @param int *numMutants, reference to integer containing number of mutants to keep in generation
 * @param int *numElites, reference to integer containing number of elites to keep in generation
 *
 * @local int i, counter
 * @local array[5], array containing integers
 * @local line[LINE_SIZE], buffer
 * @local double x, temporary variable for double
 * @local FILE *fp, file pointer
 */
void getPresets(int *numCities, int *genSize, int *numGen, int *numMutants, int *numElites) {
    int i;    
    int array[ARRAY_SIZE];
    char line[LINE_SIZE];
    double x;
    FILE *fp;

    if((fp = fopen("Presets.txt", "r")) == NULL) {
        perror("Error opening Preset.txt\n");
        exit(1);
    }
    
    for(i=0;((fgets(line, sizeof(line), fp)) != NULL);i++) {
        array[i] = atoi(line);
    }
    
        
    *numCities = array[0]-1;
    *genSize = array[1];
    *numGen = array[2];
    /* numMutants and numElites are integers in cities.txt and converted to double (percentage) */
    x = array[3];
    x = x * 0.01;    
    *numMutants = (*genSize)*x;
    x = array[4];
    x = x * 0.01;
    *numElites = (*genSize)*x;

    fclose(fp);
}
/*
 * A sorting algorithm based on bubble sort which sorts according to gaps to increase speed
 *
 * @param struct tour **tours, the array to sort
 * @param int size, size of array
 * 
 * @local const float shrink, the amount to shrink gap
 * @local struct tour *temp, temporary pointer
 * @local int i, counter
 * @local int gap, the initial gap
 * @local int swapped, boolean to determine if swapped
 */
void combSort(struct tour **tours, int size) {
    int i, gap = size; //initialize gap size
    const float shrink = 1.3f; //set shrink factor
    struct tour *temp;
    int swapped = 0;

    while((gap > 1) || swapped) {
        //gap gets smaller for each comb
        if(gap>1) { //gap must be at least 1
            gap = (int)((float)gap/shrink);
        }
        swapped = 0;
        //increment through array
        for(i=0;gap+i<size; i++) {
            if(tours[i]->total - tours[i + gap]->total > 0) {
                temp =  tours[i];
                tours[i] = tours[i + gap];
                tours[i + gap] = temp;
                swapped = 1;
            }
        }
    }
}
/*
 * Function to permute through all possible permutations in graph
 *
 * @param int n, number of cities
 * @param double graph[][GRAPH_SIZE], the graph
 *
 * @local struct timeval t1, begining time struct
 * @local struct timeval t2, ending time struct
 * @local unsigned long long i, counter
 * @local unsigned long long nfact, factorial of n
 * @local double temp1, temp2; temporary variables
 * @local double time, the time taken
 * 
 * @return double, the optimal distance between cities
 */
double perm2(int n, double graph[][GRAPH_SIZE]) {
    
    struct timeval t1, t2;
    gettimeofday(&t1, NULL); //start timer
    unsigned long long i; 
    unsigned long long nfact = 1;
    double temp1, temp2, time;
    for(i=1;i<=n;i++)
        nfact *= i;
    printf("\nRunning Brute Force Method . . .\n");
    temp1 = getDist(n,graph);   //get distance of before first permutation call
    for(i=0;i<nfact;i++) {      //permute through ever possible permutation
        perm1(n);               
        temp2 = getDist(n, graph);  //get distance
        if(temp2 < temp1) {     //if distance is smaller than temp1, temp1 becomes new distance
            temp1 = temp2;
        }
    }
    gettimeofday(&t2, NULL); //end timer
    time = (double)(t2.tv_usec - t1.tv_usec)/1000000 + (double)(t2.tv_sec - t1.tv_sec); //calculate time
    printf("\nBrute Force Method\n");
    printf("Optimal Cost = %7.3f\n", temp1);
    printf("Time taken : %7.3f seconds\n", time);

    return temp1; 
}

/*
 * Function to permute city array one time
 *
 * @param int n, the number of cities
 *
 * @local int m, k, p, q; temporary variables for swapping
 */
void perm1(int n) {
    int m, k, p, q;
    m = n-2;
    while(s[m]>s[m+1])
        m = m-1;
    k = n-1;
    while(s[m]>s[k])
        k = k-1;
    swap(m,k);
    p = m+1;
    q = n-1;
    while(p<q) {
        swap(p,q);
        p++;
        q--;
    }
}

/*
 * Function for swapping two items in global array
 *  
 * @param int i, j; array indexes for swapping
 *
 * @local int c, temporary value for swapping
 */
void swap(int i, int j) {
    int c = s[i];
    s[i] = s[j];
    s[j] = c;
}

/*
 * Function for calculating distance when going through all cities in a given permutation
 * Used for global array
 *
 * @param int n, number of cities
 * @param double graph[][GRAPH_SIZE], the graph
 *
 * @local int i, counter
 * @local double temp, the sum of the distance
 *
 * @return double, return sum of the distance
 */
double getDist(int n, double graph[][GRAPH_SIZE]) {

    int i;
    double temp;
    /* global array does not contain city 0 so city 0 is added at beginning and end when calculating distance */
    temp = graph[0][s[0]];  
    for(i=0;i<(n-1);i++)
        temp += graph[s[i]][s[i+1]];
    temp += graph[s[i]][0];

    return temp;
}

/*
 * Function for calculating distance when going through all cities in a given permutation
 * Used for local array
 *
 * @param int *arr, array containing the permutation to sum
 * @param int n, number of cities
 * @param double graph[][GRAPH_SIZE], the graph
 *
 * @local int i, counter
 * @local double temp, the sum of the distance
 *
 * @return double, return sum of the distance
 */
double getDist2(int *arr, int n, double graph[][GRAPH_SIZE]) {
    
    int i;
    double temp;
    /* cities array in tour struct does not contain city 0 so city 0 is added at beginning and end when calculating distance */
    temp = graph[0][arr[0]];
    for(i=0;i<(n-1);i++)
        temp += graph[arr[i]][arr[i+1]];
    temp += graph[arr[i]][0];

    return temp;
}

/*
 * Function to mutate a tour by swapping the first city with the last city
 *
 * @param struct tour *t, tour to mutate
 * @param int numCities, the number of cities
 * @param double graph[][GRAPH_SIZE], the graph
 *
 * @local int temp, temporary variable for swapping
 */
void mutate(struct tour *t, int numCities, double graph[][GRAPH_SIZE]) {
    
    int temp;
    temp = t->cities[1];
    t->cities[1] = t->cities[numCities-1];
    t->cities[numCities-1] = temp;
    t->total = getDist2(t->cities, numCities, graph);    
}

/*
 * Function to mutate a tour by swapping the two middle cities
 *
 * @param struct tour *t, tour to mutate
 * @param int numCities, the number of cities
 * @param double graph[][GRAPH_SIZE], the graph
 *
 * @local int temp, temporary variable for swapping
 * @local int x, variable for the middle of the array
 */
void mutate2(struct tour *t, int numCities, double graph[][GRAPH_SIZE]) {
    int temp;
    int x = numCities/2;
    temp = t->cities[x];
    t->cities[x] = t->cities[x+1];
    t->cities[x+1] = temp;
    t->total = getDist2(t->cities, numCities, graph);
}

/*
 * Function to run the genetic algorithm
 * 
 * @param double graph[][GRAPH_SIZE], the graph
 * @param int numCities, the number of cities
 * @param int genSize, the size of a generation
 * @param int numGen, the number of generations to run
 * @param int numMutants, the number of mutants in each generation
 * @param int numElites, the number of elites in each generation
 * @param double opt, the optimal number from the brute force algorithm
 *
 * @local struct gen *gen1, the generation struct
 * @local double newOpt, percentage of optimal for final tour
 * @local int c, i; counters
 * @local struct timeval t1, the beginning time struct
 * @local struct timeval t2, the ending time struct
 * @local double time, the time 
 */
void genetic(double graph[][GRAPH_SIZE], int numCities, int genSize, int numGen, int numMutants, int numElites, double opt) {

    struct gen *gen1 = newGen(numCities, genSize, numElites);
    gen1->lowest = opt*genSize; //initialize lowest
    double newOpt;
    int c, i;
    struct timeval t1, t2;
    double time;
    printf("\nRunning Genetic Algorithm . . .\n");
    gettimeofday(&t1, NULL); //start timer
    for(c=1;c<(GRAPH_SIZE-1);c++) //reinitialize global array to begining value (1, 2, 3, 4, ..., 19) 
        s[c-1] = c;
    firstGen(gen1, graph, numCities, genSize, numMutants, numElites); //create first generation
    for(c=1;c<numGen;c++) {        
        for(i=(numMutants+numElites);i<genSize;i++) { 
            generateTour(gen1->sortArray[i], numCities, graph); //generate new tours
            if(gen1->sortArray[i]->total < gen1->lowest)
                gen1->lowest = gen1->sortArray[i]->total;
        }
        combSort(gen1->sortArray, genSize); //sort array
        geneSplice(gen1, graph, numCities, genSize, numMutants, numElites); //call function to mutate
    }
    gettimeofday(&t2, NULL); // end timer
    newOpt = (gen1->lowest/opt)*100; // find percent of optimal
    time = (double)(t2.tv_usec - t1.tv_usec)/1000000 + (double)(t2.tv_sec - t1.tv_sec); //calculate time
    printf("\nGenetic Method\n");
    printf("Cost = %7.3f\n", gen1->lowest);
    printf("Percent of Optimal = %7.3f\n", newOpt);
    printf("Time taken : %7.3f seconds\n", time);
    freeGen(gen1, genSize); // free generation struct and allocate tours
}

/*
 * Function for finding elites, mutants of elites, and mutants
 *
 * @param struct gen *g, the generation struct
 * @param double graph[][GRAPH_SIZE]
 * @param int numCities, the number of cities
 * @param int genSize, size of generation
 * @param int numMutants, the number of mutants
 * @param int numElites, the number of elites
 * 
 * @local int i, j; counters
 */
void geneSplice(struct gen *g, double graph[][GRAPH_SIZE], int numCities, int genSize, int numMutants, int numElites) {
    
    int i, j;

    for(i=0;i<numElites;i++) {          // assign pointer in elites array to bottom of sortArray (temp pointers)
        g->elites[i] = g->sortArray[i];
    }
    
    j=(numMutants+numElites-1);
    i=0;
    while(i < numElites && j != (numElites-1)) {            // clone tours as clones of elites
        copyTour(g->elites[i], g->sortArray[j], numCities); // and mutate them using mutate() function
        mutate(g->sortArray[j], numCities, graph);
        if(g->sortArray[j]->total < g->lowest)
            g->lowest = g->sortArray[j]->total;
        j--;
        i++;
    }
    i=0;
    while(i < numElites && j != (numElites-1)) {            // clone tours as clones of elites
        copyTour(g->elites[i], g->sortArray[j], numCities); // and mutate them using mutate2() function
        mutate2(g->sortArray[j], numCities, graph);
        if(g->sortArray[j]->total < g->lowest)
            g->lowest = g->sortArray[j]->total;
        j--;
        i++;
    }
    
    while(j >= numElites) {         // mutate rest of mutants
        if(j%2 == 0)
            mutate(g->sortArray[j], numCities, graph);
        else
            mutate2(g->sortArray[j], numCities, graph);
        if(g->sortArray[j]->total < g->lowest)
            g->lowest = g->sortArray[j]->total;
        j--;
    }
}

/*
 * Function for creating the first generation
 * 
 * @param struct gen *g, the generation struct
 * @param double graph[][GRAPH_SIZE], the graph
 * @param int numCities, the number of cities
 * @param int genSize, the size of a generation
 * @param int numMutants, the number of mutants
 * @param int numElites, the number of elites
 *
 * @local int i, counter
 */
void firstGen(struct gen *g, double graph[][GRAPH_SIZE], int numCities, int genSize, int numMutants, int numElites) {
    
    int i;

    for(i=0;i<genSize;i++) {    // fill sortArray functions with new tours
        generateTour(g->sortArray[i], numCities, graph);
        if(g->sortArray[i]->total < g->lowest)
            g->lowest = g->sortArray[i]->total;        
    }
    combSort(g->sortArray, genSize);    //sort array
    geneSplice(g, graph, numCities, genSize, numMutants, numElites);    //call function to mutate
}

/*
 * Function to make a tour an exact replica of another tour
 *
 * @param struct tour *t, the tour to copy
 * @param struct tour *r, the tour to change
 * @param int numCities, the number of cities
 * 
 * @local int i, counter
 */
void copyTour(struct tour *t, struct tour *r, int numCities) {
    int i;
    for(i=0;i<numCities;i++) {
        r->cities[i] = t->cities[i];
    }
    r->total = t->total;
}

/*
 * Function to assign to a tour the current city permutation in global array
 * then permute again, then sum total edges
 *
 * @param struct tour *t, reference to a tour struct
 * @param int numCities, the number of cities
 * @param double graph[][GRAPH_SIZE], the graph
 *
 * @local int i, counter
 */
void generateTour(struct tour *t, int numCities, double graph[][GRAPH_SIZE]) {
    int i;
    for(i=0;i<numCities;i++) {
        t->cities[i] = s[i];
    }
    perm1(numCities);
    t->total = getDist(numCities, graph);
}
/*
 * Function to allocate a generation struct
 *
 * @param int numCities, the number of cities
 * @param int genSize, size of generation
 * @param int numElites, the number of elites
 *
 * @local struct gen *gen1, the new generation struct
 * @local int i, counter
 *
 * @return struct gen *, the new generation struct
 */
struct gen *newGen(int numCities, int genSize, int numElites) {
    
    struct gen *gen1 = malloc(sizeof(struct gen));
    if(gen1 == NULL) {
        perror("Error allocating generation struct\n");
        exit(1);
    }
    int i;
    gen1->elites = malloc(sizeof(struct tour *) * numElites);
    if(gen1->elites == NULL) {
        perror("Error allocating gen1->elites array\n");
        exit(1);
    }
    gen1->sortArray = malloc(sizeof(struct tour *) * genSize);
    if(gen1->sortArray == NULL) {
        perror("Error allocating gen1->sortArray array\n");
        exit(1);
    }
    for(i=0;i<genSize;i++)      // allocate tour struct to sortArray
        gen1->sortArray[i] = newTour(numCities);

    return gen1; 
}

/*
 * Function to free a generation struct
 *
 * @param struct gen *g, the generation struct to free
 * @param int genSize, size of generation
 *
 * @local int i, counter
 */
void freeGen(struct gen *g, int genSize) {
    int i;    
    free(g->elites);    
    for(i=0;i<genSize;i++) {    //free allocated tours
        freeTour(g->sortArray[i]);
    }
    free(g->sortArray);
    free(g);
}
