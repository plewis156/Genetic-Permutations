/***************************************************************
  Paul Lewis
  File Name: tour.c
  Genetic Algorithms and Permutations

  Contains functions for allocating and freeing a tour struct.
***************************************************************/

#include "tour.h"

/*
 * Function for allocating a new tour
 *
 * @param int numCities, the number of cities
 *
 * @local struct tour *newTour, the new tour
 *
 * @return struct tour *, the new tour
 */
struct tour *newTour(int numCities) {
    struct tour *newTour;
    
    if((newTour = malloc(sizeof(struct tour))) == NULL) {
        perror("Malloc error newTour\n");
        exit(1);
    }

    newTour->total = 0.0;

    newTour->cities = malloc(sizeof(int) * numCities);
    if(newTour->cities == NULL) {
        perror("Malloc error newTour->cities\n");
        exit(1);
    }

    return newTour;
}

/*
 * Function to free a tour
 *
 * @param struct tour *t, the tour to free
 */
void freeTour(struct tour *t) {
    free(t->cities);
    free(t);
}
