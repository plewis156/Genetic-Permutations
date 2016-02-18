/***************************************************************
  Paul Lewis
  File Name: tour.h
  Genetic Algorithms and Permutations

  Contains struct definitions, libraries, and function prototypes for tour.c
***************************************************************/

#ifndef TOUR_H
#define TOUR_H

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

/*
 * A tour struct
 *
 * @field int *cities, the cities on this tour
 * @field double total, the total of the edges on tour
 */
struct tour {
    int *cities;
    double total;
};

/*
 * Function for allocating a new tour
 *
 * @param int numCities, the number of cities
 *
 * @return struct tour *, the new tour
 */
struct tour *newTour(int numCities);
/*
 * Function to free a tour
 *
 * @param struct tour *t, the tour to free
 */
void freeTour(struct tour *t);

#endif
