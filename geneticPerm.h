/***************************************************************
  Paul Lewis
  File Name: geneticPerm.h
  Genetic Algorithms and Permutations

  Contains struct definitions, libraries, #defines, and function prototypes for geneticPerm.c
***************************************************************/

#ifndef GENETICPERM_H
#define GENETICPERM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <strings.h>
#include <time.h>
#include <sys/time.h>
#include "tour.h"

/* size of graph */
#define GRAPH_SIZE 20
/* size of buffer */
#define LINE_SIZE 8
/* size of array to hold preset numbers */
#define ARRAY_SIZE 5

/*
 * A generation struct
 *
 * @field double lowest, the current best path
 * @field struct tour **elites, array to keep elites
 * @field struct tour **sortArray, array to allocate new tours with previous mutants and elites and sort them
 */
struct gen {
    double lowest;
    struct tour **elites;
    struct tour **sortArray;
};

/* 
 * Function to populate graph with values from cities.txt
 *
 * @param double graph[][GRAPH_SIZE], the graph
 */
void populateGraph(double graph[][GRAPH_SIZE]);
/*
 * Function to parse preset information from presets.txt
 *
 * @param int *numCities, reference to integer containing number of cities
 * @param int *genSize, reference to integer containing size of generation
 * @param int *numGen, reference to integer containing number of generations to run
 * @param int *numMutants, reference to integer containing number of mutants to keep in generation
 * @param int *numElites, reference to integer containing number of elites to keep in generation
 */
void getPresets(int *numCities, int *genSize, int *numGen, int *numMutants, int *numElites);
/*
 * Function to permute through all possible permutations in graph
 *
 * @param int n, number of cities
 * @param double graph[][GRAPH_SIZE], the graph
 * 
 * @return double, the optimal distance between cities
 */
double perm2(int n, double graph[][GRAPH_SIZE]);
/*
 * Function to permute city array one time
 *
 * @param int n, the number of cities
 */
void perm1(int n);
/*
 * Function for swapping two items in global array
 *  y
 * @param int i, j; array indexes for swapping
 */
void swap(int i, int j);
/*
 * Function for calculating distance when going through all cities in a given permutation
 * Used for global array
 *
 * @param int n, number of cities
 * @param double graph[][GRAPH_SIZE], the graph
 *
 * @return double, return sum of the distance
 */
double getDist(int n, double graph[][GRAPH_SIZE]);
/*
 * Function for calculating distance when going through all cities in a given permutation
 * Used for local array
 *
 * @param int *arr, array containing the permutation to sum
 * @param int n, number of cities
 * @param double graph[][GRAPH_SIZE], the graph
 *
 * @return double, return sum of the distance
 */
double getDist2(int *arr, int n, double graph[][GRAPH_SIZE]);
/*
 * A sorting algorithm based on bubble sort which sorts according to gaps to increase speed and minimize turtles
 *
 * @param struct tour **tours, the array to sort
 * @param int size, size of array
 */
void combSort(struct tour **tours, int size);
/*
 * Function to mutate a tour by swapping the first city with the last city and the two middle cities
 *
 * @param struct tour *t, tour to mutate
 * @param int numCities, the number of cities
 * @param double graph[][GRAPH_SIZE], the graph
 */
void mutate(struct tour *t, int numCities, double graph[][GRAPH_SIZE]);
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
 */
void genetic(double graph[][GRAPH_SIZE], int numCities, int genSize, int numGen, int numMutants, int numElites, double opt);
/*
 * Function for creating the first generation
 * 
 * @param struct gen *g, the generation struct
 * @param double graph[][GRAPH_SIZE], the graph
 * @param int numCities, the number of cities
 * @param int genSize, the size of a generation
 * @param int numMutants, the number of mutants
 * @param int numElites, the number of elites
 */
void firstGen(struct gen *g, double graph[][GRAPH_SIZE], int numCities, int genSize, int numMutants, int numElites);
/*
 * Function to make a tour an exact replica of another tour
 *
 * @param struct tour *t, the tour to copy
 * @param struct tour *r, the tour to change
 * @param int numCities, the number of cities
 */
void copyTour(struct tour *t, struct tour *r, int numCities);
/*
 * Function to allocate a generation struct
 *
 * @param int numCities, the number of cities
 * @param int genSize, size of generation
 * @param int numElites, the number of elites
 * @param double opt, the optimal path from brute force alorithm
 *
 * @return struct gen *, the new generation struct
 */
struct gen *newGen(int numCities, int genSize, int numElites);
/*
 * Function to free a generation struct
 *
 * @param struct gen *g, the generation struct to free
 * @param int genSize, the size of the generation
 */
void freeGen(struct gen *g, int genSize);
/*
 * Function for finding elites, mutants of elites, and mutants
 *
 * @param struct gen *g, the generation struct
 * @param double graph[][GRAPH_SIZE]
 * @param int numCities, the number of cities
 * @param int genSize, size of generation
 * @param int numMutants, the number of mutants
 * @param int numElites, the number of elites
 */
void geneSplice(struct gen *g, double graph[][GRAPH_SIZE], int numCities, int genSize, int numMutants, int numElites);
/*
 * Function to mutate a tour by swapping the two middle cities
 *
 * @param struct tour *t, tour to mutate
 * @param int numCities, the number of cities
 * @param double graph[][GRAPH_SIZE], the graph
 */
void mutate2(struct tour *t, int numCities, double graph[][GRAPH_SIZE]);
/*
 * Function to assign to a tour the current city permutation in global array
 * then permute again, then sum total edges
 *
 * @param struct tour *t, reference to a tour struct
 * @param int numCities, the number of cities
 * @param double graph[][GRAPH_SIZE], the graph
 */
void generateTour(struct tour *t, int numCities, double graph[][GRAPH_SIZE]);

#endif
