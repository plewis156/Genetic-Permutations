# Genetic-Permutations
This project is an exploration of genetic algorithms. It attempts to come to a solution to the traveling salesman problem.

This project uses a genetic algorithm to attempt to find a 'good' if not strictly optimal solution to an NP-complete problem. It first computes the problem using a brute force method. It then does the genetic algorithm by creating a number of tours, finding the best of this set of tours and keeping them, taking the next best tours of this set and mutating them and then keeping them, then replacing the rest of the tours with new tours. It does this a given number of times then it compares the best tour with the optimal solution of the brute force method and also compares the time taken to the brute force method.
