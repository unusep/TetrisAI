package main.brain.learner.genetic.crossover;

import java.util.ArrayList;

import main.brain.learner.genetic.Gene;

public interface ICrossoverOperator {
    public abstract ArrayList<Gene> crossover(ArrayList<Gene> genes); 
}
