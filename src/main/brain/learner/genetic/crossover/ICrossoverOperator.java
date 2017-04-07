package main.brain.learner.genetic.crossover;

import java.util.ArrayList;

import main.brain.learner.genetic.Gene;

public interface ICrossoverOperator<E> {
    /**
     * Given a list of genes selected for breeding, perform a crossover and 
     * return the resulting genes in an arraylist 
     * @param genes
     * @return arraylist of baby genes
     */
    public abstract ArrayList<Gene<E>> crossover(ArrayList<Gene<E>> genes); 
}
