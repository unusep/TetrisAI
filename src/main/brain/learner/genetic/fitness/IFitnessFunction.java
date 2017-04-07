package main.brain.learner.genetic.fitness;

import main.brain.learner.genetic.Gene;

public interface IFitnessFunction<E> {
    /**
     * Given a gene, evaluate the gene's fitness
     * @param gene to evaluate
     * @return fitness of gene
     */
    public abstract double evaluateFitness(Gene<E> gene);
}
