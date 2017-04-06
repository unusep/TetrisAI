package main.brain.learner.genetic.fitness;

import main.brain.learner.genetic.Gene;

public interface IFitnessFunction {
    public abstract double evaluateFitness(Gene gene);
}
