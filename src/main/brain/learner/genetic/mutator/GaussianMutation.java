package main.brain.learner.genetic.mutator;

import java.util.ArrayList;
import java.util.Random;

import main.brain.learner.genetic.Gene;

public class GaussianMutation<E> implements IMutationOperator<E> {
    private static final double MUTATION_PROBABILITY = 0.05; 
    private static Random random = new Random();
    
    @Override
    public void mutate(Gene<E> gene) {
        ArrayList<Double> weights = gene.getChromosomeWeights();
        for (int i = 0; i < weights.size(); i++){
            if (random.nextDouble() < MUTATION_PROBABILITY){
                double weight = weights.get(i);
                weight += random.nextGaussian();
                if (weight > 1) weight = 1;
                if (weight < -1) weight = -1;
                weights.set(i, weight);
            }
        }
    }

}