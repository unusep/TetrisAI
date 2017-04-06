package main.brain.learner;
import java.util.ArrayList;

import main.brain.learner.genetic.Gene;
import main.brain.learner.genetic.Population;
import main.tetris.heuristics.Heuristic;

public class GeneticLearner extends Learner {
    Population population;
    
    public GeneticLearner(ArrayList<Heuristic> heuristics, int size) {
        super(heuristics);
        Gene gene = new Gene(heuristics);
        this.population = new Population(gene, size);
    }

    @Override
    public void trainLearner(int iterations) {
        for (int i = 0; i < iterations; i++){
            population.nextGeneration();
        }
    }

    @Override
    public ArrayList<Double> getWeights() {
        // TODO Auto-generated method stub
        return null;
    }


}
