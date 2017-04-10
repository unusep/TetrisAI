package main.brain.learner.genetic.crossover;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import javafx.util.Pair;
import main.brain.learner.genetic.Gene;
import main.tetris.heuristics.IHeuristic;

public class SinglePointCrossover implements ICrossoverOperator<IHeuristic>{
    @Override
    public ArrayList<Gene<IHeuristic>> crossover(ArrayList<Gene<IHeuristic>> genes) {
        if (genes.size() < 2) return new ArrayList<Gene<IHeuristic>>(genes);
        ArrayList<Gene<IHeuristic>> babies = new ArrayList<Gene<IHeuristic>>();
        ArrayList<IHeuristic> heuristics = genes.get(0).getChromsomes();
        
        int size = genes.size();
        for (int i = 1; i <= size/2; i++){
            ArrayList<Double> fatherChromosomeWeights = genes.get(i).getChromosomeWeights();
            ArrayList<Double> motherChromosomeWeights = genes.get(size - i).getChromosomeWeights();
            Pair<ArrayList<Double>, ArrayList<Double>> children = cross(fatherChromosomeWeights, motherChromosomeWeights);
            Gene<IHeuristic> baby1 = new Gene<IHeuristic>(heuristics, children.getKey());
            Gene<IHeuristic> baby2 = new Gene<IHeuristic>(heuristics, children.getValue());
            babies.add(baby1);
            babies.add(baby2);
        }
        
        return babies;
    }

    private Pair<ArrayList<Double>, ArrayList<Double>> cross(
            ArrayList<Double> fatherChromosomeWeights, ArrayList<Double> motherChromosomeWeights) {
        int size = fatherChromosomeWeights.size();
        Random random = new Random();
        int crossoverPoint = random.nextInt(size);
        ArrayList<Double> weights1 = new ArrayList<Double>(size);
        ArrayList<Double> weights2 = new ArrayList<Double>(size);
        
        for (int i = 0; i < crossoverPoint; i++){
            double fChromosome = fatherChromosomeWeights.get(i);
            double mChromosome = motherChromosomeWeights.get(i);  
            weights1.add(fChromosome);
            weights2.add(mChromosome);
        }
        for (int i = crossoverPoint; i < size; i++){
            double fChromosome = fatherChromosomeWeights.get(i);
            double mChromosome = motherChromosomeWeights.get(i);  
            weights2.add(fChromosome);
            weights1.add(mChromosome);
        }
        return new Pair<ArrayList<Double>, ArrayList<Double>>(weights1, weights2);
    }

}
