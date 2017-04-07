package main.brain.learner.genetic.crossover;

import java.util.ArrayList;
import java.util.Collections;

import javafx.util.Pair;
import main.brain.learner.genetic.Gene;
import main.tetris.heuristics.IHeuristic;

public class SimpleCrossover implements ICrossoverOperator<IHeuristic>{
    private static final double AMOUNT_TO_SWAP = 0.3; 

    @Override
    public ArrayList<Gene<IHeuristic>> crossover(ArrayList<Gene<IHeuristic>> genes) {
        ArrayList<Gene<IHeuristic>> babies = new ArrayList<Gene<IHeuristic>>();
        ArrayList<IHeuristic> heuristics = genes.get(0).getChromsomes();
        
        int size = genes.size();
        for (int i = 0; i < size/2; i++){
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
        int swapSize = (int) (size * AMOUNT_TO_SWAP);
        ArrayList<Integer> numbersToSwap = new ArrayList<Integer>(size); 
        for (int i = 0; i < size; i++){
            numbersToSwap.add(i);
        }
        ArrayList<Double> weights1 = new ArrayList<Double>(fatherChromosomeWeights);
        ArrayList<Double> weights2 = new ArrayList<Double>(motherChromosomeWeights); 
        Collections.shuffle(numbersToSwap);
        for (int i = 0; i < swapSize; i++){
            swap(i, weights1, weights2);
        }
        return new Pair<ArrayList<Double>, ArrayList<Double>>(weights1, weights2);
    }

    private void swap(int i, ArrayList<Double> weights1, ArrayList<Double> weights2) {
        Double save = weights1.get(i);
        weights1.set(i, weights2.get(i));
        weights2.set(i, save);
    }

}
