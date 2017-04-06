package main.brain.learner.genetic;

import java.util.ArrayList;

import main.brain.learner.genetic.crossover.*;
import main.brain.learner.genetic.fitness.*;
import main.brain.learner.genetic.mutator.*;
import main.brain.learner.genetic.selector.*;

public class Population {
    ArrayList<Gene> genePool;
    ICrossoverOperator crossOverOperator;
    IFitnessFunction fitnessFunction;
    IMutationOperator mutationOperator;
    IPopulationSelector populationSelector;
    
    public Population(Gene gene, int size) {
        this.genePool = instantiateGenePool(gene, size);
    }
    
    public Population(String filepath){
        this.genePool = instantiateGenePool(filepath);
    }

    private ArrayList<Gene> instantiateGenePool(String filepath) {
        return null;
    }

    public ArrayList<Gene> instantiateGenePool(Gene gene, int size){
        return null;
    }
    
    public void nextGeneration() {
        evaluateFitness(genePool);
        ArrayList<Gene> fittest = selectFittest(genePool);
        ArrayList<Gene> babies = crossOver(fittest);
        killWeakest(genePool);
        genePool.addAll(babies);
        save(genePool);
    }

    /**
     * Write genes to file
     * @param genes
     */
    private void save(ArrayList<Gene> genes) {
        // TODO Auto-generated method stub
        
    }

    private void killWeakest(ArrayList<Gene> genes) {
        // TODO Auto-generated method stub
        
    }

    private ArrayList<Gene> crossOver(ArrayList<Gene> genes) {
        return crossOverOperator.crossover(genes);
    }

    private ArrayList<Gene> selectFittest(ArrayList<Gene> genes) {
        // TODO Auto-generated method stub
        return null;
    }

    private void evaluateFitness(ArrayList<Gene> genes) {
        for (Gene gene : genePool){
            double fitness = fitnessFunction.evaluateFitness(gene);
            gene.setFitness(fitness);
        }
    }

}
