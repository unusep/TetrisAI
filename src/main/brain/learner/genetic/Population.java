package main.brain.learner.genetic;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;

import main.brain.learner.genetic.crossover.*;
import main.brain.learner.genetic.fitness.*;
import main.brain.learner.genetic.mutator.*;
import main.brain.learner.genetic.selector.*;

/**
 * Stores a population of genes with chromosomes of type E
 * Genetic algorithm is executed to generate nextGeneration()
 * TODO: Abstract file functions to facilitate SLAP
 */
public class Population<E> {
    private final double PERCENTAGE_TO_KILL = 0.2;
    
    private ArrayList<Gene<E>> genePool;
    private ICrossoverOperator<E> crossOverOperator;
    private IFitnessFunction<E> fitnessFunction;
    private IMutationOperator<E> mutationOperator;
    private IPopulationSelector<E> populationSelector;
    private String savePath;
    
    public Population(String filepath, ArrayList<E> chromosomes, int populationSize){
        this.savePath = filepath;
        this.genePool = instantiateGenePool(filepath, chromosomes, populationSize);
    }

    /**
     * Write genes to file
     * @param genes
     */
    public void saveToDisk(){
        File file = new File(savePath);
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(file));
            
            if (!genePool.isEmpty()){
                bw.write(genePool.get(0).getHeaders());
                bw.newLine();
            }
            
            for (Gene<E> gene : genePool){
                bw.write(gene.toString());
                bw.newLine();
            }
            
            bw.flush();
            bw.close();
        } catch (IOException e) {
            System.out.println("File failed to save");
        }
        System.out.println("Population saved to disk");
    }
    
    /**
     * Returns the fittest gene among the gene pool 
     * @return fittest gene
     */
    public Gene<E> getFittest() {
        return populationSelector.selectElite(genePool, 1).get(0);
    }
    
    /**
     * Generates the next generation by running the genetic algorithm:
     * 1) get fitness of genepool
     * 2) choose genes for breeding
     * 3) crossover to generate babies 
     * 4) mutate babies
     * 5) kill to make room
     * 6) add new babies to genepool
     */
    public void nextGeneration() {
        evaluateFitness(genePool);
//        normaliseFitness(genePool);
        int numBabies = (int) (PERCENTAGE_TO_KILL * (double) genePool.size()); 
        ArrayList<Gene<E>> parents = selectElite(genePool, numBabies);
        ArrayList<Gene<E>> babies = crossOver(parents);
        mutate(babies);
        cull(genePool, numBabies);
        genePool.addAll(babies);
        printBestGene();
    }
    
    int printCount = 0;
    private void printBestGene() {
        Gene<E> best = getFittest();
        printCount++;
        System.out.println("Turn " + printCount + " Best Fitness: " + best.getFitness() + " Weights: " + best.getChromosomeWeights());
    }

    public void setCrossOverOperator(ICrossoverOperator<E> crossOverOperator){
        this.crossOverOperator = crossOverOperator;
    }

    public void setFitnessFunction(IFitnessFunction<E> fitnessFunction){
        this.fitnessFunction = fitnessFunction;
    }

    public void setmutationOperator(IMutationOperator<E> mutationOperator){
        this.mutationOperator = mutationOperator;
    }
    
    public void setpopulationSelector(IPopulationSelector<E> populationSelector){
        this.populationSelector = populationSelector;
    }
    
    /**
     * Instantiate gene pool of given size using a text file of stored genes
     * @param gene
     * @param size
     * @return
     */
    private ArrayList<Gene<E>> instantiateGenePool(String filepath, ArrayList<E> chromosomes, int populationSize) {
        ArrayList<Gene<E>> genes = new ArrayList<Gene<E>>(populationSize);
        
        try {
            File f = new File(filepath);
            Scanner sc = new Scanner(f);
            sc.nextLine();
            
            // read from file and generate genes
            while (sc.hasNext()) {
                if (populationSize <= 0) break;
                
                String line = sc.nextLine();
                String[] result = line.split(" ");
                
                if (result.length -1 != chromosomes.size()){
                    System.out.println("Warning: chromosomes and weights are not of the same length.");
                    System.out.println("Warning: will continue by generating a new random population");
                    break;
                }
                
                double fitness = Double.parseDouble(result[result.length - 1]);
                
                ArrayList<Double> weights = new ArrayList<Double>();
                for (int i = 0; i < chromosomes.size(); i++){
                    weights.add(Double.parseDouble(result[i]));
                }

                //
                Gene<E> gene = new Gene<E>(chromosomes, weights);
                gene.setFitness(fitness);
                
                populationSize--;
                genes.add(gene);
            }
            sc.close();
        } catch (FileNotFoundException e) {
            System.out.println("Pathfile not found: " + filepath +". Will try to create one.");
        }
        
        // if file does not have enough gene, generate more automatically
        if (populationSize > 0){
            genes.addAll(instantiateGenePool(chromosomes, populationSize));
        }
        
        return genes;
    }

    /**
     * Generate a randomised gene pool of given size
     * @param gene
     * @param size
     * @return
     */
    public ArrayList<Gene<E>> instantiateGenePool(ArrayList<E> chromosomes, int populationSize){
        int fitness = 0;
        ArrayList<Gene<E>> genes = new ArrayList<Gene<E>>(populationSize); 
        for (int i = 0; i < populationSize; i++){
            ArrayList<Double> weights = new ArrayList<Double>(chromosomes.size());
            for (int j = 0; j < chromosomes.size(); j++){
                weights.add(randomDouble(-1, 1));
            }
            Gene<E> gene = new Gene<E>(chromosomes, weights);
            gene.setFitness(fitness);
            genes.add(gene);
            System.out.println("Added gene: " + gene);
        }
        return genes;
    }
    
    private double randomDouble(double min, double max){
        Random randomGenerator = new Random(); 
        return min + (randomGenerator.nextDouble() * (max - min)); 
    }
    
    
    private void mutate(ArrayList<Gene<E>> genes) {
        for (Gene<E> gene : genes){
            mutationOperator.mutate(gene);
        }
    }

    private ArrayList<Gene<E>> crossOver(ArrayList<Gene<E>> genes) {
        return crossOverOperator.crossover(genes);
    }

    /**
     * Select the genes that should be killed using the population selector
     * @param genePool
     * @param num number to select
     * @return the genes that should be removed
     */
    private void cull(ArrayList<Gene<E>> genePool, int num) {
        populationSelector.cull(genePool, num);
    }
    
    /**
     * Select the genes that are allowed to breed using the population selector
     * @param genePool
     * @param num number to kill
     * @return the genes that are allowed to breed
     */
    private ArrayList<Gene<E>> selectElite(ArrayList<Gene<E>> genes, int num) {
        return populationSelector.selectElite(genes, num);
    }

    /**
     * Compute and update the fitness of each gene in the genepool using the fitness function
     * @param genes
     */
    private void evaluateFitness(ArrayList<Gene<E>> genes) {
        for (Gene<E> gene : genePool){
            double fitness = fitnessFunction.evaluateFitness(gene);
            gene.setFitness(fitness);
        }
    }

    /**
     * Normalization means dividing the fitness value of each individual 
     * by the sum of all fitness values, so that the sum of all resulting 
     * fitness values equals 1.
     * @param genePool
     */
    private void normaliseFitness(ArrayList<Gene<E>> genePool) {
        double totalFitness = 0.0;
        for (Gene<E> gene : genePool){
            totalFitness += gene.getFitness();
        }
        for (Gene<E> gene : genePool){
            gene.setFitness(gene.getFitness()/totalFitness);
        }
    }

}
