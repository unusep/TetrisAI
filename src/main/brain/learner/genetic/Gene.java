package main.brain.learner.genetic;

import java.text.DecimalFormat;
import java.util.ArrayList;

/**
 * Gene stores chromosomes of type E and the weights of each of its chromosome
 */
public class Gene<E> implements Comparable<Gene<E>> {
    private ArrayList<Double> chromosomesWeights;
    private ArrayList<E> chromosomes;
    private double fitness;
    
    public Gene(ArrayList<E> chromosomes, ArrayList<Double> weights) {
        this.chromosomes = chromosomes; 
        this.chromosomesWeights = weights;
        this.fitness = Double.MIN_VALUE;
    }

    public double getFitness(){
        return this.fitness;
    }
    
    public void setFitness(double fitness) {
        this.fitness = fitness;
    }

    @Override
    public int compareTo(Gene<E> o) {
        return (int) (this.fitness - o.fitness);
    }
    
    public ArrayList<E> getChromsomes(){
        return chromosomes;
    }

    public ArrayList<Double> getChromosomeWeights() {
        return chromosomesWeights;
    }
    
    /**
     * Returns a space separated string representing the names of each chromosome
     */
    public String getHeaders(){
        String res = "";
        for (E n : chromosomes){
            res += n.toString() + " ";
        }
        res += "Fitness";
        return res;
    }
    
    /**
     * Returns a space separated string representation the weights of each heuristic followed by its fitness
     */
    @Override
    public String toString(){
        DecimalFormat decimalFormat = new DecimalFormat("#.0000");
        String res = "";
        for (Double weight : chromosomesWeights){
            res += decimalFormat.format(weight) + " ";
        }
        res += decimalFormat.format(fitness);
        return res;
    }
    
}
