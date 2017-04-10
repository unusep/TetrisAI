package main.brain.learner.genetic;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Comparator;

/**
 * Gene stores chromosomes of type E and the weights of each of its chromosome
 */
public class Gene<E> implements Comparable<Gene<E>> {
    public static final double INITIAL_FITNESS = -Double.MAX_VALUE;
    private ArrayList<Double> chromosomesWeights;
    private ArrayList<E> chromosomes;
    private double fitness;
    
    public Gene(ArrayList<E> chromosomes, ArrayList<Double> weights) {
        this.chromosomes = chromosomes; 
        this.chromosomesWeights = weights;
        this.fitness = INITIAL_FITNESS;
    }

    public double getFitness(){
        return this.fitness;
    }
    
    public void setFitness(double fitness) {
        this.fitness = fitness;
    }

    @Override
    public int compareTo(Gene<E> o) {
        double res = this.fitness - o.fitness;
        if (res == 0){
            return 0;
        } else if (res < 0){
            return -1;
        } else {
            return 1;
        }
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
            res += new String(decimalFormat.format(weight)) + " ";
        }
        res += decimalFormat.format(fitness);
        return res;
    }
}
