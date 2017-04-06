package main.brain.learner.genetic;

import java.util.ArrayList;

import main.tetris.heuristics.Heuristic;

public class Gene implements Comparable<Gene> {
    ArrayList<Heuristic> heuristics;
    double fitness;
    
    public Gene(ArrayList<Heuristic> heuristics) {
        this.heuristics = heuristics;
    }

    public void setFitness(double fitness) {
        this.fitness = fitness;
    }

    @Override
    public int compareTo(Gene o) {
        return (int) (this.fitness - o.fitness);
    }
}
