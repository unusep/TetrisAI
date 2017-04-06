package main.brain.learner;

import java.util.ArrayList;

import main.tetris.heuristics.Heuristic;

public abstract class Learner {
    ArrayList<Heuristic> heuristics;
    
    public Learner(ArrayList<Heuristic> heuristics){
        this.heuristics = heuristics;
    }
    
    public abstract ArrayList<Double> getWeights();
    
    public abstract void trainLearner(int iterations);
}