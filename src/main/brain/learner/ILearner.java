package main.brain.learner;

import main.brain.move.picker.IMovePicker;

public interface ILearner {
    
    public abstract IMovePicker getMovePicker();
    
    public abstract void trainLearner(int iterations);
}