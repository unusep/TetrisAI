package main.tetris.heuristics;

import main.tetris.engine.State;

public interface IHeuristic {
    public abstract double getValue(State state);
}