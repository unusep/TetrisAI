package main.tetris.heuristics;

import main.tetris.engine.State;

public interface IHeuristic {

    public double getValue(int[] move, State s);
}