package main.tetris.heuristics;

import main.tetris.engine.TetrisSimulator;

public interface IHeuristic {
    public abstract double getValue(TetrisSimulator state);
}