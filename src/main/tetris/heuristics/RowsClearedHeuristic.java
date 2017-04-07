package main.tetris.heuristics;

import main.tetris.engine.TetrisSimulator;

public class RowsClearedHeuristic implements IHeuristic {

    @Override
    public double getValue(TetrisSimulator state) {
        return state.getRowsCleared();
    }

}
