package main.tetris.heuristics;

import main.tetris.engine.State;

public class FilledLinesHeuristic implements IHeuristic {

    @Override
    public double getValue(State state) {
        return state.getRowsCleared();
    }

}
