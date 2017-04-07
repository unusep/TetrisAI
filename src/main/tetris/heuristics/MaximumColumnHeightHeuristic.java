package main.tetris.heuristics;

import main.tetris.engine.State;

public class MaximumColumnHeightHeuristic implements IHeuristic {

    @Override
    public double getValue(State state) {
        double max = 0;
        for (int col = 0; col < State.COLS; col++){
            int[] topRows = state.getTop();
            if (topRows[col] > max) max = topRows[col];
        }
        return max;
    }

}
