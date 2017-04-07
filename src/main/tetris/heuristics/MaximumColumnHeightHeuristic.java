package main.tetris.heuristics;

import main.tetris.engine.TetrisSimulator;

public class MaximumColumnHeightHeuristic implements IHeuristic {

    @Override
    public double getValue(TetrisSimulator state) {
        double max = 0;
        for (int col = 0; col < TetrisSimulator.COLS; col++){
            int[] topRows = state.getTop();
            if (topRows[col] > max) max = topRows[col];
        }
        return max;
    }

}
