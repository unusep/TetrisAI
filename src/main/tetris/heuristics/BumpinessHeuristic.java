package main.tetris.heuristics;

import main.tetris.engine.TetrisSimulator;

public class BumpinessHeuristic implements IHeuristic {

    @Override
    public double getValue(TetrisSimulator state) {
        double count = 0;
        for (int col = 1; col < TetrisSimulator.COLS; col++){
            int[] topRows = state.getTop();
            count += Math.abs(topRows[col - 1] - topRows[col]);
        }
        return count;
    }

}
