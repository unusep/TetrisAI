package main.tetris.heuristics;

import main.tetris.engine.State;
import main.tetris.engine.TetrisSimulator;

public class SumOfHeightOfColumnsHeuristic implements IHeuristic {

    public double getValue(int[][] board, int[] top) {
        double count = 0.0;
        for (int c = 0; c < top.length; c++){
            count += top[c];
        }
        return count;
    }

    @Override
    public String toString(){
        return "SumOfHeightOfColumnsHeuristic";
    }

    @Override
    public double getValue(int[] move, State s) {
        TetrisSimulator simulator = new TetrisSimulator(s);
        simulator.makeMove(move);
        return getValue(simulator.getField(), simulator.getTop());
    }
}
