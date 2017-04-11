package main.tetris.heuristics;

import main.tetris.engine.State;
import main.tetris.engine.TetrisSimulator;

public class MaximumColumnHeightHeuristic implements IHeuristic {

    public double getValue(int[][] board, int[] top) {
        double max = 0;
        for (int col = 0; col < board[0].length; col++){
            if (top[col] > max) max = top[col];
        }
        return max;
    }
    
    @Override
    public String toString(){
        return "MaximumColumnHeightHeuristic";
    }

    @Override
    public double getValue(int[] move, State s) {
        TetrisSimulator simulator = new TetrisSimulator(s);
        simulator.makeMove(move);
        return getValue(simulator.getField(), simulator.getTop());
    }

}
