package main.tetris.heuristics;

import main.tetris.engine.State;
import main.tetris.engine.TetrisSimulator;

public class NonLinearBumpinessHeuristic implements IHeuristic {

    public double getValue(int[][] board, int[] top) {
        double count = 0;
        for (int col = 1; col < board[0].length; col++){
            count += Math.pow(top[col-1] - top[col], 2);
        }
        return count;
    }
    
    @Override
    public String toString(){
        return "BumpinessHeuristic";
    }

    @Override
    public double getValue(int[] move, State s) {
        TetrisSimulator simulator = new TetrisSimulator(s);
        simulator.makeMove(move);
        return getValue(simulator.getField(), simulator.getTop());
    }


}
