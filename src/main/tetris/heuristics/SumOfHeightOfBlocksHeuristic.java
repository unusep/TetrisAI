package main.tetris.heuristics;

import main.tetris.engine.State;
import main.tetris.engine.TetrisSimulator;

public class SumOfHeightOfBlocksHeuristic implements IHeuristic {

    @Override
    public String toString(){
        return "SumOfHeightOfBlocksHeuristic";
    }


    @Override
    public double getValue(int[] move, State s) {
        TetrisSimulator simulator = new TetrisSimulator(s);
        simulator.makeMove(move);
        double res = 0;
        int[][] board = simulator.getField();
        for (int r = 0; r < State.ROWS; r++){
            for (int c = 0; c < State.COLS; c++){
                if (board[r][c] != 0) res += r;
            }
        }
        return res;
    }

}
