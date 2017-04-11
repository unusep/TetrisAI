package main.tetris.heuristics;

import main.tetris.engine.State;
import main.tetris.engine.TetrisSimulator;

public class BlockadesHeuristic implements IHeuristic {

    public double getValue(int[][] board, int[] top){
        double res = 0;
        for (int c = 0; c < board[0].length; c++){
            boolean holeFound = false;
            for (int r = 0; r <= top[c]; r++){
                if (holeFound) {
                    res += board[r][c] == 0 ? 0 : 1;
                } else {
                    if (board[r][c] == 0) holeFound = true;
                }
            }
        }
        return res;
    }

    @Override
    public String toString(){
        return "BlockadesHeuristic";
    }

    @Override
    public double getValue(int[] move, State s) {
        TetrisSimulator simulator = new TetrisSimulator(s);
        simulator.makeMove(move);
        return getValue(simulator.getField(), simulator.getTop());
    }
}
