package main.tetris.heuristics;

import main.tetris.engine.State;
import main.tetris.engine.TetrisSimulator;

public class HolesHeuristic implements IHeuristic {


    public double getValue(int[][] board, int[] top) {
        double count = 0;
        for (int col = 0; col < board[0].length; col++){
            for (int row = 0; row < top[col]; row++){
                if (board[row][col] == 0){
                    count++;
                }
            }
        }
        return count;
    }

    @Override
    public String toString(){
        return "HolesHeuristic";
    }


    @Override
    public double getValue(int[] move, State s) {
        TetrisSimulator simulator = new TetrisSimulator(s);
        simulator.makeMove(move);
        return getValue(simulator.getField(), simulator.getTop());
    }
}
