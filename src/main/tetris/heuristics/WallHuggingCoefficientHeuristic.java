package main.tetris.heuristics;

import main.tetris.engine.State;

public class WallHuggingCoefficientHeuristic implements IHeuristic {

    @Override
    public String toString(){
        return "WallHuggingCoefficientHeuristic";
    }

    @Override
    public double getValue(int[] move, State s) {
        int orient = move[State.ORIENT];
        int slot = move[State.ORIENT];
        int nextPiece = s.getNextPiece();
        int[][][] pBottom = State.getpBottom();
        int[][] pWidth = State.getpWidth();
        int[][][] pTop = State.getpTop();
        double res = 0.0;
        for(int c = 0; c < pWidth[nextPiece][orient];c++) {
            for(int h = pBottom[nextPiece][orient][c]; h < pTop[nextPiece][orient][c]; h++) {
                if (slot + c == 0 || slot + c == State.COLS -1) res++;
            }
        }
        return res;
    }

}
