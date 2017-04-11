package main.tetris.heuristics;

import main.tetris.engine.State;

public class FloorHuggingCoefficientHeuristic implements IHeuristic {

    @Override
    public String toString(){
        return "FloorHuggingCoefficientHeuristic";
    }

    @Override
    public double getValue(int[] move, State s) {
        int orient = move[State.ORIENT];
        int slot = move[State.ORIENT];
        int nextPiece = s.getNextPiece();
        int[][][] pBottom = State.getpBottom();
        int[][] pWidth = State.getpWidth();
        int[][][] pTop = State.getpTop();
        int[] top = s.getTop();
        
        double res = 0.0;
        int height = top[slot]-pBottom[nextPiece][orient][0];
        //for each column beyond the first in the piece
        for(int c = 1; c < pWidth[nextPiece][orient];c++) {
            height = Math.max(height,top[slot+c]-pBottom[nextPiece][orient][c]);
        }
        
        //for each column in the piece - fill in the appropriate blocks
        for(int i = 0; i < pWidth[nextPiece][orient]; i++) {
            //from bottom to top of brick
            for(int h = height+pBottom[nextPiece][orient][i]; h < height+pTop[nextPiece][orient][i]; h++) {
                if (h == 0) res++;
            }
        }
        return res;
    }

}
