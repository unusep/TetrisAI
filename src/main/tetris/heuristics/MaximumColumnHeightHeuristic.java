package main.tetris.heuristics;

public class MaximumColumnHeightHeuristic implements IHeuristic {

    @Override
    public double getValue(boolean[][] board, int[] top, int rowsCleared) {
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
}
