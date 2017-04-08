package main.tetris.heuristics;

public class BumpinessHeuristic implements IHeuristic {

    @Override
    public double getValue(boolean[][] board, int[] top, int rowsCleared) {
        double count = 0;
        for (int col = 1; col < board[0].length; col++){
            count += Math.abs(top[col-1] - top[col]);
        }
        return count;
    }

}
