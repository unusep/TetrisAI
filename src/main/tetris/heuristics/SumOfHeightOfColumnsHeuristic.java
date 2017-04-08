package main.tetris.heuristics;

public class SumOfHeightOfColumnsHeuristic implements IHeuristic {

    @Override
    public double getValue(boolean[][] board, int[] top, int rowsCleared) {
        double count = 0.0;
        for (int c = 0; c < top.length; c++){
            count += top[c];
        }
        return count;
    }

}
