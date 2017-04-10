package main.tetris.heuristics;

public class SumOfHeightOfColumnsHeuristic implements IHeuristic {

    public double getValue(boolean[][] board, int[] top, int rowsCleared) {
        double count = 0.0;
        for (int c = 0; c < top.length; c++){
            count += top[c];
        }
        return count;
    }

    @Override
    public String toString(){
        return "SumOfHeightOfColumnsHeuristic";
    }

    @Override
    public double getValue(boolean[][] board, int[] top, int rowsCleared, boolean[][] oldBoard,
            int oldRowsCleared, int[][][] pTop, int[][][] pBottom, int[][] pWidth, int pieceIndex, int rotationIndex, int leftPosition){
        return getValue(board, top, rowsCleared);
    }
}
