package main.tetris.heuristics;

public class RowsClearedHeuristic implements IHeuristic {

    public double getValue(int rowsCleared) {
        return rowsCleared;
    }
    
    @Override
    public String toString(){
        return "RowsClearedHeuristic";
    }

    @Override
    public double getValue(boolean[][] board, int[] top, int rowsCleared, boolean[][] oldBoard,
            int oldRowsCleared, int[][][] pTop, int[][][] pBottom, int[][] pWidth, int pieceIndex, int rotationIndex, int leftPosition){
        return getValue(rowsCleared);
    }

}
