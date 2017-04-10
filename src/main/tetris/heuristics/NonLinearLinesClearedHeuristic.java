package main.tetris.heuristics;

public class NonLinearLinesClearedHeuristic implements IHeuristic {

    @Override
    public double getValue(boolean[][] board, int[] top, int rowsCleared, boolean[][] oldBoard,
            int oldRowsCleared, int[][][] pTop, int[][][] pBottom, int[][] pWidth, int pieceIndex, int rotationIndex, int leftPosition){
        return Math.pow(2, rowsCleared - oldRowsCleared);
    }

    @Override
    public String toString(){
        return "NonLinearLinesClearedHeuristic";
    }
}
