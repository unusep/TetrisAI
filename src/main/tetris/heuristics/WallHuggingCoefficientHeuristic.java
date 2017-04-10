package main.tetris.heuristics;

public class WallHuggingCoefficientHeuristic implements IHeuristic {

    @Override
    public double getValue(boolean[][] board, int[] top, int rowsCleared, boolean[][] oldBoard,
            int oldRowsCleared, int[][][] pTop, int[][][] pBottom, int[][] pWidth, int pieceIndex, int rotationIndex, int leftPosition) {
        double wallHuggingCoeff = 0;
        if (leftPosition == 0) {
            wallHuggingCoeff += pTop[pieceIndex][rotationIndex][0] - pBottom[pieceIndex][rotationIndex][0];
        }
        if (leftPosition + pWidth[pieceIndex][rotationIndex] - 1 == board[0].length - 1) {
            wallHuggingCoeff += pTop[pieceIndex][rotationIndex][pWidth[pieceIndex][rotationIndex] - 1] -
                    pBottom[pieceIndex][rotationIndex][pWidth[pieceIndex][rotationIndex] - 1];
        }
        return wallHuggingCoeff;
    }
    
    @Override
    public String toString(){
        return "WallHuggingCoefficientHeuristic";
    }

}
