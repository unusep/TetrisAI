package main.tetris.heuristics;

public class FloorHuggingCoefficientHeuristic implements IHeuristic {

    @Override
    public double getValue(boolean[][] board, int[] top, int rowsCleared, boolean[][] oldBoard,
            int oldRowsCleared, int[][][] pTop, int[][][] pBottom, int[][] pWidth, int pieceIndex,
            int rotationIndex, int leftPosition) {
        double floorHuggingCoeff = 0;
        //height if the first column makes contact
        int height = top[leftPosition] - pBottom[pieceIndex][rotationIndex][0];
        //for each column beyond the first in the piece
        for (int c = 0; c < pWidth[pieceIndex][rotationIndex]; c++) {
            height = Math.max(height, top[leftPosition + c] - pBottom[pieceIndex][rotationIndex][c]);
            // 2. Floor-Hugging Coefficient
            if (top[leftPosition + c] == -1 && pBottom[pieceIndex][rotationIndex][c] == 0) {
                floorHuggingCoeff++;
            }
            // 2. Floor-Hugging Coefficient
        }
        return floorHuggingCoeff;
    }
    
    @Override
    public String toString(){
        return "FloorHuggingCoefficientHeuristic";
    }

}
