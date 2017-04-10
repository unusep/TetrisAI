package main.tetris.heuristics;

import java.util.ArrayList;

public class FlatteningCoefficientHeuristic implements IHeuristic {

    @Override
    public String toString(){
        return "FlatteningCoefficientHeuristic";
    }

    @Override
    public double getValue(boolean[][] board, int[] top, int rowsCleared, boolean[][] oldBoard,
            int oldRowsCleared, int[][][] pTop, int[][][] pBottom, int[][] pWidth, int pieceIndex,
            int rotationIndex, int leftPosition) {
        ArrayList<Integer> filledCoord = new ArrayList<Integer>();
        double flatteningCoeff = 0;
        int[] rows = new int[board.length + 4];
        
        int height = top[leftPosition] - pBottom[pieceIndex][rotationIndex][0];
        
        //for each column beyond the first in the piece
        for (int c = 0; c < pWidth[pieceIndex][rotationIndex]; c++) {
            height = Math.max(height, top[leftPosition + c] - pBottom[pieceIndex][rotationIndex][c]);
        }
        for (int i = 0; i < pWidth[pieceIndex][rotationIndex]; i++) {

            //from bottom to top of brick
            for (int h = height + pBottom[pieceIndex][rotationIndex][i] + 1;
                 h <= height + pTop[pieceIndex][rotationIndex][i]; h++) {

                // Tile is being filled here.
                //board[h][i + leftPosition] = true;
                rows[h] |= (1 << (i + leftPosition));

                filledCoord.add(h * board[0].length + i + leftPosition);
            }
        }
        for (int coord : filledCoord) {
            int row = coord / board[0].length;
            int col = coord % board[0].length;

            int left = row * board[0].length + col - 1;
            int right = row * board[0].length + col + 1;
            int down = (row - 1) * board[0].length + col;

            // left side
            if (col != 0 && !filledCoord.contains(left) && ((rows[row] & (1 << (col - 1))) > 0)) {
                flatteningCoeff ++;
            }

            // right side
            if (col != board[0].length - 1 && !filledCoord.contains(right) && ((rows[row] & (1 << (col + 1))) > 0)) {
                flatteningCoeff ++;
            }

            // down side
            if (row != 0 && !filledCoord.contains(down) && ((rows[row - 1] & (1 << col)) > 0)) {
                flatteningCoeff ++;
            }
        }
        return flatteningCoeff;
    }
}
