package main.tetris.heuristics;

public class BlockadesHeuristic implements IHeuristic {

    @Override
    public double getValue(boolean[][] board, int[] top, int rowsCleared, boolean[][] oldBoard,
            int oldRowsCleared, int[][][] pTop, int[][][] pBottom, int[][] pWidth, int pieceIndex, int rotationIndex, int leftPosition){
        double res = 0;
        for (int c = 0; c < board[0].length; c++){
            boolean holeFound = false;
            for (int r = 0; r <= top[c]; r++){
                if (holeFound) {
                    res += board[r][c] ? 1 : 0;
                } else {
                    if (board[r][c] == false) holeFound = true;
                }
            }
        }
        return res;
    }

    @Override
    public String toString(){
        return "BlockadesHeuristic";
    }
}
