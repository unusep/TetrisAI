package main.tetris.heuristics;

public class SumOfHeightOfBlocksHeuristic implements IHeuristic {

    @Override
    public double getValue(boolean[][] board, int[] top, int rowsCleared, boolean[][] oldBoard,
            int oldRowsCleared, int[][][] pTop, int[][][] pBottom, int[][] pWidth, int pieceIndex, int rotationIndex, int leftPosition){
        double res = 0;
        for (int r = 0; r < board.length; r++){
            for (int c = 0; c < board[0].length; c++){
                if (board[r][c]) res += r;
            }
        }
        return res;
        
    }
    
    @Override
    public String toString(){
        return "SumOfHeightOfBlocksHeuristic";
    }

}
