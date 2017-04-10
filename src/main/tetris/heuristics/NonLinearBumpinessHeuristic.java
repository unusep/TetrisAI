package main.tetris.heuristics;

public class NonLinearBumpinessHeuristic implements IHeuristic {

    public double getValue(boolean[][] board, int[] top, int rowsCleared) {
        double count = 0;
        for (int col = 1; col < board[0].length; col++){
            count += Math.pow(top[col-1] - top[col], 2);
        }
        return count;
    }
    
    @Override
    public String toString(){
        return "BumpinessHeuristic";
    }

    @Override
    public double getValue(boolean[][] board, int[] top, int rowsCleared, boolean[][] oldBoard,
            int oldRowsCleared, int[][][] pTop, int[][][] pBottom, int[][] pWidth, int pieceIndex, int rotationIndex, int leftPosition){
        // TODO Auto-generated method stub
        return getValue(board, top, rowsCleared);
    }

}
