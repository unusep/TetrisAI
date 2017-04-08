package main.tetris.heuristics;

public class HolesHeuristic implements IHeuristic {


    @Override
    public double getValue(boolean[][] board, int[] top, int rowsCleared) {
        double count = 0;
        for (int col = 0; col < board[0].length; col++){
            for (int row = 0; row < top[col]; row++){
                if (board[row][col] == false){
                    count++;
                }
            }
        }
        return count;
    }

    @Override
    public String toString(){
        return "HolesHeuristic";
    }
}
