package main.tetris.heuristics;

public class RowsClearedHeuristic implements IHeuristic {

    @Override
    public double getValue(boolean[][] board, int[] top, int rowsCleared) {
        return rowsCleared;
    }
    
    @Override
    public String toString(){
        return "RowsClearedHeuristic";
    }

}
