package main.tetris.heuristics;

public interface IHeuristic {
    public double getValue(boolean[][] board, int[] top, int rowsCleared);
}