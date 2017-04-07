package main.tetris.heuristics;

public interface IHeuristic {
    public abstract double getValue(int[][] board);
}