package main.tetris.heuristics;

public abstract class Heuristic {
    private int[][] board;
    
    public Heuristic(int[][] board){
        this.board = board;
    }
    
    public void setBoard(int[][] board){
        this.board = board;
    }
    
    public abstract double getValue();
}