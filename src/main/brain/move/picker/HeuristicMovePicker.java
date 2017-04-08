package main.brain.move.picker;

import java.util.ArrayList;

import main.tetris.engine.State;
import main.tetris.heuristics.IHeuristic;

/**
 * HeuristicMovePicker will pick the best move in the given board configuration
 * given a heuristic and its weights
 */
public class HeuristicMovePicker implements IMovePicker {
    private ArrayList<IHeuristic> heuristics;
    private ArrayList<Double> weights;
    
    public HeuristicMovePicker(ArrayList<Double> weights, ArrayList<IHeuristic> heuristics) {
        this.weights = weights;
        this.heuristics = heuristics;
    }

    @Override
    public int[] pickBest(State s) {
        double bestValue = Double.MIN_VALUE;
        
        boolean[][] newField = new boolean[State.ROWS][State.COLS];
        int[] newTop = new int[State.COLS];
        int bestRot = 0;
        int bestPos = 0;

        int nextPiece = s.getNextPiece();
        int[][] legalMoves = s.legalMoves();
        
        for (int i = 0; i < legalMoves.length; i++) {
            double score = 0.0;
            int rot = legalMoves[i][State.ORIENT];
            int pos = legalMoves[i][State.SLOT];
            int rowsCleared = performMove(s, newField, newTop, nextPiece, rot, pos);
            score = evaluateBoard(newField, newTop, rowsCleared);
            if (score > bestValue){
                bestValue = score;
                bestRot = rot;
                bestPos = pos;
            }
        }

        return new int[] { bestRot, bestPos };
    }

    private int performMove(State s, boolean[][] newField, int[] newTop, int piece, int rot, int pos) {
        // Perform Deep Copy
        for (int i = 0; i < State.ROWS; i++) {
            for (int j = 0; j < State.COLS; j++) {
                newField[i][j] = s.getField()[i][j] != 0;
            }
        }
        for (int j = 0; j < State.COLS; j++) {
            newTop[j] = s.getTop()[j];
        }

        // height if the first column makes contact
        int height = newTop[pos] - State.getpBottom()[piece][rot][0];
        // for each column beyond the first in the piece
        for (int c = 0; c < State.getpWidth()[piece][rot]; c++) {
            height = Math.max(height, newTop[pos + c] - State.getpBottom()[piece][rot][c]);
        }

        // check if game ended
        if (height + State.getpHeight()[piece][rot] >= State.ROWS) {
            return -10;
        }

        // for each column in the piece - fill in the appropriate blocks
        for (int i = 0; i < State.getpWidth()[piece][rot]; i++) {

            // from bottom to top of brick
            for (int h = height + State.getpBottom()[piece][rot][i]; h < height + State.getpTop()[piece][rot][i]; h++) {
                newField[h][i + pos] = true;
            }
        }

        // adjust top
        for (int c = 0; c < State.getpWidth()[piece][rot]; c++) {
            newTop[pos + c] = height + State.getpTop()[piece][rot][c];
        }

        int rowsCleared = 0;

        // check for full rows - starting at the top
        for (int r = height + State.getpHeight()[piece][rot] - 1; r >= height; r--) {
            // check all columns in the row
            boolean full = true;
            for (int c = 0; c < State.COLS; c++) {
                if (newField[r][c] == false) {
                    full = false;
                    break;
                }
            }
            // if the row was full - remove it and slide above stuff down
            if (full) {
                rowsCleared++;
                // for each column
                for (int c = 0; c < State.COLS; c++) {
                    // slide down all bricks
                    for (int i = r; i < newTop[c]; i++) {
                        newField[i][c] = newField[i + 1][c];
                    }
                    // lower the top
                    newTop[c]--;
                    while (newTop[c] >= 1 && newField[newTop[c] - 1][c] == false)
                        newTop[c]--;
                }
            }
        }
        return rowsCleared;
    }

    /**
     * Evaluate the value of the state given a board, rows cleared, and top rows
     * using the weights of each heuristic
     * @param board
     * @param rowsCleared 
     * @param top 
     * @return value of board
     */
    private double evaluateBoard(boolean[][] board, int[] top, int rowsCleared) {
        double score = 0.0;
        for (int i = 0; i < Math.min(weights.size(), heuristics.size()); i++){
            score += weights.get(i) * heuristics.get(i).getValue(board, top, rowsCleared);
        }
        return score;
    }

}
