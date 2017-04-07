package main.brain.move.picker;

import java.util.ArrayList;

import main.tetris.engine.State;
import main.tetris.engine.TetrisSimulator;
import main.tetris.heuristics.IHeuristic;

/**
 * HeuristicMovePicker will pick the best move in the given board configuration
 * given a heuristic and its weights
 */
public class HeuristicMovePicker implements IMovePicker {
    private ArrayList<IHeuristic> heuristics;
    private ArrayList<Double> weights;
    private final double LOST_SCORE = Double.MIN_VALUE;
    
    public HeuristicMovePicker(ArrayList<Double> weights, ArrayList<IHeuristic> heuristics) {
        this.weights = weights;
        this.heuristics = heuristics;
    }

    @Override
    public int[] pickBest(TetrisSimulator simulator, int[][] legalMoves) {
        double bestScore = 0.0;
        int[] bestMove = legalMoves[0]; // TODO: may cause null pointer exception if there are no legal moves
        simulator.saveState();
        for (int[] move : legalMoves){
            simulator.makeMove(move);
            double score;
            if (simulator.hasLost()){
                score = LOST_SCORE;
            } else {
                score = evaluateBoard(simulator);
            }
            // update score and move if it's better than current maximum
            if (score > bestScore){
                bestScore = score;
                bestMove = move;
            }
            simulator.undoMove();
        }
        return bestMove;
    }

    /**
     * Evaluate the value of the board using the weights of each heuristic
     * @param weights
     * @param board
     * @return value of board
     */
    private double evaluateBoard(TetrisSimulator state) {
        double score = 0.0;
        for (int i = 0; i < Math.min(weights.size(), heuristics.size()); i++){
            score += weights.get(i) * heuristics.get(i).getValue(state);
        }
        return score;
    }

}
