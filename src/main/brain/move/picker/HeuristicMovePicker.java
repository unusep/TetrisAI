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
        double bestScore = -Double.MAX_VALUE;
        int[][] legalMoves = s.legalMoves();
        int[] bestMove = null;
        for (int[] move : legalMoves){
            double score = evaluateMove(move, s);
            if (score > bestScore){
                bestScore = score;
                bestMove = move;
            }
        }
        return bestMove;
    }

    /**
     * evaluates the value of using a particular move in a particular state
     * using the heuristics
     * @param move
     * @param state s
     * @return score of move
     */
    private double evaluateMove(int[] move, State s) {
        double score = 0.0;
        for (int i = 0; i < Math.min(weights.size(),  heuristics.size()); i++){
            score += weights.get(i) * heuristics.get(i).getValue(move, s);
        }
        return score;
    }

}
