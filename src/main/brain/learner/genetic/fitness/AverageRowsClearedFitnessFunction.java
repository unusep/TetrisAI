package main.brain.learner.genetic.fitness;

import main.brain.learner.genetic.Gene;
import main.brain.move.picker.HeuristicMovePicker;
import main.brain.move.picker.IMovePicker;
import main.tetris.engine.TetrisSimulator;
import main.tetris.heuristics.IHeuristic;

public class AverageRowsClearedFitnessFunction implements IFitnessFunction<IHeuristic> {
    private int numPieces;
    private int numGames;
    
    /**
     * @param the number of games to play (we will take the average fitness score)
     * @param cutPieces the number of pieces to play before stopping the game
     */
    public AverageRowsClearedFitnessFunction(int numPieces, int numGames) {
        this.numGames = numGames;
        this.numPieces = numPieces;
    }
    
    /**
     * evaluate the fitness of the gene using the movePicker to pick the best move at each step
     */
    @Override
    public double evaluateFitness(Gene<IHeuristic> gene) {
        double totalFitness = 0.0;
        for (int i = 0; i < numGames; i++){
            TetrisSimulator simulator = new TetrisSimulator();
            IMovePicker movePicker = new HeuristicMovePicker(gene.getChromosomeWeights(), gene.getChromsomes());
            for (int j = 0; j < numPieces; j++){
                if (simulator.hasLost()) break;
                int[] bestMove = movePicker.pickBest(simulator, simulator.legalMoves());
                simulator.makeMove(bestMove);
            }
            totalFitness += simulator.getRowsCleared();
        }
        return totalFitness/numGames;
    }

}
