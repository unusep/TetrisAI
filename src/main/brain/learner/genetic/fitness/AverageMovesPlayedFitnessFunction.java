package main.brain.learner.genetic.fitness;

import main.brain.learner.genetic.Gene;
import main.brain.move.picker.HeuristicMovePicker;
import main.brain.move.picker.IMovePicker;
import main.tetris.engine.State;
import main.tetris.heuristics.IHeuristic;

public class AverageMovesPlayedFitnessFunction implements IFitnessFunction<IHeuristic> {
    private int numPieces;
    private int numGames;
    
    /**
     * @param numPieces the number of pieces to play before stopping the game
     * @param numGames the number of games to play (we will take the average fitness score)
     */
    public AverageMovesPlayedFitnessFunction(int numPieces, int numGames) {
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
            State simulator = new State();
            
            IMovePicker movePicker = new HeuristicMovePicker(gene.getChromosomeWeights(), gene.getChromsomes());
            int j = 0;
            for (j = 0; j < numPieces; j++){
                if (simulator.hasLost()) break;
                int[] bestMove = movePicker.pickBest(simulator);
                simulator.makeMove(bestMove);
            }
            totalFitness += j;
        }
        return totalFitness/ (double) numGames;
    }

}
