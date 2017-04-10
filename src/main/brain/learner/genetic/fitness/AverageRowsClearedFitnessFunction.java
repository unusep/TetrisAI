package main.brain.learner.genetic.fitness;

import main.brain.learner.genetic.Gene;
import main.brain.move.picker.HeuristicMovePicker;
import main.brain.move.picker.IMovePicker;
import main.tetris.engine.State;
import main.tetris.engine.TFrame;
import main.tetris.heuristics.IHeuristic;

public class AverageRowsClearedFitnessFunction implements IFitnessFunction<IHeuristic> {
    private int numPieces;
    private int numGames;
    
    /**
     * @param numPieces the number of pieces to play before stopping the game
     * @param numGames the number of games to play (we will take the average fitness score)
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
            State simulator = new State();

//            new TFrame(simulator);
            
            IMovePicker movePicker = new HeuristicMovePicker(gene.getChromosomeWeights(), gene.getChromsomes());
            
            for (int j = 0; j < numPieces; j++){
                if (simulator.hasLost()) break;
                int[] bestMove = movePicker.pickBest(simulator);
                simulator.makeMove(bestMove);
//                
//                simulator.draw();
//                simulator.drawNext(0,0);
//                try {
//                    Thread.sleep(100);
//                } catch (InterruptedException e) {
//                    e.printStackTrace();
//                }
            }
            totalFitness += simulator.getRowsCleared();
        }
        return totalFitness/numGames;
    }

}
