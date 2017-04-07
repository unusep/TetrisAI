package main.player;
import java.util.ArrayList;

import main.brain.learner.HeuristicGeneticLearner;
import main.brain.learner.ILearner;
import main.brain.learner.genetic.crossover.*;
import main.brain.learner.genetic.fitness.*;
import main.brain.learner.genetic.mutator.*;
import main.brain.learner.genetic.selector.*;
import main.tetris.engine.*;
import main.tetris.heuristics.IHeuristic;

public class PlayerSkeleton {
    public static final String POPULATION_FILEPATH = "population.txt";
    private ILearner brain;
    
	public int[] pickMove(TetrisSimulator simulator, int[][] legalMoves) {
		return brain.getMovePicker().pickBest(simulator, legalMoves);
	}
	
	/**
	 * Pass as arguments 1000 50 to train a population of 1000 people for 50 generations
	 * @param args
	 */
	public static void main(String[] args) {
		TetrisSimulator s = new TetrisSimulator();
		new TFrame(s);
		PlayerSkeleton p = new PlayerSkeleton();
		
		// set up genetic learner variables
		int populationSize = Integer.parseInt(args[0]);
		int trainingIterations = Integer.parseInt(args[1]);
		
        int trainingNumPieces = 10000; // the number of pieces to play before we cut off the fitness function
        int trainingNumGames = 1; // the number of games to play when evaluating the fitness function (will take the average fitness level) 
		
		ArrayList<IHeuristic> heuristics = new ArrayList<IHeuristic>();
		ICrossoverOperator<IHeuristic> crossOverOperator;
		IFitnessFunction<IHeuristic> fitnessFunction 
		    = new AverageRowsClearedFitnessFunction(trainingNumPieces, trainingNumGames);
		IMutationOperator<IHeuristic> mutationOperator;
		IPopulationSelector<IHeuristic> populationSelector;
		
        p.brain = new HeuristicGeneticLearner(
		        populationSize, 
		        POPULATION_FILEPATH, 
		        heuristics, 
		        crossOverOperator, 
		        fitnessFunction, 
		        mutationOperator, 
		        populationSelector
		        );
        
        p.brain.trainLearner(trainingIterations);
        
		while(!s.hasLost()) {
			s.makeMove(p.pickMove(s,s.legalMoves()));
			s.draw();
			s.drawNext(0,0);
			try {
				Thread.sleep(300);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		System.out.println("You have completed "+s.getRowsCleared()+" rows.");
	}
	
}
