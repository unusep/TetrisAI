package main.player;
import java.util.ArrayList;

import main.brain.learner.HeuristicGeneticLearner;
import main.brain.learner.ILearner;
import main.brain.learner.genetic.Gene;
import main.brain.learner.genetic.crossover.*;
import main.brain.learner.genetic.fitness.*;
import main.brain.learner.genetic.mutator.*;
import main.brain.learner.genetic.selector.*;
import main.brain.move.picker.HeuristicMovePicker;
import main.brain.move.picker.IMovePicker;
import main.tetris.engine.*;
import main.tetris.heuristics.IHeuristic;

public class PlayerSkeleton {
    public static final String POPULATION_FILEPATH = "population.txt";
    
    private ILearner<Gene<IHeuristic>> brain;
    
	//implement this function to have a working system
	public int[] pickMove(TetrisSimulator simulator, int[][] legalMoves) {
		return brain.getMovePicker().pickBest(simulator, legalMoves);
	}
	
	public static void main(String[] args) {
		TetrisSimulator s = new TetrisSimulator();
		new TFrame(s);
		PlayerSkeleton p = new PlayerSkeleton();
		
		// set up genetic learner variables
		String pathToPopulation = POPULATION_FILEPATH;
		int populationSize = Integer.parseInt(args[0]);
        int trainingNumPieces = Integer.parseInt(args[1]);
        int trainingNumGames = 1;
		
		ArrayList<IHeuristic> heuristics = new ArrayList<IHeuristic>();
		ICrossoverOperator<IHeuristic> crossOverOperator;
		IFitnessFunction<IHeuristic> fitnessFunction 
		    = new AverageRowsClearedFitnessFunction(trainingNumPieces, trainingNumGames);
		IMutationOperator<IHeuristic> mutationOperator;
		IPopulationSelector<IHeuristic> populationSelector;
		
        p.brain = new HeuristicGeneticLearner(
		        populationSize, 
		        pathToPopulation, 
		        heuristics, 
		        crossOverOperator, 
		        fitnessFunction, 
		        mutationOperator, 
		        populationSelector);
        
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
