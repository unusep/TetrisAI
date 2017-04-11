package main.player;
import java.util.ArrayList;

import main.brain.learner.HeuristicGeneticLearner;
import main.brain.learner.ILearner;
import main.brain.learner.genetic.crossover.*;
import main.brain.learner.genetic.fitness.*;
import main.brain.learner.genetic.mutator.*;
import main.brain.learner.genetic.selector.*;
import main.tetris.engine.*;
import main.tetris.heuristics.NonLinearBumpinessHeuristic;
import main.tetris.heuristics.NonLinearLinesClearedHeuristic;
import main.tetris.heuristics.BlockadesHeuristic;
import main.tetris.heuristics.FlatteningCoefficientHeuristic;
import main.tetris.heuristics.FloorHuggingCoefficientHeuristic;
import main.tetris.heuristics.HolesHeuristic;
import main.tetris.heuristics.IHeuristic;
import main.tetris.heuristics.SumOfHeightOfBlocksHeuristic;
import main.tetris.heuristics.SumOfHeightOfColumnsHeuristic;
import main.tetris.heuristics.WallHuggingCoefficientHeuristic;

public class PlayerSkeleton {
    public static final String POPULATION_FILEPATH = "population.txt";
    private ILearner brain;
    
	public int[] pickMove(State simulator, int[][] legalMoves) {
		return brain.getMovePicker().pickBest(simulator);
	}
	
	/**
	 * Pass as arguments 1000 50 to train a population of 1000 people for 50 generations
	 * @param args
	 */
	public static void main(String[] args) {
		State s = new State();
		new TFrame(s);
		PlayerSkeleton p = new PlayerSkeleton();
		
		// set up genetic learner variables
		int TRAINING_MAX_PIECES = 10000;
		int FITNESS_BEST_OF_NUM_GAMES = 1;
		
		int POPULATION_SIZE = 1000;
		int NUMBER_OF_GENERATIONS = 100;
		double PERCENTAGE_TO_CULL = 0.3; // 0.0-0.49 what percentage of the population to cull
		double CHROMOSOME_MUTATION_PROBABILITY = 0.2; // 0.0-1.0 how likely it is for a particular chromosome in a gene to mutate 
		double GENE_MUTATION_PROBABILITY = 0.2; // how like it is for a particular gene in the population to be selected for mutation
		
		// end of set up genetic learner variables 
		
		ArrayList<IHeuristic> heuristics = new ArrayList<IHeuristic>();
        heuristics.add(new NonLinearLinesClearedHeuristic());
        heuristics.add(new SumOfHeightOfBlocksHeuristic());
		heuristics.add(new SumOfHeightOfColumnsHeuristic());
		heuristics.add(new NonLinearBumpinessHeuristic());
		heuristics.add(new HolesHeuristic());
        heuristics.add(new BlockadesHeuristic());
        heuristics.add(new WallHuggingCoefficientHeuristic());
        heuristics.add(new FloorHuggingCoefficientHeuristic());
        heuristics.add(new FlatteningCoefficientHeuristic());
		ICrossoverOperator<IHeuristic> crossOverOperator = new SinglePointCrossover();
		IFitnessFunction<IHeuristic> fitnessFunction 
		    = new AverageRowsClearedFitnessFunction(TRAINING_MAX_PIECES, FITNESS_BEST_OF_NUM_GAMES);
		IMutationOperator<IHeuristic> mutationOperator = new UniformMutation<IHeuristic>(CHROMOSOME_MUTATION_PROBABILITY);
		IPopulationSelector<IHeuristic> populationSelector = new TruncationFitnessSelector<IHeuristic>();
		
        p.brain = new HeuristicGeneticLearner(
                POPULATION_SIZE, 
		        POPULATION_FILEPATH, 
		        heuristics, 
		        crossOverOperator, 
		        fitnessFunction, 
		        mutationOperator, 
		        populationSelector,
		        PERCENTAGE_TO_CULL,
		        GENE_MUTATION_PROBABILITY);
        
        p.brain.trainLearner(NUMBER_OF_GENERATIONS);
        
		while(!s.hasLost()) {
			s.makeMove(p.pickMove(s, s.legalMoves()));
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
