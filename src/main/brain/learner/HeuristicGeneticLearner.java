package main.brain.learner;

import java.util.ArrayList;

import main.brain.learner.genetic.Gene;
import main.brain.learner.genetic.Population;
import main.brain.learner.genetic.crossover.ICrossoverOperator;
import main.brain.learner.genetic.fitness.IFitnessFunction;
import main.brain.learner.genetic.mutator.IMutationOperator;
import main.brain.learner.genetic.selector.IPopulationSelector;
import main.brain.move.picker.HeuristicMovePicker;
import main.brain.move.picker.IMovePicker;
import main.tetris.heuristics.IHeuristic;

/**
 * A genetic algorithm implementation that learns the best tetris move
 * given a gene of heuristics
 *
 * @param <IHeuristic>
 */
public class HeuristicGeneticLearner implements ILearner {
    private Population<IHeuristic> population;
    
    public HeuristicGeneticLearner(
            int populationSize, 
            String pathToPopulation, 
            ArrayList<IHeuristic> heuristics,
            ICrossoverOperator<IHeuristic> crossOverOperator, 
            IFitnessFunction<IHeuristic> fitnessFunction, 
            IMutationOperator<IHeuristic> mutationOperator, 
            IPopulationSelector<IHeuristic> populationSelector,
            double PERCENTAGE_TO_CULL,
            double GENE_MUTATION_PROBABILITY
            ){
        this.population = new Population<IHeuristic>(pathToPopulation, heuristics, populationSize, PERCENTAGE_TO_CULL, GENE_MUTATION_PROBABILITY);
        population.setCrossOverOperator(crossOverOperator);
        population.setFitnessFunction(fitnessFunction);
        population.setmutationOperator(mutationOperator);
        population.setpopulationSelector(populationSelector);
    }

    @Override
    public void trainLearner(int iterations) {
        for (int i = 0; i < iterations; i++){
            population.nextGeneration();
            population.saveToDisk();
        }
    }

    /**
     * Simple hack to return a move picker assuming that we are using genes of heuristics
     */
    @Override
    public IMovePicker getMovePicker() {
        Gene<IHeuristic> fittestGene = population.getFittest();
        return new HeuristicMovePicker(fittestGene.getChromosomeWeights(), fittestGene.getChromsomes());
    }


}
