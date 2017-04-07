package main.brain.learner.genetic.selector;

import java.util.ArrayList;

import main.brain.learner.genetic.Gene;

public interface IPopulationSelector<E> {

    /**
     * Given a genePool, select the genes to breed
     * @param genePool
     * @return genes selected for breed
     */
    public ArrayList<Gene<E>> selectToBreed(ArrayList<Gene<E>> genePool);
    
    /**
     * Given a genePool, select the genes to kill
     * @param genePool
     * @return genes selected for killing
     */
    public ArrayList<Gene<E>> selectToKill(ArrayList<Gene<E>> genePool);
}
