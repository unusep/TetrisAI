package main.brain.learner.genetic.selector;

import java.util.ArrayList;

import main.brain.learner.genetic.Gene;

public interface IPopulationSelector<E> {

    /**
     * Given a genePool, select the genes to breed
     * @param genePool
     * @param num number of genes to select
     * @return genes selected for breeding
     */
    public ArrayList<Gene<E>> selectElite(ArrayList<Gene<E>> genePool, int num);
    
    /**
     * Given a genePool, select the genes to kill
     * @param genePool
     * @param num number of genes to kill
     * @return genepool after culling
     */
    public void cull(ArrayList<Gene<E>> genePool, int num);
}
