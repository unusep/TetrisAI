package main.brain.learner.genetic.selector;

import java.util.ArrayList;
import java.util.Collections;

import main.brain.learner.genetic.Gene;

public class PopulationFitnessSelector<E> implements IPopulationSelector<E> {

    @Override
    public ArrayList<Gene<E>> selectElite(ArrayList<Gene<E>> genePool, int num) {
        Collections.sort(genePool);
        int size = genePool.size();
        ArrayList<Gene<E>> result = new ArrayList<Gene<E>>();
        for (int i = 1; i <= num; i++){
            Gene<E> gene = genePool.get(size - i);
            if (gene.getFitness() <= 0) break;
            result.add(genePool.get(size - i));
        }
        return result;
    }

    @Override
    public void cull(ArrayList<Gene<E>> genePool, int num) {
        Collections.sort(genePool, Collections.reverseOrder());
        int size = genePool.size();
        for (int i = 1; i <= num; i++){
            genePool.remove(size - i);
        }
    }
    
}
