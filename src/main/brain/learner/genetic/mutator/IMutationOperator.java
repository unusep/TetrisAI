package main.brain.learner.genetic.mutator;

import main.brain.learner.genetic.Gene;

public interface IMutationOperator<E> {

    /**
     * Performs a mutation on a gene of type E
     * @param gene to mutate
     */
    void mutate(Gene<E> gene);

}
