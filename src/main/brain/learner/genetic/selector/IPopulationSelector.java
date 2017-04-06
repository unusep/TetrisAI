package main.brain.learner.genetic.selector;

import java.util.ArrayList;

import main.brain.learner.genetic.Gene;

public interface IPopulationSelector {
    public ArrayList<Gene> selectFittest(ArrayList<Gene> genes);
}
