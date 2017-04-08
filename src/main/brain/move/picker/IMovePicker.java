package main.brain.move.picker;

import main.tetris.engine.State;

/**
 * An instance of IMovePicker is an object that 
 * tells us the best move to make given a board
 */
public interface IMovePicker {

    /**
     * Given a state return the best move to make
     * @param state
     * @return move
     */
    public int[] pickBest(State s);

}
