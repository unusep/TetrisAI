package main.brain.move.picker;

import main.tetris.engine.TetrisSimulator;

/**
 * An instance of IMovePicker is an object that 
 * tells us the best move to make given a board
 */
public interface IMovePicker {

    /**
     * Given a move, 
     * @param simulator
     * @param legalMoves
     * @return
     */
    public int[] pickBest(TetrisSimulator simulator, int[][] legalMoves);

}
