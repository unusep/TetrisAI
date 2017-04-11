package main.tetris.heuristics;

import main.tetris.engine.State;
import main.tetris.engine.TetrisSimulator;

public class NonLinearLinesClearedHeuristic implements IHeuristic {

    @Override
    public String toString(){
        return "NonLinearLinesClearedHeuristic";
    }

    @Override
    public double getValue(int[] move, State s) {
        TetrisSimulator simulator = new TetrisSimulator(s);
        int before = simulator.getRowsCleared();
        simulator.makeMove(move);
        int after = simulator.getRowsCleared();
        return Math.pow(2, after - before);
    }
}
