package main.tetris.heuristics;

import main.tetris.engine.TetrisSimulator;

public class HolesHeuristic implements IHeuristic {

    @Override
    public double getValue(TetrisSimulator state) {
        double count = 0;
        
        for(int col = 0; col < TetrisSimulator.COLS; col++){
            int topRow = state.getTop()[col];
            
            //Check the rows from bottom up
            //Increment when an empty space that is below the highest block is found
            for(int row = 0; row < topRow; row++){
                if(state.getField()[row][col] == 0){
                    count++;
                }
            }
        }
        return count;
    }

}
