package main.player;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Container;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.RenderingHints;
import java.awt.Toolkit;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelListener;
import java.awt.geom.Arc2D;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.DataBuffer;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.Scanner;

import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;

import javafx.util.Pair;
import main.brain.learner.ILearner;
import main.brain.learner.genetic.Gene;
import main.brain.learner.genetic.Population;
import main.brain.learner.genetic.crossover.ICrossoverOperator;
import main.brain.learner.genetic.fitness.IFitnessFunction;
import main.brain.learner.genetic.mutator.IMutationOperator;
import main.brain.learner.genetic.selector.IPopulationSelector;
import main.brain.move.picker.HeuristicMovePicker;
import main.brain.move.picker.IMovePicker;
import main.tetris.engine.State;
import main.tetris.engine.TFrame;
import main.tetris.engine.TLabel;
import main.tetris.engine.TetrisSimulator;
import main.tetris.heuristics.IHeuristic;

public class PlayerSkeleton {
    /**
     * computes the number of edges of the new piece touching the walls 
     */
    public class WallHuggingCoefficientHeuristic implements IHeuristic {

        @Override
        public String toString(){
            return "WallHuggingCoefficientHeuristic";
        }

        @Override
        public double getValue(int[] move, State s) {
            int orient = move[State.ORIENT];
            int slot = move[State.ORIENT];
            int nextPiece = s.getNextPiece();
            int[][][] pBottom = State.getpBottom();
            int[][] pWidth = State.getpWidth();
            int[][][] pTop = State.getpTop();
            double res = 0.0;
            for(int c = 0; c < pWidth[nextPiece][orient];c++) {
                for(int h = pBottom[nextPiece][orient][c]; h < pTop[nextPiece][orient][c]; h++) {
                    if (slot + c == 0 || slot + c == State.COLS -1) res++;
                }
            }
            return res;
        }

    }

    /**
     * computes the sum of heights of each column on the board
     */
    public class SumOfHeightOfColumnsHeuristic implements IHeuristic {

        public double getValue(int[][] board, int[] top) {
            double count = 0.0;
            for (int c = 0; c < top.length; c++){
                count += top[c];
            }
            return count;
        }

        @Override
        public String toString(){
            return "SumOfHeightOfColumnsHeuristic";
        }

        @Override
        public double getValue(int[] move, State s) {
            TetrisSimulator simulator = new TetrisSimulator(s);
            simulator.makeMove(move);
            return getValue(simulator.getField(), simulator.getTop());
        }
    }

    /**
     * computes the sum of heights of each block in the board
     *
     */
    public class SumOfHeightOfBlocksHeuristic implements IHeuristic {

        @Override
        public String toString(){
            return "SumOfHeightOfBlocksHeuristic";
        }


        @Override
        public double getValue(int[] move, State s) {
            TetrisSimulator simulator = new TetrisSimulator(s);
            simulator.makeMove(move);
            double res = 0;
            int[][] board = simulator.getField();
            for (int r = 0; r < State.ROWS; r++){
                for (int c = 0; c < State.COLS; c++){
                    if (board[r][c] != 0) res += r;
                }
            }
            return res;
        }

    }

    /**
     * computes the number of lines a given move clears
     */
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
    /**
     * computes how "bumpy" a board is after placing a piece 
     * bumpiness defined as the difference in heights between adjacent columns
     */
    public class NonLinearBumpinessHeuristic implements IHeuristic {
    
        public double getValue(int[][] board, int[] top) {
            double count = 0;
            for (int col = 1; col < board[0].length; col++){
                count += Math.pow(top[col-1] - top[col], 2);
            }
            return count;
        }
        
        @Override
        public String toString(){
            return "BumpinessHeuristic";
        }
    
        @Override
        public double getValue(int[] move, State s) {
            TetrisSimulator simulator = new TetrisSimulator(s);
            simulator.makeMove(move);
            return getValue(simulator.getField(), simulator.getTop());
        }
    
    
    }

    /**
     * computes the number of holes there will be on the field after placing a move
     */
    public class HolesHeuristic implements IHeuristic {


        public double getValue(int[][] board, int[] top) {
            double count = 0;
            for (int col = 0; col < board[0].length; col++){
                for (int row = 0; row < top[col]; row++){
                    if (board[row][col] == 0){
                        count++;
                    }
                }
            }
            return count;
        }

        @Override
        public String toString(){
            return "HolesHeuristic";
        }


        @Override
        public double getValue(int[] move, State s) {
            TetrisSimulator simulator = new TetrisSimulator(s);
            simulator.makeMove(move);
            return getValue(simulator.getField(), simulator.getTop());
        }
    }
    /**
     * computes how many edges of the new piece touches the floor
     */
    public class FloorHuggingCoefficientHeuristic implements IHeuristic {

        @Override
        public String toString(){
            return "FloorHuggingCoefficientHeuristic";
        }

        @Override
        public double getValue(int[] move, State s) {
            int orient = move[State.ORIENT];
            int slot = move[State.ORIENT];
            int nextPiece = s.getNextPiece();
            int[][][] pBottom = State.getpBottom();
            int[][] pWidth = State.getpWidth();
            int[][][] pTop = State.getpTop();
            int[] top = s.getTop();
            
            double res = 0.0;
            int height = top[slot]-pBottom[nextPiece][orient][0];
            //for each column beyond the first in the piece
            for(int c = 1; c < pWidth[nextPiece][orient];c++) {
                height = Math.max(height,top[slot+c]-pBottom[nextPiece][orient][c]);
            }
            
            //for each column in the piece - fill in the appropriate blocks
            for(int i = 0; i < pWidth[nextPiece][orient]; i++) {
                //from bottom to top of brick
                for(int h = height+pBottom[nextPiece][orient][i]; h < height+pTop[nextPiece][orient][i]; h++) {
                    if (h == 0) res++;
                }
            }
            return res;
        }

    }
    /**
     * computes how many edges of the new piece will touch an existing piece
     *
     */
    public class FlatteningCoefficientHeuristic implements IHeuristic {

        @Override
        public String toString(){
            return "FlatteningCoefficientHeuristic";
        }

        public double getValue(int[][] board, int[] top, int[][][] pTop, int[][][] pBottom, int[][] pWidth, int pieceIndex,
                int rotationIndex, int leftPosition) {
            ArrayList<Integer> filledCoord = new ArrayList<Integer>();
            double flatteningCoeff = 0;
            int[] rows = new int[board.length + 4];
            
            int height = top[leftPosition] - pBottom[pieceIndex][rotationIndex][0];
            
            //for each column beyond the first in the piece
            for (int c = 0; c < pWidth[pieceIndex][rotationIndex]; c++) {
                height = Math.max(height, top[leftPosition + c] - pBottom[pieceIndex][rotationIndex][c]);
            }
            for (int i = 0; i < pWidth[pieceIndex][rotationIndex]; i++) {

                //from bottom to top of brick
                for (int h = height + pBottom[pieceIndex][rotationIndex][i] + 1;
                     h <= height + pTop[pieceIndex][rotationIndex][i]; h++) {

                    // Tile is being filled here.
                    //board[h][i + leftPosition] = true;
                    rows[h] |= (1 << (i + leftPosition));

                    filledCoord.add(h * board[0].length + i + leftPosition);
                }
            }
            for (int coord : filledCoord) {
                int row = coord / board[0].length;
                int col = coord % board[0].length;

                int left = row * board[0].length + col - 1;
                int right = row * board[0].length + col + 1;
                int down = (row - 1) * board[0].length + col;

                // left side
                if (col != 0 && !filledCoord.contains(left) && ((rows[row] & (1 << (col - 1))) > 0)) {
                    flatteningCoeff ++;
                }

                // right side
                if (col != board[0].length - 1 && !filledCoord.contains(right) && ((rows[row] & (1 << (col + 1))) > 0)) {
                    flatteningCoeff ++;
                }

                // down side
                if (row != 0 && !filledCoord.contains(down) && ((rows[row - 1] & (1 << col)) > 0)) {
                    flatteningCoeff ++;
                }
            }
            return flatteningCoeff;
        }

        @Override
        public double getValue(int[] move, State s) {
            return getValue(s.getField(), s.getTop(), State.getpTop(), State.getpBottom(), State.getpWidth(), s.getNextPiece(), move[State.ORIENT], move[State.SLOT]);
        }
    }


    /**
     * computes how many blocks there are above a hole
     *
     */
    public class BlockadesHeuristic implements IHeuristic {
    
        public double getValue(int[][] board, int[] top){
            double res = 0;
            for (int c = 0; c < board[0].length; c++){
                boolean holeFound = false;
                for (int r = 0; r <= top[c]; r++){
                    if (holeFound) {
                        res += board[r][c] == 0 ? 0 : 1;
                    } else {
                        if (board[r][c] == 0) holeFound = true;
                    }
                }
            }
            return res;
        }
    
        @Override
        public String toString(){
            return "BlockadesHeuristic";
        }
    
        @Override
        public double getValue(int[] move, State s) {
            TetrisSimulator simulator = new TetrisSimulator(s);
            simulator.makeMove(move);
            return getValue(simulator.getField(), simulator.getTop());
        }
    }

    /**
     * IHeuristic interface supports the ability to get a value of a move given a specific game state
     *
     */
    public interface IHeuristic {

        public double getValue(int[] move, State s);
    }
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

    /**
     * HeuristicMovePicker will pick the best move in the given board configuration
     * given a heuristic and its weights
     */
    public class HeuristicMovePicker implements IMovePicker {
        private ArrayList<IHeuristic> heuristics;
        private ArrayList<Double> weights;
        
        public HeuristicMovePicker(ArrayList<Double> weights, ArrayList<IHeuristic> heuristics) {
            this.weights = weights;
            this.heuristics = heuristics;
        }

        @Override
        public int[] pickBest(State s) {
            double bestScore = -Double.MAX_VALUE;
            int[][] legalMoves = s.legalMoves();
            int[] bestMove = null;
            for (int[] move : legalMoves){
                double score = evaluateMove(move, s);
                if (score > bestScore){
                    bestScore = score;
                    bestMove = move;
                }
            }
            return bestMove;
        }

        /**
         * evaluates the value of using a particular move in a particular state
         * using the heuristics
         * @param move
         * @param state s
         * @return score of move
         */
        private double evaluateMove(int[] move, State s) {
            double score = 0.0;
            for (int i = 0; i < Math.min(weights.size(),  heuristics.size()); i++){
                score += weights.get(i) * heuristics.get(i).getValue(move, s);
            }
            return score;
        }

    }

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
    public class TruncationFitnessSelector<E> implements IPopulationSelector<E> {

        @Override
        public ArrayList<Gene<E>> selectElite(ArrayList<Gene<E>> genePool, int num) {
            Collections.sort(genePool, Collections.reverseOrder());
            ArrayList<Gene<E>> result = new ArrayList<Gene<E>>();
            for (int i = 0; i < num; i++){
                Gene<E> gene = genePool.get(i);
                if (gene.getFitness() <= 0) break;
                result.add(genePool.get(i));
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

    public interface IMutationOperator<E> {

        /**
         * Performs a mutation on a gene of type E
         * @param gene to mutate
         */
        void mutate(Gene<E> gene);

    }
    public class UniformMutation<E> implements IMutationOperator<E> {
        private double MUTATION_PROBABILITY = 0.20; 
        
        public UniformMutation(double mutation){
            MUTATION_PROBABILITY = mutation;
        }
        
        @Override
        public void mutate(Gene<E> gene) {
            Random random = new Random();
            ArrayList<Double> weights = gene.getChromosomeWeights();
            for (int i = 0; i < weights.size(); i++){
                if (random.nextDouble() < MUTATION_PROBABILITY){
                    double weight = randomDouble(-1, 1);
                    weights.set(i, weight);
                }
            }
        }
        
    
        private double randomDouble(double min, double max){
            Random randomGenerator = new Random(); 
            return min + (randomGenerator.nextDouble() * (max - min)); 
        }
    
    }
    public interface IFitnessFunction<E> {
        /**
         * Given a gene, evaluate the gene's fitness
         * @param gene to evaluate
         * @return fitness of gene
         */
        public abstract double evaluateFitness(Gene<E> gene);
    }

    public class AverageRowsClearedFitnessFunction implements IFitnessFunction<IHeuristic> {
        private int numPieces;
        private int numGames;
        
        /**
         * @param numPieces the number of pieces to play before stopping the game
         * @param numGames the number of games to play (we will take the average fitness score)
         */
        public AverageRowsClearedFitnessFunction(int numPieces, int numGames) {
            this.numGames = numGames;
            this.numPieces = numPieces;
        }
        
        /**
         * evaluate the fitness of the gene using the movePicker to pick the best move at each step
         */
        @Override
        public double evaluateFitness(Gene<IHeuristic> gene) {
            IMovePicker movePicker = new HeuristicMovePicker(gene.getChromosomeWeights(), gene.getChromsomes());
            double totalFitness = 0.0;
            for (int i = 0; i < numGames; i++){
                State simulator = new State();
                
                for (int j = 0; j < numPieces; j++){
                    if (simulator.hasLost()) break;
                    int[] bestMove = movePicker.pickBest(simulator);
                    simulator.makeMove(bestMove);
                    
                }
                totalFitness += simulator.getRowsCleared();
            }
            return totalFitness/numGames;
        }

    }

    /**
     * A single point crossover implements the crossoveroperator interface 
     * randomly selects a point at which to crossover chromosomes in both parents
     */
    public class SinglePointCrossover implements ICrossoverOperator<IHeuristic>{
        @Override
        public ArrayList<Gene<IHeuristic>> crossover(ArrayList<Gene<IHeuristic>> genes) {
            if (genes.size() < 2) return new ArrayList<Gene<IHeuristic>>(genes);
            ArrayList<Gene<IHeuristic>> babies = new ArrayList<Gene<IHeuristic>>();
            ArrayList<IHeuristic> heuristics = genes.get(0).getChromsomes();
            
            int size = genes.size();
            for (int i = 1; i <= size/2; i++){
                ArrayList<Double> fatherChromosomeWeights = genes.get(i).getChromosomeWeights();
                ArrayList<Double> motherChromosomeWeights = genes.get(size - i).getChromosomeWeights();
                Pair<ArrayList<Double>, ArrayList<Double>> children = cross(fatherChromosomeWeights, motherChromosomeWeights);
                Gene<IHeuristic> baby1 = new Gene<IHeuristic>(heuristics, children.getKey());
                Gene<IHeuristic> baby2 = new Gene<IHeuristic>(heuristics, children.getValue());
                babies.add(baby1);
                babies.add(baby2);
            }
            
            return babies;
        }

        private Pair<ArrayList<Double>, ArrayList<Double>> cross(
                ArrayList<Double> fatherChromosomeWeights, ArrayList<Double> motherChromosomeWeights) {
            int size = fatherChromosomeWeights.size();
            Random random = new Random();
            int crossoverPoint = random.nextInt(size);
            ArrayList<Double> weights1 = new ArrayList<Double>(size);
            ArrayList<Double> weights2 = new ArrayList<Double>(size);
            
            for (int i = 0; i < crossoverPoint; i++){
                double fChromosome = fatherChromosomeWeights.get(i);
                double mChromosome = motherChromosomeWeights.get(i);  
                weights1.add(fChromosome);
                weights2.add(mChromosome);
            }
            for (int i = crossoverPoint; i < size; i++){
                double fChromosome = fatherChromosomeWeights.get(i);
                double mChromosome = motherChromosomeWeights.get(i);  
                weights2.add(fChromosome);
                weights1.add(mChromosome);
            }
            return new Pair<ArrayList<Double>, ArrayList<Double>>(weights1, weights2);
        }

    }
    public interface ICrossoverOperator<E> {
        /**
         * Given a list of genes selected for breeding, perform a crossover and 
         * return the resulting genes in an arraylist 
         * @param genes
         * @return arraylist of baby genes
         */
        public abstract ArrayList<Gene<E>> crossover(ArrayList<Gene<E>> genes); 
    }
    /**
     * Gene stores chromosomes of type E and the weights of each of its chromosome
     */
    public class Gene<E> implements Comparable<Gene<E>> {
        public static final double INITIAL_FITNESS = 0;
        private ArrayList<Double> chromosomesWeights;
        private ArrayList<E> chromosomes;
        private double fitness;
        public boolean evaluated;
        
        public Gene(ArrayList<E> chromosomes, ArrayList<Double> weights) {
            this.chromosomes = chromosomes; 
            this.chromosomesWeights = weights;
            this.fitness = INITIAL_FITNESS;
            this.evaluated = false;
        }
        
        public double getFitness(){
            return this.fitness;
        }
        
        public void setFitness(double fitness) {
            this.fitness = fitness;
        }
        
        @Override
        public boolean equals(Object obj) {
            if (obj == null) {
                return false;
            }
            if (!Gene.class.isAssignableFrom(obj.getClass())) {
                return false;
            }
            final Gene<E> other = (Gene<E>) obj;
            for (int i = 0; i < chromosomesWeights.size(); i++){
                if (this.chromosomesWeights.get(i) != other.chromosomesWeights.get(i)) return false;
            }
            return true;
        }

        @Override
        public int compareTo(Gene<E> o) {
            double res = this.fitness - o.fitness;
            if (res == 0){
                return 0;
            } else if (res < 0){
                return -1;
            } else {
                return 1;
            }
        }
        
        public ArrayList<E> getChromsomes(){
            return chromosomes;
        }

        public ArrayList<Double> getChromosomeWeights() {
            return chromosomesWeights;
        }
        
        /**
         * Returns a space separated string representing the names of each chromosome
         */
        public String getHeaders(){
            String res = "";
            for (E n : chromosomes){
                res += n.toString() + " ";
            }
            res += "Fitness";
            return res;
        }
        
        /**
         * Returns a space separated string representation the weights of each heuristic followed by its fitness
         */
        @Override
        public String toString(){
            DecimalFormat decimalFormat = new DecimalFormat("#.0000");
            String res = "";
            for (Double weight : chromosomesWeights){
                res += new String(decimalFormat.format(weight)) + " ";
            }
            res += decimalFormat.format(fitness);
            return res;
        }
    }
    /**
     * Stores a population of genes with chromosomes of type E
     * Genetic algorithm is executed to generate nextGeneration()
     */
    public class Population<E> {
        private double PERCENTAGE_TO_KILL = 0.2;
        private double MUTATION_PROBABILITY = 0.2;
        
        private ArrayList<Gene<E>> genePool;
        private ICrossoverOperator<E> crossOverOperator;
        private IFitnessFunction<E> fitnessFunction;
        private IMutationOperator<E> mutationOperator;
        private IPopulationSelector<E> populationSelector;
        private String savePath;
        
        public Population(String filepath, ArrayList<E> chromosomes, int populationSize, double percToCull, double mutationProbability){
            this.savePath = filepath;
            this.genePool = instantiateGenePool(filepath, chromosomes, populationSize);
            PERCENTAGE_TO_KILL = percToCull;
            MUTATION_PROBABILITY = mutationProbability;
        }

        /**
         * Write genes to file
         * @param genes
         */
        public void saveToDisk(){
            File file = new File(savePath);
            try {
                BufferedWriter bw = new BufferedWriter(new FileWriter(file));
                
                if (!genePool.isEmpty()){
                    bw.write(genePool.get(0).getHeaders());
                    bw.newLine();
                }
                
                for (Gene<E> gene : genePool){
                    bw.write(gene.toString());
                    bw.newLine();
                }
                
                bw.flush();
                bw.close();
            } catch (IOException e) {
                System.out.println("File failed to save");
            }
            System.out.println("Population saved to disk");
        }
        
        /**
         * Returns the fittest gene among the gene pool 
         * @return fittest gene
         */
        public Gene<E> getFittest() {
            Collections.sort(genePool);
            return genePool.get(genePool.size() - 1);
        }
        
        /**
         * Generates the next generation by running the genetic algorithm:
         * 1) get fitness of genepool
         * 2) choose genes for breeding
         * 3) crossover to generate babies 
         * 4) kill to make room
         * 5) add new babies to genepool
         * 6) mutate population
         */
        public void nextGeneration() {
            evaluateFitness(genePool);
//            normaliseFitness(genePool);
            int numBabies = (int) (PERCENTAGE_TO_KILL * (double) genePool.size()); 
            ArrayList<Gene<E>> parents = selectElite(genePool, numBabies);
            ArrayList<Gene<E>> babies = crossOver(parents);
            cull(genePool, babies.size());
            mutate(babies);
            genePool.addAll(babies);
            printBestGene();
        }
        
        int printCount = 0;
        private void printBestGene() {
            Gene<E> best = getFittest();
            printCount++;
            System.out.println("Turn " + printCount + " Best Fitness: " + best.getFitness() + " Weights: " + best.getChromosomeWeights());
        }

        public void setCrossOverOperator(ICrossoverOperator<E> crossOverOperator){
            this.crossOverOperator = crossOverOperator;
        }

        public void setFitnessFunction(IFitnessFunction<E> fitnessFunction){
            this.fitnessFunction = fitnessFunction;
        }

        public void setmutationOperator(IMutationOperator<E> mutationOperator){
            this.mutationOperator = mutationOperator;
        }
        
        public void setpopulationSelector(IPopulationSelector<E> populationSelector){
            this.populationSelector = populationSelector;
        }
        
        /**
         * Instantiate gene pool of given size using a text file of stored genes
         * @param gene
         * @param size
         * @return
         */
        private ArrayList<Gene<E>> instantiateGenePool(String filepath, ArrayList<E> chromosomes, int populationSize) {
            ArrayList<Gene<E>> genes = new ArrayList<Gene<E>>(populationSize);
            
            try {
                File f = new File(filepath);
                Scanner sc = new Scanner(f);
                sc.nextLine();
                
                // read from file and generate genes
                while (sc.hasNext()) {
                    if (populationSize <= 0) break;
                    
                    String line = sc.nextLine();
                    String[] result = line.split(" ");
                    
                    if (result.length -1 != chromosomes.size()){
                        System.out.println("Warning: chromosomes and weights are not of the same length.");
                        System.out.println("Warning: will continue by generating a new random population");
                        break;
                    }
                    
                    double fitness = Double.parseDouble(result[result.length - 1]);
                    
                    ArrayList<Double> weights = new ArrayList<Double>();
                    for (int i = 0; i < chromosomes.size(); i++){
                        weights.add(Double.parseDouble(result[i]));
                    }

                    //
                    Gene<E> gene = new Gene<E>(chromosomes, weights);
                    gene.setFitness(fitness);
                    
                    populationSize--;
                    genes.add(gene);
                }
                sc.close();
            } catch (FileNotFoundException e) {
                System.out.println("Pathfile not found: " + filepath +". Will try to create one.");
            }
            
            // if file does not have enough gene, generate more automatically
            if (populationSize > 0){
                genes.addAll(instantiateGenePool(chromosomes, populationSize));
            }
            
            return genes;
        }

        /**
         * Generate a randomised gene pool of given size
         * @param gene
         * @param size
         * @return
         */
        public ArrayList<Gene<E>> instantiateGenePool(ArrayList<E> chromosomes, int populationSize){
            int fitness = 0;
            ArrayList<Gene<E>> genes = new ArrayList<Gene<E>>(populationSize); 
            for (int i = 0; i < populationSize; i++){
                ArrayList<Double> weights = new ArrayList<Double>(chromosomes.size());
                for (int j = 0; j < chromosomes.size(); j++){
                    weights.add(randomDouble(-1, 1));
                }
                Gene<E> gene = new Gene<E>(chromosomes, weights);
                gene.setFitness(fitness);
                genes.add(gene);
                System.out.println("Added gene: " + gene);
            }
            return genes;
        }
        
        private double randomDouble(double min, double max){
            Random randomGenerator = new Random(); 
            return min + (randomGenerator.nextDouble() * (max - min)); 
        }
        
        
        private void mutate(ArrayList<Gene<E>> genes) {
            Random ran = new Random();
            for (Gene<E> gene : genes){
                if (ran.nextDouble() < MUTATION_PROBABILITY) {
                    mutationOperator.mutate(gene);
                }
            }
        }

        private ArrayList<Gene<E>> crossOver(ArrayList<Gene<E>> genes) {
            return crossOverOperator.crossover(genes);
        }

        /**
         * Select the genes that should be killed using the population selector
         * @param genePool
         * @param num number to select
         * @return the genes that should be removed
         */
        private void cull(ArrayList<Gene<E>> genePool, int num) {
            populationSelector.cull(genePool, num);
        }
        
        /**
         * Select the genes that are allowed to breed using the population selector
         * @param genePool
         * @param num the maximum number of elites to select
         * @return the genes that are allowed to breed
         */
        private ArrayList<Gene<E>> selectElite(ArrayList<Gene<E>> genes, int num) {
            return populationSelector.selectElite(genes, num);
        }

        /**
         * Compute and update the fitness of each gene in the genepool using the fitness function
         * @param genes
         */
        private void evaluateFitness(ArrayList<Gene<E>> genes) {
            for (Gene<E> gene : genes){
                double fitness = fitnessFunction.evaluateFitness(gene);
                gene.setFitness(fitness);
            }
        }

        /**
         * Normalization means dividing the fitness value of each individual 
         * by the sum of all fitness values, so that the sum of all resulting 
         * fitness values equals 1.
         * @param genePool
         */
        private void normaliseFitness(ArrayList<Gene<E>> genePool) {
            double totalFitness = 0.0;
            for (Gene<E> gene : genePool){
                totalFitness += gene.getFitness();
            }
            for (Gene<E> gene : genePool){
                gene.setFitness(gene.getFitness()/totalFitness);
            }
        }

    }

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

        /**
         * runs the learner through /iterations/ number of generations
         */
        @Override
        public void trainLearner(int iterations) {
            for (int i = 0; i < iterations; i++){
                population.nextGeneration();
                population.saveToDisk();
            }
        }

        /**
         * returns a move picker based on heuristics
         */
        @Override
        public IMovePicker getMovePicker() {
            Gene<IHeuristic> fittestGene = population.getFittest();
            return new HeuristicMovePicker(fittestGene.getChromosomeWeights(), fittestGene.getChromsomes());
        }


    }
    /**
     * A learner interface which supports only two functions. 
     * Our genetic algorithm heuristic learner will support only the ability to learn 
     * and to generate a tetris move picker
     */
    public interface ILearner {
        
        public abstract IMovePicker getMovePicker();
        
        public abstract void trainLearner(int iterations);
    }
    
    public static class TLabel{
        private static final long serialVersionUID = 1L;
        
        public JLabel draw;
        
        // pre-defined colors
        public static final Color BLACK      = Color.BLACK;
        public static final Color BLUE       = Color.BLUE;
        public static final Color CYAN       = Color.CYAN;
        public static final Color DARK_GRAY  = Color.DARK_GRAY;
        public static final Color GRAY       = Color.GRAY;
        public static final Color GREEN      = Color.GREEN;
        public static final Color LIGHT_GRAY = Color.LIGHT_GRAY;
        public static final Color MAGENTA    = Color.MAGENTA;
        public static final Color ORANGE     = Color.ORANGE;
        public static final Color PINK       = Color.PINK;
        public static final Color RED        = Color.RED;
        public static final Color WHITE      = Color.WHITE;
        public static final Color YELLOW     = Color.YELLOW;
        public static final Color NICEGREEN = new Color(0,153,0);
        // default colors
        public static final Color DEFAULT_PEN_COLOR   = BLACK;
        public static final Color DEFAULT_CLEAR_COLOR = WHITE;

        // current pen color
        private static Color penColor;

        // default canvas size is SIZE-by-SIZE
        static final int SIZE = 512;
        public int width  = SIZE;
        private int height = SIZE;

        // default pen radius
        private static final double DEFAULT_PEN_RADIUS = 0.002;

        // current pen radius
        private static double penRadius;

        // boundary of drawing canvas, 0% border
        public double BORDER = 0.00;
        private static final double DEFAULT_XMIN = 0.0;
        private static final double DEFAULT_XMAX = 1.0;
        private static final double DEFAULT_YMIN = 0.0;
        private static final double DEFAULT_YMAX = 1.0;
        public double xmin, ymin, xmax, ymax;

        // default font
        private final Font DEFAULT_FONT = new Font("Serif", Font.PLAIN, 16);

        // current font
        private static Font font;

        // double buffered graphics
        private BufferedImage offscreenImage, onscreenImage;
        protected Graphics2D offscreen, onscreen;
        
        // change the user coordinate system
        public void setXscale() { setXscale(DEFAULT_XMIN, DEFAULT_XMAX); }
        public void setYscale() { setYscale(DEFAULT_YMIN, DEFAULT_YMAX); }
        public void setXscale(double min, double max) {
            double size = max - min;
            xmin = min - BORDER * size;
            xmax = max + BORDER * size;
        }

        public void setYscale(double min, double max) {
            double size = max - min;
            ymin = min - BORDER * size;
            ymax = max + BORDER * size;
        }


        // helper functions that scale from user coordinates to screen coordinates and back
        protected double scaleX (double x) { return width  * (x - xmin) / (xmax - xmin); }
        protected double scaleY (double y) { return height * (ymax - y) / (ymax - ymin); }
        public double factorX(double w) { return w * width  / Math.abs(xmax - xmin);  }
        public double factorY(double h) { return h * height / Math.abs(ymax - ymin);  }
        public double userX  (double x) { return xmin + x * (xmax - xmin) / width;    }
        public double userY  (double y) { return ymax - y * (ymax - ymin) / height;   }


        //create a frame, insert self in frame, then show self
        public void showInFrame() {
            JFrame j = new JFrame();
            j.setTitle("Configuration");
            j.setContentPane(this.draw);
            j.setVisible(true);
            j.pack();
            show();
        }
        
        
        // clear the screen with given color
        public void clear() { clear(DEFAULT_CLEAR_COLOR); }
        public void clear(Color color) {
            offscreen.setColor(color);
            offscreen.fillRect(0, 0, width, height);
            offscreen.setColor(penColor);
        }

        public void clear(double x1, double x2, double y1, double y2)
        { clear(DEFAULT_CLEAR_COLOR, x1, x2, y1, y2);  }
        public void clear(Color color, double x1, double x2, double y1, double y2) {

            int ix1, ix2, iy1, iy2;
            ix1 = (int) scaleX(x1);
            ix2 = (int) scaleX(x2);
            iy1 = (int) scaleY(y1);
            iy2 = (int) scaleY(y2);

            offscreen.setColor(color);
            offscreen.fillRect(ix1, iy1, ix2, iy2);
            offscreen.setColor(penColor);
            //show();
        }

        // set the pen size
        public void setPenRadius() { setPenRadius(DEFAULT_PEN_RADIUS); }
        public void setPenRadius(double r) {
            penRadius = r * SIZE;
            BasicStroke stroke = new BasicStroke((float) penRadius);
            offscreen.setStroke(stroke);
        }



        // set the pen color
        public void setPenColor() { setPenColor(DEFAULT_PEN_COLOR); }
        public void setPenColor(Color color) {
            penColor = color;
            offscreen.setColor(penColor);
        }

        // write the given string in the current font
        public void setFont() { setFont(DEFAULT_FONT); }
        public void setFont(Font f) { 
            Toolkit toolkit = java . awt . Toolkit . getDefaultToolkit ();
            double x = toolkit.getScreenSize().getWidth();
            double y = toolkit.getScreenSize().getHeight();
            double xscale = x/1400.0;
            double yscale = y/1050.0;
            double scale = Math.sqrt((xscale*xscale+yscale*yscale)/2);

            font = f.deriveFont((float) (f.getSize()*scale));
        }
        
        public TLabel(int w, int h){
            width = w;
            height = h;
            offscreenImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
            onscreenImage  = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
            offscreen = offscreenImage.createGraphics();
            onscreen  = onscreenImage.createGraphics();
            setXscale();
            setYscale();
            offscreen.setColor(DEFAULT_CLEAR_COLOR);
            offscreen.fillRect(0, 0, width, height);
            setPenColor();
            setPenRadius();
            setFont();
            clear();

            // add antialiasing
            RenderingHints hints = new RenderingHints(RenderingHints.KEY_ANTIALIASING,
                    RenderingHints.VALUE_ANTIALIAS_ON);
            hints.put(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
            offscreen.addRenderingHints(hints);

            // frame stuff
            ImageIcon icon = new ImageIcon(onscreenImage);
            draw = new JLabel(icon);
        }
        
        public void add(Container frame, String spot){
            frame.add(draw, spot);
        }
        

        
        public void addML(MouseListener frame){
            draw.addMouseListener(frame);
        }
        
        public void addMML(MouseMotionListener frame){
            draw.addMouseMotionListener(frame);
        }
        
        public void addKL(KeyListener frame){
            draw.addKeyListener(frame);
        }
        
        public void addMWL(MouseWheelListener frame) {
            draw.addMouseWheelListener(frame);
        }
        
        public void remML(MouseListener frame){
            draw.removeMouseListener(frame);
        }
        
        public void remMML(MouseMotionListener frame){
            draw.removeMouseMotionListener(frame);
        }
        
        public void remKL(KeyListener frame){
            draw.removeKeyListener(frame);
        }
        
        public void remMWL(MouseWheelListener frame){
            draw.removeMouseWheelListener(frame);
        }
        
        // draw a line from (x0, y0) to (x1, y1)
        public void line(double x0, double y0, double x1, double y1) {
//          System.out.println("drawing a line from " + new Point(x0, y0).toString()+ " to " + new Point(x1,y1).toString());
            offscreen.draw(new Line2D.Double(scaleX(x0), scaleY(y0), scaleX(x1), scaleY(y1)));
        }



        // draw one pixel at (x, y)
        public void pixel(double x, double y) {
            offscreen.fillRect((int) Math.round(scaleX(x)), (int) Math.round(scaleY(y)), 1, 1);
        }
        
        // draw one pixel at (x, y)
        public void pixelP(double x, double y) {
            offscreen.fillRect((int) Math.round(x), (int) Math.round(y), 1, 1);
        }
        
        // draw one pixel at (x, y)
        public void pixelP(double x, double y, Color c) {
            setPenColor(c);
            offscreen.fillRect((int) Math.round(x), (int) Math.round(y), 1, 1);
        }

        // draw point at (x, y)
        public void point(double x, double y) {
            double xs = scaleX(x);
            double ys = scaleY(y);
            double r = penRadius;
            // double ws = factorX(2*r);
            // double hs = factorY(2*r);
            // if (ws <= 1 && hs <= 1) pixel(x, y);
            if (r <= 1) pixel(x, y);
            else offscreen.fill(new Ellipse2D.Double(xs - r/2, ys - r/2, r, r));
        }

        public void arc(double x, double y, double r, double startAngle, double arcRange) {
            double xs = scaleX(x);
            double ys = scaleY(y);
            double ws = factorX(2*r);
            double hs = factorY(2*r);
            if (ws <= 1 && hs <= 1) pixel(x, y);
            else offscreen.draw(new Arc2D.Double(xs - ws/2, ys - hs/2, ws, hs,startAngle, arcRange, Arc2D.OPEN));

        }
        
        // draw circle of radius r, centered on (x, y); degenerate to pixel if small
        public void circle(double x, double y, double r) {
            double xs = scaleX(x);
            double ys = scaleY(y);
            double ws = factorX(2*r);
            double hs = factorY(2*r);
            if (ws <= 1 && hs <= 1) pixel(x, y);
            else offscreen.draw(new Ellipse2D.Double(xs - ws/2, ys - hs/2, ws, hs));
        }
        
        public void circleP(double x, double y, double r, Color col) {
            setPenColor(col);
            circleP(x,y,r);
        }
        
        public void circleP(double x, double y, double r) {
            double ws = 2*r;
            double hs = 2*r;
            double xs = scaleX(x);
            double ys = scaleY(y);
            if (ws <= 1 && hs <= 1) pixel(x, y);
            else offscreen.draw(new Ellipse2D.Double(xs - ws/2, ys - hs/2, ws, hs));
        }


        
        public void circle(double x, double y, double r, Color color) {
            setPenColor(color);
            circle(x,y,r);
        }
        
        
        // draw filled circle of radius r, centered on (x, y); degenerate to pixel if small
        public void filledCircle(double x, double y, double r) {
            double xs = scaleX(x);
            double ys = scaleY(y);
            double ws = factorX(2*r);
            double hs = factorY(2*r);
            if (ws <= 1 && hs <= 1) pixel(x, y);
            else offscreen.fill(new Ellipse2D.Double(xs - ws/2, ys - hs/2, ws, hs));
        }


        
        public void filledCircleP(double x, double y, double r) {
            double ws = 2*r;
            double hs = 2*r;
            if (ws <= 1 && hs <= 1) pixel(x, y);
            else offscreen.fill(new Ellipse2D.Double(x - ws/2, y - hs/2, ws, hs));
        }

        // draw squared of side length 2r, centered on (x, y); degenerate to pixel if small
        public void square(double x, double y, double r) {
            // screen coordinates
            double xs = scaleX(x);
            double ys = scaleY(y);
            double ws = factorX(2*r);
            double hs = factorY(2*r);
            if (ws <= 1 && hs <= 1) pixel(x, y);
            else offscreen.draw(new Rectangle2D.Double(xs - ws/2, ys - hs/2, ws, hs));
        }

        // draw squared of side length 2r, centered on (x, y); degenerate to pixel if small
        public void filledSquare(double x, double y, double r) {
            // screen coordinates
            double xs = scaleX(x);
            double ys = scaleY(y);
            double ws = factorX(2*r);
            double hs = factorY(2*r);
            if (ws <= 1 && hs <= 1) pixel(x, y);
            else offscreen.fill(new Rectangle2D.Double(xs - ws/2, ys - hs/2, ws, hs));
        }


        //Draw an arrow of the appropriate scale, color, position
        public void Arrow(double x, double y, double w, double h, double Scale, Color Color) {
            rectangle(x, y, w, h, Color);
            double[] xarray = {x+w/2-w/10+Scale*1/10*Math.sqrt(h*h*25+w*w)*Math.sqrt(3)/3,x+w/2-w/10,x+w/2-w/10};
            double[] yarray = {y,y+Scale*1/10*Math.sqrt(h*h*25+w*w)*Math.sqrt(2)/2,y-Scale*1/10*Math.sqrt(h*h*25+w*w)*Math.sqrt(2)/2};
            setPenColor(Color);
            filledPolygon(xarray, yarray);
            setPenColor();
        }

        // draw a polygon with the given (x[i], y[i]) coordinates
        public void polygon(double[] x, double[] y) {
            int N = x.length;
            GeneralPath path = new GeneralPath();
            path.moveTo((float) scaleX(x[0]), (float) scaleY(y[0]));
            for (int i = 0; i < N; i++)
                path.lineTo((float) scaleX(x[i]), (float) scaleY(y[i]));
            path.closePath();
            offscreen.draw(path);
        }



    //  draw a polygon with the given (x[i], y[i]) coordinates
        public void polygonP(double[] x, double[] y) {
            int N = x.length;
            GeneralPath path = new GeneralPath();
            path.moveTo((float) x[0], (float) y[0]);
            for (int i = 0; i < N; i++)
                path.lineTo((float) x[i], (float) y[i]);
            path.closePath();
            offscreen.draw(path);
        }

        // draw a filled polygon with the given (x[i], y[i]) coordinates
        public void filledPolygon(double[] x, double[] y) {
            int N = x.length;
            GeneralPath path = new GeneralPath();
            path.moveTo((float) scaleX(x[0]), (float) scaleY(y[0]));
            for (int i = 0; i < N; i++)
                path.lineTo((float) scaleX(x[i]), (float) scaleY(y[i]));
            path.closePath();
            offscreen.fill(path);
        }

    //  draw a filled polygon with the given (x[i], y[i]) coordinates
        public void filledPolygonP(double[] x, double[] y) {
            int N = x.length;
            GeneralPath path = new GeneralPath();
            path.moveTo((float) x[0], (float) y[0]);
            for (int i = 0; i < N; i++)
                path.lineTo((float) x[i], (float) y[i]);
            path.closePath();
            offscreen.fill(path);
        }

        //Draw rectangle at the given coordinates
        public void rectangle(double x, double y, double w, double h) {
            double[] xarray = {x-w/2,x-w/2,x+w/2,x+w/2};
            double[] yarray = {y-h/2,y+h/2,y+h/2,y-h/2};
            polygon(xarray,yarray);
        }
        
        public void rectangleLL(double x, double y, double w, double h) {
            double[] xarray = {x,x,x+w,x+w};
            double[] yarray = {y,y+h,y+h,y};
            polygon(xarray,yarray);
        }
        
        public void rectangleP(double x, double y, double w, double h) {
            double[] xarray = {x,x,x+w,x+w};
            double[] yarray = {y,y+h,y+h,y};
            polygonP(xarray,yarray);
        }

        public void rectangle(double x, double y, double w, double h, Color c) {
            double[] xarray = {x-w/2,x-w/2,x+w/2,x+w/2};
            double[] yarray = {y-h/2,y+h/2,y+h/2,y-h/2};
            setPenColor(c);
            filledPolygon(xarray,yarray);
            setPenColor(DEFAULT_PEN_COLOR);
        }
        
        public void rectangleC(double x, double y, double w, double h, Color c) {
            double[] xarray = {x,x,x+w,x+w};
            double[] yarray = {y,y+h,y+h,y};
            setPenColor(c);
            filledPolygon(xarray,yarray);
            setPenColor(DEFAULT_PEN_COLOR);
        }

        public void filledRectangleP(double x, double y, double w, double h, Color c) {
            double[] xarray = {x,x,x+w,x+w};
            double[] yarray = {y,y+h,y+h,y};
            setPenColor(c);
            filledPolygonP(xarray,yarray);
            setPenColor(DEFAULT_PEN_COLOR);
        }
        
        public void filledRectangleLL(double x, double y, double w, double h, Color c) {
            double[] xarray = {x,x,x+w,x+w};
            double[] yarray = {y,y+h,y+h,y};
            setPenColor(c);
            filledPolygon(xarray,yarray);
            setPenColor(DEFAULT_PEN_COLOR);
        }

        public void rectangle(double x, double y, double w, double h, Color c,boolean Border, Color BorderColor) {

            double[] xarray = {x-w/2,x-w/2,x+w/2,x+w/2};
            double[] yarray = {y-h/2,y+h/2,y+h/2,y-h/2};
            if (c!=null) {
                setPenColor(c);
                filledPolygon(xarray,yarray);
            }
            setPenColor(BorderColor);
            if (Border) polygon(xarray, yarray);
            setPenColor();
        }
        

        
        // draw picture (gif, jpg, or png) upperLeft on (x, y), rescaled to w-by-h
        public void image(double x, double y, Image image, double w, double h) {

            double xs = scaleX(x);
            double ys = scaleY(y);
            double ws = factorX(w);
            double hs = factorY(h);
            if (ws <= 1 && hs <= 1) pixel(x, y);
            else {
                offscreen.drawImage(image, (int) Math.round(xs),
                        (int) Math.round(ys),
                        (int) Math.round(ws),
                        (int) Math.round(hs), null);
            }
        }
        
        // draw picture (gif, jpg, or png) upperLeft on (x, y), rescaled to w-by-h
        public void imageP(double x, double y, Image image) {

            //if (ws <= 1 && hs <= 1) pixel(x, y);

            offscreen.drawImage(image, (int) Math.round(x),(int) Math.round(y), null);
            
        }
        


        //Invert an image
        public BufferedImage invert(Image image) {
            BufferedImage b1 = new BufferedImage(image.getWidth(null),image.getHeight(null),BufferedImage.TYPE_INT_RGB);
            Graphics bg = b1.getGraphics();
            bg.drawImage(image, 0, 0, null);
            bg.dispose();

            BufferedImage b2 = new BufferedImage(image.getWidth(null),image.getHeight(null),BufferedImage.TYPE_INT_RGB);
            DataBuffer db1 = b1.getRaster().getDataBuffer();
            DataBuffer db2 = b2.getRaster().getDataBuffer();

            for (int i = db1.getSize() - 1, j = 0; i >= 0; --i, j++) {
                db2.setElem(j, db1.getElem(i));

            }
            for (int i = db1.getSize() - 1, j = 0; i >= 0; --i, j++) {
                db1.setElem(i, db1.getElem(i));

            }
            return b2;    
        }




        // write the given text string in the current font, center on (x, y)
        public void text(double x, double y, String s) {
            offscreen.setFont(font);
            FontMetrics metrics = offscreen.getFontMetrics();
            double xs = scaleX(x);
            double ys = scaleY(y);
            int ws = metrics.stringWidth(s);
            int hs = metrics.getDescent();
            offscreen.drawString(s, (float) (xs - ws/2.0), (float) (ys + hs));
        }
        
        public void textTop(double x, double y, String s) {
            offscreen.setFont(font);
            FontMetrics metrics = offscreen.getFontMetrics();
            double xs = scaleX(x);
            double ys = scaleY(y);
            int ws = metrics.stringWidth(s);
            int hs = metrics.getDescent();
            offscreen.drawString(s, (float) (xs - ws/2.0), (float) (ys + hs*3));
        }

        public void textLeft(double x, double y, String s,Color c) {
            setPenColor(c);
            offscreen.setFont(font);
            FontMetrics metrics = offscreen.getFontMetrics();
            double xs = scaleX(x);
            double ys = scaleY(y);
            //int ws = metrics.stringWidth(s);
            int hs = metrics.getDescent();
            offscreen.drawString(s, (float) (xs), (float) (ys + hs));
            setPenColor();
        }


        //Draw text at the appropriate point and color
        public void text(double x, double y, String s,Color c) {
            setPenColor(c);
            offscreen.setFont(font);
            FontMetrics metrics = offscreen.getFontMetrics();
            double xs = scaleX(x);
            double ys = scaleY(y);
            int ws = metrics.stringWidth(s);
            int hs = metrics.getDescent();
            offscreen.drawString(s, (float) (xs - ws/2.0), (float) (ys + hs));
            setPenColor();
        }

    //  write the given text string in the current font, center on (x, y) sized to w, h
        public void text(double x, double y, String s, double w, double h) {
            offscreen.setFont(font);
            //FontMetrics metrics = offscreen.getFontMetrics();
            double xs = scaleX(x);
            double ys = scaleY(y);
            double ws = factorX(w);
            double hs = factorY(h);
            offscreen.drawString(s, (float) (xs - ws/2.0), (float) (ys + hs));
        }

        public void absText(String s, int x, int y) {
            offscreen.drawString(s, (float) (x), (float) (y));
        }
        
        public void textinvert(double x, double y, String s) {
            offscreen.setFont(font);
            FontMetrics metrics = offscreen.getFontMetrics();
            double xs = scaleX(x);
            double ys = scaleY(y);
            int ws = metrics.stringWidth(s);
            int hs = metrics.getDescent();
            BufferedImage bimage = new BufferedImage(ws, hs,BufferedImage.TYPE_INT_ARGB);
            Graphics2D bimagegraphics = bimage.createGraphics();
            bimagegraphics.drawString(s, (float) (xs - ws/2.0), (float) (ys + hs));
            BufferedImage bimage2 = invert(bimage);
            offscreen.drawImage(bimage2, (int) Math.round(xs - ws/2.0),
                    (int) Math.round(ys + hs), null);
        }
        
        

        
        // view on-screen, creating new frame if necessary
        public void show() {
            onscreen.drawImage(offscreenImage, 0, 0, null);
            try{
                draw.repaint();
                //frame.paint(frame.getGraphics());
            }
            catch(NullPointerException e){
                System.out.println("Null Pointer Exception in showatonce");
            }

        }

        


        

        
        

        
    }

    public static class TFrame extends JFrame implements KeyListener{
        private static final long serialVersionUID = 1L;
        public TLabel label = new TLabel(300,700);
        public State s;
        
        public int orient, slot;
        
        public static final int MANUAL = 0;
        public static final int NONE = 1;
        
        public int mode = MANUAL;
        
        //constructor
        public TFrame (State s){
            this.s = s;
            s.label = label;
            setResizable(false);
            setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);            // closes all windows when this is closed
            setTitle("Tetris BKW");
            setContentPane(label.draw);
            pack();
            label.BORDER = .05;
            label.setXscale(0, State.COLS);
            label.setYscale(0, State.ROWS+5);
            this.addKeyListener(this);  //may be unnecessary (not certain)
            setVisible(true);
        }
        
        //switches which state is attached to this TFrame
        public void bindState(State s) {
            if(s!= null)    s.label = null;
            this.s = s;
            s.label = label;
        }
        
        ///
        /// ADDED BY DON (AKA Pimp Masta) 1/22/09
        ///
        public TFrame (){
            s.label = label;
            setResizable(false);
            setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);            // closes all windows when this is closed
            setTitle("Eric Whitman's Tetris Simulator");
            setContentPane(label.draw);
            pack();
            label.BORDER = .05;
            label.setXscale(0, State.COLS);
            label.setYscale(0, State.ROWS+5);
            this.addKeyListener(this);  //may be unnecessary (not certain)
            setVisible(true);
        }

        public void keyPressed(KeyEvent e) {
            switch(mode) {
                case(MANUAL): {
                    switch(e.getKeyCode()) {
                        case(KeyEvent.VK_RIGHT): {
                            if(slot < State.COLS-State.pWidth[s.nextPiece][orient]) slot++;
                            s.clearNext();
                            s.drawNext(slot, orient);
                            break;
                        }
                        case(KeyEvent.VK_LEFT): {
                            if(slot > 0)    slot--;
                            s.clearNext();
                            s.drawNext(slot, orient);
                            break;
                        }
                        case(KeyEvent.VK_UP): {
                            orient++;
                            if(orient%State.pOrients[s.nextPiece]==0)   orient = 0;
                            if(slot > State.COLS-State.pWidth[s.nextPiece][orient])
                                slot = State.COLS-State.pWidth[s.nextPiece][orient];
                            s.clearNext();
                            s.drawNext(slot, orient);
                            break;
                        }
                        case(KeyEvent.VK_DOWN): {
                            if(!s.makeMove(orient, slot))   mode = NONE;
                            if(orient >= State.pOrients[s.nextPiece])   orient = 0;
                            if(slot > State.COLS-State.pWidth[s.nextPiece][orient])
                                slot = State.COLS-State.pWidth[s.nextPiece][orient];
                            
                            s.draw();
                            if(mode == NONE)    {
                                label.text(State.COLS/2.0, State.ROWS/2.0, "You Lose");
                            }
                            s.clearNext();
                            s.drawNext(slot, orient);
                            break;
                        }
                        default:
                            break;
                    }
                }
                case(NONE): break;
                default:
                    System.out.println("unknown mode");
                    break;
            }
            
            
            
            
        }


        public void keyReleased(KeyEvent e) {
        }


        public void keyTyped(KeyEvent e) {

        }
        
        public void save(String filename) {
            File file = new File(filename);
            String suffix = filename.substring(filename.lastIndexOf('.') + 1);


            BufferedImage bImage = new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_INT_ARGB);
            Graphics2D graphic = (Graphics2D)bImage.getGraphics();
            paint(graphic);
            graphic.drawImage(bImage, 0, 0, null);
//             png files
            if (suffix.toLowerCase().equals("png")) {
                try { ImageIO.write(bImage, suffix, file); }
                catch (IOException e) { e.printStackTrace(); }
            }
            else System.out.println("unknown extension");
        }
        
        public static void main(String[] args) {
            State s = new State();
            TFrame t = new TFrame(s);
            s.draw();
            s.drawNext(0,0);
            //t.save("picture.png");
            
        }
    }

    /**
     * A simulator that allows us to make moves without affecting actual game state
     * mimics the actual game State without side effects
     */
    public static class TetrisSimulator {

        public static int COLS = 10;
        public static int ROWS = 21;
        public static int N_PIECES = 7;



        public boolean lost = false;




        public TLabel label;

        //current turn
        private int turn = 0;
        private int cleared = 0;

        //each square in the grid - int means empty - other values mean the turn it was placed
        private int[][] field = new int[ROWS][COLS];
        //top row+1 of each column
        //0 means empty
        private int[] top = new int[COLS];


        //number of next piece
        protected int nextPiece;



        //all legal moves - first index is piece type - then a list of 2-length arrays
        protected static int[][][] legalMoves = new int[N_PIECES][][];

        //indices for legalMoves
        public static final int ORIENT = 0;
        public static final int SLOT = 1;

        //possible orientations for a given piece type
        protected static int[] pOrients = {1,2,4,4,4,2,2};

        //the next several arrays define the piece vocabulary in detail
        //width of the pieces [piece ID][orientation]
        protected static int[][] pWidth = {
                {2},
                {1,4},
                {2,3,2,3},
                {2,3,2,3},
                {2,3,2,3},
                {3,2},
                {3,2}
        };
        //height of the pieces [piece ID][orientation]
        private static int[][] pHeight = {
                {2},
                {4,1},
                {3,2,3,2},
                {3,2,3,2},
                {3,2,3,2},
                {2,3},
                {2,3}
        };
        private static int[][][] pBottom = {
                {{0,0}},
                {{0},{0,0,0,0}},
                {{0,0},{0,1,1},{2,0},{0,0,0}},
                {{0,0},{0,0,0},{0,2},{1,1,0}},
                {{0,1},{1,0,1},{1,0},{0,0,0}},
                {{0,0,1},{1,0}},
                {{1,0,0},{0,1}}
        };
        private static int[][][] pTop = {
                {{2,2}},
                {{4},{1,1,1,1}},
                {{3,1},{2,2,2},{3,3},{1,1,2}},
                {{1,3},{2,1,1},{3,3},{2,2,2}},
                {{3,2},{2,2,2},{2,3},{1,2,1}},
                {{1,2,2},{3,2}},
                {{2,2,1},{2,3}}
        };


        //constructor
        public TetrisSimulator(State s) {
            this.turn = s.getTurnNumber();
            this.nextPiece = s.getNextPiece();
            this.cleared = s.getRowsCleared();
            this.lost = s.hasLost();
            COLS = State.COLS;
            ROWS = State.ROWS;
            N_PIECES = State.N_PIECES;
            this.field = new int[ROWS][COLS];
            int[][] sField = s.getField();
            for (int r = 0; r < ROWS; r++){
                for (int c = 0; c < COLS; c++){
                    this.field[r][c] = sField[r][c];
                }
            }
            this.top = new int[COLS];
            for (int i = 0; i < COLS; i++){
                this.top[i] = s.getTop()[i];
            }
            pOrients = State.getpOrients();
            pWidth = State.getpWidth();
            pHeight = State.getpHeight();
            pBottom = State.getpBottom();
            pTop = State.getpTop();
            for(int i = 0; i < N_PIECES; i++) {
                //figure number of legal moves
                int n = 0;
                for(int j = 0; j < pOrients[i]; j++) {
                    //number of locations in this orientation
                    n += COLS+1-pWidth[i][j];
                }
                //allocate space
                legalMoves[i] = new int[n][2];
                //for each orientation
                n = 0;
                for(int j = 0; j < pOrients[i]; j++) {
                    //for each slot
                    for(int k = 0; k < COLS+1-pWidth[i][j];k++) {
                        legalMoves[i][n][ORIENT] = j;
                        legalMoves[i][n][SLOT] = k;
                        n++;
                    }
                }
            }
        }
        //constructor
        public TetrisSimulator(TetrisSimulator s) {
            this.turn = s.getTurnNumber();
            this.nextPiece = s.getNextPiece();
            this.cleared = s.getRowsCleared();
            this.lost = s.hasLost();
            COLS = State.COLS;
            ROWS = State.ROWS;
            N_PIECES = State.N_PIECES;
            this.field = new int[ROWS][COLS];
            for (int r = 0; r < ROWS; r++){
                for (int c = 0; c < COLS; c++){
                    this.field[r][c] = s.getField()[r][c];
                }
            }
            this.top = new int[COLS];
            for (int i = 0; i < COLS; i++){
                this.top[i] = s.getTop()[i];
            }
            pOrients = State.getpOrients();
            pWidth = State.getpWidth();
            pHeight = State.getpHeight();
            pBottom = State.getpBottom();
            pTop = State.getpTop();
            for(int i = 0; i < N_PIECES; i++) {
                //figure number of legal moves
                int n = 0;
                for(int j = 0; j < pOrients[i]; j++) {
                    //number of locations in this orientation
                    n += COLS+1-pWidth[i][j];
                }
                //allocate space
                legalMoves[i] = new int[n][2];
                //for each orientation
                n = 0;
                for(int j = 0; j < pOrients[i]; j++) {
                    //for each slot
                    for(int k = 0; k < COLS+1-pWidth[i][j];k++) {
                        legalMoves[i][n][ORIENT] = j;
                        legalMoves[i][n][SLOT] = k;
                        n++;
                    }
                }
            }
        }
        
        //initialize legalMoves
        {
            //for each piece type
            for(int i = 0; i < N_PIECES; i++) {
                //figure number of legal moves
                int n = 0;
                for(int j = 0; j < pOrients[i]; j++) {
                    //number of locations in this orientation
                    n += COLS+1-pWidth[i][j];
                }
                //allocate space
                legalMoves[i] = new int[n][2];
                //for each orientation
                n = 0;
                for(int j = 0; j < pOrients[i]; j++) {
                    //for each slot
                    for(int k = 0; k < COLS+1-pWidth[i][j];k++) {
                        legalMoves[i][n][ORIENT] = j;
                        legalMoves[i][n][SLOT] = k;
                        n++;
                    }
                }
            }

        }


        public int[][] getField() {
            return field;
        }

        public int[] getTop() {
            return top;
        }

        public static int[] getpOrients() {
            return pOrients;
        }

        public static int[][] getpWidth() {
            return pWidth;
        }

        public static int[][] getpHeight() {
            return pHeight;
        }

        public static int[][][] getpBottom() {
            return pBottom;
        }

        public static int[][][] getpTop() {
            return pTop;
        }


        public int getNextPiece() {
            return nextPiece;
        }

        public boolean hasLost() {
            return lost;
        }

        public int getRowsCleared() {
            return cleared;
        }

        public int getTurnNumber() {
            return turn;
        }




        //random integer, returns 0-6
        private int randomPiece() {
            return (int)(Math.random()*N_PIECES);
        }




        //gives legal moves for 
        public int[][] legalMoves() {
            return legalMoves[nextPiece];
        }

        //make a move based on the move index - its order in the legalMoves list
        public void makeMove(int move) {
            makeMove(legalMoves[nextPiece][move]);
        }

        //make a move based on an array of orient and slot
        public void makeMove(int[] move) {
            makeMove(move[ORIENT],move[SLOT]);
        }

        //returns false if you lose - true otherwise
        public boolean makeMove(int orient, int slot) {
            turn++;
            //height if the first column makes contact
            int height = top[slot]-pBottom[nextPiece][orient][0];
            //for each column beyond the first in the piece
            for(int c = 1; c < pWidth[nextPiece][orient];c++) {
                height = Math.max(height,top[slot+c]-pBottom[nextPiece][orient][c]);
            }

            //check if game ended
            if(height+pHeight[nextPiece][orient] >= ROWS) {
                lost = true;
                return false;
            }


            //for each column in the piece - fill in the appropriate blocks
            for(int i = 0; i < pWidth[nextPiece][orient]; i++) {

                //from bottom to top of brick
                for(int h = height+pBottom[nextPiece][orient][i]; h < height+pTop[nextPiece][orient][i]; h++) {
                    field[h][i+slot] = turn;
                }
            }

            //adjust top
            for(int c = 0; c < pWidth[nextPiece][orient]; c++) {
                top[slot+c]=height+pTop[nextPiece][orient][c];
            }

            int rowsCleared = 0;

            //check for full rows - starting at the top
            for(int r = height+pHeight[nextPiece][orient]-1; r >= height; r--) {
                //check all columns in the row
                boolean full = true;
                for(int c = 0; c < COLS; c++) {
                    if(field[r][c] == 0) {
                        full = false;
                        break;
                    }
                }
                //if the row was full - remove it and slide above stuff down
                if(full) {
                    rowsCleared++;
                    cleared++;
                    //for each column
                    for(int c = 0; c < COLS; c++) {

                        //slide down all bricks
                        for(int i = r; i < top[c]; i++) {
                            field[i][c] = field[i+1][c];
                        }
                        //lower the top
                        top[c]--;
                        while(top[c]>=1 && field[top[c]-1][c]==0)   top[c]--;
                    }
                }
            }


            //pick a new piece
            nextPiece = randomPiece();



            return true;
        }

        public void draw() {
            label.clear();
            label.setPenRadius();
            //outline board
            label.line(0, 0, 0, ROWS+5);
            label.line(COLS, 0, COLS, ROWS+5);
            label.line(0, 0, COLS, 0);
            label.line(0, ROWS-1, COLS, ROWS-1);

            //show bricks

            for(int c = 0; c < COLS; c++) {
                for(int r = 0; r < top[c]; r++) {
                    if(field[r][c] != 0) {
                        drawBrick(c,r);
                    }
                }
            }

            for(int i = 0; i < COLS; i++) {
                label.setPenColor(Color.red);
                label.line(i, top[i], i+1, top[i]);
                label.setPenColor();
            }

            label.show();


        }

        public static final Color brickCol = Color.gray; 

        private void drawBrick(int c, int r) {
            label.filledRectangleLL(c, r, 1, 1, brickCol);
            label.rectangleLL(c, r, 1, 1);
        }

        public void drawNext(int slot, int orient) {
            for(int i = 0; i < pWidth[nextPiece][orient]; i++) {
                for(int j = pBottom[nextPiece][orient][i]; j <pTop[nextPiece][orient][i]; j++) {
                    drawBrick(i+slot, j+ROWS+1);
                }
            }
            label.show();
        }

        //visualization
        //clears the area where the next piece is shown (top)
        public void clearNext() {
            label.filledRectangleLL(0, ROWS+.9, COLS, 4.2, TLabel.DEFAULT_CLEAR_COLOR);
            label.line(0, 0, 0, ROWS+5);
            label.line(COLS, 0, COLS, ROWS+5);
        }


    }

    public static class State {
        public static final int COLS = 10;
        public static final int ROWS = 21;
        public static final int N_PIECES = 7;

        

        public boolean lost = false;
        
        
        

        
        public TLabel label;
        
        //current turn
        private int turn = 0;
        private int cleared = 0;
        
        //each square in the grid - int means empty - other values mean the turn it was placed
        private int[][] field = new int[ROWS][COLS];
        //top row+1 of each column
        //0 means empty
        private int[] top = new int[COLS];
        
        
        //number of next piece
        protected int nextPiece;
        
        
        
        //all legal moves - first index is piece type - then a list of 2-length arrays
        protected static int[][][] legalMoves = new int[N_PIECES][][];
        
        //indices for legalMoves
        public static final int ORIENT = 0;
        public static final int SLOT = 1;
        
        //possible orientations for a given piece type
        protected static int[] pOrients = {1,2,4,4,4,2,2};
        
        //the next several arrays define the piece vocabulary in detail
        //width of the pieces [piece ID][orientation]
        protected static int[][] pWidth = {
                {2},
                {1,4},
                {2,3,2,3},
                {2,3,2,3},
                {2,3,2,3},
                {3,2},
                {3,2}
        };
        //height of the pieces [piece ID][orientation]
        private static int[][] pHeight = {
                {2},
                {4,1},
                {3,2,3,2},
                {3,2,3,2},
                {3,2,3,2},
                {2,3},
                {2,3}
        };
        private static int[][][] pBottom = {
            {{0,0}},
            {{0},{0,0,0,0}},
            {{0,0},{0,1,1},{2,0},{0,0,0}},
            {{0,0},{0,0,0},{0,2},{1,1,0}},
            {{0,1},{1,0,1},{1,0},{0,0,0}},
            {{0,0,1},{1,0}},
            {{1,0,0},{0,1}}
        };
        private static int[][][] pTop = {
            {{2,2}},
            {{4},{1,1,1,1}},
            {{3,1},{2,2,2},{3,3},{1,1,2}},
            {{1,3},{2,1,1},{3,3},{2,2,2}},
            {{3,2},{2,2,2},{2,3},{1,2,1}},
            {{1,2,2},{3,2}},
            {{2,2,1},{2,3}}
        };
        
        //initialize legalMoves
        {
            //for each piece type
            for(int i = 0; i < N_PIECES; i++) {
                //figure number of legal moves
                int n = 0;
                for(int j = 0; j < pOrients[i]; j++) {
                    //number of locations in this orientation
                    n += COLS+1-pWidth[i][j];
                }
                //allocate space
                legalMoves[i] = new int[n][2];
                //for each orientation
                n = 0;
                for(int j = 0; j < pOrients[i]; j++) {
                    //for each slot
                    for(int k = 0; k < COLS+1-pWidth[i][j];k++) {
                        legalMoves[i][n][ORIENT] = j;
                        legalMoves[i][n][SLOT] = k;
                        n++;
                    }
                }
            }
        
        }
        
        
        public int[][] getField() {
            return field;
        }

        public int[] getTop() {
            return top;
        }

        public static int[] getpOrients() {
            return pOrients;
        }
        
        public static int[][] getpWidth() {
            return pWidth;
        }

        public static int[][] getpHeight() {
            return pHeight;
        }

        public static int[][][] getpBottom() {
            return pBottom;
        }

        public static int[][][] getpTop() {
            return pTop;
        }


        public int getNextPiece() {
            return nextPiece;
        }
        
        public boolean hasLost() {
            return lost;
        }
        
        public int getRowsCleared() {
            return cleared;
        }
        
        public int getTurnNumber() {
            return turn;
        }
        
        
        
        //constructor
        public State() {
            nextPiece = randomPiece();

        }
        
        //random integer, returns 0-6
        private int randomPiece() {
            return (int)(Math.random()*N_PIECES);
        }
        


        
        //gives legal moves for 
        public int[][] legalMoves() {
            return legalMoves[nextPiece];
        }
        
        //make a move based on the move index - its order in the legalMoves list
        public void makeMove(int move) {
            makeMove(legalMoves[nextPiece][move]);
        }
        
        //make a move based on an array of orient and slot
        public void makeMove(int[] move) {
            makeMove(move[ORIENT],move[SLOT]);
        }
        
        //returns false if you lose - true otherwise
        public boolean makeMove(int orient, int slot) {
            turn++;
            //height if the first column makes contact
            int height = top[slot]-pBottom[nextPiece][orient][0];
            //for each column beyond the first in the piece
            for(int c = 1; c < pWidth[nextPiece][orient];c++) {
                height = Math.max(height,top[slot+c]-pBottom[nextPiece][orient][c]);
            }
            
            //check if game ended
            if(height+pHeight[nextPiece][orient] >= ROWS) {
                lost = true;
                return false;
            }

            
            //for each column in the piece - fill in the appropriate blocks
            for(int i = 0; i < pWidth[nextPiece][orient]; i++) {
                
                //from bottom to top of brick
                for(int h = height+pBottom[nextPiece][orient][i]; h < height+pTop[nextPiece][orient][i]; h++) {
                    field[h][i+slot] = turn;
                }
            }
            
            //adjust top
            for(int c = 0; c < pWidth[nextPiece][orient]; c++) {
                top[slot+c]=height+pTop[nextPiece][orient][c];
            }
            
            int rowsCleared = 0;
            
            //check for full rows - starting at the top
            for(int r = height+pHeight[nextPiece][orient]-1; r >= height; r--) {
                //check all columns in the row
                boolean full = true;
                for(int c = 0; c < COLS; c++) {
                    if(field[r][c] == 0) {
                        full = false;
                        break;
                    }
                }
                //if the row was full - remove it and slide above stuff down
                if(full) {
                    rowsCleared++;
                    cleared++;
                    //for each column
                    for(int c = 0; c < COLS; c++) {

                        //slide down all bricks
                        for(int i = r; i < top[c]; i++) {
                            field[i][c] = field[i+1][c];
                        }
                        //lower the top
                        top[c]--;
                        while(top[c]>=1 && field[top[c]-1][c]==0)   top[c]--;
                    }
                }
            }
        

            //pick a new piece
            nextPiece = randomPiece();
            

            
            return true;
        }
        
        public void draw() {
            label.clear();
            label.setPenRadius();
            //outline board
            label.line(0, 0, 0, ROWS+5);
            label.line(COLS, 0, COLS, ROWS+5);
            label.line(0, 0, COLS, 0);
            label.line(0, ROWS-1, COLS, ROWS-1);
            
            //show bricks
                    
            for(int c = 0; c < COLS; c++) {
                for(int r = 0; r < top[c]; r++) {
                    if(field[r][c] != 0) {
                        drawBrick(c,r);
                    }
                }
            }
            
            for(int i = 0; i < COLS; i++) {
                label.setPenColor(Color.red);
                label.line(i, top[i], i+1, top[i]);
                label.setPenColor();
            }
            
            label.show();
            
            
        }
        
        public static final Color brickCol = Color.gray; 
        
        private void drawBrick(int c, int r) {
            label.filledRectangleLL(c, r, 1, 1, brickCol);
            label.rectangleLL(c, r, 1, 1);
        }
        
        public void drawNext(int slot, int orient) {
            for(int i = 0; i < pWidth[nextPiece][orient]; i++) {
                for(int j = pBottom[nextPiece][orient][i]; j <pTop[nextPiece][orient][i]; j++) {
                    drawBrick(i+slot, j+ROWS+1);
                }
            }
            label.show();
        }
        
        //visualization
        //clears the area where the next piece is shown (top)
        public void clearNext() {
            label.filledRectangleLL(0, ROWS+.9, COLS, 4.2, TLabel.DEFAULT_CLEAR_COLOR);
            label.line(0, 0, 0, ROWS+5);
            label.line(COLS, 0, COLS, ROWS+5);
        }

    }
    
    public static final String POPULATION_FILEPATH = "population.txt";
    private ILearner brain;
    
	public int[] pickMove(State simulator, int[][] legalMoves) {
		return brain.getMovePicker().pickBest(simulator);
	}
	
	/**
	 * Pass as arguments 1000 50 to train a population of 1000 people for 50 generations
	 * @param args
	 */
	public static void main(String[] args) {
		State s = new State();
//		new TFrame(s);
		PlayerSkeleton p = new PlayerSkeleton();
		
		// set up genetic learner variables
		int TRAINING_MAX_PIECES = 10000;
		int FITNESS_BEST_OF_NUM_GAMES = 5;
		
		int POPULATION_SIZE = 100;
		int NUMBER_OF_GENERATIONS = 1000;
		double PERCENTAGE_TO_CULL = 0.3; // 0.0-0.49 what percentage of the population to cull
		double CHROMOSOME_MUTATION_PROBABILITY = 0.2; // 0.0-1.0 how likely it is for a particular chromosome in a gene to mutate 
		double GENE_MUTATION_PROBABILITY = 0.3; // how like it is for a particular gene in the population to be selected for mutation
		
		// end of set up genetic learner variables 
		
		ArrayList<IHeuristic> heuristics = new ArrayList<IHeuristic>();
        heuristics.add(p.new NonLinearLinesClearedHeuristic());
        heuristics.add(p.new SumOfHeightOfBlocksHeuristic());
		heuristics.add(p.new SumOfHeightOfColumnsHeuristic());
		heuristics.add(p.new NonLinearBumpinessHeuristic());
		heuristics.add(p.new HolesHeuristic());
        heuristics.add(p.new BlockadesHeuristic());
        heuristics.add(p.new WallHuggingCoefficientHeuristic());
        heuristics.add(p.new FloorHuggingCoefficientHeuristic());
        heuristics.add(p.new FlatteningCoefficientHeuristic());
		ICrossoverOperator<IHeuristic> crossOverOperator = p.new SinglePointCrossover();
		IFitnessFunction<IHeuristic> fitnessFunction 
		    = p.new AverageRowsClearedFitnessFunction(TRAINING_MAX_PIECES, FITNESS_BEST_OF_NUM_GAMES);
		IMutationOperator<IHeuristic> mutationOperator = p.new UniformMutation<IHeuristic>(CHROMOSOME_MUTATION_PROBABILITY);
		IPopulationSelector<IHeuristic> populationSelector = p.new TruncationFitnessSelector<IHeuristic>();
		
        p.brain = p.new HeuristicGeneticLearner(
                POPULATION_SIZE, 
		        POPULATION_FILEPATH, 
		        heuristics, 
		        crossOverOperator, 
		        fitnessFunction, 
		        mutationOperator, 
		        populationSelector,
		        PERCENTAGE_TO_CULL,
		        GENE_MUTATION_PROBABILITY);
        
        p.brain.trainLearner(NUMBER_OF_GENERATIONS);
        
		while(!s.hasLost()) {
			s.makeMove(p.pickMove(s, s.legalMoves()));
//			s.draw();
//			s.drawNext(0,0);
//			try {
//				Thread.sleep(300);
//			} catch (InterruptedException e) {
//				e.printStackTrace();
//			}
		}
		System.out.println("You have completed "+s.getRowsCleared()+" rows.");
	}
	
}
