/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package GeneticAlgorithm;

import GUI.ProgressFrame;
import java.awt.Toolkit;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.concurrent.ExecutionException;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JOptionPane;
import javax.swing.SwingWorker;
import javax.xml.stream.XMLStreamException;
import org.sbml.jsbml.validator.ModelOverdeterminedException;

/**
 *
 * @author chuanfuyap
 */

public class GAlgorithmGUI extends SwingWorker<double[], Void>  {
    private int popsize;          //Population size
    private int maxgens;       //Number of generations
    private int nparam;                 //number of parameters
    private int plateaulimit;
    private int plagueCount;
    
    private int numCore;
    
    private double[] params;            //array to hold the final estimated parameter

    private Population population;
    private Population newgeneration;                    //arraylist of new generation produced by crossing over
    private SystemToSolve bioSystem;
    public FitnessEvaluation fitness;
    private boolean cancellation = false;
    
    private ArrayList<Individual> elite = new ArrayList<>();                            //arraylist of elite individuals, which is eessentially a top % population according to fitness                   
    private ArrayList<Integer> unwantedclones= new ArrayList<>();                       //for some odd reason
    private ArrayList<Individual> clonearmy = new ArrayList<>();                        //creating a population clone to allow for elite arraylist to pull out top 10% of fit individual
                                                                                        //this clone would prevent the reordering of the population
    private ArrayList fitnessScore = new ArrayList<>();
    
    public GAlgorithmGUI(SystemToSolve bioSystem, Boolean SS, int popsize, int maxgen, int plateau, int plague, int numCore, double[][] general_kp_range) {
        this.bioSystem = bioSystem;                                 //generate systemtosolve object from the sbml file and inputdata txt file
        this.popsize=popsize;
        this.maxgens=maxgen;
        this.plateaulimit=plateau;
        this.plagueCount=plague;
        
        this.numCore=numCore;
       
        nparam = bioSystem.getParametersCount();                    //get the number of parameters from systemtosolve object
        population = new Population(nparam, bioSystem, general_kp_range);    //generate new population
        population.initializePopulation(bioSystem, popsize);                          //initialize the new population    
        
        fitness = new FitnessEvaluation(bioSystem);      //object for fitness evaluation 
    }
    
    @Override
    public double[] doInBackground() {
        
        double[] parameters = new double[nparam];
        
        try {
            
            double best_ever = 1e-35;
            Individual best = new Individual(nparam);                               //placeholder for fittest individual
            ProgressFrame.estimationProgressOutput.append("starting size: "+popsize+"\n");
            fitness.InitialEvaluation(population,numCore);                                  //initial evaluation of the entire population
            ProgressFrame.estimationProgressOutput.append("INITIAL EVALUATION DONE\n");
            int miss = fitness.getMissing();                                         //num would be the number of individuals that are deleted as a result of f=1e-100
            System.out.println(miss +"\t bad individuals");
            
            while(population.pop.size()<popsize && !cancellation){                                   //determine if they are individuals that produce f=1e-100
                Population replacement = new Population(nparam);
                replacement.initializePopulation(bioSystem, popsize);
                fitness.InitialEvaluation(replacement,numCore);
                
                population.pop.addAll(replacement.pop);
                System.out.println(replacement.pop.size()+"\t of replacement added\t|new pop count:\t"+population.pop.size());
            }
            
            for(int j =0; j<population.pop.size();j++){                             //determining the fittest individual in all the population
                if(population.pop.get(j).getF()>best_ever){
                    best_ever=population.pop.get(j).getF();
                    best = population.pop.get(j);
                }
            }
            
            population.setBest(best);
            fitnessScore.add(best_ever);
            
            ProgressFrame.estimationProgressOutput.append("BEST IS:\t" + best.getF()+"\n");
            
            for (int i = 0; i < maxgens && !cancellation; i++) {                                     //optimization with GA begins
                int counter = i+1;
                
                removeClones();
            
                if((counter)%plagueCount==0){                                       //plague to kill off the weakest individuals, refreshing the population
                    if(population.pop.size()>population.getInitialCount() ){
                    population.plague();
                    }
                }

                elite = searchforElites((int) (population.pop.size()*0.1));

                newgeneration = spawnNextgen();

                population.pop.addAll(newgeneration.pop);

                removeClones();

                best_ever=1e-50;
                for(int j =0; j<population.pop.size();j++){                         //reconfirm the fittest individual
                    if(population.pop.get(j).getF()>best_ever){
                        best_ever=population.pop.get(j).getF();
                        best = population.pop.get(j);
                    }
                }
                if (best.getF()>population.getBest().getF()){
                    population.setBest(best);
                }
                fitnessScore.add(best_ever);
                
                ProgressFrame.estimationProgressOutput.append("Generation:\t"+counter+"\tBest:\t"+population.getBest().getF()+"\n");
                
                if ((counter) % 3 == 0) {
                    parameters = new double[nparam];
                    System.arraycopy(population.getBest().x, 0, parameters, 0, nparam);
                    if ((best_ever == 1 ) /*|| (best_ever>0.8)*/) {
                        break;                
                    }
                }
                
                if ((counter)>plateaulimit){
                    double difference = population.getBest().getF() - (double)fitnessScore.get(i-plateaulimit);
                    if(difference<1E-8){
                        parameters = new double[nparam];
                        System.arraycopy(population.getBest().x, 0, parameters, 0, nparam);
                        break;
                    }
                }
            }
            
        } catch (ModelOverdeterminedException | InstantiationException | IllegalAccessException | IllegalArgumentException | NoSuchMethodException | InterruptedException | XMLStreamException | IOException ex) {
            Logger.getLogger(GAlgorithmGUI.class.getName()).log(Level.SEVERE, null, ex);
            JOptionPane.showMessageDialog(null, "System ran out of Memory, please restart the process.", "Memory Error", JOptionPane.ERROR_MESSAGE);            
        }
        
       return parameters;
    }
    
    public double[] getParamaters() {        
        return params;
    }
    @Override
    public void done() {
        Toolkit.getDefaultToolkit().beep();
        try {
            params= get();
        } catch (InterruptedException | ExecutionException ex) {
            Logger.getLogger(GAlgorithmGUI.class.getName()).log(Level.SEVERE, null, ex);
        }
        ProgressFrame.estimationProgressOutput.append("Estimation Complete!\n");
        System.out.println("TASK COMPLETE");        
    }
    
    public void cancelGA(boolean cancel){
        cancellation=cancel;
    }
    
    private void removeClones(){
        for(int j=0;j<population.pop.size()-1; j++){                          //loop to remove individuals that have duplicated fitness
            unwantedclones.clear();
            double[] current = population.pop.get(j).getOutput();
            for(int k=(j+1);k<population.pop.size(); k++){                  //idea of this is that it will check everything after as indeces before this would have been scanned/checked
                
                double[] after = population.pop.get(k).getOutput();
                boolean clone = isitSame(current, after);
                if(clone){                        
                    unwantedclones.add(k);
                }
            }
            if(unwantedclones.size()>0){
                Collections.reverse(unwantedclones);

                for(Object index : unwantedclones){
                    int num = (int) index;
                    population.pop.remove(num);
                }
            }
        }
    }
    
    private boolean isitSame(double[] current, double[] after){
        boolean same=false;
        double total=0;
        
        for(int i=0; i<current.length; i++){
            double difference=current[i]-after[i];
            total+=difference*difference;
        }
        
        total = total/current.length;
        
        if(Math.sqrt(total)<1E-10){
            same=true;
        }            
        return same;
    }
    
    private Population spawnNextgen() throws ModelOverdeterminedException, InstantiationException, IllegalAccessException, IllegalArgumentException, NoSuchMethodException, InterruptedException, XMLStreamException, IOException{
        ArrayList<Individual> children = population.nextGeneration();       //produce offspring for the population +10% growth
        ArrayList<Individual> mutants = population.randomMutation();        //mutate 10% of the population 
        ArrayList<Individual> evolve = population.mutate(elite);            //mutate the elites
        Individual royalty = population.royalbirth(elite);                  //1 offspring from the elite
        Individual mutantroyalty = population.clone(royalty);
        population.mutateIndividual(mutantroyalty);                         //mutation on the elite offspring;

        Population nextgen = new Population(nparam);
        nextgen.pop.addAll(children);
        nextgen.pop.addAll(mutants);
        nextgen.pop.addAll(evolve);
        nextgen.pop.add(royalty);
        nextgen.pop.add(mutantroyalty);

        fitness.evaluatePop(nextgen,numCore);
        
        return nextgen;
    }
    
    private ArrayList<Individual> searchforElites(int size){
        ArrayList<Individual> search = new ArrayList();
        clonearmy.clear();                                                  //clear the clonearmy to be ready for next generation
        int elitecount = size;

        for(int k =0;k<population.pop.size();k++){
            Individual clone = population.clone(population.pop.get(k));
            clonearmy.add(clone);
        }

        for(int j =0;j<elitecount;j++){
            double big =1e-100;
            int track=-1;
            for(int k =0;k<clonearmy.size();k++){                           //find the fittest individual in the population
                if(clonearmy.get(k).getF()>big){
                    big = clonearmy.get(k).getF();
                    track=k;
                }
            }
            Individual clone = population.clone(clonearmy.get(track));
            search.add(clone);
            clonearmy.remove(track);                                        //remove previously fittest to rinse repeat and find the second fittest and so on
        }
        
        return search;
    }
    
}
