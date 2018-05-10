/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package GeneticAlgorithm;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import javax.xml.stream.XMLStreamException;
import org.sbml.jsbml.validator.ModelOverdeterminedException;

/**
 *
 * @author chuanfuyap
 */

public class GAlgorithm {
    private int popsize;          //Population size
    private int maxgens;       //Number of generations
    private int nparam;                 //number of parameters
    private int plateaulimit;
    private int plagueCount;
    
    private ArrayList top30params = new ArrayList();            //array to hold the final estimated parameter

    private Population Main_Population;
    private Population newgeneration;                    //arraylist of new generation produced by crossing over
    private boolean SS;                 //SS to determine data being estimated is of SteadyState or of Time Course data, TRUE IF SS
    private SystemToSolve bioSystem;
    private FitnessEvaluation fitness;
   
    private ArrayList<Individual> elite = new ArrayList<>();                            //arraylist of elite individuals, which is eessentially a top % population according to fitness                   
    private ArrayList<Integer> unwantedclones= new ArrayList<>();                       //for some odd reason
    private ArrayList<Individual> clonearmy = new ArrayList<>();                        //creating a population clone to allow for elite arraylist to pull out top 10% of fit individual
                                                                                        //this clone would prevent the reordering of the population
    private ArrayList fitnessScore = new ArrayList<>();
    
    public GAlgorithm(SystemToSolve bioSystem, Boolean SS, int popsize, int maxgen, int plateau, int plague, double[][] general_kp_range){
        this.bioSystem = bioSystem;                                 //generate systemtosolve object from the sbml file and inputdata txt file
        this.popsize=popsize;
        this.maxgens=maxgen;
        this.plateaulimit=plateau;
        this.plagueCount=plague;
        
        nparam = bioSystem.getParametersCount();                    //get the number of parameters from systemtosolve object
        Main_Population = new Population(nparam, bioSystem, general_kp_range);    //generate new population
        Main_Population.initializePopulation(bioSystem, popsize);                          //initialize the new population    
        
        fitness = new FitnessEvaluation(bioSystem);      //object for fitness evaluation 
    }
        
    public void run() throws ModelOverdeterminedException, InstantiationException, IllegalAccessException, IllegalArgumentException, NoSuchMethodException, InterruptedException, XMLStreamException, IOException{
        double best_ever = 1e-35;
         
        Individual best = new Individual(nparam);                               //placeholder for fittest individual
        System.out.println("starting size: "+popsize);
            
        fitness.InitialEvaluation(Main_Population);                                  //initial evaluation of the entire population
        System.out.println("\033[34m" + "INITIAL EVALUATION DONE; HEALTHY SIZE IS:"+"\u001B[0m"+Main_Population.pop.size());
            
        int miss = fitness.getMissing();                                         //num would be the number of individuals that are deleted as a result of f=1e-100 
        System.out.println("\033[36m" + miss +"\t bad individuals"+"\u001B[0m");
                
        while(Main_Population.pop.size()<popsize){                                   //determine if they are individuals that produce f=1e-100            
            Population replacement = new Population(nparam, bioSystem, Main_Population.getRange());
            replacement.initializePopulation(bioSystem, popsize);            
            fitness.InitialEvaluation(replacement);
            
            Main_Population.pop.addAll(replacement.pop);
            
            System.out.println(replacement.pop.size()+"\t of replacement added\t|new pop count:\t"+Main_Population.pop.size());
        }

        Main_Population.updateInitialCount(popsize);
                
        for(int j =0; j<Main_Population.pop.size();j++){                             //determining the fittest individual in all the population
            if(Main_Population.pop.get(j).getF()>best_ever){
                best_ever=Main_Population.pop.get(j).getF();
                best = Main_Population.pop.get(j);
            }
        }
        
        Main_Population.setBest(best);
        fitnessScore.add(best_ever);
        
        System.out.println("BEST IS:\t" + best.getF());
        
        for (int i = 0; i < maxgens; i++) {                                     //optimization with GA begins 
            int counter = i+1;
            
            removeClones();
            
            if((counter)%plagueCount==0){                                       //plague to kill off the weakest individuals, refreshing the population
                if(Main_Population.pop.size()>Main_Population.getInitialCount() ){
                Main_Population.plague();
                }
            }
            
            elite = searchforElites((int) (Main_Population.pop.size()*0.1));
            
            newgeneration = spawnNextgen();
            
            Main_Population.pop.addAll(newgeneration.pop);
                        
            removeClones();
            
            best_ever=1e-50;
            for(int j =0; j<Main_Population.pop.size();j++){                         //reconfirm the fittest individual
                if(Main_Population.pop.get(j).getF()>best_ever){
                    best_ever=Main_Population.pop.get(j).getF();
                    best = Main_Population.pop.get(j);
                }
            }
            if (best.getF()>Main_Population.getBest().getF()){
                Main_Population.setBest(best);
            }
            fitnessScore.add(best_ever);
            
            System.out.println("Generation:\t"+counter+"\tBest:\t"+Main_Population.getBest().getF()+"\tCurrent population size:"+Main_Population.pop.size());
            
                                                        //clear the clonearmy to be ready for next generation
            if ((counter) % 10 == 0) {
                top30params.clear();
                elite = searchforElites(30);
                for(int l =0;l<30;l++){
                    double[]params = new double[nparam];
                    System.arraycopy(elite.get(l).x, 0, params, 0, nparam);
                    top30params.add(params);
                }
                                
                if ((best_ever == 1 ) /*|| (best_ever>0.8)*/) {
                    break;
                }
            }
            
            if ((counter)>(plateaulimit+1)){
                double difference = Main_Population.getBest().getF() - (double)fitnessScore.get(i-(plateaulimit+1));
                if(difference<1E-8){
                    top30params.clear();                                               
                    elite = searchforElites(30);
                    for(int l =0;l<30;l++){
                        double[]params = new double[nparam];
                        System.arraycopy(elite.get(l).x, 0, params, 0, nparam);
                        top30params.add(params);
                    }
                    System.out.println("\nGENETIC ALGORITHM FINISHED WITH CONVERGENCE FOR "+plateaulimit+" GENERATIONS\n");
                    
                    break;
                }
            }
            
            if ((counter)>plateaulimit){
                double difference = Main_Population.getBest().getF() - (double)fitnessScore.get(i-plateaulimit);
                if(difference<1E-6){
                    //adaptive mutation to break out of plateau
                    ArrayList<Individual> adaptMutants = Main_Population.mutateAll(); 
                    Population adapt = new Population(nparam, bioSystem, Main_Population.getRange());
                    adapt.pop.addAll(adaptMutants);
                    fitness.evaluatePop(adapt);
                    Main_Population.pop.addAll(adapt.pop);
                }
            }
        }
    }
    
    public ArrayList getParameters() {
        return top30params;
    }
    
    public ArrayList getFitnessTrack(){
        return fitnessScore;
    }
    
    private void removeClones(){
        for(int j=0;j<Main_Population.pop.size()-1; j++){                          //loop to remove individuals that have duplicated fitness
            unwantedclones.clear();
            double[] current = Main_Population.pop.get(j).getOutput();
            for(int k=(j+1);k<Main_Population.pop.size(); k++){                  //idea of this is that it will check everything after as indeces before this would have been scanned/checked
                
                double[] after = Main_Population.pop.get(k).getOutput();
                boolean clone = isitSame(current, after);
                if(clone){                        
                    unwantedclones.add(k);
                }
            }
            if(unwantedclones.size()>0){
                Collections.reverse(unwantedclones);

                for(Object index : unwantedclones){
                    int num = (int) index;
                    Main_Population.pop.remove(num);
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
        ArrayList<Individual> children = Main_Population.nextGeneration();       //produce offspring for the population +10% growth
        ArrayList<Individual> mutants = Main_Population.randomMutation();        //mutate 10% of the population 
        ArrayList<Individual> evolve = Main_Population.mutate(elite);            //mutate the elites
        Individual royalty = Main_Population.royalbirth(elite);                  //1 offspring from the elite
        Individual mutantroyalty = Main_Population.clone(royalty);
        Main_Population.mutateIndividual(mutantroyalty);                         //mutation on the elite offspring;

        Population nextgen = new Population(nparam, bioSystem, Main_Population.getRange());
        nextgen.pop.addAll(children);
        nextgen.pop.addAll(mutants);
        nextgen.pop.addAll(evolve);
        nextgen.pop.add(royalty);
        nextgen.pop.add(mutantroyalty);

        fitness.evaluatePop(nextgen);
        
        return nextgen;
    }
    
    private ArrayList<Individual> searchforElites(int size){
        ArrayList<Individual> search = new ArrayList();
        clonearmy.clear();                                                  //clear the clonearmy to be ready for next generation
        int elitecount = size;

        for(int k =0;k<Main_Population.pop.size();k++){
            Individual clone = Main_Population.clone(Main_Population.pop.get(k));
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
            Individual clone = Main_Population.clone(clonearmy.get(track));
            search.add(clone);
            clonearmy.remove(track);                                        //remove previously fittest to rinse repeat and find the second fittest and so on
        }
        
        return search;
    }    
    
}
