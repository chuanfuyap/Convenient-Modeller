/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package GeneticAlgorithm;

import java.util.ArrayList;
import java.util.Random;
import org.apache.commons.math3.random.RandomDataGenerator;

/**
 *
 * @author chuanfuyap
 */
public class Population {
    
    public ArrayList<Individual> pop = new ArrayList<>();
    private int mutationrate;
    private int nparam;
    private Individual best;
    private RandomDataGenerator random = new RandomDataGenerator();     //apache maths random number generator for uniform distribution
    private Random generator = new Random();
    private int initialcount;
    
    private ArrayList<double[]> Parameter_Range = new ArrayList<double[]>();
    
    private int[] fRIndex;
    private int[] rRIndex;
    private int[] kaIndex;
    private int[] kiIndex;
    
    public Population(int nparam){
        this.nparam=nparam;        
    }
    
    public Population(int nparam,SystemToSolve bioSystem, ArrayList<double[]> All_Range){
        this.nparam=nparam;
        Parameter_Range = All_Range;        
    }
    
    public void initializePopulation(SystemToSolve bioSystem, int popsize) {
        initialcount=popsize;
        
        int[] fRIndex=bioSystem.givefRIndex();          //storing all the forward rate parameter index
        int[] rRIndex=bioSystem.giverRIndex();          //storing all the reverse rate parameter index
        int[] kiIndex=bioSystem.giverKiIndex();
        int[] kaIndex=bioSystem.giverKaIndex();
        best = new Individual(nparam);                  //fittest individual in the entire population, ergo BEST
        
        for (int i = 0; i < popsize; i++) {             
            int counter=0;
            int counter2=0;
            int counter3=0;
            int counter4=0;
            Individual brandnew = new Individual(nparam);       //instead of having each individual store the number as string of binary code, 
            for (int j = 0; j < nparam; j++) {                  //a random with a given range is generated instead, and is the index value for base 10
                if(counter<fRIndex.length && j==fRIndex[counter]) {       //for forward rate constant and backward rate constant range which is thought to have a smaller range in general
//                    double num = this.RN(3, 5.5); // yeast
                    double num = this.RN(1, 3); // leaf
                    brandnew.x[j]=Math.pow(10, num);
                    counter++;
                }else if (counter2<rRIndex.length && j==rRIndex[counter2]) {
                    double num=this.RN(-4, -1); // yeast
                    brandnew.x[j]=Math.pow(10, num);
                    counter2++;
                }else if (counter3<kaIndex.length && j==kaIndex[counter3]) {
                    double num=this.RN(0.5, 2); // yeast
                    brandnew.x[j]=Math.pow(10, num);
                    counter3++;
                }else if (counter4<kiIndex.length && j==kiIndex[counter4]) {
                    double num=this.RN(-5, 2.5); // yeast
                    brandnew.x[j]=Math.pow(10, num);
                    counter4++;
                }else{
                    double num = this.RN(-2.5, 3); //leaf
                    brandnew.x[j]=Math.pow(10, num);
                }
            }
            
            pop.add(brandnew);
        }
    }
    
    public Population(int nparam,SystemToSolve bioSystem, double[][] General_Range){
        this.nparam=nparam;
        fRIndex=bioSystem.givefRIndex();          //storing all the forward rate parameter index
        rRIndex=bioSystem.giverRIndex();          //storing all the reverse rate parameter index
        kaIndex=bioSystem.giverKaIndex();
        kiIndex=bioSystem.giverKiIndex();
        
        int counter=0;
        int counter2=0;
        int counter3=0;
        int counter4=0;
        for(int i =0; i<nparam; i++){
            if(counter<fRIndex.length && i==fRIndex[counter]) {       //for forward rate constant and backward rate constant range which is thought to have a smaller range in general
                Parameter_Range.add(General_Range[0]);
                counter++;
            }else if (counter2<rRIndex.length && i==rRIndex[counter2]) {
                Parameter_Range.add(General_Range[1]);
                counter2++;
            }else if (counter3<kaIndex.length && i==kaIndex[counter3]) {
                Parameter_Range.add(General_Range[2]);
                counter3++;
            }else if (counter4<kiIndex.length && i==kiIndex[counter4]) {
                Parameter_Range.add(General_Range[3]);
                counter4++;
            }else{
                Parameter_Range.add(General_Range[4]);
            }
        }
    }
    
    public ArrayList<Individual> nextGeneration(){
        int growth= (int) (pop.size()*0.1);        
        ArrayList<Individual> children = new ArrayList<>();
        
        for(int i =0; i<growth;i++){
            Individual parentA = tournament();
            Individual parentB = tournament();
            
            Individual child = new Individual(nparam);
            crossover(parentA, parentB, child);
            children.add(child);
        }
        
        return children;
    }
    
    public ArrayList<Individual> randomMutation(){
        
        ArrayList<Individual> mutants = new ArrayList<>();
        mutationrate = (int) (pop.size()*0.1);
        for(int i =0 ; i<mutationrate;i++){
            int mem = (int) (this.RN() * pop.size());
            Individual mutant = clone(pop.get(mem));
            mutateIndividual(mutant);
            mutants.add(mutant);
        }
        
        return mutants;
    }
    
    public ArrayList<Individual> mutate(ArrayList<Individual> elite){
        ArrayList<Individual> mutants = new ArrayList<>();
        for(int i =0 ; i<elite.size();i++){
            Individual mutant = clone(elite.get(i));
            mutateIndividual(mutant);
            mutants.add(mutant);
        }
        return mutants;
    }
    
    public Individual royalbirth(ArrayList<Individual> elite){
        Individual royalty = new Individual(nparam);
        int mem = (int) (this.RN() * elite.size());
        int mem2 = (int) (this.RN() * elite.size());
        
        crossover(elite.get(mem),elite.get(mem2), royalty);
        
        return royalty;
    }
    
    public void setBest(Individual bestindividual){
        best = bestindividual;
    }
    public Individual getBest(){
        return best;
    }
    public ArrayList<Individual> getPop(){
        return pop;
    }
    
    private double RN(double lower, double upper) {
        return random.nextUniform(lower, upper);
    }
    
    private double RN() {
        return generator.nextDouble();
    }
    
    private Individual tournament(){

        int tourn_size = pop.size()/20;//3;
        int i;
        int mem;

        int fittest = -1;
        double best_fitness = 1e-200;
        i = 0;
        while (i < tourn_size) {
            mem = (int) (this.RN() * pop.size());
            if (pop.get(mem).getF() > best_fitness) {
                best_fitness =pop.get(mem).getF();
                fittest = mem;
            }
            i++;
        }
        if (fittest == -1) {
            System.out.println("Fittest uninitialized");
        }
        
        return pop.get(fittest);

    }
    
    private void crossover(Individual p1, Individual p2, Individual ch) {

        // uniform crossover 
        for (int i = 0; i < nparam; i++) {
            if (this.RN() < 0.5) {
                ch.x[i] = p1.x[i];
            } else {
                ch.x[i] = p2.x[i];
            }
        }

    }
    
    public void mutateIndividual(Individual mem) {
        
        int AllorOne = 1;
        if(AllorOne==0){
            int increaseORdrop = (int) (0.5 + this.RN());
            Double targetparam = this.RN()*nparam;
            int targetparam1 = targetparam.intValue();
        
            if(increaseORdrop==0){
                mem.x[targetparam1]=mem.x[targetparam1]*1.13;
            }else if (increaseORdrop==1){
                mem.x[targetparam1]=mem.x[targetparam1]*0.87;
            }            
        }else{            
            for(int i = 0; i<nparam; i++){
                int increaseORdrop = (int) (0.5 + this.RN());                   //to determine if it goes up or down in value
                int change = (int) (0.5 + this.RN());                           //to determine if this parameter is going to be modified or not 
                if(change==1){
                    if(increaseORdrop==0){
                        mem.x[i]=mem.x[i]*1.13;
                    }else if (increaseORdrop==1){
                        mem.x[i]=mem.x[i]*0.87;
                    }
                }
            }            
        }

    }
    
    public void plague(){
        int deathcount=(int) (pop.size()-initialcount);
        
        for(int count =0; count<deathcount;count++){
            double small =1e100;
            int track=-1;
        
            for(int i =0;i<pop.size();i++){
                if((double)pop.get(i).getF()<small){
                    small=(double)pop.get(i).getF();
                    track=i;
                }
            }
            pop.remove(track);
        }
    }
    
    public Individual clone(Individual ind){
        Individual clone = new Individual(ind.size);
        clone.setF(ind.getF());
        for(int i=0;i<ind.size;i++){
            clone.x[i]=ind.x[i];
        }        
        return clone;
    }
    
    public void updateInitialCount(int newcount){
        initialcount = newcount;
    }
    
    public int getInitialCount(){
        return initialcount;
    }
    
    public ArrayList<double[]> getRange(){
        return Parameter_Range;
    }
    
    public ArrayList<Individual> mutateAll(){
        ArrayList<Individual> mutants = new ArrayList<>();
        
        for(int i =0 ; i<pop.size();i++){
            Individual mutant = clone(pop.get(i));
            mutateIndividual(mutant);
            mutants.add(mutant);
        }
        
        return mutants;
    }
}
