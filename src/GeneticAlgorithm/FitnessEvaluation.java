/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package GeneticAlgorithm;

import Main.ModelReaction;
import Species.Compound;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import javax.xml.stream.XMLStreamException;
import org.sbml.jsbml.validator.ModelOverdeterminedException;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

/**
 *
 * @author chuanfuyap
 */
public class FitnessEvaluation  {
    
    SystemToSolve bioSystem;
    ArrayList<Compound> Compounds;
    ArrayList<ModelReaction> Reactions;
    int missing = 0;            //to keep track of bad individuals
    ArrayList badindex = new ArrayList<>();
    ArrayList fitness = new ArrayList<>();
    ThreadPoolExecutor executor = null;
    public FitnessEvaluation(SystemToSolve bioSystem){
        this.bioSystem=bioSystem;
        this.Compounds=bioSystem.getCompounds();
        this.Reactions=bioSystem.getReactions();
    }
    
    public void InitialEvaluation(Population pop) throws ModelOverdeterminedException, InstantiationException, IllegalAccessException, IllegalArgumentException, NoSuchMethodException, InterruptedException, XMLStreamException, IOException{
    
        executor = (ThreadPoolExecutor) Executors.newFixedThreadPool(10);
        for(int i =0; i < pop.pop.size();i++){        
            PercentageError PET = new PercentageError(bioSystem, pop, i);
            executor.execute(PET);
        }
        executor.shutdown();
        while (!executor.awaitTermination(30, TimeUnit.MINUTES)) {
        }
        missing = 0;
        badindex.clear();
        for(int i =0; i < pop.pop.size();i++){
            if(pop.pop.get(i).getF()==1e-100){
                missing++;
                badindex.add(i);
            }
        }
        
        if(missing>0){
            Collections.reverse(badindex);

            for(Object index : badindex){
                int num = (int) index;
                pop.pop.remove(num);
            }
        }
    }
    
    public void evaluatePop(Population pop) throws ModelOverdeterminedException, InstantiationException, IllegalAccessException, IllegalArgumentException, NoSuchMethodException, InterruptedException, XMLStreamException, IOException{

        executor = (ThreadPoolExecutor) Executors.newFixedThreadPool(10);
        for(int i =0; i < pop.pop.size();i++){        
            PercentageError PET = new PercentageError(bioSystem, pop, i);
            executor.execute(PET);
        }
        executor.shutdown();
        while (!executor.awaitTermination(30, TimeUnit.MINUTES)) {
        }
    }
    
    public int getMissing(){
        return missing;
    }
    public void stopGA(){
        executor.shutdownNow();
    }
}
