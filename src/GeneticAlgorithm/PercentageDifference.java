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
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.xml.stream.XMLStreamException;
import org.sbml.jsbml.validator.ModelOverdeterminedException;

/**
 *
 * @author chuanfuyap
 */

public class PercentageDifference implements Runnable {
    
    HashMap expMet;             //expected metconc values
    HashMap expFlux;            //expected flux values
    SystemToSolve bioSystem;
    ArrayList<Compound> Compounds;
    ArrayList<ModelReaction> Reactions;
    int missing = 0;            //to keep track of bad individuals
    ArrayList badindex = new ArrayList<>();
    Individual Ind;
    int index;
    Population pop;
    
    public PercentageDifference (SystemToSolve bioSystem, HashMap expMet, HashMap expFlux, Population pop, int index){
        this.index=index;
        this.expMet=expMet;
        this.expFlux=expFlux;
        this.bioSystem=bioSystem;
        this.pop=pop;
        this.Compounds=bioSystem.getCompounds();
        this.Reactions=bioSystem.getReactions();
    }
    
    @Override
    public void run() {
        double summed_error = 0.0;
        double[] params = pop.pop.get(index).x;
        double f;
        double mettotal=0;
        double fluxtotal=0;
        
        try {
            bioSystem.runEst(params);
        } catch (ModelOverdeterminedException | InstantiationException |IllegalAccessException | IllegalArgumentException | NoSuchMethodException | XMLStreamException | IOException ex) {
            Logger.getLogger(PercentageDifference.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        HashMap estMet = bioSystem.getEstMetMap();      //estimated metconc values
        HashMap estFlux = bioSystem.getEstFluxMap();    //estimated flux values
        
        pop.pop.get(index).setMet(estMet);
        pop.pop.get(index).setFlux(estFlux);
        int mettrack;
        int fluxtrack;
        if(expMet!=null){
            mettrack = expMet.size();
        }else{
            mettrack = 0;
        }
        if(expFlux!=null){
            fluxtrack = expFlux.size();
        }else{
            fluxtrack = 9;
        }
        
        
        if(bioSystem.goodsolution()==false){
            pop.pop.get(index).setF(1e-100);
        }else{
            if(expMet!=null){
            for (Compound comp : Compounds){
                if(comp.getBoundaryCondition()==false){
                    String ID = comp.getID();
                    if(expMet.containsKey(ID)==true && (double) expMet.get(ID)>0){
                        double difference =(double) expMet.get(ID)-(double)estMet.get(ID);
                        double percentagedifference=0;
                        if((double) expMet.get(ID)>0){
                            percentagedifference = (Math.abs(difference)/(double) expMet.get(ID))*100;
                        }else{
                            percentagedifference = Math.abs(difference)*100;
                        }                        
                        mettotal+=Math.pow(percentagedifference, 2);
                    }else if (expMet.get(ID)==null){
                        mettrack--;
                    }else{
                        mettrack--;
                    }
                }
            }
            }else{
                mettrack=0;
            }
            
            if(expFlux!=null){
            for(ModelReaction reaction : Reactions){
                String ID = reaction.getReactionID();
                if(expFlux.containsKey(ID)==true && (double) expFlux.get(ID)>0){
                    double difference = (double) expFlux.get(ID) - (double) estFlux.get(ID);
                    double percentagedifference=0;
                    if((double) expFlux.get(ID)>0){
                        percentagedifference = (Math.abs(difference)/(double) expFlux.get(ID))*100;
                    }else{
                        percentagedifference = Math.abs(difference)*100;
                    }
                    fluxtotal+=Math.pow(percentagedifference, 2);
                }else if (expFlux.get(ID)==null){
                    fluxtrack--;
                }else{
                    fluxtrack--;
                }
            }
            }else{
                fluxtrack=0;
            }
            
            if(mettrack==0){
                double fluxdifference = fluxtotal/fluxtrack;
                summed_error=fluxdifference;
            }else if (fluxtrack==0){
                double metdifference = mettotal/mettrack;
                summed_error=metdifference;
            }else{
                double metdifference = mettotal/mettrack;
                metdifference*=0.5;
                double fluxdifference = fluxtotal/fluxtrack;
                fluxdifference*=0.5;

                summed_error+=metdifference+fluxdifference;
            }

            f = 1.0 / (1.0 + summed_error);         //bigger the error, lower the f, closer to 1 the better the fitness
            System.out.println(f);
            pop.pop.get(index).setF(f);
        }
        
    }
    
}
