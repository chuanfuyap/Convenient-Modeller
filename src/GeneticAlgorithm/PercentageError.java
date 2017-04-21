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

public class PercentageError implements Runnable {
    
    HashMap expMet;             //expected metconc values
    HashMap expFlux;            //expected flux values
    HashMap TCexpMet;             //expected metconc values
    HashMap TCexpFlux;            //expected flux values
    SystemToSolve bioSystem;
    ArrayList<Compound> Compounds;
    ArrayList<ModelReaction> Reactions;
    int missing = 0;            //to keep track of bad individuals
    ArrayList badindex = new ArrayList<>();
    Individual Ind;
    int index;
    Population pop;
    boolean SSorTC;
    
    public PercentageError (SystemToSolve bioSystem, Population pop, int index){
        this.index=index;
        this.bioSystem=bioSystem;
        this.pop=pop;
        this.Compounds=bioSystem.getCompounds();
        this.Reactions=bioSystem.getReactions();
        this.SSorTC=bioSystem.isitSSorTC();
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
            Logger.getLogger(PercentageError.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        if(SSorTC==true){
        
        HashMap estMet = bioSystem.getEstMetMap();      //estimated metconc values
        HashMap estFlux = bioSystem.getEstFluxMap();    //estimated flux values
        expMet = bioSystem.getMetMap();      //estimated metconc values
        expFlux = bioSystem.getFluxMap();    //estimated flux values
        
        int mettrack=0;
        int fluxtrack=0;

        if(bioSystem.goodsolution()==false){
            pop.pop.get(index).setF(1e-100);
        }else{
            if(expMet!=null){
            for (Compound comp : Compounds){
                if(comp.getBoundaryCondition()==false){
                    String ID = comp.getID();
                    if(expMet.containsKey(ID)==true && (double) expMet.get(ID)>=0){
                        double difference =(double) expMet.get(ID)-(double)estMet.get(ID);
                        double percentagedifference=0;
                        if((double) expMet.get(ID)>0){
                            percentagedifference = (Math.abs(difference)/(double) expMet.get(ID))*100;
                        }else{
                            percentagedifference = Math.abs(difference)*100;
                        }                        
                        mettotal+=Math.pow(percentagedifference, 2);
                        mettrack++;
                    }
                }
            }
            }
            
            if(expFlux!=null){
            for(ModelReaction reaction : Reactions){
                String ID = reaction.getReactionID();
                if(expFlux.containsKey(ID)==true && (double) expFlux.get(ID)>=0){
                    double difference = (double) expFlux.get(ID) - (double) estFlux.get(ID);
                    double percentagedifference=0;
                    if((double) expFlux.get(ID)>0){
                        percentagedifference = (Math.abs(difference)/(double) expFlux.get(ID))*100;
                    }else{
                        percentagedifference = Math.abs(difference)*100;
                    }
                    fluxtotal+=Math.pow(percentagedifference, 2);
                    fluxtrack++;
                }
            }
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
            pop.pop.get(index).setF(f);
        }
        
        }
        else{
            HashMap TCMetEst = bioSystem.getTCestmetMap();
            HashMap TCFluxEst = bioSystem.getTCestfluxMap();
            TCexpMet = bioSystem.getTCmetMap();
            TCexpFlux = bioSystem.getTCfluxMap();
            
            int mettrack=0;
            int fluxtrack=0;
            int timepoints=bioSystem.timelength();
            
            if(TCexpMet!=null){
                for (Compound comp : Compounds){
                    if(comp.getBoundaryCondition()==false){
                        String ID = comp.getID();
                        if(TCexpMet.containsKey(ID)){
                            double[] estMetTC =(double[]) TCMetEst.get(ID);
                            double[] expMetTC =(double[]) TCexpMet.get(ID);
                            double percentagedifference=0;
                            int track = 0;
                            for(int i = 0; i<timepoints; i++){
                                if(expMetTC[i]>=0){
                                    double difference =Math.abs(expMetTC[i]-estMetTC[i]);
                                    if(expMetTC[i]>0){
                                        percentagedifference += (difference/expMetTC[i])*100;
                                    }
                                    else{
                                        percentagedifference += difference*100;
                                    }
                                    track++;
                                }
                            }
                            double avgdiff = percentagedifference/track;
                            mettotal+=Math.pow(avgdiff, 2);
                            mettrack++;
                        }

                    }
                }
            }
            
            if(TCexpFlux!=null){
                for (ModelReaction reaction : Reactions){
                    String ID = reaction.getReactionID();
                    if(TCexpFlux.containsKey(ID)){
                        double[] estFluxTC =(double[]) TCFluxEst.get(ID);
                        double[] expFluxTC =(double[]) TCexpFlux.get(ID);
                        double percentagedifference=0;
                        int track = 0;
                        for(int i = 0; i<timepoints; i++){
                            if(expFluxTC[i]>=0){
                                double difference =Math.abs(expFluxTC[i]-estFluxTC[i]);
                                if(expFluxTC[i]>0){
                                    percentagedifference += (difference/expFluxTC[i])*100;
                                }
                                else{
                                    percentagedifference += difference*100;
                                }
                                track++;
                            }
                        }
                        double avgdiff = percentagedifference/track;
                        fluxtotal+=Math.pow(avgdiff, 2);
                        fluxtrack++;
                    }
                }
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
            pop.pop.get(index).setF(f);
        }
    }
    
}
