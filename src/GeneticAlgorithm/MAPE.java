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

public class MAPE implements Runnable {
    
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
    
    public MAPE (SystemToSolve bioSystem, Population pop, int index){
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
        double mettotal=0;
        double fluxtotal=0;        
        
        try {
            bioSystem.solve_ODE_model(params, Integer.toHexString(this.hashCode()));
        } catch (ModelOverdeterminedException | InstantiationException |IllegalAccessException | IllegalArgumentException | NoSuchMethodException | XMLStreamException | IOException ex) {
            Logger.getLogger(MAPE.class.getName()).log(Level.SEVERE, null, ex);
        }
        ArrayList<ConditionInfo> conditions_list = bioSystem.get_Condition_List();
        
        if(SSorTC==true){       // for STEADY STATE output comparison and scoring
            double[] f = new double[conditions_list.size()];
        for(int i =0;i<conditions_list.size();i++){
            HashMap estMet = conditions_list.get(i).get_Solved_metmap();      //estimated metconc values
            HashMap estFlux = conditions_list.get(i).get_Solved_fluxmap();    //estimated flux values
            HashMap expMet = conditions_list.get(i).get_metabolites_info();      //expected metconc values
            HashMap expFlux = conditions_list.get(i).get_fluxes_info();    //expected flux values

            int mettrack=0;
            int fluxtrack=0;

            if(bioSystem.goodsolution()==false){
                f[i] = 1e-100;
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

                f[i] = 1.0 / (1.0 + summed_error);         //bigger the error, lower the f, closer to 1 the better the fitness
                
            }
        }
            double sumf = 0;
            for(double d: f)sumf+=d;
            pop.pop.get(index).setF(sumf/f.length);
            
            HashMap estMet = conditions_list.get(0).get_Solved_metmap();      //estimated metconc values
            HashMap estFlux = conditions_list.get(0).get_Solved_fluxmap();    //estimated flux values
            
            double[] output = new double[estMet.size()+estFlux.size()];
            if (!estMet.isEmpty()){
              
            int counter=0;
            for (Compound comp : Compounds){
            if(comp.getBoundaryCondition()==false){
                String ID = comp.getID();
                double expected_met_value = (double) estMet.get(ID);
                output[counter]=(double) estMet.get(ID);
                counter++;
                }
            }
            for(ModelReaction react : Reactions ){
                String ID = react.getMyReaction().getReactionID();
                output[counter]=(double) estFlux.get(ID);
                counter++;
            }
            
            }
            else{
                for (int i =0;i<output.length;i++){
                    output[i]=0;
                }
            }
            pop.pop.get(index).setOutput(output);
        
        }
        else{       //for TIME SERIES output comparison and scoring
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
            
            double[] output = new double[TCMetEst.size()+TCMetEst.size()];
            int counter=0;
            for (Compound comp : Compounds){
            if(comp.getBoundaryCondition()==false){
                String ID = comp.getID();
                double[] estFluxTC =(double[]) TCFluxEst.get(ID);
                output[counter]=estFluxTC[timepoints-1];
                counter++;
                }
            }
            for(ModelReaction react : Reactions ){
                String ID = react.getMyReaction().getReactionID();
                double[] estFluxTC =(double[]) TCFluxEst.get(ID);
                output[counter]=estFluxTC[timepoints-1];
                counter++;
            }
            
            pop.pop.get(index).setOutput(output);
        }
    }
    
}