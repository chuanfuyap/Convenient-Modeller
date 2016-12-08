/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package GeneticAlgorithm;

import Main.ModelReaction;
import Main.SBMLFile;
import Species.Compound;
import Species.Enzyme;
import SteadyState.SolveSteadyState;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import javax.xml.stream.XMLStreamException;
import odeBridge.TCodesolve;
import odeBridge.TCodesolvePK;
import odeBridge.odesolve;
import org.sbml.jsbml.Model;
import org.sbml.jsbml.validator.ModelOverdeterminedException;

/**
 *
 * @author chuanfuyap
 */
public class SystemToSolve {
    private int totalParameters=0;
    private double[] metConc;
    private double[] flux;
    private double[] activeConc;
    private ArrayList<Compound> Compounds;
    private ArrayList<ModelReaction> Reactions;
    private ArrayList<Enzyme> Enzymes;
    private Model Model;
    private int activeCompNum;
    private double[] estimatedMet;
    private double[] estimatedFlux;
    private ArrayList<Integer> fRateCIndex = new ArrayList<>(); //rate constant index
    private ArrayList<Integer> rRateCIndex = new ArrayList<>();
    private HashMap metMap = new HashMap<>();
    private HashMap fluxMap = new HashMap<>();
    private HashMap estMetMap = new HashMap<>();
    private HashMap estFluxMap = new HashMap<>();
    private HashMap TCmetMap = new HashMap<>();
    private HashMap TCprotMap = new HashMap<>();
    private HashMap TCfluxMap = new HashMap<>();
    private HashMap TCestMetMap = new HashMap<>();
    private HashMap TCestProtMap = new HashMap<>();
    private HashMap TCestFluxMap = new HashMap<>();
    private SBMLFile modelreactions;
    private String link;
    private boolean SS;
    private ReadData data;
    private boolean goodbad=false;
    private boolean reachplateau=false;
       
    public SystemToSolve(String Link, ReadData data) throws IOException, XMLStreamException, ModelOverdeterminedException{
        this.data=data;
        SS = data.SSorTC();
        File SBMLFile = new File(Link);
        this.modelreactions = new SBMLFile(SBMLFile);
        this.Compounds=modelreactions.returncompounds();
        this.Reactions=modelreactions.returnReactionList();
        this.Enzymes=modelreactions.returnenzymes();
        this.link=Link;
        
        //to store number of active compounds that would be varying in concentration;
        activeCompNum=0;
        for (Compound comp: Compounds){
            if (comp.getBoundaryCondition()==false){
                activeCompNum++;
            }
        }
        
        HashMap ssdata;
        HashMap tcdata;
        if(SS==true){
            ssdata = data.getSSdata();
            //storing the initial concentration of the active compounds and the SS value from input data for PE;
            activeConc = new double[activeCompNum];
            metConc = new double[activeCompNum];
            int ticktock=0;
            for (Compound compound : Compounds){
                if (compound.getBoundaryCondition()==false){
                    metConc[ticktock] = (double) ssdata.get(compound.getID());
                    metMap.put(compound.getID(), ssdata.get(compound.getID()));
                    activeConc[ticktock] = metConc[ticktock];
                    ticktock++;                
                }
            }

            flux = new double[Reactions.size()];
            ticktock=0;
            for (ModelReaction reaction : Reactions){
                int numSub = reaction.getSubstrates().size();
                int numProd = reaction.getProducts().size();

                int regulation = reaction.getRegulation();
                int para =2;
                if(regulation>1){
                    para+=1;
                }
                para+=numSub;
                para+=numProd;
                fRateCIndex.add(totalParameters);
                rRateCIndex.add(totalParameters+1);
                totalParameters+=para;

                flux[ticktock] = (double) ssdata.get(reaction.getReactionID());
                fluxMap.put(reaction.getReactionID(),ssdata.get(reaction.getReactionID()));
                ticktock++;
            }
        }else{
            tcdata= data.getTCdata();
            for (Compound compound : Compounds){
                if(compound.getBoundaryCondition()==false){
                
                TCmetMap.put(compound.getID(), (double[]) tcdata.get(compound.getID()));
                
                }
            }
            for (Enzyme enzyme : Enzymes){
                TCprotMap.put(enzyme.getID(), (double[]) tcdata.get(enzyme.getID()));
            }
            for (ModelReaction reaction : Reactions){
                
                int numSub = reaction.getSubstrates().size();
                int numProd = reaction.getProducts().size();

                int regulation = reaction.getRegulation();
                int para =2;
                if(regulation>1){
                    para+=1;
                }
                para+=numSub;
                para+=numProd;
                fRateCIndex.add(totalParameters);
                rRateCIndex.add(totalParameters+1);
                totalParameters+=para;
                
                TCfluxMap.put(reaction.getReactionID(), (double[]) tcdata.get(reaction.getReactionID()));
            }
        }        
    }
    
    public void runEst(double[] parameters) throws ModelOverdeterminedException, InstantiationException, IllegalAccessException, IllegalArgumentException, NoSuchMethodException, XMLStreamException, IOException {
        estimatedMet = new double[activeCompNum];
        estimatedFlux = new double[Reactions.size()];
        
        if(SS==true){
            SolveSteadyState steadystateresults = new SolveSteadyState(Compounds, Reactions, activeConc, parameters, link);
            
            if(steadystateresults.solveSS()==true){
                estMetMap = steadystateresults.getSolvedMMap();
                estFluxMap = steadystateresults.getSolvedFMap();
                goodbad=true;
            }else{
                goodbad=false;
            }
            
        }else{
            double[] time = data.returnTime();
            double[] parameter = new double[getParametersCount()];  
            int activeCompNum=0;
            for (Compound comp: Compounds){
                if (comp.getBoundaryCondition()==false){
                    activeCompNum++;
                }
            }

            String[] paranames = modelreactions.getParameterNames();
            String[] reactionID = modelreactions.getReactionID();
            String[] variables= new String[activeCompNum];
            String[] reactionID2 = modelreactions.getReactionID2();
            String[] proteinID = modelreactions.getProteinID();

            int count=0;
            for(Compound comp : Compounds){
                if(comp.getBoundaryCondition()==false){
                    variables[count]=comp.getID();
                    count++;
                }
            }

            if(proteinID.length>0){
                TCodesolvePK ode = new TCodesolvePK(link, variables, getParametersCount(), paranames, reactionID, reactionID2, proteinID);

                double output[][] = ode.runsolver(parameter, time);

                for(int j=1;j<variables.length+1+reactionID2.length+proteinID.length;j++){
                    double[] temparray = new double[time.length];
                    for(int i =0;i<time.length;i++){                   
                        temparray[i]=output[i][j];
                    }
                }
            }else{
                TCodesolve ode = new TCodesolve(link, variables, getParametersCount(), paranames, reactionID, reactionID2);

                double output[][] = ode.runsolver(parameter, time);

                for(int j=1;j<variables.length+1+reactionID2.length+proteinID.length;j++){
                    double[] temparray = new double[time.length];
                    for(int i =0;i<time.length;i++){                   
                        temparray[i]=output[i][j];
                    }
                    if((j-1)<activeCompNum){
                        TCestMetMap.put(variables[j], temparray);
                    }else{
                        TCestFluxMap.put(reactionID2[j-activeCompNum],temparray);
                    }
                }
            }
        }
    }

    public double[] getEstMet(){
        return estimatedMet;
    }
    public HashMap getEstMetMap(){
        return estMetMap;
    }
    public double[] getEstFlux(){
        return estimatedFlux;
    }
    public HashMap getEstFluxMap(){
        return estFluxMap;
    }
    
    public HashMap getTCmetMap(){
        return TCmetMap;
    }
    
    public HashMap getTCprotMap(){
        return TCprotMap;
    }
    
    public HashMap getTCfluxMap(){
        return TCfluxMap;
    }
    
    public HashMap getTCestmetMap(){
        return TCestMetMap;
    }
    
    public HashMap getTCestprotMap(){
        return TCestProtMap;
    }
    
    public HashMap getTCestfluxMap(){
        return TCestFluxMap;
    }
    public int[] givefRIndex(){
        int[] rateindex = new int[fRateCIndex.size()];
        for(int i=0;i<rateindex.length; i++){
            rateindex[i]=fRateCIndex.get(i);
        }
        
        return rateindex;
    }
    
    public int[] giverRIndex(){
        int[] rateindex = new int[rRateCIndex.size()];
        for(int i=0;i<rateindex.length; i++){
            rateindex[i]=rRateCIndex.get(i);
        }
        
        return rateindex;
    }
    
    public int getActiveCount(){
        return activeCompNum;
    }
    
    public int getParametersCount(){
        return totalParameters;
    }
    
    public double[] getMetConc(){
        return metConc;
    }
    
    public HashMap getMetMap(){
        return metMap;
    }
    
    public double[] getFlux(){
        return flux;
    }
    
    public HashMap getFluxMap(){
        return fluxMap;
    }
    
    public ArrayList<Compound> getCompounds(){
        return Compounds;
    }
    
    public ArrayList<ModelReaction> getReactions(){
        return Reactions;
    }
    
    public boolean goodsolution(){
        return goodbad;
    }
}
