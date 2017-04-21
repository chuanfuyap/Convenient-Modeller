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
    private boolean SSorTC;
    private ReadData data;
    private boolean goodbad=false;
       
    public SystemToSolve(String Link, ReadData data) throws IOException, XMLStreamException, ModelOverdeterminedException{
        this.data=data;
        SSorTC = data.SSorTC();
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
        if(SSorTC==true){
            ssdata = data.getSSdata();
            boolean do_we_have_metabolites = is_there_metabolites(ssdata);
            boolean do_we_have_fluxes = is_there_fluxes(ssdata);
            //storing the initial concentration of the active compounds and the SS value from input data for PE;
            activeConc = new double[activeCompNum];
            metConc = new double[activeCompNum];
            int ticktock=0;
            if (do_we_have_metabolites!=false){
            for (Compound compound : Compounds){
                if (compound.getBoundaryCondition()==false){
                    if(ssdata.containsKey(compound.getID())==true){
                        metConc[ticktock] = (double) ssdata.get(compound.getID());
                        metMap.put(compound.getID(), ssdata.get(compound.getID()));
                        activeConc[ticktock] = metConc[ticktock];
                        ticktock++;
                    }else{
                        activeConc[ticktock] = compound.getConcentration();
                        ticktock++;
                    }
                }
            }
            }else{
                metMap = null;
                for (Compound compound : Compounds){
                    if (compound.getBoundaryCondition()==false){
                        activeConc[ticktock] = compound.getConcentration();
                        ticktock++;
                    }
                }
            }
            
            
            flux = new double[Reactions.size()];
            
            for(ModelReaction reaction : Reactions){
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
            }
            
            ticktock=0;
            if (do_we_have_fluxes!=false){
            for (ModelReaction reaction : Reactions){
                
                if(ssdata.containsKey(reaction.getReactionID())==true){
                    flux[ticktock] = (double) ssdata.get(reaction.getReactionID());
                    fluxMap.put(reaction.getReactionID(),ssdata.get(reaction.getReactionID()));
                    ticktock++;
                }
            }
            }else{
                fluxMap =null;
            }
        }
        
        else{
            tcdata= data.getTCdata();
            boolean do_we_have_metabolites = is_there_metabolites(tcdata);
            boolean do_we_have_fluxes = is_there_fluxes(tcdata);
            boolean do_we_have_proteins = is_there_proteins(tcdata);
            if (do_we_have_metabolites!=false){
            for (Compound compound : Compounds){
                if(compound.getBoundaryCondition()==false){
                if(tcdata.containsKey(compound.getID())==true){
                TCmetMap.put(compound.getID(), (double[]) tcdata.get(compound.getID()));
                }
                }
            }
            }else{
                TCmetMap =null;
            }
            
            if (do_we_have_proteins!=false){
            for (Enzyme enzyme : Enzymes){
                if(tcdata.containsKey(enzyme.getID())==true){
                TCprotMap.put(enzyme.getID(), (double[]) tcdata.get(enzyme.getID()));
                }
            }
            }else{
                TCprotMap =null;
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
            }
            
            if (do_we_have_fluxes!=false){
            for (ModelReaction reaction : Reactions){
                if(tcdata.containsKey(reaction.getReactionID())==true){
                    TCfluxMap.put(reaction.getReactionID(), (double[]) tcdata.get(reaction.getReactionID()));
                }
            }
            }else{
                TCfluxMap =null;
            }
        }
    }
    
    public void runEst(double[] parameters) throws ModelOverdeterminedException, InstantiationException, IllegalAccessException, IllegalArgumentException, NoSuchMethodException, XMLStreamException, IOException {
        estimatedMet = new double[activeCompNum];
        estimatedFlux = new double[Reactions.size()];
        
        if(SSorTC==true){
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
            double[] parameter = parameters;  
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
                
                for(int j=1;j<variables.length+1+reactionID2.length;j++){
                    double[] temparray = new double[time.length];
                    for(int k =0;k<time.length;k++){                   
                        temparray[k]=output[k][j];
                    }
                    if((j-1)<activeCompNum){
                        TCestMetMap.put(variables[j-1], temparray);
                    }else{
                        TCestFluxMap.put(reactionID2[j-activeCompNum-1],temparray);
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
    
    public int timelength(){
        return data.returnTime().length;
    }
    
    private boolean is_there_metabolites(HashMap data){
        boolean is_it_there = false;
        
        for (Compound compound : Compounds){
            if(compound.getBoundaryCondition()==false){
                String compound_ID = compound.getID();
                is_it_there = data.containsKey(compound_ID);
            }
        }
        
        return is_it_there;
    }
    
    private boolean is_there_fluxes(HashMap data){
        boolean is_it_there = false;
        
        for (ModelReaction reaction : Reactions){
            String reaction_ID = reaction.getReactionID();
            is_it_there = data.containsKey(reaction_ID);
        }
        
        return is_it_there;
    }
    
    private boolean is_there_proteins(HashMap data) {
        boolean is_it_there = false;
        
        for (Enzyme enzyme : Enzymes){
            String enzyme_ID = enzyme.getID();
            is_it_there = data.containsKey(enzyme_ID);
        }
        
        return is_it_there; 
    }
    
    public boolean isitSSorTC(){
        return SSorTC;
    }
}
