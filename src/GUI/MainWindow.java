/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package GUI;

import GeneticAlgorithm.ReadData;
import GeneticAlgorithm.SystemToSolve;
import Main.ModelReaction;
import Main.TsvToModel;
import Reaction.MetabolicReaction;
import SBML.ExportSBML;
import SBML.ImportSBML;
import Species.Compound;
import Species.Enzyme;
import java.awt.HeadlessException;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.KeyEvent;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.DefaultCellEditor;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableColumn;
import javax.xml.stream.XMLStreamException;
import org.sbml.jsbml.KineticLaw;
import org.sbml.jsbml.ListOf;
import org.sbml.jsbml.LocalParameter;
import org.sbml.jsbml.Model;
import org.sbml.jsbml.ModifierSpeciesReference;
import org.sbml.jsbml.Reaction;
import org.sbml.jsbml.SBMLDocument;
import org.sbml.jsbml.SBMLException;
import org.sbml.jsbml.SBMLWriter;
import org.sbml.jsbml.Species;
import org.sbml.jsbml.SpeciesReference;
import org.sbml.jsbml.validator.ModelOverdeterminedException;

/**
 *
 * @author chuanfuyap
 */
public class MainWindow extends javax.swing.JFrame {

    /**
     * Creates new form MainWindow
     */
    
    
    public MainWindow() {
        initComponents();
        SpeciesBox.addItem("N/A");
        SpeciesTable.getModel().addTableModelListener(SpeciesTableChanged);
        EnzymesTable.getModel().addTableModelListener(EnzymeTableChanged);
        ReactionTable.getModel().addTableModelListener(ReactionTableChanged);
        
    }
    
    private void ImportModel(ImportSBML model) {
        importedmodel=true;
        try{
            SpeciesBox.removeAllItems();
            EnzymeBox.removeAllItems();
            SpeciesBox.addItem("N/A");
            allthereactions = new ArrayList<>();
            allthespecies = new ArrayList<>();
            alltheenzymes = new ArrayList<>();
            Model = model.getModel();
            listofspecies = model.getSpecies();
            listofreactions = model.getReactions();
            modelName = model.getModelName();
            if (modelName == null) {
                modelName = Model.getId();
            }
            if (Model.getId()== null || Model.getId().isEmpty()){
                modelID=Model.getName().replaceAll("\\s+","");
            }
            
            ModelNameTextField.setText(modelName);
            
            addSpeciesToSpeciesList(listofspecies);
            addReactionsToReactionList(listofreactions);
            
            
        } catch (NullPointerException e) {
            JOptionPane.showMessageDialog(null, "Incompatible File format.", "Import Error", JOptionPane.ERROR_MESSAGE);
        } catch (Exception ee) {
            JOptionPane.showMessageDialog(null, "Please check SBML document.", "Import Error", JOptionPane.ERROR_MESSAGE);
        }
    }
    
    private void addReactionsToReactionList(ListOf<Reaction> listofreactions){
        
        for (int i = 0; i < listofreactions.size(); i++){           //reads the model's list of reactions.
            Reaction reaction = listofreactions.get(i);             //extracts each reaction found in the list, go through them 1 by 1 in the For loop.
            MetabolicReaction mReaction = new MetabolicReaction(reaction);
            ListOf modifiers = (ListOf) reaction.getListOfModifiers();
            ListOf reactionsubstrates = (ListOf) reaction.getListOfReactants();
            ListOf reactionproducts = (ListOf) reaction.getListOfProducts();
            Enzyme e = null;
            ArrayList<LocalParameter> parameters = new ArrayList<LocalParameter>();
            int regulation =1;
            ArrayList<Compound> mod = new ArrayList();
            
            for (int j = 0; j <reactionsubstrates.size(); j++){
                SpeciesReference sf = (SpeciesReference) reactionsubstrates.get(j);
                Species substratespecies = sf.getSpeciesInstance();
                Compound c = new Compound(substratespecies);
                mReaction.addSubstrate(c, (int) sf.getStoichiometry());
            }

            for (int k = 0; k <reactionproducts.size(); k++){
                SpeciesReference sf = (SpeciesReference) reactionproducts.get(k);
                Species productspecies = sf.getSpeciesInstance();
                Compound c = new Compound(productspecies);
                mReaction.addProduct(c, (int) sf.getStoichiometry());
            }

            KineticLaw KL=reaction.getKineticLaw();
            LLP=KL.getListOfLocalParameters();

            for (int count = 0; count < LLP.size(); count++){
                parameters.add(LLP.get(count));
            }

            if (modifiers.size()==0){
                Species enzymeSpecies = new Species("EnZ"+ reaction.getId());
                enzymeSpecies.setName("EnZ"+reaction.getName());
                enzymeSpecies.setInitialConcentration(1);
                e = new Enzyme(enzymeSpecies);
                addEnzymesFromReaction(e);
            }else{
                //sbo 20/206/207/597 are inhibitors
                //sbo 461/535/534/533 are activators;
                //sbo 460 is enzyme catalyst
                boolean old = true;
                for ( int l=0;l<modifiers.size();l++){
                    int sbo = modifiers.get(l).getSBOTerm();

                    if(sbo!=-1){
                        old=false;
                    }
                }
                if(old==true){
                    e = oldSBML( mReaction,  reaction,  modifiers, e);
                }else{
                    for ( int l=0;l<modifiers.size();l++){
                        int sbo = modifiers.get(l).getSBOTerm();
                        if(sbo==460){
                            ModifierSpeciesReference msr = (ModifierSpeciesReference) modifiers.get(l);
                            Species s = msr.getSpeciesInstance();
                            e = new Enzyme(s);
                        }else if (sbo==461||sbo==535||sbo==534||sbo==533){
                            regulation = 2; 
                            ModifierSpeciesReference msr = (ModifierSpeciesReference) modifiers.get(l);
                            Species s = msr.getSpeciesInstance();
                            mod.add( new Compound(s) );
                        }else if (sbo==20||sbo==206||sbo==207||sbo==597){
                            regulation = 3; 
                            ModifierSpeciesReference msr = (ModifierSpeciesReference) modifiers.get(l);
                            Species s = msr.getSpeciesInstance();
                            mod.add( new Compound(s) );
                        }
                    }
                }
            }

            if (regulation>1){
                mReaction.setEnzyme(e);
                mReaction.setModifier(mod);
                if(mod.size()==2){
                    mReaction.setRegulation(4);
                }else if (mod.size()==1){
                    mReaction.setRegulation(regulation);
                }
                ModelReaction im = new ModelReaction(e, mReaction, parameters);
                allthereactions.add(im);
            }else{
                mReaction.setEnzyme(e);
                ModelReaction im = new ModelReaction(e, mReaction);
                allthereactions.add(im);
            }
            
        }
        Object[][] reactioninfo = new Object[allthereactions.size()][5];
        for(int i =0;i< allthereactions.size();i++){
            reactioninfo[i][0]=i+1;
            reactioninfo[i][1]=allthereactions.get(i).getEnzyme().getName();
            if(allthereactions.get(i).getRegulation()>1){
                switch (allthereactions.get(i).getRegulation()) {
                    case 2:
                        reactioninfo[i][2]=allthereactions.get(i).getModifier().get(0).getName();
                        reactioninfo[i][3]="N/A";
                        break;
                    case 3:
                        reactioninfo[i][2]="N/A";
                        reactioninfo[i][3]=allthereactions.get(i).getModifier().get(0).getName();
                        break;
                    case 4:
                        reactioninfo[i][2]=allthereactions.get(i).getModifier().get(0).getName();
                        reactioninfo[i][3]=allthereactions.get(i).getModifier().get(1).getName();
                        break;
                    default:
                        break;
                }
            }else {
                reactioninfo[i][2]="N/A";
                reactioninfo[i][3]="N/A";
            }
            
            reactioninfo[i][4]=allthereactions.get(i).toString();
            
        }
        setUpReactionTable(reactioninfo);
        ///for building the reactiontable
        numberofreactions.setText(Integer.toString(allthereactions.size()));
    }
    
    private void addReactionsToReactionList_v2(ArrayList<ModelReaction> allthereactions){
        
        Object[][] reactioninfo = new Object[allthereactions.size()][5];
        for(int i =0;i< allthereactions.size();i++){
            reactioninfo[i][0]=i+1;
            reactioninfo[i][1]=allthereactions.get(i).getEnzyme().getName();
            if(allthereactions.get(i).getRegulation()>1){
                switch (allthereactions.get(i).getRegulation()) {
                    case 2:
                        reactioninfo[i][2]=allthereactions.get(i).getModifier().get(0).getName();
                        reactioninfo[i][3]="N/A";
                        break;
                    case 3:
                        reactioninfo[i][2]="N/A";
                        reactioninfo[i][3]=allthereactions.get(i).getModifier().get(0).getName();
                        break;
                    case 4:
                        reactioninfo[i][2]=allthereactions.get(i).getModifier().get(0).getName();
                        reactioninfo[i][3]=allthereactions.get(i).getModifier().get(1).getName();
                        break;
                    default:
                        break;
                }
            }else {
                reactioninfo[i][2]="N/A";
                reactioninfo[i][3]="N/A";
            }
            
            reactioninfo[i][4]=allthereactions.get(i).toString();
            
        }
        setUpReactionTable(reactioninfo);
        ///for building the reactiontable
        numberofreactions.setText(Integer.toString(allthereactions.size()));
    }
    
    private Enzyme oldSBML(MetabolicReaction mReaction, Reaction reaction, ListOf modifiers, Enzyme e){
        
        System.out.println("there is more than 1 modifier in reaction: "+ reaction.getId());    
                
        String currentreaction = mReaction.toString();
        String message = "More than 1 modifier found in this reaction. \n"+currentreaction+"\nPlease select 1 enzyme for this reaction.";

        ArrayList<String> modifierlist = new ArrayList<String>();
        ArrayList<Species> modifierlist2 = new ArrayList<Species>();
        for (int l = 0; l <modifiers.size();l++){
            ModifierSpeciesReference msr = (ModifierSpeciesReference) modifiers.get(l);
            Species s = msr.getSpeciesInstance();
            e = new Enzyme(s);
            modifierlist.add(e.getName());
            modifierlist2.add(s);
        }
        
        String message2 = "More than 1 modifier found in this reaction. \n"+currentreaction+"\nIs this modifier ("+ modifierlist.get(0)+") an enzyme?";
        String modifierId;
        
        if(modifierlist.size()==1){
            Object[] options = {"Yes", "No"}; 
            int input;
            input = (int)JOptionPane.showOptionDialog(null, message2,"Please pick an enzyme.",JOptionPane.DEFAULT_OPTION, JOptionPane.INFORMATION_MESSAGE,null,options,options[0]);
            if(input==1){
               Species enzymeSpecies = new Species("EnZ"+ reaction.getId());
               enzymeSpecies.setName("EnZ"+reaction.getName());
               enzymeSpecies.setInitialConcentration(1);
               e = new Enzyme(enzymeSpecies);
               addEnzymesFromReaction(e);
            }else{
               e = new Enzyme(modifierlist2.get(input));
           }
        }

        if (modifierlist.size()==2){

           Object[] options = {modifierlist.get(0),modifierlist.get(1), "None of the them."}; 

           int input;
           input = (int)JOptionPane.showOptionDialog(null, message,"Please pick an enzyme.",JOptionPane.DEFAULT_OPTION, JOptionPane.INFORMATION_MESSAGE,null,options,options[0]);

           if (input == 2 || input == -1){
               Species enzymeSpecies = new Species("EnZ"+ reaction.getId());
               enzymeSpecies.setName("EnZ"+reaction.getName());
               enzymeSpecies.setInitialConcentration(1);
               e = new Enzyme(enzymeSpecies);
               addEnzymesFromReaction(e);
           }else{
               e = new Enzyme(modifierlist2.get(input));
           }

        }else if (modifierlist.size()==3){
            Object[] options = {modifierlist.get(0),modifierlist.get(1),modifierlist.get(2), "None of the them."}; 
            int input;
           input = (int)JOptionPane.showOptionDialog(null, message,"Please pick an enzyme.",JOptionPane.DEFAULT_OPTION, JOptionPane.INFORMATION_MESSAGE,null,options,options[0]);

           if (input == 3 || input == -1){
               Species enzymeSpecies = new Species("EnZ"+ reaction.getId());
               enzymeSpecies.setName("EnZ"+reaction.getName());
               enzymeSpecies.setInitialConcentration(1);
               e = new Enzyme(enzymeSpecies);
               addEnzymesFromReaction(e);
           }else{
               e = new Enzyme(modifierlist2.get(input));
           }
        }else if (modifierlist.size()==4){
            Object[] options = {modifierlist.get(0),modifierlist.get(1),modifierlist.get(2),modifierlist.get(3), "None of the them."}; 
            int input;
           input = (int)JOptionPane.showOptionDialog(null, message,"Please pick an enzyme.",JOptionPane.DEFAULT_OPTION, JOptionPane.INFORMATION_MESSAGE,null,options,options[0]);

           if (input == 4 || input == -1){
               Species enzymeSpecies = new Species("EnZ"+ reaction.getId());
               enzymeSpecies.setName("EnZ"+reaction.getName());
               enzymeSpecies.setInitialConcentration(1);
               e = new Enzyme(enzymeSpecies);
               addEnzymesFromReaction(e);
           }else{
               e = new Enzyme(modifierlist2.get(input));
           }
        }else if (modifierlist.size()==5){
            Object[] options = {modifierlist.get(0),modifierlist.get(1),modifierlist.get(2),modifierlist.get(3),modifierlist.get(4), "None of the them."}; 
            int input;
           input = (int)JOptionPane.showOptionDialog(null, message,"Please pick an enzyme.",JOptionPane.DEFAULT_OPTION, JOptionPane.INFORMATION_MESSAGE,null,options,options[0]);

           if (input == 5 || input == -1){
               Species enzymeSpecies = new Species("EnZ"+ reaction.getId());
               enzymeSpecies.setName("EnZ"+reaction.getName());
               enzymeSpecies.setInitialConcentration(1);
               e = new Enzyme(enzymeSpecies);
               addEnzymesFromReaction(e);
           }else{
               e = new Enzyme(modifierlist2.get(input));
           }
        }
        return e;
    }
    private void addEnzymesFromReaction(Enzyme enz){
        
        if(isinEnzymeList(enz) == false ){
            int i = 0;
            alltheenzymes.add(enz);
            EnzymeBox.addItem(enz.getName());
        }
        Object[][] enzymeinfo = new Object[alltheenzymes.size()][3];
        for(int i =0; i<alltheenzymes.size(); i++){
            enzymeinfo[i][0]=i+1;
            enzymeinfo[i][1]=alltheenzymes.get(i).getName();
            enzymeinfo[i][2]=alltheenzymes.get(i).getConcentration();
        }
        setUpEnzymeTable(enzymeinfo);
        numberofenzymes.setText(Integer.toString(alltheenzymes.size()));
    }
    
    public boolean isinEnzymeList(Enzyme e){
        boolean isin = false;
        
        for (int i = 0; i <alltheenzymes.size(); i++){
            if (e.equals(alltheenzymes.get(i)) || e.getID().matches(alltheenzymes.get(i).getID())){
                isin =true;
            }
        }
        return isin;
    }
    
    private void addSpeciesToSpeciesList(ListOf<Species> listofspecies){
        for (int i = 0; i < listofspecies.size(); i++) {
            String name = listofspecies.get(i).getName();
            String ID = listofspecies.get(i).getId();
            int SBO = listofspecies.get(i).getSBOTerm();
            Compound newCompound = new Compound(listofspecies.get(i));            
            String sub = ID.substring(0,1);
            String sub2="";
            if(ID.length()>1){
                sub2 = ID.substring(1,2);
            }
            if(SBO==252 || (sub.equals("E") && sub2.matches("\\d"))){
                EnzymeBox.addItem(name);
                alltheenzymes.add(new Enzyme(listofspecies.get(i)));
            }else if (SBO ==250){
            }else{
                SpeciesBox.addItem(name);
                allthespecies.add(newCompound);
            }
        }
        Object[][] speciesinfo = new Object[allthespecies.size()][4];
        Object[][] enzymeinfo = new Object[alltheenzymes.size()][3];
        
        for(int i =0; i<allthespecies.size(); i++){
            speciesinfo[i][0]=i+1;
            speciesinfo[i][1]=allthespecies.get(i).getName();
            speciesinfo[i][2]=allthespecies.get(i).getConcentration();
            speciesinfo[i][3]=allthespecies.get(i).getBoundaryCondition();
        }
        setUpSpeciesTable(speciesinfo);
        for(int i =0; i<alltheenzymes.size(); i++){
            enzymeinfo[i][0]=i+1;
            enzymeinfo[i][1]=alltheenzymes.get(i).getName();
            enzymeinfo[i][2]=alltheenzymes.get(i).getConcentration();
        }
        setUpEnzymeTable(enzymeinfo);
        numberofspecies.setText(Integer.toString(allthespecies.size()));
        numberofenzymes.setText(Integer.toString(alltheenzymes.size()));
    }
    
    private void addSpeciesToSpeciesList_v2(ArrayList<Compound> allthespecies, ArrayList<Enzyme> alltheenzymes){
        
        Object[][] speciesinfo = new Object[allthespecies.size()][4];
        Object[][] enzymeinfo = new Object[alltheenzymes.size()][3];
        
        for(int i =0; i<allthespecies.size(); i++){
            speciesinfo[i][0]=i+1;
            speciesinfo[i][1]=allthespecies.get(i).getName();
            speciesinfo[i][2]=allthespecies.get(i).getConcentration();
            speciesinfo[i][3]=allthespecies.get(i).getBoundaryCondition();
        }
        setUpSpeciesTable(speciesinfo);
        for(int i =0; i<alltheenzymes.size(); i++){
            enzymeinfo[i][0]=i+1;
            enzymeinfo[i][1]=alltheenzymes.get(i).getName();
            enzymeinfo[i][2]=alltheenzymes.get(i).getConcentration();
        }
        setUpEnzymeTable(enzymeinfo);
        numberofspecies.setText(Integer.toString(allthespecies.size()));
        numberofenzymes.setText(Integer.toString(alltheenzymes.size()));
    }
    private void setUpSpeciesTable(Object[][] Info){
        
        SpeciesTable.setModel(new javax.swing.table.DefaultTableModel(Info,SpeciesTableHeader)
    {
        Class[] types = new Class [] {
            java.lang.Integer.class, java.lang.String.class, java.lang.Double.class, java.lang.Boolean.class
        };
        boolean[] canEdit = new boolean [] {
            false, true, true, true
        };

        public Class getColumnClass(int columnIndex) {
            return types [columnIndex];
        }

        public boolean isCellEditable(int rowIndex, int columnIndex) {
            return canEdit [columnIndex];
        }
            });

        SpeciesTable.getTableHeader().setReorderingAllowed(false);

        SpeciesScrollPane.setViewportView(SpeciesTable);
        if (SpeciesTable.getColumnModel().getColumnCount() > 0) {
                    SpeciesTable.getColumnModel().getColumn(0).setMinWidth(28);
                    SpeciesTable.getColumnModel().getColumn(0).setPreferredWidth(28);
                    SpeciesTable.getColumnModel().getColumn(0).setMaxWidth(28);
                    SpeciesTable.getColumnModel().getColumn(3).setMinWidth(65);            
                    SpeciesTable.getColumnModel().getColumn(3).setPreferredWidth(65);
                    SpeciesTable.getColumnModel().getColumn(3).setMaxWidth(65);
                }
        SpeciesTable.getModel().addTableModelListener(SpeciesTableChanged);
    }
    
    private void setUpEnzymeTable(Object[][] Info){
        EnzymesTable.setModel(new javax.swing.table.DefaultTableModel(
            Info,
            EnzymeTableHeader
                ) {
            Class[] types = new Class [] {
                java.lang.Integer.class, java.lang.String.class, java.lang.Double.class
            };
            boolean[] canEdit = new boolean [] {
                false, true, true
            };

            public Class getColumnClass(int columnIndex) {
                return types [columnIndex];
            }

            public boolean isCellEditable(int rowIndex, int columnIndex) {
                return canEdit [columnIndex];
            }
                });
                EnzymesTable.getColumnModel().getColumn(0).setMinWidth(28);
                    EnzymesTable.getColumnModel().getColumn(0).setPreferredWidth(28);
                    EnzymesTable.getColumnModel().getColumn(0).setMaxWidth(28);


                EnzymeScrollPane.setViewportView(EnzymesTable);
                EnzymesTable.getModel().addTableModelListener(EnzymeTableChanged);
    }
    
    private void setUpReactionTable(Object[][] Info){
        ReactionTable.setModel(new javax.swing.table.DefaultTableModel(
            Info,
            ReactionTableHeader
        ) {
            Class[] types = new Class [] {
                java.lang.Integer.class, java.lang.Object.class, java.lang.Object.class, java.lang.Object.class, java.lang.Object.class
            };
            boolean[] canEdit = new boolean [] {
                false, true, true, true, false
            };

            public Class getColumnClass(int columnIndex) {
                return types [columnIndex];
            }

            public boolean isCellEditable(int rowIndex, int columnIndex) {
                return canEdit [columnIndex];
            }
            });

            ReactionTable.getColumnModel().getColumn(0).setMinWidth(28);
            ReactionTable.getColumnModel().getColumn(0).setPreferredWidth(28);
            ReactionTable.getColumnModel().getColumn(0).setMaxWidth(28);
            ReactionTable.getColumnModel().getColumn(2).setMaxWidth(115);
            ReactionTable.getColumnModel().getColumn(1).setMaxWidth(115);
            ReactionTable.getColumnModel().getColumn(3).setMaxWidth(115);
            ReactionTable.getColumnModel().getColumn(2).setMinWidth(115);
            ReactionTable.getColumnModel().getColumn(1).setMinWidth(115);
            ReactionTable.getColumnModel().getColumn(3).setMinWidth(115);
            ReactionScrollPane.setViewportView(ReactionTable);
            ReactionTable.getModel().addTableModelListener(ReactionTableChanged);
    }
    private void setUpEnzymeColumn(JTable table, TableColumn enzymeColumn) {
        enzymeColumn.setCellEditor(new DefaultCellEditor(EnzymeBox));

        //Set up tool tips for the sport cells.
        DefaultTableCellRenderer renderer =
                new DefaultTableCellRenderer();
        renderer.setToolTipText("Click for combo box");
        enzymeColumn.setCellRenderer(renderer);
    }
    private void setUpActivatorColumn(JTable table, TableColumn enzymeColumn) {
        enzymeColumn.setCellEditor(new DefaultCellEditor(SpeciesBox));

        //Set up tool tips for the sport cells.
        DefaultTableCellRenderer renderer =
                new DefaultTableCellRenderer();
        renderer.setToolTipText("Click for combo box");
        enzymeColumn.setCellRenderer(renderer);
    }
    private void setUpInhibitorColumn(JTable table, TableColumn enzymeColumn) {
        enzymeColumn.setCellEditor(new DefaultCellEditor(SpeciesBox));

        //Set up tool tips for the sport cells.
        DefaultTableCellRenderer renderer =
                new DefaultTableCellRenderer();
        renderer.setToolTipText("Click for combo box");
        enzymeColumn.setCellRenderer(renderer);
    }
    
    private boolean check_for_repeat(String sps_name){
        boolean species_exist=false;
        
        for (int i =0; i<allthespecies.size();i++){
            String name = allthespecies.get(i).getName();
                if(name.equals(sps_name)){
                    species_exist=true;
                }
            }
        
        return species_exist;
    }
    
    private TableModelListener SpeciesTableChanged = new TableModelListener() {

        public void tableChanged(TableModelEvent e) {
            
            if(delete_metabolite){
                delete_metabolite=false;
            }else{
            
            int row = e.getFirstRow();
            int column = e.getColumn();
            
            if(!SpeciesTable.getValueAt(row,1).equals("")){
            menuSaveModel.setEnabled(true);  
            saveModelButton.setEnabled(true);
            String newname = (String) SpeciesTable.getValueAt(row,1);
            if(column==1){
                int marker=-1;
                double conc = 0.1;
                if(SpeciesTable.getValueAt(row,2)!=null){
                    conc= (double)SpeciesTable.getValueAt(row,2);
                }
                for (int i =0; i<allthespecies.size();i++){
                    String name = allthespecies.get(i).getName();
                    if(name.equals(tempname)){
                        marker=i;
                    }
                }
                if(marker==-1){
                    if(check_for_repeat(newname)){
                        
                    }else{
                        Compound newspecies = new Compound(newname, conc);
                        newspecies.setBoundaryCondition((boolean) SpeciesTable.getValueAt(row,3));
                        allthespecies.add(newspecies);
                        SpeciesBox.addItem(newname);
                    }
                }else{
                    allthespecies.get(marker).setName(newname);
                    SpeciesBox.removeAllItems();
                    SpeciesBox.addItem("N/A");
                    for(Compound eh : allthespecies){
                        String meh =eh.getName();
                        SpeciesBox.addItem(meh);
                    }
                }
            }
            numberofspecies.setText(Integer.toString(allthespecies.size()));
            
            if(column==2){
                double conc = (double) SpeciesTable.getValueAt(row,column);
                if(newname!=null){
                    for (int i =0; i<allthespecies.size();i++){
                        String name = allthespecies.get(i).getName();
                        if(name.equals(newname)){
                            allthespecies.get(i).setConcentration(conc);
                        }
                        
                    }
                }
            }
            
            if(column==3){
                boolean yesno = (boolean) SpeciesTable.getValueAt(row,column);
                
                if(newname!=null){
                    for (int i =0; i<allthespecies.size();i++){
                        String name = allthespecies.get(i).getName();
                        if(name.equals(newname)){
                            allthespecies.get(i).setBoundaryCondition(yesno);
                        }
                        
                    }
                }
            }
            }
            }
        }
    };
    
    private TableModelListener EnzymeTableChanged = new TableModelListener() {

        public void tableChanged(TableModelEvent e) {
            
            if(delete_enzyme){
                delete_enzyme=false;
            }else{
            
            int row = e.getFirstRow();
            int column = e.getColumn();
            
            if(!EnzymesTable.getValueAt(row,1).equals("")){
            String newname = (String) EnzymesTable.getValueAt(row,1);
            if(column==1){
                int marker=-1;
                
                double conc = 1;
                if(EnzymesTable.getValueAt(row,2)!=null){
                    conc= (double)EnzymesTable.getValueAt(row,2);
                }
                for (int i =0; i<alltheenzymes.size();i++){
                    String name = alltheenzymes.get(i).getName();
                    if(name.equals(tempname)){
                        marker=i;
                    }
                }
                if(marker==-1){
                    Enzyme newspecies = new Enzyme(newname, conc);
                    alltheenzymes.add(newspecies);
                    EnzymeBox.addItem(newname);
                }else{
                    alltheenzymes.get(marker).setName(newname);
                    EnzymeBox.removeAllItems();
                    for(Enzyme eh : alltheenzymes){
                        String meh = eh.getName();
                        EnzymeBox.addItem(meh);
                    }
                }
            }
            numberofenzymes.setText(Integer.toString(alltheenzymes.size()));
            
            if(column==2){
                double conc = (double) EnzymesTable.getValueAt(row,column);
                if(newname!=null){
                    for (int i =0; i<alltheenzymes.size();i++){
                        String name = alltheenzymes.get(i).getName();
                        if(name.equals(newname)){
                            alltheenzymes.get(i).setConcentration(conc);
                        }
                    }
                }
            }
            }
            }
        }
    };
    
    private TableModelListener ReactionTableChanged = new TableModelListener() {

        public void tableChanged(TableModelEvent e) {
            
            if(delete_reaction){
                delete_reaction=false;
            }else{
            
            int row = e.getFirstRow();
            int column = e.getColumn();
            String reaction = (String) ReactionTable.getValueAt(row,4);
            if(reaction!=null){
                int editreaction=-1;
                String[] tokens = reaction.split("\\W");
                for (int i =0; i <allthereactions.size(); i++){
                    String id = allthereactions.get(i).getReactionID();
                    if(id.matches(tokens[0])){
                        editreaction=i;
                    }
                }
                if(column==1){
                    String enzyme = (String) ReactionTable.getValueAt(row,1);
                    Enzyme newenzyme=null;
                    for(int i =0; i<alltheenzymes.size();i++){
                        if(enzyme.equals(alltheenzymes.get(i).getName())){
                            newenzyme=alltheenzymes.get(i);
                        }
                    }
                    allthereactions.get(editreaction).setEnzyme(newenzyme);
                }
                                
                if(column==2||column==3){
                String activator = (String) ReactionTable.getValueAt(row, 2);
                String inhibitor = (String) ReactionTable.getValueAt(row, 3);
                int regulation;
                ArrayList<Compound> modifier= new ArrayList();
                if(!activator.equals("N/A")&&!inhibitor.equals("N/A")){
                    regulation=4;
                    if(activator.equals(inhibitor)){
                        JOptionPane.showMessageDialog( null, "Having the same activator and inhibitor is not advised.", "Not a good idea", JOptionPane.ERROR_MESSAGE);
                    }
                    for ( int i =0; i<allthespecies.size();i++){
                        if(activator.equals(allthespecies.get(i).getName())){
                            modifier.add(allthespecies.get(i));
                        }
                    }
                    for ( int i =0; i<allthespecies.size();i++){
                        if(inhibitor.equals(allthespecies.get(i).getName())){
                            modifier.add(allthespecies.get(i));
                        }
                    }
                } else{
                    if(activator.equals("N/A")){
                        if(inhibitor.equals("N/A")){
                            regulation=1;
                        }else{
                            regulation=3;
                            for ( int i =0; i<allthespecies.size();i++){
                                if(inhibitor.equals(allthespecies.get(i).getName())){
                                    modifier.add(allthespecies.get(i));
                                }
                            }
                        }
                    }else {
                        regulation=2;
                        for ( int i =0; i<allthespecies.size();i++){
                            if(activator.equals(allthespecies.get(i).getName())){
                                modifier.add(allthespecies.get(i));
                            }
                        }
                    }
                }
                    allthereactions.get(editreaction).setRegulation(regulation);
                    allthereactions.get(editreaction).setModifier(modifier);
                }
            }
            
        }
        }
    };
    
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        ReactionBuilder = new javax.swing.JPanel();
        jPanel1 = new javax.swing.JPanel();
        prod1spin = new javax.swing.JSpinner();
        prod1box = new javax.swing.JComboBox();
        prod2box = new javax.swing.JComboBox();
        prod2spin = new javax.swing.JSpinner();
        prod3box = new javax.swing.JComboBox();
        prod3spin = new javax.swing.JSpinner();
        jLabel2 = new javax.swing.JLabel();
        jLabel4 = new javax.swing.JLabel();
        jPanel2 = new javax.swing.JPanel();
        sub1spin = new javax.swing.JSpinner();
        sub1box = new javax.swing.JComboBox();
        sub2box = new javax.swing.JComboBox();
        sub2spin = new javax.swing.JSpinner();
        sub3box = new javax.swing.JComboBox();
        sub3spin = new javax.swing.JSpinner();
        jLabel1 = new javax.swing.JLabel();
        jLabel3 = new javax.swing.JLabel();
        GAparameters = new javax.swing.JPanel();
        PopulationLabel = new javax.swing.JLabel();
        GenerationLabel = new javax.swing.JLabel();
        PlateauLabel = new javax.swing.JLabel();
        PlagueLabel = new javax.swing.JLabel();
        populationText = new javax.swing.JTextField();
        generationText = new javax.swing.JTextField();
        plateauText = new javax.swing.JTextField();
        plagueText = new javax.swing.JTextField();
        PlateauLabel1 = new javax.swing.JLabel();
        PlateauLabel2 = new javax.swing.JLabel();
        PlateauLabel3 = new javax.swing.JLabel();
        PlateauLabel4 = new javax.swing.JLabel();
        PlateauLabel5 = new javax.swing.JLabel();
        number_of_core = new javax.swing.JTextField();
        vf_low = new javax.swing.JTextField();
        PlateauLabel6 = new javax.swing.JLabel();
        vf_upp = new javax.swing.JTextField();
        PlateauLabel7 = new javax.swing.JLabel();
        PlateauLabel8 = new javax.swing.JLabel();
        vr_low = new javax.swing.JTextField();
        PlateauLabel9 = new javax.swing.JLabel();
        vr_upp = new javax.swing.JTextField();
        PlateauLabel10 = new javax.swing.JLabel();
        ka_low = new javax.swing.JTextField();
        PlateauLabel11 = new javax.swing.JLabel();
        ka_upp = new javax.swing.JTextField();
        PlateauLabel12 = new javax.swing.JLabel();
        ki_low = new javax.swing.JTextField();
        PlateauLabel13 = new javax.swing.JLabel();
        ki_upp = new javax.swing.JTextField();
        PlateauLabel14 = new javax.swing.JLabel();
        PlateauLabel15 = new javax.swing.JLabel();
        km_low = new javax.swing.JTextField();
        PlateauLabel16 = new javax.swing.JLabel();
        km_upp = new javax.swing.JTextField();
        missingSpeciesPanel = new javax.swing.JPanel();
        jPanel3 = new javax.swing.JPanel();
        jLabel7 = new javax.swing.JLabel();
        jLabel6 = new javax.swing.JLabel();
        jLabel5 = new javax.swing.JLabel();
        jScrollPane3 = new javax.swing.JScrollPane();
        missingSpeciesText = new javax.swing.JTextArea();
        jScrollPane4 = new javax.swing.JScrollPane();
        missingEnzymeText = new javax.swing.JTextArea();
        RunParameterEstimation = new javax.swing.JButton();
        saveModelButton = new javax.swing.JButton();
        SpeciesPanel = new javax.swing.JPanel();
        SpeciesLabel = new javax.swing.JLabel();
        SpeciesScrollPane = new javax.swing.JScrollPane();
        SpeciesTable = new javax.swing.JTable();
        SpeciesCounter = new javax.swing.JLabel();
        AddSpecies = new javax.swing.JButton();
        DeleteSpecies = new javax.swing.JButton();
        numberofspecies = new javax.swing.JLabel();
        EnzymePanel = new javax.swing.JPanel();
        EnzymeLabel = new javax.swing.JLabel();
        EnzymeScrollPane = new javax.swing.JScrollPane();
        EnzymesTable = new javax.swing.JTable();
        EnzymeCounter = new javax.swing.JLabel();
        AddEnzyme = new javax.swing.JButton();
        DeleteEnzyme = new javax.swing.JButton();
        numberofenzymes = new javax.swing.JLabel();
        ReactionPanel = new javax.swing.JPanel();
        ReactionLabel = new javax.swing.JLabel();
        ReactionCounter = new javax.swing.JLabel();
        AddReaction = new javax.swing.JButton();
        DeleteReaction = new javax.swing.JButton();
        jSeparator1 = new javax.swing.JSeparator();
        numberofreactions = new javax.swing.JLabel();
        ReactionScrollPane = new javax.swing.JScrollPane();
        ReactionTable = new javax.swing.JTable();
        jLabel8 = new javax.swing.JLabel();
        generateDataHeader = new javax.swing.JButton();
        ModelNameLabel = new javax.swing.JLabel();
        ModelNameTextField = new HintTextField("Insert Model Name");
        InsertFittingData = new javax.swing.JButton();
        GRaPeLogo = new javax.swing.JLabel();
        jLabel9 = new javax.swing.JLabel();
        MenuBar = new javax.swing.JMenuBar();
        jMenu1 = new javax.swing.JMenu();
        menuNewModel = new javax.swing.JMenuItem();
        tsv_to_model = new javax.swing.JMenuItem();
        menuSaveModel = new javax.swing.JMenuItem();
        menuOpenModel = new javax.swing.JMenuItem();
        menuSeparator = new javax.swing.JPopupMenu.Separator();
        menuExit = new javax.swing.JMenuItem();

        jPanel1.setBorder(javax.swing.BorderFactory.createEtchedBorder());

        prod1spin.setValue(1);

        prod2spin.setValue(1);

        prod3spin.setRequestFocusEnabled(false);
        prod3spin.setValue(1);

        jLabel2.setText("Product(s)");

        jLabel4.setText("Stoichiometry");

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(jLabel2)
                    .addComponent(prod1box, 0, 80, Short.MAX_VALUE)
                    .addComponent(prod2box, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(prod3box, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addGap(18, 18, 18)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel4)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addGap(51, 51, 51)
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(prod1spin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(prod3spin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(prod2spin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))))
                .addContainerGap(23, Short.MAX_VALUE))
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap(14, Short.MAX_VALUE)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel2)
                    .addComponent(jLabel4))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(prod1box, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(prod1spin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(prod2box, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(prod2spin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(prod3box, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(prod3spin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(14, Short.MAX_VALUE))
        );

        jPanel2.setBorder(javax.swing.BorderFactory.createEtchedBorder());

        sub1spin.setValue(1);

        sub2spin.setValue(1);

        sub3spin.setValue(1);

        jLabel1.setText("Substrate(s)");

        jLabel3.setText("Stoichiometry");

        javax.swing.GroupLayout jPanel2Layout = new javax.swing.GroupLayout(jPanel2);
        jPanel2.setLayout(jPanel2Layout);
        jPanel2Layout.setHorizontalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel2Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(jLabel1)
                    .addComponent(sub1box, 0, 80, Short.MAX_VALUE)
                    .addComponent(sub2box, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(sub3box, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addGap(18, 18, 18)
                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel3)
                    .addGroup(jPanel2Layout.createSequentialGroup()
                        .addGap(51, 51, 51)
                        .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(sub1spin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(sub3spin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(sub2spin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))))
                .addContainerGap(16, Short.MAX_VALUE))
        );
        jPanel2Layout.setVerticalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel2Layout.createSequentialGroup()
                .addContainerGap(14, Short.MAX_VALUE)
                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel1)
                    .addComponent(jLabel3))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(sub1box, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(sub1spin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(sub2box, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(sub2spin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(sub3box, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(sub3spin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(14, Short.MAX_VALUE))
        );

        javax.swing.GroupLayout ReactionBuilderLayout = new javax.swing.GroupLayout(ReactionBuilder);
        ReactionBuilder.setLayout(ReactionBuilderLayout);
        ReactionBuilderLayout.setHorizontalGroup(
            ReactionBuilderLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, ReactionBuilderLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jPanel2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(jPanel1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(12, 12, 12))
        );
        ReactionBuilderLayout.setVerticalGroup(
            ReactionBuilderLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(ReactionBuilderLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(ReactionBuilderLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.CENTER)
                    .addComponent(jPanel2, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jPanel1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        PopulationLabel.setText("Starting Population: ");

        GenerationLabel.setText("Max Generations:");

        PlateauLabel.setText("Plateau Limit:");

        PlagueLabel.setText("Plague at every N Generation: ");

        populationText.setText("50");

        generationText.setText("500");

        plateauText.setText("25");

        plagueText.setText("5");

        PlateauLabel1.setText("Number of CPU");
        PlateauLabel1.setMaximumSize(new java.awt.Dimension(175, 14));
        PlateauLabel1.setMinimumSize(new java.awt.Dimension(175, 14));
        PlateauLabel1.setPreferredSize(new java.awt.Dimension(175, 14));

        PlateauLabel2.setText("Vf parameter range");
        PlateauLabel2.setMaximumSize(new java.awt.Dimension(175, 14));
        PlateauLabel2.setMinimumSize(new java.awt.Dimension(175, 14));
        PlateauLabel2.setPreferredSize(new java.awt.Dimension(175, 14));

        PlateauLabel3.setText("Ka parameter range");
        PlateauLabel3.setMaximumSize(new java.awt.Dimension(175, 14));
        PlateauLabel3.setMinimumSize(new java.awt.Dimension(175, 14));
        PlateauLabel3.setPreferredSize(new java.awt.Dimension(175, 14));

        PlateauLabel4.setText("Vr parameter range");
        PlateauLabel4.setMaximumSize(new java.awt.Dimension(175, 14));
        PlateauLabel4.setMinimumSize(new java.awt.Dimension(175, 14));
        PlateauLabel4.setPreferredSize(new java.awt.Dimension(175, 14));

        PlateauLabel5.setText("Ki parameter range");
        PlateauLabel5.setMaximumSize(new java.awt.Dimension(175, 14));
        PlateauLabel5.setMinimumSize(new java.awt.Dimension(175, 14));
        PlateauLabel5.setPreferredSize(new java.awt.Dimension(175, 14));

        number_of_core.setText("3");

        vf_low.setText("3");

        PlateauLabel6.setText("10^");
        PlateauLabel6.setMaximumSize(new java.awt.Dimension(175, 14));
        PlateauLabel6.setMinimumSize(new java.awt.Dimension(175, 14));
        PlateauLabel6.setPreferredSize(new java.awt.Dimension(175, 14));

        vf_upp.setText("5");

        PlateauLabel7.setText("10^");
        PlateauLabel7.setMaximumSize(new java.awt.Dimension(175, 14));
        PlateauLabel7.setMinimumSize(new java.awt.Dimension(175, 14));
        PlateauLabel7.setPreferredSize(new java.awt.Dimension(175, 14));

        PlateauLabel8.setText("10^");
        PlateauLabel8.setMaximumSize(new java.awt.Dimension(175, 14));
        PlateauLabel8.setMinimumSize(new java.awt.Dimension(175, 14));
        PlateauLabel8.setPreferredSize(new java.awt.Dimension(175, 14));

        vr_low.setText("-3");

        PlateauLabel9.setText("10^");
        PlateauLabel9.setMaximumSize(new java.awt.Dimension(175, 14));
        PlateauLabel9.setMinimumSize(new java.awt.Dimension(175, 14));
        PlateauLabel9.setPreferredSize(new java.awt.Dimension(175, 14));

        vr_upp.setText("1");

        PlateauLabel10.setText("10^");
        PlateauLabel10.setMaximumSize(new java.awt.Dimension(175, 14));
        PlateauLabel10.setMinimumSize(new java.awt.Dimension(175, 14));
        PlateauLabel10.setPreferredSize(new java.awt.Dimension(175, 14));

        ka_low.setText("0.5");

        PlateauLabel11.setText("10^");
        PlateauLabel11.setMaximumSize(new java.awt.Dimension(175, 14));
        PlateauLabel11.setMinimumSize(new java.awt.Dimension(175, 14));
        PlateauLabel11.setPreferredSize(new java.awt.Dimension(175, 14));

        ka_upp.setText("3");

        PlateauLabel12.setText("10^");
        PlateauLabel12.setMaximumSize(new java.awt.Dimension(175, 14));
        PlateauLabel12.setMinimumSize(new java.awt.Dimension(175, 14));
        PlateauLabel12.setPreferredSize(new java.awt.Dimension(175, 14));

        ki_low.setText("-3");

        PlateauLabel13.setText("10^");
        PlateauLabel13.setMaximumSize(new java.awt.Dimension(175, 14));
        PlateauLabel13.setMinimumSize(new java.awt.Dimension(175, 14));
        PlateauLabel13.setPreferredSize(new java.awt.Dimension(175, 14));

        ki_upp.setText("2");

        PlateauLabel14.setText("Km parameter range");
        PlateauLabel14.setMaximumSize(new java.awt.Dimension(175, 14));
        PlateauLabel14.setMinimumSize(new java.awt.Dimension(175, 14));
        PlateauLabel14.setPreferredSize(new java.awt.Dimension(175, 14));

        PlateauLabel15.setText("10^");
        PlateauLabel15.setMaximumSize(new java.awt.Dimension(175, 14));
        PlateauLabel15.setMinimumSize(new java.awt.Dimension(175, 14));
        PlateauLabel15.setPreferredSize(new java.awt.Dimension(175, 14));

        km_low.setText("-3");

        PlateauLabel16.setText("10^");
        PlateauLabel16.setMaximumSize(new java.awt.Dimension(175, 14));
        PlateauLabel16.setMinimumSize(new java.awt.Dimension(175, 14));
        PlateauLabel16.setPreferredSize(new java.awt.Dimension(175, 14));

        km_upp.setText("2");

        javax.swing.GroupLayout GAparametersLayout = new javax.swing.GroupLayout(GAparameters);
        GAparameters.setLayout(GAparametersLayout);
        GAparametersLayout.setHorizontalGroup(
            GAparametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(GAparametersLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(GAparametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(GAparametersLayout.createSequentialGroup()
                        .addGroup(GAparametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                            .addComponent(PlateauLabel1, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addGroup(GAparametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addGroup(GAparametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(PlagueLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, GAparametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                        .addComponent(GenerationLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .addComponent(PopulationLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                                .addComponent(PlateauLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(GAparametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(populationText)
                            .addComponent(generationText)
                            .addComponent(plateauText)
                            .addComponent(plagueText)
                            .addComponent(number_of_core, javax.swing.GroupLayout.Alignment.TRAILING)))
                    .addGroup(GAparametersLayout.createSequentialGroup()
                        .addComponent(PlateauLabel4, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(PlateauLabel8, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(vr_low, javax.swing.GroupLayout.PREFERRED_SIZE, 41, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(PlateauLabel9, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(vr_upp, javax.swing.GroupLayout.PREFERRED_SIZE, 41, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(GAparametersLayout.createSequentialGroup()
                        .addComponent(PlateauLabel3, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(PlateauLabel10, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(ka_low, javax.swing.GroupLayout.PREFERRED_SIZE, 41, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(PlateauLabel11, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(ka_upp, javax.swing.GroupLayout.PREFERRED_SIZE, 41, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(GAparametersLayout.createSequentialGroup()
                        .addComponent(PlateauLabel5, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(PlateauLabel12, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(ki_low, javax.swing.GroupLayout.PREFERRED_SIZE, 41, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(PlateauLabel13, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(ki_upp, javax.swing.GroupLayout.PREFERRED_SIZE, 41, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(GAparametersLayout.createSequentialGroup()
                        .addComponent(PlateauLabel2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(12, 12, 12)
                        .addComponent(PlateauLabel6, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(vf_low, javax.swing.GroupLayout.PREFERRED_SIZE, 41, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(PlateauLabel7, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(vf_upp, javax.swing.GroupLayout.PREFERRED_SIZE, 41, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, GAparametersLayout.createSequentialGroup()
                        .addComponent(PlateauLabel14, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(PlateauLabel15, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(km_low, javax.swing.GroupLayout.PREFERRED_SIZE, 41, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(PlateauLabel16, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(km_upp, javax.swing.GroupLayout.PREFERRED_SIZE, 41, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap())
        );

        GAparametersLayout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {GenerationLabel, PlagueLabel, PlateauLabel, PopulationLabel});

        GAparametersLayout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {PlateauLabel1, PlateauLabel2, PlateauLabel3, PlateauLabel4, PlateauLabel5});

        GAparametersLayout.setVerticalGroup(
            GAparametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(GAparametersLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(GAparametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(PopulationLabel)
                    .addComponent(populationText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(GAparametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(GenerationLabel)
                    .addComponent(generationText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(GAparametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(plateauText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(PlateauLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(GAparametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(PlagueLabel)
                    .addComponent(plagueText, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(GAparametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(PlateauLabel1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(number_of_core, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(GAparametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(PlateauLabel2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(vf_low, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(PlateauLabel6, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(vf_upp, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(PlateauLabel7, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(GAparametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(PlateauLabel4, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(vr_low, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(PlateauLabel8, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(vr_upp, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(PlateauLabel9, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(GAparametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(PlateauLabel3, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(ka_low, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(PlateauLabel10, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(ka_upp, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(PlateauLabel11, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(GAparametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(PlateauLabel5, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(ki_low, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(PlateauLabel12, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(ki_upp, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(PlateauLabel13, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(GAparametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(PlateauLabel14, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(km_low, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(PlateauLabel15, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(km_upp, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(PlateauLabel16, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap())
        );

        GAparametersLayout.linkSize(javax.swing.SwingConstants.VERTICAL, new java.awt.Component[] {generationText, plagueText, plateauText, populationText});

        jLabel7.setText("Specie(s) and/or Enzyme(s) found in reactions but");

        jLabel6.setText("not created in model.");

        jLabel5.setText("Unable to Export Model!");

        javax.swing.GroupLayout jPanel3Layout = new javax.swing.GroupLayout(jPanel3);
        jPanel3.setLayout(jPanel3Layout);
        jPanel3Layout.setHorizontalGroup(
            jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel3Layout.createSequentialGroup()
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel7)
                    .addComponent(jLabel6)
                    .addComponent(jLabel5))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        jPanel3Layout.setVerticalGroup(
            jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel3Layout.createSequentialGroup()
                .addComponent(jLabel5)
                .addGap(9, 9, 9)
                .addComponent(jLabel7)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jLabel6))
        );

        missingSpeciesText.setColumns(20);
        missingSpeciesText.setLineWrap(true);
        missingSpeciesText.setRows(5);
        jScrollPane3.setViewportView(missingSpeciesText);

        missingEnzymeText.setColumns(20);
        missingEnzymeText.setLineWrap(true);
        missingEnzymeText.setRows(5);
        jScrollPane4.setViewportView(missingEnzymeText);

        javax.swing.GroupLayout missingSpeciesPanelLayout = new javax.swing.GroupLayout(missingSpeciesPanel);
        missingSpeciesPanel.setLayout(missingSpeciesPanelLayout);
        missingSpeciesPanelLayout.setHorizontalGroup(
            missingSpeciesPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(missingSpeciesPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(missingSpeciesPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(missingSpeciesPanelLayout.createSequentialGroup()
                        .addComponent(jScrollPane3, javax.swing.GroupLayout.PREFERRED_SIZE, 0, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jScrollPane4, javax.swing.GroupLayout.PREFERRED_SIZE, 0, Short.MAX_VALUE)
                        .addContainerGap())
                    .addComponent(jPanel3, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
        );
        missingSpeciesPanelLayout.setVerticalGroup(
            missingSpeciesPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(missingSpeciesPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jPanel3, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(missingSpeciesPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jScrollPane3, javax.swing.GroupLayout.DEFAULT_SIZE, 189, Short.MAX_VALUE)
                    .addComponent(jScrollPane4))
                .addContainerGap())
        );

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setTitle("GRaPe 2.0");

        RunParameterEstimation.setText("Run Parameter Estimation");
        RunParameterEstimation.setEnabled(false);
        RunParameterEstimation.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                RunParameterEstimationActionPerformed(evt);
            }
        });

        saveModelButton.setText("Save");
        saveModelButton.setEnabled(false);
        saveModelButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                saveModelButtonActionPerformed(evt);
            }
        });

        SpeciesPanel.setMaximumSize(new java.awt.Dimension(32767, 250));
        SpeciesPanel.setPreferredSize(new java.awt.Dimension(335, 194));

        SpeciesLabel.setText("Reacting Species");

        SpeciesScrollPane.setToolTipText("");

        SpeciesTable.setModel(new javax.swing.table.DefaultTableModel(SpeciesTableInfo,SpeciesTableHeader)
            {
                Class[] types = new Class [] {
                    java.lang.Integer.class, java.lang.String.class, java.lang.Double.class, java.lang.Boolean.class
                };
                boolean[] canEdit = new boolean [] {
                    false, true, true, true
                };

                public Class getColumnClass(int columnIndex) {
                    return types [columnIndex];
                }

                public boolean isCellEditable(int rowIndex, int columnIndex) {
                    return canEdit [columnIndex];
                }
            });
            SpeciesTable.getTableHeader().setReorderingAllowed(false);
            SpeciesTable.addMouseListener(new java.awt.event.MouseAdapter() {
                public void mouseClicked(java.awt.event.MouseEvent evt) {
                    SpeciesTableMouseClicked(evt);
                }
            });
            SpeciesTable.addKeyListener(new java.awt.event.KeyAdapter() {
                public void keyReleased(java.awt.event.KeyEvent evt) {
                    SpeciesTableKeyReleased(evt);
                }
            });
            SpeciesScrollPane.setViewportView(SpeciesTable);
            if (SpeciesTable.getColumnModel().getColumnCount() > 0) {
                SpeciesTable.getColumnModel().getColumn(0).setMinWidth(28);
                SpeciesTable.getColumnModel().getColumn(0).setPreferredWidth(28);
                SpeciesTable.getColumnModel().getColumn(0).setMaxWidth(28);
                SpeciesTable.getColumnModel().getColumn(3).setMinWidth(65);
                SpeciesTable.getColumnModel().getColumn(3).setPreferredWidth(65);
                SpeciesTable.getColumnModel().getColumn(3).setMaxWidth(65);
            }

            SpeciesCounter.setText("Number of Species:");
            SpeciesCounter.setBorder(javax.swing.BorderFactory.createCompoundBorder());

            AddSpecies.setText("Add");
            AddSpecies.addActionListener(new java.awt.event.ActionListener() {
                public void actionPerformed(java.awt.event.ActionEvent evt) {
                    AddSpeciesActionPerformed(evt);
                }
            });

            DeleteSpecies.setText("Delete");
            DeleteSpecies.addActionListener(new java.awt.event.ActionListener() {
                public void actionPerformed(java.awt.event.ActionEvent evt) {
                    DeleteSpeciesActionPerformed(evt);
                }
            });

            numberofspecies.setText("0");

            javax.swing.GroupLayout SpeciesPanelLayout = new javax.swing.GroupLayout(SpeciesPanel);
            SpeciesPanel.setLayout(SpeciesPanelLayout);
            SpeciesPanelLayout.setHorizontalGroup(
                SpeciesPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(SpeciesPanelLayout.createSequentialGroup()
                    .addContainerGap()
                    .addGroup(SpeciesPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                        .addGroup(SpeciesPanelLayout.createSequentialGroup()
                            .addGroup(SpeciesPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addGroup(SpeciesPanelLayout.createSequentialGroup()
                                    .addComponent(SpeciesLabel)
                                    .addGap(0, 0, Short.MAX_VALUE))
                                .addGroup(SpeciesPanelLayout.createSequentialGroup()
                                    .addComponent(SpeciesCounter)
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                    .addComponent(numberofspecies, javax.swing.GroupLayout.PREFERRED_SIZE, 43, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 201, Short.MAX_VALUE)
                                    .addComponent(AddSpecies)
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                    .addComponent(DeleteSpecies)))
                            .addContainerGap())
                        .addComponent(SpeciesScrollPane, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, 0, Short.MAX_VALUE)))
            );

            SpeciesPanelLayout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {AddSpecies, DeleteSpecies});

            SpeciesPanelLayout.setVerticalGroup(
                SpeciesPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, SpeciesPanelLayout.createSequentialGroup()
                    .addComponent(SpeciesLabel)
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                    .addComponent(SpeciesScrollPane, javax.swing.GroupLayout.DEFAULT_SIZE, 213, Short.MAX_VALUE)
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                    .addGroup(SpeciesPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(AddSpecies)
                        .addComponent(DeleteSpecies)
                        .addComponent(SpeciesCounter, javax.swing.GroupLayout.PREFERRED_SIZE, 26, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(numberofspecies))
                    .addContainerGap())
            );

            SpeciesPanelLayout.linkSize(javax.swing.SwingConstants.VERTICAL, new java.awt.Component[] {AddSpecies, DeleteSpecies});

            EnzymePanel.setMaximumSize(new java.awt.Dimension(32767, 250));
            EnzymePanel.setPreferredSize(new java.awt.Dimension(335, 205));

            EnzymeLabel.setText("Enzymes");

            EnzymesTable.setModel(new javax.swing.table.DefaultTableModel(
                EnzymeTableInfo,
                EnzymeTableHeader
            ) {
                Class[] types = new Class [] {
                    java.lang.Integer.class, java.lang.String.class, java.lang.Double.class
                };
                boolean[] canEdit = new boolean [] {
                    false, true, true
                };

                public Class getColumnClass(int columnIndex) {
                    return types [columnIndex];
                }

                public boolean isCellEditable(int rowIndex, int columnIndex) {
                    return canEdit [columnIndex];
                }
            });
            EnzymesTable.getColumnModel().getColumn(0).setMinWidth(28);
            EnzymesTable.getColumnModel().getColumn(0).setPreferredWidth(28);
            EnzymesTable.getColumnModel().getColumn(0).setMaxWidth(28);
            EnzymesTable.addMouseListener(new java.awt.event.MouseAdapter() {
                public void mouseClicked(java.awt.event.MouseEvent evt) {
                    EnzymesTableMouseClicked(evt);
                }
            });
            EnzymesTable.addKeyListener(new java.awt.event.KeyAdapter() {
                public void keyReleased(java.awt.event.KeyEvent evt) {
                    EnzymesTableKeyReleased(evt);
                }
            });
            EnzymeScrollPane.setViewportView(EnzymesTable);

            EnzymeCounter.setText("Number of Enzymes:");

            AddEnzyme.setText("Add");
            AddEnzyme.addActionListener(new java.awt.event.ActionListener() {
                public void actionPerformed(java.awt.event.ActionEvent evt) {
                    AddEnzymeActionPerformed(evt);
                }
            });

            DeleteEnzyme.setText("Delete");
            DeleteEnzyme.addActionListener(new java.awt.event.ActionListener() {
                public void actionPerformed(java.awt.event.ActionEvent evt) {
                    DeleteEnzymeActionPerformed(evt);
                }
            });

            numberofenzymes.setText("0");

            javax.swing.GroupLayout EnzymePanelLayout = new javax.swing.GroupLayout(EnzymePanel);
            EnzymePanel.setLayout(EnzymePanelLayout);
            EnzymePanelLayout.setHorizontalGroup(
                EnzymePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addComponent(EnzymeScrollPane, javax.swing.GroupLayout.PREFERRED_SIZE, 0, Short.MAX_VALUE)
                .addGroup(EnzymePanelLayout.createSequentialGroup()
                    .addContainerGap()
                    .addComponent(EnzymeCounter)
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                    .addComponent(numberofenzymes, javax.swing.GroupLayout.PREFERRED_SIZE, 44, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 196, Short.MAX_VALUE)
                    .addComponent(AddEnzyme)
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                    .addComponent(DeleteEnzyme)
                    .addContainerGap())
                .addGroup(EnzymePanelLayout.createSequentialGroup()
                    .addComponent(EnzymeLabel)
                    .addGap(0, 0, Short.MAX_VALUE))
            );

            EnzymePanelLayout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {AddEnzyme, DeleteEnzyme});

            EnzymePanelLayout.setVerticalGroup(
                EnzymePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(EnzymePanelLayout.createSequentialGroup()
                    .addComponent(EnzymeLabel)
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                    .addComponent(EnzymeScrollPane, javax.swing.GroupLayout.DEFAULT_SIZE, 213, Short.MAX_VALUE)
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                    .addGroup(EnzymePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(DeleteEnzyme)
                        .addComponent(AddEnzyme)
                        .addComponent(EnzymeCounter, javax.swing.GroupLayout.PREFERRED_SIZE, 26, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(numberofenzymes))
                    .addGap(12, 12, 12))
            );

            EnzymePanelLayout.linkSize(javax.swing.SwingConstants.VERTICAL, new java.awt.Component[] {AddEnzyme, DeleteEnzyme});

            ReactionPanel.setPreferredSize(new java.awt.Dimension(690, 590));

            ReactionLabel.setText("Reactions");

            ReactionCounter.setText("Number of Reactions:");

            AddReaction.setText("Add");
            AddReaction.addActionListener(new java.awt.event.ActionListener() {
                public void actionPerformed(java.awt.event.ActionEvent evt) {
                    AddReactionActionPerformed(evt);
                }
            });

            DeleteReaction.setText("Delete");
            DeleteReaction.addActionListener(new java.awt.event.ActionListener() {
                public void actionPerformed(java.awt.event.ActionEvent evt) {
                    DeleteReactionActionPerformed(evt);
                }
            });

            numberofreactions.setText("0");

            ReactionTable.setModel(new javax.swing.table.DefaultTableModel(
                ReactionTableInfo,
                ReactionTableHeader
            ) {
                Class[] types = new Class [] {
                    java.lang.Integer.class, java.lang.Object.class, java.lang.Object.class, java.lang.Object.class, java.lang.Object.class
                };
                boolean[] canEdit = new boolean [] {
                    false, true, true, true, false
                };

                public Class getColumnClass(int columnIndex) {
                    return types [columnIndex];
                }

                public boolean isCellEditable(int rowIndex, int columnIndex) {
                    return canEdit [columnIndex];
                }
            });
            ReactionTable.addMouseListener(new java.awt.event.MouseAdapter() {
                public void mouseClicked(java.awt.event.MouseEvent evt) {
                    ReactionTableMouseClicked(evt);
                }
            });
            ReactionTable.getColumnModel().getColumn(0).setMinWidth(28);
            ReactionTable.getColumnModel().getColumn(0).setPreferredWidth(28);
            ReactionTable.getColumnModel().getColumn(0).setMaxWidth(28);
            ReactionTable.getColumnModel().getColumn(2).setMaxWidth(115);
            ReactionTable.getColumnModel().getColumn(1).setMaxWidth(115);
            ReactionTable.getColumnModel().getColumn(3).setMaxWidth(115);
            ReactionTable.getColumnModel().getColumn(2).setMinWidth(115);
            ReactionTable.getColumnModel().getColumn(1).setMinWidth(115);
            ReactionTable.getColumnModel().getColumn(3).setMinWidth(115);
            ReactionScrollPane.setViewportView(ReactionTable);

            jLabel8.setText("Model Fitting");

            javax.swing.GroupLayout ReactionPanelLayout = new javax.swing.GroupLayout(ReactionPanel);
            ReactionPanel.setLayout(ReactionPanelLayout);
            ReactionPanelLayout.setHorizontalGroup(
                ReactionPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(ReactionPanelLayout.createSequentialGroup()
                    .addContainerGap()
                    .addGroup(ReactionPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                        .addGroup(ReactionPanelLayout.createSequentialGroup()
                            .addComponent(jLabel8)
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                            .addComponent(jSeparator1))
                        .addGroup(ReactionPanelLayout.createSequentialGroup()
                            .addComponent(ReactionLabel)
                            .addGap(0, 0, Short.MAX_VALUE))
                        .addGroup(ReactionPanelLayout.createSequentialGroup()
                            .addGroup(ReactionPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addGroup(ReactionPanelLayout.createSequentialGroup()
                                    .addComponent(ReactionCounter)
                                    .addGap(29, 29, 29)
                                    .addComponent(numberofreactions, javax.swing.GroupLayout.PREFERRED_SIZE, 69, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .addComponent(AddReaction)
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                    .addComponent(DeleteReaction))
                                .addComponent(ReactionScrollPane))
                            .addContainerGap())))
            );

            ReactionPanelLayout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {AddReaction, DeleteReaction});

            ReactionPanelLayout.setVerticalGroup(
                ReactionPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, ReactionPanelLayout.createSequentialGroup()
                    .addComponent(ReactionLabel)
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                    .addComponent(ReactionScrollPane, javax.swing.GroupLayout.DEFAULT_SIZE, 336, Short.MAX_VALUE)
                    .addGap(12, 12, 12)
                    .addGroup(ReactionPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(ReactionCounter)
                        .addComponent(numberofreactions)
                        .addComponent(DeleteReaction)
                        .addComponent(AddReaction))
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                    .addGroup(ReactionPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                        .addComponent(jSeparator1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(jLabel8))
                    .addGap(0, 0, 0))
            );

            ReactionPanelLayout.linkSize(javax.swing.SwingConstants.VERTICAL, new java.awt.Component[] {AddReaction, DeleteReaction});

            generateDataHeader.setText("Generate Header File");
            generateDataHeader.setToolTipText("Generate Header File for Parameter Estimation");
            generateDataHeader.setActionCommand("generateHeaderFile");
            generateDataHeader.addActionListener(new java.awt.event.ActionListener() {
                public void actionPerformed(java.awt.event.ActionEvent evt) {
                    generateDataHeaderActionPerformed(evt);
                }
            });

            ModelNameLabel.setText("Model Name: ");

            ModelNameTextField.addActionListener(new java.awt.event.ActionListener() {
                public void actionPerformed(java.awt.event.ActionEvent evt) {
                    ModelNameTextFieldActionPerformed(evt);
                }
            });

            InsertFittingData.setText("Insert Fitting Data");
            InsertFittingData.addActionListener(new java.awt.event.ActionListener() {
                public void actionPerformed(java.awt.event.ActionEvent evt) {
                    InsertFittingDataActionPerformed(evt);
                }
            });

            GRaPeLogo.setIcon(new javax.swing.ImageIcon(getClass().getResource("/GUI/grapeLogo2.png"))); // NOI18N

            jLabel9.setIcon(new javax.swing.ImageIcon(getClass().getResource("/GUI/sbml.png"))); // NOI18N

            jMenu1.setText("File");

            menuNewModel.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_N, java.awt.event.InputEvent.CTRL_MASK));
            menuNewModel.setText("New Model");
            menuNewModel.addActionListener(new java.awt.event.ActionListener() {
                public void actionPerformed(java.awt.event.ActionEvent evt) {
                    menuNewModelActionPerformed(evt);
                }
            });
            jMenu1.add(menuNewModel);

            tsv_to_model.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_I, java.awt.event.InputEvent.CTRL_MASK));
            tsv_to_model.setText("Tsv to Model");
            tsv_to_model.addActionListener(new java.awt.event.ActionListener() {
                public void actionPerformed(java.awt.event.ActionEvent evt) {
                    tsv_to_modelActionPerformed(evt);
                }
            });
            jMenu1.add(tsv_to_model);

            menuSaveModel.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_S, java.awt.event.InputEvent.CTRL_MASK));
            menuSaveModel.setIcon(new javax.swing.ImageIcon(getClass().getResource("/GUI/save.png"))); // NOI18N
            menuSaveModel.setText("Export Model");
            menuSaveModel.setEnabled(false);
            menuSaveModel.addActionListener(new java.awt.event.ActionListener() {
                public void actionPerformed(java.awt.event.ActionEvent evt) {
                    menuSaveModelsaveModelAsSBMLFile(evt);
                }
            });
            jMenu1.add(menuSaveModel);

            menuOpenModel.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_O, java.awt.event.InputEvent.CTRL_MASK));
            menuOpenModel.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_O, java.awt.event.InputEvent.CTRL_MASK));
            menuOpenModel.setIcon(new javax.swing.ImageIcon(getClass().getResource("/GUI/load.png"))); // NOI18N
            menuOpenModel.setText("Import Model");
            menuOpenModel.addActionListener(new java.awt.event.ActionListener() {
                public void actionPerformed(java.awt.event.ActionEvent evt) {
                    menuOpenModelopenExistingModel(evt);
                }
            });
            jMenu1.add(menuOpenModel);
            jMenu1.add(menuSeparator);

            menuExit.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_W, java.awt.event.InputEvent.CTRL_MASK));
            menuExit.setText("Exit");
            menuExit.addActionListener(new java.awt.event.ActionListener() {
                public void actionPerformed(java.awt.event.ActionEvent evt) {
                    menuExitActionPerformed(evt);
                }
            });
            jMenu1.add(menuExit);

            MenuBar.add(jMenu1);

            setJMenuBar(MenuBar);

            javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
            getContentPane().setLayout(layout);
            layout.setHorizontalGroup(
                layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addComponent(ReactionPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 1081, Short.MAX_VALUE)
                .addGroup(layout.createSequentialGroup()
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                        .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                            .addGap(0, 0, Short.MAX_VALUE)
                            .addComponent(generateDataHeader)
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                            .addComponent(InsertFittingData)
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                            .addComponent(RunParameterEstimation)
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                            .addComponent(saveModelButton, javax.swing.GroupLayout.PREFERRED_SIZE, 76, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                            .addComponent(SpeciesPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 531, Short.MAX_VALUE)
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                            .addComponent(EnzymePanel, javax.swing.GroupLayout.DEFAULT_SIZE, 532, Short.MAX_VALUE))
                        .addGroup(layout.createSequentialGroup()
                            .addContainerGap()
                            .addComponent(ModelNameLabel)
                            .addGap(4, 4, 4)
                            .addComponent(ModelNameTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 165, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jLabel9, javax.swing.GroupLayout.PREFERRED_SIZE, 68, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                            .addComponent(GRaPeLogo)))
                    .addContainerGap())
            );
            layout.setVerticalGroup(
                layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(layout.createSequentialGroup()
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                        .addGroup(layout.createSequentialGroup()
                            .addGap(2, 2, 2)
                            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                .addComponent(ModelNameLabel)
                                .addComponent(ModelNameTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                        .addComponent(GRaPeLogo)
                        .addComponent(jLabel9))
                    .addGap(8, 8, 8)
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                        .addComponent(SpeciesPanel, javax.swing.GroupLayout.PREFERRED_SIZE, 283, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(EnzymePanel, javax.swing.GroupLayout.PREFERRED_SIZE, 283, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                    .addComponent(ReactionPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 411, Short.MAX_VALUE)
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(RunParameterEstimation)
                        .addComponent(saveModelButton)
                        .addComponent(generateDataHeader)
                        .addComponent(InsertFittingData))
                    .addContainerGap())
            );

            layout.linkSize(javax.swing.SwingConstants.VERTICAL, new java.awt.Component[] {EnzymePanel, SpeciesPanel});

            pack();
        }// </editor-fold>//GEN-END:initComponents

    private void menuOpenModelopenExistingModel(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_menuOpenModelopenExistingModel
        // TODO add your handling code here:
        try {
            boolean openFile = true;
            JFileChooser fc = new JFileChooser();
            fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
            int returnVal = fc.showOpenDialog(this);
            String fileName;
            if (returnVal == JFileChooser.APPROVE_OPTION) {
                File file = fc.getSelectedFile();
                fileName = file.getPath();
                //This is where a real application would open the chosenfile.
                if (fileName.contains(".xml") || fileName.contains(".sbml")) {

                    System.out.println("Importing SBML File: " + file.getName());
                    ImportSBML sbmlmodel = new ImportSBML(file);
                    ImportModel(sbmlmodel);
                    menuSaveModel.setEnabled(true);
                    saveModelButton.setEnabled(true);

                } else {
                    JOptionPane.showMessageDialog(null, "File Format NOT supported", "File Open Error", JOptionPane.ERROR_MESSAGE);
                }

            }

        } catch (Exception e) {

            e.printStackTrace();
            JOptionPane.showMessageDialog(null,
                                            "Unable to open model. Check SBML version and validity.",
                                            "SBML Error",
                                            JOptionPane.WARNING_MESSAGE);

        }
    }//GEN-LAST:event_menuOpenModelopenExistingModel

    private void menuSaveModelsaveModelAsSBMLFile(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_menuSaveModelsaveModelAsSBMLFile
        // TODO add your handling code here:
        try {
            if(checkModel()){
            modelName = ModelNameTextField.getText();
            
            JFileChooser fc = new JFileChooser();
            int returnVal = fc.showSaveDialog(this);
            
            File file;
            if (returnVal == JFileChooser.APPROVE_OPTION) {
                
                ExportSBML modelConvertor;
                if(GAdone){
                    modelConvertor = new ExportSBML(allthereactions, allthespecies, alltheenzymes, modelName, parameters);
                }else{
                    modelConvertor = new ExportSBML(allthereactions, allthespecies, alltheenzymes, modelName);
                }
                SBMLDocument sbml2doc = new SBMLDocument(modelConvertor.getsbmldoc());
                            
                if (importedmodel == true){
                    sbml2doc.getModel().setId(modelID);
                }else{
                    sbml2doc.getModel().setId(modelName.replaceAll("\\s+",""));
                }
                
                file = fc.getSelectedFile();

                if (file.exists()) {
                    int response = JOptionPane.showConfirmDialog(null,
                            "Overwrite existing file?", "Confirm Overwrite",
                            JOptionPane.OK_CANCEL_OPTION,
                            JOptionPane.QUESTION_MESSAGE);
                    if (response == JOptionPane.OK_OPTION) {
                        SBMLWriter.write(sbml2doc, file.getPath(), file.getPath(), "X");    //(sbmldoc, file name, program name, version)
                    }
                } else {
                    if (file.getPath().contains(".xml"))
                {
                    SBMLWriter.write(sbml2doc, file.getPath(), file.getPath(), "X");
                    
                }else {
                    SBMLWriter.write(sbml2doc, file.getPath()+".xml", file.getPath(), "X");
                }
                }
                
                System.out.println("Saving model: " + file.getPath());
                importedmodel=true;
            }else{
                System.out.println("Exporting SBML File command cancelled by user");
            }
            }else{
                System.out.println("something wong");
                if(missingSpecies.size()>0){
                    missingSpeciesText.append("Missing Species: \n");
                    for(String missing : missingSpecies){
                        missingSpeciesText.append(missing+"\n");
                    }
                }else{
                    missingSpeciesText.append("No Species Missing!");
                }
                if(missingEnzymes.size()>0){
                    missingEnzymeText.append("Missing Enzymes: \n");
                    for(String missing : missingEnzymes){
                        missingEnzymeText.append(missing+"\n");
                    }
                }else{
                    missingEnzymeText.append("No Enzymes Missing!");
                }
                JOptionPane.showMessageDialog(null, missingSpeciesPanel, "EXPORT ERROR", JOptionPane.ERROR_MESSAGE);
                missingEnzymeText.setText(null);
                missingSpeciesText.setText(null);
            }
            
        } catch (HeadlessException | XMLStreamException | FileNotFoundException | SBMLException e) {
            System.err.println(e.getMessage());
            System.out.println("Cant Save Model");
        }
    }//GEN-LAST:event_menuSaveModelsaveModelAsSBMLFile

    private void DeleteSpeciesActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_DeleteSpeciesActionPerformed
        // TODO add your handling code here:
        delete_metabolite = true;
        int row = SpeciesTable.getSelectedRow();
        if(row !=-1){
        String localname = (String) SpeciesTable.getValueAt(row, 1);
        for (int i =0; i < allthespecies.size(); i++){
            String name = allthespecies.get(i).getName();
            if (name.equals(localname)){
                allthespecies.remove(i);
                numberofspecies.setText(Integer.toString(allthespecies.size()));
            }
        }
        SpeciesBox.removeItem(SpeciesTable.getValueAt(row, 1));
        }
        if(row!=-1){
        DefaultTableModel apple = (DefaultTableModel) SpeciesTable.getModel();
        apple.removeRow(row);
        
        int totalrow = SpeciesTable.getRowCount();
        Object[][] temp = new Object[totalrow][4];
        
        for (int i =0; i < SpeciesTable.getRowCount(); i++){
            for (int column = 1; column<4; column++){
                temp[i][column]=SpeciesTable.getModel().getValueAt(i, column);
            }
        }
        
        for (int i =0; i < totalrow; i++){
            temp[i][0]=i+1;
        }
        
        SpeciesTableInfo = temp;
        setUpSpeciesTable(temp);
        }
        
    }//GEN-LAST:event_DeleteSpeciesActionPerformed
        
    private void DeleteReactionActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_DeleteReactionActionPerformed
        // TODO add your handling code here:
        int row = ReactionTable.getSelectedRow();
        delete_reaction = true;
        
        if(row !=-1){
        String localname = (String) ReactionTable.getValueAt(row, 4);
        for (int i =0; i < allthereactions.size(); i++){
            String name = allthereactions.get(i).toString();
            if (name.equals(localname)){
                allthereactions.remove(i);
                numberofreactions.setText(Integer.toString(allthereactions.size()));
            }
        }
        }
        
        if(row!=-1){
            DefaultTableModel apple = (DefaultTableModel) ReactionTable.getModel();
            apple.removeRow(row);
        
            int totalrow = ReactionTable.getRowCount();
            Object[][] temp = new Object[totalrow][5];

            for (int i =0; i < totalrow; i++){
                temp[i][0]=i+1;
            }
            for (int i =0; i < ReactionTable.getRowCount(); i++){

                for (int column = 1; column<5; column++){
                        temp[i][column]=ReactionTable.getModel().getValueAt(i, column);
                }
            }
            ReactionTableInfo=temp;
            setUpReactionTable(temp);
        
        }
    }//GEN-LAST:event_DeleteReactionActionPerformed

    private void DeleteEnzymeActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_DeleteEnzymeActionPerformed
        // TODO add your handling code here:
        delete_enzyme = true;
        int row = EnzymesTable.getSelectedRow();
        if(row!=-1){
        String localname = (String) EnzymesTable.getValueAt(row, 1);
        for (int i =0; i < alltheenzymes.size(); i++){
            String name = alltheenzymes.get(i).getName();
            if (name.equals(localname)){
                alltheenzymes.remove(i);
                numberofenzymes.setText(Integer.toString(alltheenzymes.size()));
            }
        }
        EnzymeBox.removeItem(EnzymesTable.getValueAt(row, 1));
        }
        if(row !=-1){
        DefaultTableModel apple = (DefaultTableModel) EnzymesTable.getModel();
        apple.removeRow(row);
        int totalrow = EnzymesTable.getRowCount();
        Object[][] temp = new Object[totalrow][3];
        
        for (int i =0; i < totalrow; i++){
            temp[i][0]=i+1;
        }
        for (int i =0; i < EnzymesTable.getRowCount(); i++){
           
            for (int column = 1; column<3; column++){
                    temp[i][column]=EnzymesTable.getModel().getValueAt(i, column);
            }
        }
        EnzymeTableInfo=temp;
        setUpEnzymeTable(temp);
        }
    }//GEN-LAST:event_DeleteEnzymeActionPerformed

    private void saveModelButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_saveModelButtonActionPerformed
       
        try {
            
            if(checkModel()){
            modelName = ModelNameTextField.getText();
            
            JFileChooser fc = new JFileChooser();
            int returnVal = fc.showSaveDialog(this);
            
            File file;
            if (returnVal == JFileChooser.APPROVE_OPTION) {
                
                ExportSBML modelConvertor;
                if(frame==null){
                    GAdone=false;
                }
                else if(frame.GeneticAlgorithm.isDone()==true){
                    parameters = frame.getParameters();
                    GAdone=true;
                }
                if(GAdone){                    
                    modelConvertor = new ExportSBML(allthereactions, allthespecies, alltheenzymes, modelName, parameters);
                }else{
                    modelConvertor = new ExportSBML(allthereactions, allthespecies, alltheenzymes, modelName);
                }
                SBMLDocument sbml2doc = new SBMLDocument(modelConvertor.getsbmldoc());
                            
                if (importedmodel == true){
                    sbml2doc.getModel().setId(modelID);
                }else{
                    sbml2doc.getModel().setId(modelName.replaceAll("\\s+",""));
                }
                
                file = fc.getSelectedFile();

                if (file.exists()) {
                    int response = JOptionPane.showConfirmDialog(null,
                            "Overwrite existing file?", "Confirm Overwrite",
                            JOptionPane.OK_CANCEL_OPTION,
                            JOptionPane.QUESTION_MESSAGE);
                    if (response == JOptionPane.OK_OPTION) {
                        SBMLWriter.write(sbml2doc, file.getPath(), file.getPath(), "X");    //(sbmldoc, file name, program name, version)                       
                    }
                } else {
                    if (file.getPath().contains(".xml"))
                {
                    SBMLWriter.write(sbml2doc, file.getPath(), file.getPath(), "X");
                    
                }else {
                    SBMLWriter.write(sbml2doc, file.getPath()+".xml", file.getPath(), "X");
                }
                }                
                System.out.println("Saving model: " + file.getPath());
                importedmodel=true;
            }else{
                System.out.println("Exporting SBML File command cancelled by user");
            }
            }else{
                System.out.println("something wrong");
                if(missingSpecies.size()>0){
                    missingSpeciesText.append("Missing Species: \n");
                    for(String missing : missingSpecies){
                        missingSpeciesText.append(missing+"\n");
                    }
                }else{
                    missingSpeciesText.append("No Species Missing!");
                }
                if(missingEnzymes.size()>0){
                    missingEnzymeText.append("Missing Enzymes: \n");
                    for(String missing : missingEnzymes){
                        missingEnzymeText.append(missing+"\n");
                    }
                }else{
                    missingEnzymeText.append("No Enzymes Missing!");
                }
                JOptionPane.showMessageDialog(null, missingSpeciesPanel, "EXPORT ERROR", JOptionPane.ERROR_MESSAGE);
                missingEnzymeText.setText(null);
                missingSpeciesText.setText(null);
            }
            
        } catch (HeadlessException | XMLStreamException | FileNotFoundException | SBMLException e) {
            System.err.println(e.getMessage());
            System.out.println("Cant Save Model");
        }
    }//GEN-LAST:event_saveModelButtonActionPerformed

    private void AddEnzymeActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_AddEnzymeActionPerformed
        // TODO add your handling code here:
    int totalrow = EnzymesTable.getRowCount()+1;
        Object[][] temp = new Object[totalrow][3];
        
        for (int i =0; i < totalrow; i++){
            temp[i][0]=i+1;
        }
    for (int i =0; i < EnzymesTable.getRowCount(); i++){
           
            for (int column = 1; column<3; column++){
                    temp[i][column]=EnzymesTable.getModel().getValueAt(i, column);
            }
        }
    EnzymeTableInfo=temp;
    setUpEnzymeTable(temp);
    }//GEN-LAST:event_AddEnzymeActionPerformed

    private void AddSpeciesActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_AddSpeciesActionPerformed
        // TODO add your handling code here:
        int totalrow = SpeciesTable.getRowCount()+1;
        Object[][] temp = new Object[totalrow][4];
        for (int i =0; i < SpeciesTable.getRowCount(); i++){
            for (int column = 1; column<4; column++){
                temp[i][column]=SpeciesTable.getModel().getValueAt(i, column);
            }
        }
        
        for (int i =0; i < totalrow; i++){
            temp[i][0]=i+1;
        }
        temp[totalrow-1][3]=false;
        
        SpeciesTableInfo = temp;
        setUpSpeciesTable(temp);
    }//GEN-LAST:event_AddSpeciesActionPerformed

    private void ReactionTableMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_ReactionTableMouseClicked
        // TODO add your handling code here:
        if (evt.getClickCount() == 1) {
          JTable target = (JTable)evt.getSource();
          int column = target.getSelectedColumn();
          if (column == 1){
            TableColumn enzymeColumn = ReactionTable.getColumnModel().getColumn(1);
            setUpEnzymeColumn(ReactionTable, enzymeColumn);
          }
          if (column == 2){
            TableColumn enzymeColumn = ReactionTable.getColumnModel().getColumn(2);
            setUpActivatorColumn(ReactionTable, enzymeColumn);
          }
          if (column == 3){
            TableColumn enzymeColumn = ReactionTable.getColumnModel().getColumn(3);
            setUpInhibitorColumn(ReactionTable, enzymeColumn);
          }
          
        }
        if (evt.getClickCount() == 2) {
            JTable target = (JTable)evt.getSource();
            int column = target.getSelectedColumn();
            int row = target.getSelectedRow();
        
            if (column == 4){
                String activator = (String) ReactionTable.getValueAt(row, 2);
                String inhibitor = (String) ReactionTable.getValueAt(row, 3);
                ArrayList<Compound> modifier = new ArrayList();
                String enz = (String) ReactionTable.getValueAt(row, 1);
                int editreaction =-1;
                boolean edit=false;
                String[] subs= {"N/A","N/A","N/A"}, prods = {"N/A","N/A","N/A"};
                int[] s = {1,1,1} ,p = {1,1,1};
                int regulation = 1;
                
                sub1box.removeAllItems();
                sub2box.removeAllItems();
                sub3box.removeAllItems();
                prod1box.removeAllItems();
                prod2box.removeAllItems();
                prod3box.removeAllItems();
                for(int i =0; i <SpeciesBox.getItemCount();i++){
                    sub1box.addItem(SpeciesBox.getItemAt(i));
                    sub2box.addItem(SpeciesBox.getItemAt(i));
                    sub3box.addItem(SpeciesBox.getItemAt(i));
                    prod1box.addItem(SpeciesBox.getItemAt(i));
                    prod2box.addItem(SpeciesBox.getItemAt(i));
                    prod3box.addItem(SpeciesBox.getItemAt(i));
                }
                
                if(ReactionTable.getValueAt(row, column)!=null){
                    String line = (String) ReactionTable.getValueAt(row, column);
                    String[] tokens = line.split("\\W");
                    ModelReaction empty = null;
                    for (int i =0; i <allthereactions.size(); i++){
                        String id = allthereactions.get(i).getReactionID();
                        if(id.matches(tokens[0])){
                            empty=allthereactions.get(i);
                            editreaction=i;
                        }
                    }
                    ArrayList<Compound> substrates = empty.getSubstrates();
                    ArrayList<Compound> products = empty.getProducts();
                    ArrayList substoichio = empty.getSubstratesStoichio();
                    ArrayList prodstoichio = empty.getProductsStoichio();
                    
                    for ( int i=0;i<substrates.size();i++){
                        subs[i]=substrates.get(i).getName();
                    }
                    for ( int i=0;i<products.size();i++){
                        prods[i]=products.get(i).getName();
                    }
                    for ( int i=0;i<substoichio.size();i++){
                        s[i]=(int) substoichio.get(i);
                    }
                    for ( int i=0;i<prodstoichio.size();i++){
                        p[i]=(int) prodstoichio.get(i);
                    }
                    
                    sub1box.setSelectedItem(subs[0]);
                    sub2box.setSelectedItem(subs[1]);
                    sub3box.setSelectedItem(subs[2]);
                    prod1box.setSelectedItem(prods[0]);
                    prod2box.setSelectedItem(prods[1]);
                    prod3box.setSelectedItem(prods[2]);
                    sub1spin.setValue(s[0]);
                    sub2spin.setValue(s[1]);
                    sub3spin.setValue(s[2]);
                    prod1spin.setValue(p[0]);
                    prod2spin.setValue(p[1]);
                    prod3spin.setValue(p[2]);
                    
                    edit=true;
                }else{
                    sub1spin.setValue(1);
                    sub2spin.setValue(1);
                    sub3spin.setValue(1);
                    prod1spin.setValue(1);
                    prod2spin.setValue(1);
                    prod3spin.setValue(1);
                }
                
                int result = JOptionPane.showConfirmDialog(null, ReactionBuilder, "Please build your reaction", JOptionPane.OK_CANCEL_OPTION);
                
                if (result == JOptionPane.OK_OPTION) {
                    
                    boolean error=false;
                    subs[0]=(String) sub1box.getSelectedItem();
                    subs[1]=(String) sub2box.getSelectedItem();
                    subs[2]=(String) sub3box.getSelectedItem();
                    prods[0]=(String) prod1box.getSelectedItem();
                    prods[1]=(String) prod2box.getSelectedItem();
                    prods[2]=(String) prod3box.getSelectedItem();
                    s[0] = (int) sub1spin.getValue();
                    s[1] = (int) sub2spin.getValue();
                    s[2] = (int) sub3spin.getValue();
                    p[0] = (int) prod1spin.getValue();
                    p[1] = (int) prod2spin.getValue();
                    p[2] = (int) prod3spin.getValue();
                    
                    ArrayList substrates = new ArrayList<>();
                    ArrayList products = new ArrayList<>();
                    ArrayList substoi = new ArrayList<>();
                    ArrayList prodstoi = new ArrayList<>();
                    for (int i=0; i <3; i++) {
                        if(subs[i].equals("N/A")){
                            
                        }else{
                            substrates.add(subs[i]);
                            substoi.add(s[i]);
                        }
                        if(prods[i].equals("N/A")){
                            
                        }else{
                            products.add(prods[i]);
                            prodstoi.add(p[i]);
                        }
                    }
                    if (substrates.size()==0){
                        JOptionPane.showMessageDialog(null, "Please select AT LEAST one substrate", "Reaction Error", JOptionPane.ERROR_MESSAGE);
                        error=true;
                    }
                    if (products.size()==0){
                        JOptionPane.showMessageDialog(null, "Please select AT LEAST one product", "Reaction Error", JOptionPane.ERROR_MESSAGE);
                        error=true;
                    }
                    
                    if(error==false){
                        if(enz!=null){
                            if(enz.equals("N/A")){
                                JOptionPane.showMessageDialog( null, EnzymeBox, "Please pick an enzyme", JOptionPane.QUESTION_MESSAGE);
                                enz = (String) EnzymeBox.getSelectedItem();
                                if(enz==null){
                                    JOptionPane.showMessageDialog( null, "Cannot Build a Reaction Without an Enzyme!", "Please pick an enzyme", JOptionPane.ERROR_MESSAGE);
                                }
                            }
                        }else{
                            JOptionPane.showMessageDialog( null, EnzymeBox, "Please pick an enzyme", JOptionPane.QUESTION_MESSAGE);
                            enz = (String) EnzymeBox.getSelectedItem();                            
                        }
                        
                        if(!activator.equals("N/A")&&!inhibitor.equals("N/A")){
                            regulation=4;
                            if(activator.equals(inhibitor)){
                                JOptionPane.showMessageDialog( null, "Having the same activator and inhibitor is not advised.", "Not a good idea", JOptionPane.ERROR_MESSAGE);
                            }
                            for (int i =0; i<allthespecies.size();i++){
                                if(activator.equals(allthespecies.get(i).getName())){
                                    modifier.add(allthespecies.get(i));
                                }
                            }
                            for ( int i =0; i<allthespecies.size();i++){
                                if(inhibitor.equals(allthespecies.get(i).getName())){
                                    modifier.add(allthespecies.get(i));
                                }
                            }
                        } else{
                            if(activator.equals("N/A")){
                                if(inhibitor.equals("N/A")){
                                    regulation=1;
                                }else{
                                    regulation=3;
                                    for ( int i =0; i<allthespecies.size();i++){
                                        if(inhibitor.equals(allthespecies.get(i).getName())){
                                            modifier.add(allthespecies.get(i));
                                        }
                                    }
                                }
                            }else {
                                regulation=2;
                                for ( int i =0; i<allthespecies.size();i++){
                                    if(activator.equals(allthespecies.get(i).getName())){
                                        modifier.add(allthespecies.get(i));
                                    }
                                }
                            }
                        }
                        
                        ArrayList<Compound> subS = new ArrayList<>();
                        ArrayList<Compound> prodS = new ArrayList<>();
                        
                        for(int i=0; i < substrates.size();i++){
                            for(int j=0; j < allthespecies.size(); j++){
                                String compname = allthespecies.get(j).getName();
                                if(substrates.get(i).equals(compname)){
                                    subS.add(allthespecies.get(j));
                                }
                            }
                        }
                        for(int i=0; i < products.size();i++){
                            for(int j=0; j < allthespecies.size(); j++){
                                String compname = allthespecies.get(j).getName();
                                if(products.get(i).equals(compname)){
                                    prodS.add(allthespecies.get(j));
                                }
                            }
                        }
                        Enzyme newenzyme = null;
                        for(int i=0; i < alltheenzymes.size();i++){
                            if (enz.equals(alltheenzymes.get(i).getName())){
                                newenzyme=alltheenzymes.get(i);
                            }
                        }
                        
                        if(edit==true){
                            String id = allthereactions.get(editreaction).getReactionID();
                            String name = allthereactions.get(editreaction).getName();
                            MetabolicReaction editmetabolicreaction = new MetabolicReaction(enz);
                            editmetabolicreaction.setID(id);
                            editmetabolicreaction.setName(name);
                            editmetabolicreaction.addSubstrates(subS, substoi);
                            editmetabolicreaction.addProducts(prodS, prodstoi);
                            editmetabolicreaction.setRegulation(regulation);
                            if(regulation>1){
                                editmetabolicreaction.setModifier(modifier);
                            }
                            ModelReaction editedreaction = new ModelReaction(newenzyme, editmetabolicreaction);
                            allthereactions.remove(editreaction);
                            allthereactions.add(editreaction, editedreaction);
                        }else{
                            if(newenzyme!=null){
                            MetabolicReaction newreaction = new MetabolicReaction(newenzyme.getName());
                            for(int i =0; i<subS.size();i++){
                                newreaction.addSubstrate(subS.get(i), (int) substoi.get(i));
                            }
                            for(int i =0; i<prodS.size();i++){
                                newreaction.addProduct(prodS.get(i), (int) prodstoi.get(i));
                            }
                            newreaction.setRegulation(regulation);
                            if(regulation>1){
                                newreaction.setModifier(modifier);
                            }
                            ModelReaction reactiontoadd = new ModelReaction(newenzyme, newreaction);
                            allthereactions.add(reactiontoadd);
                            }
                        }
                        
                    }
                    ////
                    if(edit==true){
                        ReactionTable.setValueAt(allthereactions.get(editreaction).toString(), row,4 );
                    }else{
                        ReactionTable.setValueAt(allthereactions.get(allthereactions.size()-1).toString(), row,4 );
                        ReactionTable.setValueAt(allthereactions.get(allthereactions.size()-1).getEnzyme().getName(), row,1 );
                    }
                    numberofreactions.setText(Integer.toString(allthereactions.size()));
                }
            }
        }
    }//GEN-LAST:event_ReactionTableMouseClicked
    
    private void AddReactionActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_AddReactionActionPerformed
        // TODO add your handling code here:
        int totalrow = ReactionTable.getRowCount()+1;
        Object[][] temp = new Object[totalrow][5];
        
        for (int i =0; i < totalrow; i++){
            temp[i][0]=i+1;
        }
        temp[totalrow-1][2]="N/A";
        temp[totalrow-1][3]="N/A";
        for (int i =0; i < ReactionTable.getRowCount(); i++){
            for (int column = 1; column<5; column++){
                    temp[i][column]=ReactionTable.getModel().getValueAt(i, column);
            }
        }
        ReactionTableInfo=temp;
        setUpReactionTable(temp);
    }//GEN-LAST:event_AddReactionActionPerformed
        
    private void SpeciesTableKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_SpeciesTableKeyReleased
        // TODO add your handling code here:
        if (evt.getKeyCode ()== KeyEvent.VK_ENTER ||evt.getKeyCode ()== KeyEvent.VK_UP||evt.getKeyCode ()== KeyEvent.VK_DOWN||evt.getKeyCode ()== KeyEvent.VK_LEFT||evt.getKeyCode ()== KeyEvent.VK_RIGHT){
            JTable target = (JTable)evt.getSource();
            int column = target.getSelectedColumn();
            int row = target.getSelectedRow();
            if(column==1){                
                tempname = (String) SpeciesTable.getValueAt(row,column);
            }
        }
    }//GEN-LAST:event_SpeciesTableKeyReleased

    private void SpeciesTableMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_SpeciesTableMouseClicked
        // TODO add your handling code here:
        JTable target = (JTable)evt.getSource();
        int column = target.getSelectedColumn();
        int row = target.getSelectedRow();
        if(column==1){
            tempname = (String) SpeciesTable.getValueAt(row,column);
        }
    }//GEN-LAST:event_SpeciesTableMouseClicked

    private void EnzymesTableKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_EnzymesTableKeyReleased
        // TODO add your handling code here:
        if (evt.getKeyCode ()== KeyEvent.VK_ENTER ||evt.getKeyCode ()== KeyEvent.VK_UP||evt.getKeyCode ()== KeyEvent.VK_DOWN||evt.getKeyCode ()== KeyEvent.VK_LEFT||evt.getKeyCode ()== KeyEvent.VK_RIGHT){
            JTable target = (JTable)evt.getSource();
            int column = target.getSelectedColumn();
            int row = target.getSelectedRow();
            if(column==1){
                tempname = (String) EnzymesTable.getValueAt(row,column);
            }
        }
    }//GEN-LAST:event_EnzymesTableKeyReleased

    private void EnzymesTableMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_EnzymesTableMouseClicked
        // TODO add your handling code here:
        JTable target = (JTable)evt.getSource();
        int column = target.getSelectedColumn();
        int row = target.getSelectedRow();
        if(column==1){
            tempname = (String) EnzymesTable.getValueAt(row,column);
        }
    }//GEN-LAST:event_EnzymesTableMouseClicked

    private void generateDataHeaderActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_generateDataHeaderActionPerformed
        try {
            Object[] options = { "Steady State", "Time Course" };
            int SSorTC = JOptionPane.showOptionDialog(null, "What Type of Data is it?", 
                    "Please Choose One", JOptionPane.YES_NO_CANCEL_OPTION, 
                    JOptionPane.QUESTION_MESSAGE,null, options, options[0]);
            if(SSorTC==JOptionPane.YES_OPTION){//SS type data
                String headerFileID = "";
                String headerFileName = "";
                for (int c = 0; c < allthespecies.size(); c++) {
                    headerFileName += allthespecies.get(c).getName() + "\t";
                    headerFileName += allthespecies.get(c).getID() + "\n";
                }
                for (int r = 0; r < allthereactions.size(); r++) {
                    headerFileName += allthereactions.get(r).getName() + "\t";
                    headerFileName += allthereactions.get(r).getID() + "\n";
                }

                JFileChooser fc = new JFileChooser();
                int returnVal = fc.showSaveDialog(this);
                File file = fc.getSelectedFile();
                if (returnVal == JFileChooser.APPROVE_OPTION) {
                    if (file.exists()) {
                        int response = JOptionPane.showConfirmDialog(null,
                                "Overwrite existing file?", "Confirm Overwrite",
                                JOptionPane.OK_CANCEL_OPTION,
                                JOptionPane.QUESTION_MESSAGE);
                        if (response == JOptionPane.OK_OPTION) {
                            writeHeader(file.getPath(), headerFileName, headerFileID);
                        }
                    } else {
                        writeHeader(file.getPath()+".txt", headerFileName, headerFileID);
                    }
                }
            }else if (SSorTC==JOptionPane.NO_OPTION){
                String headerFileID = "Time\t";
                String headerFileName = "N/A\t";
                for (int c = 0; c < allthespecies.size(); c++) {
                    headerFileName += allthespecies.get(c).getName() + "\t";
                    headerFileID += allthespecies.get(c).getID() + "\t";
                }
                for (int e = 0; e < alltheenzymes.size(); e++) {
                    headerFileName += alltheenzymes.get(e).getName() + "\t";
                    headerFileID += alltheenzymes.get(e).getID() + "\t";
                }
                for (int r = 0; r < allthereactions.size(); r++) {
                    headerFileName += allthereactions.get(r).getName() + "\t";
                    headerFileID += allthereactions.get(r).getID() + "\t";
                }

                JFileChooser fc = new JFileChooser();
                int returnVal = fc.showSaveDialog(this);
                File file = fc.getSelectedFile();
                if (returnVal == JFileChooser.APPROVE_OPTION) {
                    if (file.exists()) {
                        int response = JOptionPane.showConfirmDialog(null,
                                "Overwrite existing file?", "Confirm Overwrite",
                                JOptionPane.OK_CANCEL_OPTION,
                                JOptionPane.QUESTION_MESSAGE);
                        if (response == JOptionPane.OK_OPTION) {
                            writeHeader(file.getPath(), headerFileName, headerFileID);
                        }
                    } else {
                        writeHeader(file.getPath()+".txt", headerFileName, headerFileID);
                    }
                }
            }
        } catch (Exception e) {
            System.err.println(e.getMessage());
        }
    }//GEN-LAST:event_generateDataHeaderActionPerformed

    private void RunParameterEstimationActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_RunParameterEstimationActionPerformed
        Object[] options = {"Run","Cancel"};
        populationText.setToolTipText("Number of models to be solved (high values would lead to longer parameter estimation time)");
        generationText.setToolTipText("Maximum number of generations to go for before stopping even if optimal parameters aren’t found. (high values would lead to longer parameter estimation time)");
        plateauText.setToolTipText("Number of generations where there is no further increase in fitness score. After reaching the given number, the process would stop. ");
        plagueText.setToolTipText("Process of removing unfit individuals and maintaining only the initial population size at given number of generation.");
        
        
        int result = JOptionPane.showOptionDialog(null, GAparameters, "Parameters for Genetic Algorithm", JOptionPane.OK_CANCEL_OPTION, JOptionPane.NO_OPTION, null, options, options[0]);
        
        if (result == JOptionPane.OK_OPTION) {
            int popsize = Integer.parseInt(populationText.getText());
            int maxgen = Integer.parseInt(generationText.getText());
            int plateau = Integer.parseInt(plateauText.getText());
            int plague = Integer.parseInt(plagueText.getText());
            int numCore= Integer.parseInt(number_of_core.getText());
                        
            double[][] general_kp_range=new double[5][2];
            double[] fr={Double.parseDouble(vf_low.getText()), Double.parseDouble(vf_upp.getText())};
            general_kp_range[0]=fr;
            double[] rr={Double.parseDouble(vr_low.getText()), Double.parseDouble(vr_upp.getText())};
            general_kp_range[1]=rr;
            double[] ka={Double.parseDouble(ka_low.getText()), Double.parseDouble(ka_upp.getText())};
            general_kp_range[2]=ka;
            double[] ki={Double.parseDouble(ki_low.getText()), Double.parseDouble(ki_upp.getText())};
            general_kp_range[3]=ki;
            double[] km={Double.parseDouble(km_low.getText()), Double.parseDouble(km_upp.getText())};
            general_kp_range[4]=km;
            
            if(FittingData!=null){
                try {
                    boolean SS = FittingData.SSorTC();
                    String tmplink = "/tmp/tmp.xml";
                    File tmpfile = new File(tmplink);
                    ExportSBML modelConvertor = new ExportSBML(allthereactions, allthespecies, alltheenzymes, modelName);
                    SBMLDocument sbml2doc = new SBMLDocument(modelConvertor.getsbmldoc());
                    SBMLWriter.write(sbml2doc, tmpfile.getPath(), tmpfile.getPath(), "X");
                    SystemToSolve bioSystem = new SystemToSolve(tmplink, FittingData);                    
                    
                    frame = new ProgressFrame();
                    frame.setVisible(true);
                    frame.setLocation(300, 180);
                    frame.runGA(bioSystem, SS, popsize, maxgen, plateau, plague, numCore, general_kp_range);
                                        
                } catch (XMLStreamException | SBMLException | ModelOverdeterminedException | IllegalArgumentException ex) {
                    Logger.getLogger(MainWindow.class.getName()).log(Level.SEVERE, null, ex);
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(MainWindow.class.getName()).log(Level.SEVERE, null, ex);
                } catch (IOException ex) {
                    Logger.getLogger(MainWindow.class.getName()).log(Level.SEVERE, null, ex);
                }
            }else{
                JOptionPane.showMessageDialog(null, "Please Insert Fitting Data", "Missing Fitting Data",JOptionPane.ERROR_MESSAGE);
            }
        }
    }//GEN-LAST:event_RunParameterEstimationActionPerformed

    private void InsertFittingDataActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_InsertFittingDataActionPerformed
        try {

            JFileChooser fc = new JFileChooser();
            fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);

            int returnVal = fc.showOpenDialog(this);

            File datafile = null;
            String fileName = null;

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                datafile = fc.getSelectedFile();
                fileName = datafile.getPath();
                //This is where a real application would open the file.
                if (fileName.contains(".txt")){
                    FittingData = new ReadData(datafile);                
                    RunParameterEstimation.setEnabled(true);
                    
                    System.out.println("Selected Data File: " + datafile.getName() + "\n");
                }else{
                    JOptionPane.showMessageDialog(null, "File Format NOT supported", "File Open Error", JOptionPane.ERROR_MESSAGE);
                }

            } else {
                System.out.println("Importing Data File command cancelled by user");
            }

        } catch (Exception e) {
            System.err.println("File not found");
        }
    }//GEN-LAST:event_InsertFittingDataActionPerformed

    private void menuExitActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_menuExitActionPerformed
        // TODO add your handling code here:
        System.exit(0);
    }//GEN-LAST:event_menuExitActionPerformed

    private void menuNewModelActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_menuNewModelActionPerformed
        importedmodel=false;
        
        SpeciesBox.removeAllItems();
        EnzymeBox.removeAllItems();
        SpeciesBox.addItem("N/A");
        allthereactions = new ArrayList<>();
        allthespecies = new ArrayList<>();
        alltheenzymes = new ArrayList<>();
        if(Model!=null){
        Model.clearUserObjects();
        listofspecies.clear();
        listofreactions.clear();
        modelName=null;
        }        
        
        ModelNameTextField.setText("");
        
        numberofspecies.setText(Integer.toString(0));
        numberofreactions.setText(Integer.toString(0));
        numberofenzymes.setText(Integer.toString(0));
        
        setUpSpeciesTable(null);
        setUpEnzymeTable(null);
        setUpReactionTable(null);
        menuSaveModel.setEnabled(false);
        saveModelButton.setEnabled(false);
    }//GEN-LAST:event_menuNewModelActionPerformed

    private void ModelNameTextFieldActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ModelNameTextFieldActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_ModelNameTextFieldActionPerformed

    private void tsv_to_modelActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_tsv_to_modelActionPerformed
        // TODO add your handling code here:
        try {
            TsvToModel tsm = new TsvToModel();
            
            JFileChooser fc = new JFileChooser();
            fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
            int returnVal = fc.showOpenDialog(this);
            String fileName;
            if (returnVal == JFileChooser.APPROVE_OPTION) {
                File file = fc.getSelectedFile();
                fileName = file.getPath();
                //This is where a real application would open the chosenfile.
                if (fileName.contains(".txt")) {

                    System.out.println("Converting Txt File " + file.getName() +" to Model.");
                    tsm.build_model_from_tsv(fileName);
                    build_model_using_tsv(tsm.get_metabolites(),  tsm.get_enzymes(), tsm.get_reactions());
                    
                    menuSaveModel.setEnabled(true);
                    saveModelButton.setEnabled(true);

                } else {
                    JOptionPane.showMessageDialog(null, "File Format NOT supported", "File Open Error", JOptionPane.ERROR_MESSAGE);
                }

            }

        } catch (Exception e) {

            e.printStackTrace();
            JOptionPane.showMessageDialog(null,
                                            "File Format NOT supported.",
                                            "File Error",
                                            JOptionPane.WARNING_MESSAGE);

        }
    }//GEN-LAST:event_tsv_to_modelActionPerformed

    private void build_model_using_tsv(ArrayList metabolites, ArrayList enzymes, ArrayList reactions) {
        try{
            SpeciesBox.removeAllItems();
            EnzymeBox.removeAllItems();
            SpeciesBox.addItem("N/A");
            allthereactions = reactions;
            allthespecies = metabolites;
            alltheenzymes = enzymes;
            addReactionsToReactionList_v2(allthereactions);
            addSpeciesToSpeciesList_v2(allthespecies, alltheenzymes);


        } catch (NullPointerException e) {
             e.printStackTrace();
            JOptionPane.showMessageDialog(null, "Incompatible File format.", "Import Error", JOptionPane.ERROR_MESSAGE);
            
        }
    }
    
    private boolean checkModel() {
        boolean goodModel = false;
        missingSpecies.clear();
        missingEnzymes.clear();
        ArrayList inReactions = new ArrayList<>();
        ArrayList inSpecies = new ArrayList<>();
        ArrayList inEnzymes = new ArrayList<>();
        
        HashSet speciesTempSet = new HashSet<>();
        HashSet enzymeTempSet = new HashSet<>();
        
        for(ModelReaction reaction : allthereactions){            
            for(Compound comp : reaction.getSubstrates()){
                speciesTempSet.add(comp.getName());
            }
            
            for(Compound comp : reaction.getProducts()){
                speciesTempSet.add(comp.getName());
            }
            
            enzymeTempSet.add(reaction.getEnzyme().getName());        
        }
        
        for(Object name : speciesTempSet){
            inReactions.add(name);
        }
        for(Compound comp : allthespecies ){
            inSpecies.add(comp.getName());
        }
        for(Enzyme enz : alltheenzymes ){
            inEnzymes.add(enz.getName());
        }
        enzymeTempSet.removeAll(inEnzymes);
        if(enzymeTempSet.size()>0){missingEnzymes.addAll(enzymeTempSet);}
        
        inReactions.removeAll(inSpecies);
        if(inReactions.size()>0){missingSpecies.addAll(inReactions);}
        
        if(inReactions.size()>0 || enzymeTempSet.size()>0){
            goodModel=false;
        }else{
            goodModel=true;
        }
        return goodModel;
    }
    class HintTextField extends JTextField implements FocusListener {

    private final String hint;
    private boolean showingHint;

    public HintTextField(final String hint) {
        super(hint);
        this.hint = hint;
        if(importedmodel){
            this.showingHint = false;
        }else{
            this.showingHint = true;
        }
        super.addFocusListener(this);
    }

    @Override
    public void focusGained(FocusEvent e) {
        if(this.getText().isEmpty()) {
            if(importedmodel){
                super.setText(modelName);
            }else{
                super.setText("");
            }
            showingHint = false;
        }
    }
    @Override
    public void focusLost(FocusEvent e) {
        if(this.getText().isEmpty()) {
            super.setText(hint);
            showingHint = true;
        }
    }

    @Override
    public String getText() {
        return showingHint ? "" : super.getText();
    }
}
    private void writeHeader(String outputFile, String headerName, String headerID) {
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
            out.write("Copy to Excel. \nRemove the Names. \nKeep only the IDs."
                    + "\nPaste data back into txt file to be inserted into GRaPe, make sure to avoid empty lines.\n\n");
            out.write(headerName + " \n");
            out.write(headerID+" \n\n");
            
            out.close();
        } catch (Exception e) {
            //ignore
        }
    }
    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(MainWindow.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(MainWindow.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(MainWindow.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(MainWindow.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
        //</editor-fold>
        //</editor-fold>
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new MainWindow().setVisible(true);
            }
        });
    }
    
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton AddEnzyme;
    private javax.swing.JButton AddReaction;
    private javax.swing.JButton AddSpecies;
    private javax.swing.JButton DeleteEnzyme;
    private javax.swing.JButton DeleteReaction;
    private javax.swing.JButton DeleteSpecies;
    private javax.swing.JLabel EnzymeCounter;
    private javax.swing.JLabel EnzymeLabel;
    private javax.swing.JPanel EnzymePanel;
    private javax.swing.JScrollPane EnzymeScrollPane;
    private javax.swing.JTable EnzymesTable;
    private javax.swing.JPanel GAparameters;
    private javax.swing.JLabel GRaPeLogo;
    private javax.swing.JLabel GenerationLabel;
    private javax.swing.JButton InsertFittingData;
    private javax.swing.JMenuBar MenuBar;
    private javax.swing.JLabel ModelNameLabel;
    private javax.swing.JTextField ModelNameTextField;
    private javax.swing.JLabel PlagueLabel;
    private javax.swing.JLabel PlateauLabel;
    private javax.swing.JLabel PlateauLabel1;
    private javax.swing.JLabel PlateauLabel10;
    private javax.swing.JLabel PlateauLabel11;
    private javax.swing.JLabel PlateauLabel12;
    private javax.swing.JLabel PlateauLabel13;
    private javax.swing.JLabel PlateauLabel14;
    private javax.swing.JLabel PlateauLabel15;
    private javax.swing.JLabel PlateauLabel16;
    private javax.swing.JLabel PlateauLabel2;
    private javax.swing.JLabel PlateauLabel3;
    private javax.swing.JLabel PlateauLabel4;
    private javax.swing.JLabel PlateauLabel5;
    private javax.swing.JLabel PlateauLabel6;
    private javax.swing.JLabel PlateauLabel7;
    private javax.swing.JLabel PlateauLabel8;
    private javax.swing.JLabel PlateauLabel9;
    private javax.swing.JLabel PopulationLabel;
    private javax.swing.JPanel ReactionBuilder;
    private javax.swing.JLabel ReactionCounter;
    private javax.swing.JLabel ReactionLabel;
    private javax.swing.JPanel ReactionPanel;
    private javax.swing.JScrollPane ReactionScrollPane;
    private javax.swing.JTable ReactionTable;
    private javax.swing.JButton RunParameterEstimation;
    private javax.swing.JLabel SpeciesCounter;
    private javax.swing.JLabel SpeciesLabel;
    private javax.swing.JPanel SpeciesPanel;
    private javax.swing.JScrollPane SpeciesScrollPane;
    private javax.swing.JTable SpeciesTable;
    private javax.swing.JButton generateDataHeader;
    private javax.swing.JTextField generationText;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel6;
    private javax.swing.JLabel jLabel7;
    private javax.swing.JLabel jLabel8;
    private javax.swing.JLabel jLabel9;
    private javax.swing.JMenu jMenu1;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel jPanel3;
    private javax.swing.JScrollPane jScrollPane3;
    private javax.swing.JScrollPane jScrollPane4;
    private javax.swing.JSeparator jSeparator1;
    private javax.swing.JTextField ka_low;
    private javax.swing.JTextField ka_upp;
    private javax.swing.JTextField ki_low;
    private javax.swing.JTextField ki_upp;
    private javax.swing.JTextField km_low;
    private javax.swing.JTextField km_upp;
    private javax.swing.JMenuItem menuExit;
    private javax.swing.JMenuItem menuNewModel;
    private javax.swing.JMenuItem menuOpenModel;
    private javax.swing.JMenuItem menuSaveModel;
    private javax.swing.JPopupMenu.Separator menuSeparator;
    private javax.swing.JTextArea missingEnzymeText;
    private javax.swing.JPanel missingSpeciesPanel;
    private javax.swing.JTextArea missingSpeciesText;
    private javax.swing.JTextField number_of_core;
    private javax.swing.JLabel numberofenzymes;
    private javax.swing.JLabel numberofreactions;
    private javax.swing.JLabel numberofspecies;
    private javax.swing.JTextField plagueText;
    private javax.swing.JTextField plateauText;
    private javax.swing.JTextField populationText;
    private javax.swing.JComboBox prod1box;
    private javax.swing.JSpinner prod1spin;
    private javax.swing.JComboBox prod2box;
    private javax.swing.JSpinner prod2spin;
    private javax.swing.JComboBox prod3box;
    private javax.swing.JSpinner prod3spin;
    public static javax.swing.JButton saveModelButton;
    private javax.swing.JComboBox sub1box;
    private javax.swing.JSpinner sub1spin;
    private javax.swing.JComboBox sub2box;
    private javax.swing.JSpinner sub2spin;
    private javax.swing.JComboBox sub3box;
    private javax.swing.JSpinner sub3spin;
    private javax.swing.JMenuItem tsv_to_model;
    private javax.swing.JTextField vf_low;
    private javax.swing.JTextField vf_upp;
    private javax.swing.JTextField vr_low;
    private javax.swing.JTextField vr_upp;
    // End of variables declaration//GEN-END:variables
    private Object[][] SpeciesTableInfo = new Object [][] {
                {1, "",  null, false}
    };
    private String[] SpeciesTableHeader = new String [] {
        "No.", "Species", "Concentration", "Boundary"
    };
    
    private Object[][] EnzymeTableInfo = new Object [][] {
                {1, "",  null}
            };
    private String[] EnzymeTableHeader = new String [] {
        "No.", "Enzyme", "Concentration"
    };
    
    private Object[][] ReactionTableInfo = new Object [][] {
                {1, "N/A", "N/A", "N/A", null}
            };
    private String[] ReactionTableHeader = new String [] {
        "No.", "Enzyme", "Activator", "Inhibitor", "Reaction"};
    
    private ListOf<Species> listofspecies;    
    private ListOf<Reaction> listofreactions;
    private Model Model;
    private String modelName;
    private String modelID;
    private javax.swing.JComboBox EnzymeBox = new JComboBox(); 
    private javax.swing.JComboBox SpeciesBox = new JComboBox();
    private boolean importedmodel=false;
    private ArrayList<Compound> allthespecies = new ArrayList<>();
    private ArrayList<Enzyme> alltheenzymes = new ArrayList<>();
    private ListOf<LocalParameter> LLP;
    private ArrayList<ModelReaction> allthereactions = new ArrayList<>();
    private ReadData FittingData = null;
    private String tempname = null;
    private ArrayList<String> missingSpecies = new ArrayList<>();
    private ArrayList<String> missingEnzymes = new ArrayList<>();
    public boolean GAdone = false;
    public double[] parameters;
    private ProgressFrame frame = null;
    private boolean delete_reaction = false;
    private boolean delete_enzyme = false;
    private boolean delete_metabolite = false;
}
