/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Main;

import GUI.MainWindow;
import java.io.IOException;
import javax.swing.JFrame;
import javax.xml.stream.XMLStreamException;
import org.sbml.jsbml.validator.ModelOverdeterminedException;

/**
 *
 * @author chuanfuyap
 */
public class Main {

    /**
     * @param args the command line arguments
     * @throws javax.xml.stream.XMLStreamException
     * @throws java.io.IOException
     * @throws org.sbml.jsbml.validator.ModelOverdeterminedException
     * @throws java.lang.InstantiationException
     * @throws java.lang.IllegalAccessException
     * @throws java.lang.NoSuchMethodException
     */
    public static void main(String[] args) throws XMLStreamException, IOException, ModelOverdeterminedException, InstantiationException, IllegalAccessException, IllegalArgumentException, NoSuchMethodException, InterruptedException {
        int popsize=500;
        int maxgen=20;
        int plateau=1;
        int plague=5;

        double[][] general_range = new double[5][2];
        double[] fr={3, 4.2};
        general_range[0]=fr;
        double[] rr={-3, -1};
        general_range[1]=rr;
        double[] ka={0.5, 2};
        general_range[2]=ka;
        double[] ki={-3.5, 2.5};
        general_range[3]=ki;
        double[] km={-5, 3};
        general_range[4]=km;
        
        if (args[0].toLowerCase().equals("gui")){
            JFrame f = new MainWindow();       
            f.pack();
            f.setVisible(true);
            f.setLocation(600, 100);                                 // determines where in the comp screen it shows up       
        }
        if (args[0].toLowerCase().equals("makemodel")){
            String tsvfile = args[1];
            String outputlink = args[2];
            
            TsvToModel test = new TsvToModel(tsvfile, outputlink);
        }
        else if (args[0].toLowerCase().equals("runGA")){
            
            
            String modelfile = args[1];
            String datafile = args[2];
            //optional
            String outputlink = args[3];
           
            //String rangefile = args[3]; 
            
            Script runScript = new Script();

            runScript.addModel(modelfile);
            runScript.addData(datafile);
            runScript.build();
            runScript.runGA(outputlink, popsize, maxgen, plateau, plague, general_range);

        }
        
    }
  
}