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
        int numCore=7;

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
        try {
        switch (args[0].toLowerCase()) {
            case "gui":
                {
                    JFrame f = new MainWindow();
                    f.pack();
                    f.setVisible(true);
                    f.setLocation(600, 100);                                 // determines where in the comp screen it shows up       
                    break;
                }
            case "makemodel":
                {
                    String tsvfile = args[1];
                    String outputlink = args[2];
                    TsvToModel test = new TsvToModel(tsvfile, outputlink);
                    break;
                }
            case "runGA":
                {
                    String modelfile = args[1];
                    String datafile = args[2];
                    String outputlink = args[3];
                    Script runScript = new Script();
                    runScript.addModel(modelfile);
                    runScript.addData(datafile);
                    runScript.build();
                    runScript.runGA(outputlink, popsize, maxgen, plateau, plague, numCore, general_range);
                    break;
                }
            case "runGA_advanced":
                {
                    String modelfile = args[1];
                    String datafile = args[2];
                    String outputlink = args[3];
                    Script runScript = new Script();
                    
                    read_general_range reader = new read_general_range(args[4]);
                    general_range = reader.get_ranges();
                    
                    popsize = Integer.parseInt(args[5]);
                    maxgen = Integer.parseInt(args[6]);
                    plateau = Integer.parseInt(args[7]);
                    plague = Integer.parseInt(args[8]);
                    numCore = Integer.parseInt(args[9]);
                                           
                    runScript.addModel(modelfile);
                    runScript.addData(datafile);
                    runScript.build();
                    runScript.runGA(outputlink, popsize, maxgen, plateau, plague, numCore, general_range);
                    break;
                }
            default:
                {
                    System.out.println("Please Enter One of the Following Options\n"
                            + "GUI       - to use the graphical user interface of the tool.\n"
                            + "makemodel - to convert a tab separated file to SBML file\n"
                            + "            e.g. makemodel model_info.txt output_sbml.xml\n"
                            + "runGA     - to run parameter estimation when given model and fitting data\n"
                            + "            e.g. runGA model.xml fitting_data.txt output_sbml.xml");
                    break;
                }
        }
        
        }catch(ArrayIndexOutOfBoundsException e){
            System.out.println("Please Enter One of the Following Options\n"
                            + "GUI       - to use the graphical user interface of the tool.\n"
                            + "makemodel - to convert a tab separated file to SBML file\n"
                            + "            e.g. makemodel model_info.txt output_sbml.xml\n"
                            + "runGA     - to run parameter estimation when given model and fitting data\n"
                            + "            e.g. runGA model.xml fitting_data.txt output_sbml.xml");
        }
        
    }
  
}