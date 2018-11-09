/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Main;

import java.io.File;
import java.io.IOException;
import org.apache.commons.math3.random.RandomDataGenerator;
import java.util.ArrayList;
import javax.xml.stream.XMLStreamException;
import org.sbml.jsbml.validator.ModelOverdeterminedException;

/**
// *
 * @author chuanfuyap
 */
public class Testing {
    
    public static void main(String[] args) throws XMLStreamException, IOException, ModelOverdeterminedException, InstantiationException, IllegalAccessException, IllegalArgumentException, NoSuchMethodException, InterruptedException {
        int popsize=250;
        int maxgen=500;
        int plateau=11;
        int plague=5;
        
        double[][] general_range = new double[5][2];
        double[] fr={-3.5,3.5};
        general_range[0]=fr;
        double[] rr={-3.5,3.5};
        general_range[1]=rr;
        double[] ka={-2, 2};
        general_range[2]=ka;
        double[] ki={-2, 2};
        general_range[3]=ki;
        double[] km={-3.5,3.5};
        general_range[4]=km;
        String sslink = "/home/chuanfuyap/ToyModel/Toy_Model.xml";
        
        String link = "/home/chuanfuyap/ToyModel/Toy_Model.xml";
        
        String datalink = "/home/chuanfuyap/ToyModel/temp_test.txt";
        
        String outputlink = "/home/chuanfuyap/ToyModel/Toy_Output";
                
        Script runScript = new Script();
        
        int todo = 1;
                
        switch (todo) {
            case 1://Steady State Solution
                runScript.addModel(sslink);
                runScript.addData(datalink);
//                runScript.build();
//                runScript.SSsolution();
                break;
            
            default:
                break;
        }
            
    }
    private static double RN(double lower, double upper) {
        RandomDataGenerator random = new RandomDataGenerator();
        return random.nextUniform(lower, upper);
    }

}