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
        
        JFrame f = new MainWindow();       
        f.pack();
        f.setVisible(true);
        f.setLocation(600, 100);                                 // determines where in the comp screen it shows up       
        
    }
  
}
