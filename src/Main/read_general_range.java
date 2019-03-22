/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Main;

import Reaction.MetabolicReaction;
import SBML.ExportSBML;
import Species.Compound;
import Species.Enzyme;
import java.io.File;
import java.io.IOException;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import javax.xml.stream.XMLStreamException;
import org.sbml.jsbml.SBMLDocument;
import org.sbml.jsbml.SBMLWriter;

/**
 *
 * @author chuanfuyap
 */
public class read_general_range {
    
    private double[][] general_range = new double[5][2];
    
    public read_general_range(){}
    
    public read_general_range(String filelink) throws IOException, XMLStreamException{

        Path path = FileSystems.getDefault().getPath(filelink);
        List<String> datas = Files.readAllLines(path);
        
        datas.stream().map(s -> s.trim()).filter(s -> !s.isEmpty()).forEach((line) -> {
            String[] temp = line.split("\\s");
            switch (temp[0]) {
                case "fr":
                    double[] fr={Double.parseDouble(temp[1]), Double.parseDouble(temp[2])};
                    general_range[0]=fr;
                    break;
                case "rr":
                    double[] rr={Double.parseDouble(temp[1]), Double.parseDouble(temp[2])};
                    general_range[1]=rr;
                    break;
                case "ki":
                    double[] ki={Double.parseDouble(temp[1]), Double.parseDouble(temp[2])};
                    general_range[2]=ki;
                    break;
                case "ka":
                    double[] ka={Double.parseDouble(temp[1]), Double.parseDouble(temp[2])};
                    general_range[3]=ka;
                    break;
                case "km":
                    double[] km={Double.parseDouble(temp[1]), Double.parseDouble(temp[2])};
                    general_range[4]=km;
                    break;
                default:
                    break;
            }
            
        });
        
    }
    
    public double[][] get_ranges(){
        return general_range;
    }
    
}