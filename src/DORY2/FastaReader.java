/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package DORY2;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Locale;

/**
 *
 * @author Mutt
 */
public class FastaReader {

    private String proteinName = null;
    private BufferedReader reader = null;
    public final long fileLenth;
    public long currentPosition;
    public int process; //0-100 indicate the current progress through the file
    private final StringBuilder sb;
    public FastaReader(String FileName) throws IOException {
        File file = new File(FileName);
        reader = new BufferedReader(new FileReader(file), 1048576);
        fileLenth = file.length();
        sb=new StringBuilder();
    }
    
    public void closeFile() throws IOException{
        reader.close();
    }

    public Protein readNext() throws IOException {
        String strLine;
        Protein protein=new Protein();
        sb.setLength(0);

        while ((strLine = reader.readLine()) != null) {
            if (!strLine.isEmpty()) {
                currentPosition += strLine.length() + 1;
                if (strLine.charAt(0) == '>') {
                    if (proteinName == null) {
                        proteinName = strLine;
                    } else {
                        protein.name= proteinName;
                        proteinName = strLine;
                        break;
                    }
                } else if (proteinName != null) {
                    strLine = strLine.toUpperCase(Locale.ENGLISH);
                    sb.append(strLine.replaceAll("[^\\*ARNDCQEGHILKMFPSTWYVBZX]", ""));
                }
            } else {
                currentPosition++;
            }
            process = (int)(100 * currentPosition / fileLenth);
        }
        
        protein.sequence=sb.toString();
        if (strLine == null) {
            process=100;
            if (sb.length()>0) {
                protein.name = proteinName;
            } else if (protein.name == null) {
                protein = null;
            }
        }
        return protein;
    }
}
