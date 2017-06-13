/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package DORY2;

/**
 *
 * @author Mutt
 */
public class KASHTailHunter {

    private final int tmdFrameLength;// the transmembrane domain frame length;
    private final int tmdHydrophobicThreshold;//the Threshold to be a TMD;
    private final int maxKASHTailLength; // the maximum of the KASH domain length;
    private final int minKASHTailLength; // the minimum of the KASH domain length; 
    private final String regex;

    public KASHTailHunter(int TmdFrameLength,
            int TmdHydrophobicThreshold,
            int MaxKASHTailLength,
            int MinKASHTailLength,
            String Regex) {
        tmdFrameLength = TmdFrameLength;
        tmdHydrophobicThreshold = TmdHydrophobicThreshold;
        maxKASHTailLength = MaxKASHTailLength;
        minKASHTailLength = MinKASHTailLength;
        regex = Regex;
    }

    public String searchKASH(String proteinSeq) {//return the KASH tail sequence if found, or return "";
        String tailSeq;
        int tmdPosition;
        int tailLength;
        tmdPosition = searchTMD(proteinSeq);
        if (tmdPosition != -1) {
            tailLength = proteinSeq.length() - tmdPosition - tmdFrameLength;
            if (tailLength >= minKASHTailLength && tailLength <= maxKASHTailLength) {
                tailSeq = proteinSeq.substring(tmdPosition + tmdFrameLength);
                if (tailSeq.matches(regex)) {
                    return tailSeq;
                }
            }

        }
        return "";
    }

    private double getHydrophobicValue(char c) {
        /*hydrophobicity value is from 
         *A simple method for displaying the hydropathic character of a protein. Kyte J, Doolittle RF. J Mol Biol. 1982 May 5;157(1):105-32.
         *http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html#anote
         *a 20 amino acid fragment is normally considered as a TMD if the hydrophobicity value of this fragment is >20*1.6
         *which is the same method used in http://www.vivo.colostate.edu/molkit/hydropathy/index.html Protein Hydrophobicity Plots        
         */
        double value = 0;
        switch (c) {
            case 'I':
                value = 4.5;
                break;
            case 'V':
                value = 4.2;
                break;
            case 'L':
                value = 3.8;
                break;
            case 'F':
                value = 2.8;
                break;
            case 'C':
                value = 2.5;
                break;
            case 'M':
                value = 1.9;
                break;
            case 'A':
                value = 1.8;
                break;
            case 'G':
                value = -0.4;
                break;
            case 'T':
                value = -0.7;
                break;
            case 'S':
                value = -0.8;
                break;
            case 'W':
                value = -0.9;
                break;
            case 'Y':
                value = -1.3;
                break;
            case 'P':
                value = -1.6;
                break;
            case 'H':
                value = -3.2;
                break;
            case 'E':
                value = -3.5;
                break;
            case 'Q':
                value = -3.5;
                break;
            case 'D':
                value = -3.5;
                break;
            case 'N':
                value = -3.5;
                break;
            case 'K':
                value = -3.9;
                break;
            case 'R':
                value = -4.5;
                break;
        }
        return value;
    }

    private int searchTMD(String Sequence) {
        //if only 1 TMD is found in the sequence, return the first char position of a TMD; the first amino acid of the sequence is position 0;
        //otherwise return -1;
        char[] cs = Sequence.toCharArray();

        double maxRatio = 0; //maxium hydrophobic ratio of a TMD frame;
        double ratio=0; //current hydrophobic ratio of a TMD frame;
        int maxRatioIndex = -1; //the beginning position of the maximum hydrophobic frame;

        for (int i = 0; i < tmdFrameLength; i++) {
            ratio+=getHydrophobicValue(cs[i]);
        }    
        for (int i = tmdFrameLength; i < cs.length; i++) {
            ratio+=(getHydrophobicValue(cs[i])-getHydrophobicValue(cs[i-tmdFrameLength]));//the first amino acid M are not considered to be part of any TMD;
            if (ratio >= tmdHydrophobicThreshold) {
                if (maxRatioIndex != -1 && i - maxRatioIndex > 10) {
                    //maxRatioIndex!=-1 means this is the second high ratio found; 
                    //i-maxRatioIndex>10 means the two high ratio value is more than 10 amino acids away indicating multiple TMD or a long (>30 amino acids) hygrophobic seqence;
                    return -1;
                }
                if (maxRatio <= ratio) {//mark the highest hydrophobic position;
                    maxRatio = ratio;
                    maxRatioIndex = i;
                }
            }
        }
        return maxRatioIndex-tmdFrameLength;

    }

}
