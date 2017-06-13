
package DORY2;

/**
 *
 * @author Xiao Zhou
 * Protein class
 * 
 */
public class Protein {
    public String name;
    public String sequence;
    
    public int length(){
        return sequence.length();
    }
    public Protein(String ProteinName, String ProteinSequence){
        name=ProteinName;
        sequence=ProteinSequence;
    }
    public Protein(){
        super();
    }

}
