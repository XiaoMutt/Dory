/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package DORY2;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.swing.JTextArea;
import javax.swing.SwingWorker;

/**
 *
 * @author Mutt
 */
public class KASHTailFilter extends SwingWorker<String, String> {

    private final int tmdFrameLength;// the transmembrane domain frame length;
    private final int tmdHydrophobicThreshold;//the Threshold to be a TMD;
    private final int maxKASHTailLength; // the maximum of the KASH domain length;
    private final int minKASHTailLength; // the minimum of the KASH domain length; 
    private final int proteinLengthFrom;//the minimum protein length considered for KASHFilter searching;
    private final int proteinLengthTo;//the maximum protein length considered for KASHFilter searching;
    private final String inputFileName;//input file name;
    private final String workingFolder; //the folder path of the Input File
    private final String speciesFileName; //the file name of the species names
    private final boolean outputKASHTail;//whether to output KASH tail into a file;
    private final boolean leftPadKASHTail;//whether to left pad KASH tail;
    private final String regex;//the Regex to use;
    private final boolean queryTax;//whether to query the NCBI Taxonomy Brower
    private final String eukaryoteCheckText; //the text used for check whether the species is eukaryotic;
    private final JTextArea logList;//the log list in the UI;
    private final PrintWriter logTextPrinter;//output log file named "DORYLog.txt"; 

    private ArrayList<String> speciesNames = null;//the list to hold species names;  
    private long db_length = 0; //total database length;
    private long db_num_seqs = 0; //total number of sequences;  
    private long candidateCount;
    private final String kashFilterResult;
    private final KASHTailHunter kashTailHunter;//class to search for KASH Tail;

    //private static final String seqW = "*ARNDCQEGHILKMFPSTWYVBZX"; // Amino acide order in the BLAST's scoring matrix (e.g.,Blosum62). 
    public KASHTailFilter(int TmdFrameLength,
            int TmdHydrophobicThreshold,
            int MaxKASHTailLength,
            int MinKASHTailLength,
            int ProteinLengthFrom,
            int ProteinLengthTo,
            String InputFileName,
            String SpeciesFileName,
            boolean OutputKASHTail,
            boolean LeftPadKASHTail,
            String Regex,
            boolean QueryTax,
            String EukaryoteCheckText,
            String WorkingFolder,
            JTextArea LogList,
            PrintWriter LogTextPrinter) throws IOException {
        tmdFrameLength = TmdFrameLength;
        tmdHydrophobicThreshold = TmdHydrophobicThreshold;
        maxKASHTailLength = MaxKASHTailLength;
        minKASHTailLength = MinKASHTailLength;
        proteinLengthFrom = ProteinLengthFrom;
        proteinLengthTo = ProteinLengthTo;
        inputFileName = InputFileName;
        speciesFileName = SpeciesFileName;
        outputKASHTail = OutputKASHTail;
        leftPadKASHTail = LeftPadKASHTail;
        regex = Regex;
        queryTax = QueryTax;
        eukaryoteCheckText = EukaryoteCheckText;
        logList = LogList;
        workingFolder = WorkingFolder;
        logTextPrinter = LogTextPrinter;
        candidateCount = 0;
        kashFilterResult = workingFolder + File.separator + "KASHFilterResult.txt";
        kashTailHunter = new KASHTailHunter(tmdFrameLength,
                tmdHydrophobicThreshold,
                maxKASHTailLength,
                minKASHTailLength,
                regex);

    }

    private class ProteinLoader implements Runnable {

        private final FastaReader proteinReader;
        private final BlockingQueue<Protein> proteinQ;

        public ProteinLoader(String InputFileName, BlockingQueue ProteinQ) throws IOException {
            proteinReader = new FastaReader(InputFileName);
            proteinQ = ProteinQ;
        }

        @Override
        public void run() {
            Protein protein;
            try {
                while ((protein = proteinReader.readNext()) != null && !isCancelled()) {
                    proteinQ.put(protein);
                    db_length += protein.length();//count database length;
                    db_num_seqs++;//count protein sequence;
                    setProgress(proteinReader.process);
                }
                proteinQ.put(new Protein("end", ""));
                proteinReader.closeFile();
            } catch (IOException | InterruptedException ex) {
                publish("Error :" + ex.getMessage());
            }
        }
    }

    private class ProteinChecker implements Runnable {

        private final BlockingQueue<Protein> inputQ;
        private final BlockingQueue<Protein> outputQ;
        private final BlockingQueue<Protein> kashTailQ;
        private final KASHTailHunter kashTailHunter;

        public ProteinChecker(BlockingQueue InputQ, BlockingQueue OutputQ, BlockingQueue KashTailQ) throws IOException {
            inputQ = InputQ;
            outputQ = OutputQ;
            kashTailQ = KashTailQ;
            kashTailHunter = new KASHTailHunter(tmdFrameLength,
                    tmdHydrophobicThreshold,
                    maxKASHTailLength,
                    minKASHTailLength,
                    regex);
        }

        @Override
        public void run() {
            Protein protein;
            String KASHSequence;
            try {
                while (!isCancelled()) {
                    protein = inputQ.take();
                    if ("end".equals(protein.name)) {
                        inputQ.put(protein);
                        outputQ.put(protein);
                        break;
                    }

                    if (checkSpecies(protein.name) && (protein.length() >= proteinLengthFrom && protein.length() <= proteinLengthTo)) {
                        KASHSequence = kashTailHunter.searchKASH(protein.sequence);
                        if (!KASHSequence.isEmpty()) {
                            //out put result
                            if (!queryTax || queryNCBITaxonomyBrowser(protein.name)) {
                                outputQ.put(protein);

                                if (outputKASHTail) {
                                    if (leftPadKASHTail) {
                                        kashTailQ.put(new Protein(protein.name, leftPad(KASHSequence)));
                                    } else {
                                        kashTailQ.put(new Protein(protein.name, KASHSequence));
                                    }
                                }
                            }
                        }
                    }
                }

            } catch (InterruptedException ex) {
                publish("Error :" + ex.getMessage());
            }
        }
    }

    @Override//return the filename of KASHFilterResult
    public String doInBackground() throws InterruptedException, IOException {
        if (!speciesFileName.isEmpty()) {
            loadSpecies();
            publish("   The protein names should contain the following speices:)");
            for (String species : speciesNames) {
                publish("       " + species);
            }
        } else {
            publish("   No species name checking during search.");
        }
        publish("---KASHTailFilter started---");
        filter();
        return kashFilterResult;
    }

    @Override
    public void done() {
        //output log
        publish(">>>KASHFilter searched " + db_num_seqs + " protein(s) (total length " + db_length + " amino acids).");
        publish(">>>" + candidateCount + " candidate(s) were found.");
        publish(">>>Please see results in KASHFilterResult.txt.");
        if (isCancelled()) {
            publish("***KASHFilter cancelled***");
        } else {
            publish("***KASHFilter finished***");
        }

    }

    @Override
    public void process(List<String> logs) {
        for (String log : logs) {
            logList.append(log);
            logList.append("\n");
            logTextPrinter.println(log);
        }
    }

    private void filter() throws IOException, InterruptedException {
        Protein protein, kashTail;
        PrintWriter positivePrintWriter = new PrintWriter(new FileWriter(kashFilterResult)); //Output the postiive proteins to a file
        PrintWriter KASHTailOut; //Output The KASH domain to a file
        int capacity = 1048576;
        BlockingQueue<Protein> inputQ = new LinkedBlockingQueue(capacity);
        BlockingQueue<Protein> outputQ = new LinkedBlockingQueue(capacity);
        BlockingQueue<Protein> kashTailQ = null;

        if (outputKASHTail) {
            KASHTailOut = new PrintWriter(new FileWriter(workingFolder + File.separator + "KASHTail.txt"));
            kashTailQ = new LinkedBlockingQueue(capacity);
        } else {
            KASHTailOut = null;
        }

        setProgress(0);
        //load proteins using a thread;
        Thread proteinLoader = new Thread(new ProteinLoader(inputFileName, inputQ));
        proteinLoader.start();

        int cores = Runtime.getRuntime().availableProcessors();//number of threads = number of cpu cores; 
        Thread[] proteinChecker = new Thread[cores];
        for (int i = 0; i < cores; i++) {
            proteinChecker[i] = new Thread(new ProteinChecker(inputQ, outputQ, kashTailQ));
            proteinChecker[i].start();
        }

        while (!isCancelled() && cores > 0) {
            protein = outputQ.take();
            if ("end".equals(protein.name)) {//indicates proteinChecker threads end;
                cores--;
            } else {
                positivePrintWriter.println(protein.name);
                positivePrintWriter.println(protein.sequence);
                if (kashTailQ != null && KASHTailOut != null && (kashTail = kashTailQ.take()) != null) {
                    KASHTailOut.println(kashTail.name);
                    KASHTailOut.println(kashTail.sequence);
                }
                candidateCount++;
            }
        }
        setProgress(100);
        positivePrintWriter.close();
        if (KASHTailOut != null) {
            KASHTailOut.close();
        }
    }

    private boolean checkSpecies(String inputName) {
        if (speciesNames == null) {
            return true;
        } else {
            for (String name : speciesNames) {
                if (inputName.contains(name)) {
                    return true;
                }
            }
            return false;
        }
    }

    private void loadSpecies() {
        BufferedReader fileIn;
        String strLine;
        speciesNames = new ArrayList<>();
        try {
            fileIn = new BufferedReader(new FileReader(speciesFileName));
            while ((strLine = fileIn.readLine()) != null && !strLine.isEmpty()) {
                speciesNames.add(strLine);
            }
            fileIn.close();
        } catch (IOException e) {
            publish("Error during loading species names: " + e.getMessage());
        }
    }

    //return ture if the species is eukaryote or if the species name is not found or error happens.
    //return false if the speices is not found.
    private boolean queryNCBITaxonomyBrowser(String proteinName) {
        String speciesName;
        Pattern pattern = Pattern.compile("\\[(.+?)\\]");
        Matcher matcher = pattern.matcher(proteinName);
        if (matcher.find()) {
            speciesName = matcher.group(1);
        } else {
            return true;
        }
        try {
            String query = "name=" + URLEncoder.encode(speciesName, "UTF-8");
            URLConnection connection = new URL("http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi" + "?" + query).openConnection();
            try (BufferedReader in = new BufferedReader(new InputStreamReader(
                    connection.getInputStream()))) {
                String inputLine;
                while ((inputLine = in.readLine()) != null) {
                    if (inputLine.contains(eukaryoteCheckText)) {
                        in.close();
                        return true;
                    }
                }
            }
        } catch (IOException E) {
            publish("Error: " + E.getMessage());
            return true;
        }
        return false;
    }

    private String leftPad(String Str) {
        int num;
        num = maxKASHTailLength - Str.length();
        char[] filling = new char[num];
        for (int i = 0; i < filling.length; i++) {
            filling[i] = '-';
        }
        return new String(filling) + Str;
    }

}
