/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package DORY2;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.concurrent.ExecutionException;
import javax.swing.JTextArea;
import javax.swing.SwingWorker;

/**
 *
 * @author Mutt
 */
public class DoryEngine extends SwingWorker<Void, String> {

    private final int runType;//the run type;    
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

    private final double eValueCutOff;

    private final boolean finalCheckNR;
    private final String nrDB;
    private final JTextArea logList;//the log list in the UI;
    private final PrintWriter logTextPrinter;//output log file named "DORYLog.txt"; 
    private final String searchTag;//the random unique tag for this run;
    private KASHTailFilter ktf = null;
    private HomologyFilter hf = null;
    private FinalNRCheck nrf = null;

    public DoryEngine(int RunType,
            int TmdFrameLength,
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
            double EValueCutOff,
            boolean FinalCheckNR,
            String NRDB,
            String WorkingFolder,
            JTextArea LogList) throws IOException {
        runType = RunType;
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
        eValueCutOff = EValueCutOff;
        finalCheckNR = FinalCheckNR;
        nrDB = NRDB;
        searchTag = java.util.UUID.randomUUID().toString();
        workingFolder = WorkingFolder + File.separator + "DORY_" + searchTag;
        //making working folder;
        new File(workingFolder).mkdir();
        //initiate logTextPrinter;
        logTextPrinter = new PrintWriter(new FileWriter(workingFolder + File.separator + "DORYLog" + ".txt"));
    }

    private void createFinalNRCheckWorker() throws InterruptedException, ExecutionException {
        nrf = new FinalNRCheck(hf.get(),
                workingFolder + File.separator + "homology_filter_result",
                eValueCutOff,
                nrDB,
                tmdFrameLength,
                tmdHydrophobicThreshold,
                maxKASHTailLength,
                minKASHTailLength,
                regex,
                logList,
                logTextPrinter);
        nrf.addPropertyChangeListener(new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                if ("progress".equals(evt.getPropertyName())) {
                    setProgress((Integer) evt.getNewValue());
                }
            }

        });
        nrf.execute();
        nrf.get();

    }

    private void createHomologyFilterWorker() throws IOException, InterruptedException, ExecutionException {
        hf = new HomologyFilter(ktf.get(),
                workingFolder,
                eValueCutOff,
                logList,
                logTextPrinter);

        hf.addPropertyChangeListener(new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                if ("progress".equals(evt.getPropertyName())) {
                    setProgress((Integer) evt.getNewValue());
                }
            }
        });
        hf.execute();
        hf.get();

    }

    private void createKASHTailFilterWorker() throws IOException, InterruptedException, ExecutionException {
        ktf = new KASHTailFilter(tmdFrameLength,
                tmdHydrophobicThreshold,
                maxKASHTailLength,
                minKASHTailLength,
                proteinLengthFrom,
                proteinLengthTo,
                inputFileName,
                speciesFileName,
                outputKASHTail,
                leftPadKASHTail,
                regex,
                queryTax,
                eukaryoteCheckText,
                workingFolder,
                logList,
                logTextPrinter
        );
        ktf.addPropertyChangeListener(new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                if ("progress".equals(evt.getPropertyName())) {
                    setProgress((Integer) evt.getNewValue());
                }
            }
        });
        ktf.execute();
        ktf.get();
    }

    @Override

    protected Void doInBackground() {

        //run search
        try {
            publish("");
            publish("***DORY2 Started with the following settings:***");
            publish("   Search tag: " + searchTag);
            publish("   Please see results in the folder: " + workingFolder);
            publish("   input fasta file: " + inputFileName);
            publish("   tmdFrameLength = " + tmdFrameLength);
            publish("   tmdHydrophobicThreshold = " + tmdHydrophobicThreshold);
            publish("   maxKASHTailLength = " + maxKASHTailLength);
            publish("   minKASHTailLength = " + minKASHTailLength);
            publish("   proteins with a length from " + proteinLengthFrom + " to " + proteinLengthTo + " amino acids are taken into consideration");
            publish("   output KASH in a file = " + outputKASHTail + (leftPadKASHTail ? " with" : " without") + " leftpadding to the maxKASHTailLength");
            publish("   Regex for the KASH tail =\"" + regex + "\"");
            publish("   query NCBI Taxonomy Browser = " + queryTax + ", and positive when return text contains \"" + eukaryoteCheckText + "\"");
            publish("   eValueCutOff for homologous proteins" + eValueCutOff);

            if (finalCheckNR) {
                publish("   Perform a final BLASTP for identified KASH candidates agains nr database which is at " + nrDB);
            } else {
                publish(" Do not perform a final BLASTP for identified KASH candidates agains nr database");
            }

            switch (runType) {
                case 0:
                    publish("---run type = \"Run Full Search\"---");
                    createKASHTailFilterWorker();

                    createHomologyFilterWorker();
                    if (finalCheckNR) {
                        createFinalNRCheckWorker();
                    }
                    break;
                case 1:
                    publish("---run type = \"Run KASHFilter only\"---");
                    createKASHTailFilterWorker();
                    break;
                case 2:
                    publish("---run type = \"Run HomologyFilter only\"---");
                    createHomologyFilterWorker();
                    if (finalCheckNR) {
                        createFinalNRCheckWorker();
                    }
                    break;
                default:
                    break;
            }
        } catch (IOException | InterruptedException | ExecutionException E) {
            publish("Error during search:" + E.getMessage());
        }

        return null;

    }

    @Override
    public void process(List<String> logs) {
        for (String log : logs) {
            logList.append(log);
            logList.append("\n");
            logTextPrinter.println(log);
        }
        if (isDone() || isCancelled()) {
            logTextPrinter.close();
        }
    }

    @Override
    public void done() {
        if (isCancelled()) {
            if (ktf != null) {
                publish("---Waiting for KASHTailFilter to cancel---");
                ktf.cancel(true);
            }
            if (hf != null) {
                publish("---Waiting for HomologyFilter to cancel---");
                hf.cancel(true);
            }
            if (nrf != null) {
                publish("---Waiting for final nr database BLAST to cancel---");
                nrf.cancel(true);
            }
            publish("***DORY Search cancelled***");
        } else {
            publish("***DORY Search finished***");
        }
        publish("");

    }

}
