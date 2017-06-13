package DORY2;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.List;
import javax.swing.JTextArea;
import javax.swing.SwingWorker;

/**
 *
 * @author Mutt
 */
public class HomologyFilter extends SwingWorker<Long, String> {

    private final String inputFileName;
    private final String workingFolder;
    private final double eValueCutOff;
    private final JTextArea logList;//the log list in the UI;
    private final PrintWriter logTextPrinter;//output log file named "DORYLog.txt"; 
    private long homologyGroupCount;
    private Process proc;

    public HomologyFilter(String InputFileName,
            String WorkingFolder,
            double EValueCutOff,
            JTextArea LogList,
            PrintWriter LogTextPrinter) throws IOException, InterruptedException {
        inputFileName = InputFileName;
        File tempfile = new File(inputFileName);
        workingFolder = tempfile.getParent() + File.separator + "homology_filter_result";
        new File(workingFolder).mkdir();
        eValueCutOff = EValueCutOff;
        logList = LogList;
        logTextPrinter = LogTextPrinter;
        homologyGroupCount = 0;
        proc = null;
    }

    private boolean testBlastp() throws IOException, InterruptedException {
        proc = Runtime.getRuntime().exec("blastp -version");
        BufferedReader reader
                = new BufferedReader(new InputStreamReader(proc.getInputStream()));
        proc.waitFor();
        String str;
        boolean res = false;
        publish("   Installed blastp version:");
        while ((str = reader.readLine()) != null) {
            publish("\t" + str);
            if (str.contains("build")) {
                res = true;
            }
        }
        return res;
    }
    private boolean testMakeBlastDB() throws IOException, InterruptedException {
        proc = Runtime.getRuntime().exec("makeblastdb -version");
        BufferedReader reader
                = new BufferedReader(new InputStreamReader(proc.getInputStream()));
        proc.waitFor();
        String str;
        boolean res = false;
        publish("   Installed makeblastdb version:");
        while ((str = reader.readLine()) != null) {
            publish("\t" + str);
            if (str.contains("build")) {
                res = true;
            }
        }
        return res;
    }    

    private void prepareBlastDB() throws IOException, InterruptedException {
        //get blast results to sequence IDs;
        String command = "makeblastdb -in " + inputFileName + " -parse_seqids -dbtype prot";
        proc = Runtime.getRuntime().exec(command);
        proc.waitFor();
    }

    private boolean checkIfAlreadyGrouped(Protein protein) throws IOException {
        FastaReader readProtein;
        Protein tempProtein;
        for (int i = 0; i < homologyGroupCount; i++) {
            //search duplicates;
            readProtein = new FastaReader(workingFolder + File.separator + i + ".fa");
            while ((tempProtein = readProtein.readNext()) != null) {
                if (protein.sequence.equals(tempProtein.sequence)) {
                    readProtein.closeFile();
                    return true;
                }
            }
            readProtein.closeFile();
        }
        return false;
    }

    private void filter() throws IOException, InterruptedException {
        Protein protein;
        PrintWriter tempfile;//file writer to save postitives to a file
        File fileHandler;//file handler to manipulate files such as renaming and deletion.
        String queryFile;
        String blastResultFile;
        String blastTempFile;

        int cores=Math.max(Runtime.getRuntime().availableProcessors()-1,2);
        FastaReader readProtein = new FastaReader(inputFileName);

        //set progressbar;
        setProgress(0);

        while ((protein = readProtein.readNext()) != null && !isCancelled()) {

            if (checkIfAlreadyGrouped(protein)) {
                continue;
            }

            //prepare query file for blast
            queryFile = workingFolder + File.separator + "query.fa";
            tempfile = new PrintWriter(new FileWriter(queryFile));
            tempfile.println(protein.name);
            tempfile.println(protein.sequence);
            tempfile.close();

            //get blast results to sequence IDs;
            blastTempFile = workingFolder + File.separator + homologyGroupCount + ".ids";
            String command = "blastp -db " + inputFileName + " -query " + queryFile + " -outfmt \"6 sseqid\" -evalue " + eValueCutOff + " -out " + blastTempFile+" -num_threads "+cores;
            proc = Runtime.getRuntime().exec(command);
            proc.waitFor();
            //delete temporary query file;
            fileHandler = new File(queryFile);
            fileHandler.delete();

            //get rid of any duplicates, not necessary for nr database;
            //stripDuplicates(blastTempFile);
            //get sequences from 
            fileHandler = new File(blastTempFile);
            if (fileHandler.length() != 0) {
                blastResultFile = workingFolder + File.separator + homologyGroupCount + ".fa";
                command = "blastdbcmd -db " + inputFileName + " -entry_batch " + blastTempFile + " -outfmt \"%f\" -out " + blastResultFile;
                proc = Runtime.getRuntime().exec(command);
                proc.waitFor();

                homologyGroupCount++;
            }
            fileHandler.delete();
            setProgress(readProtein.process);
        }
        setProgress(100);
        readProtein.closeFile();
    }

    @Override//return the count of homology groups
    public Long doInBackground() throws InterruptedException, IOException {
        if (testMakeBlastDB()&&testBlastp()) {
            setProgress(0);
            publish("---HomologyFilter started. Please wait patiently---");
            prepareBlastDB();
            filter();
        } else {
            publish("Error: blast+ not installed");
        }
        return homologyGroupCount;
    }

    @Override
    public void done() {

        if (isCancelled()) {
            if (proc != null) {
                proc.destroy();
            }
            publish("***HomoloyFilter cancelled***");
        } else {
            publish(">>>Please see results in homology_filter_result folder. Each group is numbered");
            publish("***HomologyFilter finished***");
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
}
