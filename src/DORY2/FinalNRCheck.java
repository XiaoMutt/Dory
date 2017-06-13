/*
 *  require ncbi blast program to run
 *  require ncbi blast database (not in fasta format to run
 *  can be downloaded at ftp://ftp.ncbi.nlm.nih.gov/blast/db/
 *  or use given update_blastdb.pl script in the installed blast+ folder
 *  see ftp://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html for detail
 *
 */
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
public class FinalNRCheck extends SwingWorker<Void, String> {

    private final long homologyGroupCount;
    private final String readingFolder;
    private final String workingFolder;
    private final double eValueCutOff;
    private final String nrDB;

    private final int tmdFrameLength;// the transmembrane domain frame length;
    private final int tmdHydrophobicThreshold;//the Threshold to be a TMD;
    private final int maxKASHTailLength; // the maximum of the KASH domain length;
    private final int minKASHTailLength; // the minimum of the KASH domain length; 
    private final String regex;

    private final JTextArea logList;
    private final PrintWriter logTextPrinter;
    private long blastGroupCount;
    private Process proc;

    public FinalNRCheck(long HomologyGroupCount,
            String ReadingFolder,
            double EValueCutOff,
            String NRDB,
            int TmdFrameLength,
            int TmdHydrophobicThreshold,
            int MaxKASHTailLength,
            int MinKASHTailLength,
            String Regex,
            JTextArea LogList,
            PrintWriter LogTextPrinter) {
        homologyGroupCount = HomologyGroupCount;
        readingFolder = ReadingFolder;
        workingFolder = new File(readingFolder).getParent() + File.separator + "nr_blastp_homology_group";
        new File(workingFolder).mkdir();
        eValueCutOff = EValueCutOff;
        nrDB = new File(NRDB).getParent() + File.separator + "nr";
        tmdFrameLength = TmdFrameLength;
        tmdHydrophobicThreshold = TmdHydrophobicThreshold;
        maxKASHTailLength = MaxKASHTailLength;
        minKASHTailLength = MinKASHTailLength;
        regex = Regex;
        logList = LogList;
        logTextPrinter = LogTextPrinter;
        blastGroupCount = 0;
        proc = null;
    }

    @Override
    protected Void doInBackground() throws Exception {
        if (testBlastp()) {
            setProgress(0);
            publish("****************************************************");
            publish("   Final nr BLASTP started. Please wait patiently");
            publish("****************************************************");
            runBlastp();
            publish("***BLASTP finished, generating final report***");
            generateReport();
        } else {
            publish("Error: blast+ not installed");
        }

        return null;
    }

    private boolean testBlastp() throws IOException, InterruptedException {
        proc = Runtime.getRuntime().exec("blastp -version");
        BufferedReader reader
                = new BufferedReader(new InputStreamReader(proc.getInputStream()));
        proc.waitFor();
        String str;
        boolean res = false;

        publish("   Intalled blastp version:");
        while ((str = reader.readLine()) != null) {
            publish("\t" + str);
            if (str.contains("build")) {
                res = true;
            }
        }
        return res;
    }

    private void runBlastp() throws IOException, InterruptedException {
        Protein protein;
        FastaReader readProtein;
        PrintWriter tempfile;//file writer to save postitives to a file
        File fileHandler;//file handler to manipulate files such as renaming and deletion.
        String queryFile;
        String blastResultFile;
        String blastTempFile;
        int cores=Math.max(Runtime.getRuntime().availableProcessors()-1,2);

        for (long i = 0; i < homologyGroupCount && !isCancelled(); i++) {
            setProgress((int) (i * 100 / homologyGroupCount));
            readProtein = new FastaReader(readingFolder + File.separator + i + ".fa");

            if ((protein = readProtein.readNext()) != null) {

                //prepare query file for blast
                queryFile = workingFolder + File.separator + "query.fa";
                tempfile = new PrintWriter(new FileWriter(queryFile));
                tempfile.println(protein.name);
                tempfile.println(protein.sequence);
                tempfile.close();

                //get blast results to sequence IDs;
                blastTempFile = workingFolder + File.separator + i + ".ids";
                String command = "blastp -db " + nrDB + " -query " + queryFile + " -outfmt \"6 sseqid\" -evalue " + eValueCutOff + " -out " + blastTempFile+" -num_threads "+cores;
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
                    blastResultFile = workingFolder + File.separator + blastGroupCount + ".txt";
                    command = "blastdbcmd -db " + nrDB + " -entry_batch " + blastTempFile + " -outfmt \"%f\" -out " + blastResultFile;
                    proc = Runtime.getRuntime().exec(command);
                    proc.waitFor();
                }
                fileHandler.delete();
                blastGroupCount++;
            }
            readProtein.closeFile();
        }
        setProgress(100);
    }

    private void generateReport() throws IOException {
        PrintWriter report = new PrintWriter(new FileWriter(workingFolder + File.separator + "report.html"));
        report.println("<html><body><table>");
        long[] temp;
        report.println("<tr>");
        report.println("<td>|   File#</td>");
        report.println("<td>|   Number of Total Protein</td>");
        report.println("<td>|   Number of Protein with Positive KASH Tail</td>");
        report.println("<td>|   Positive%</td>");
        report.println("</tr>");
        for (int i = 0; i < blastGroupCount; i++) {
            temp = countKASH(workingFolder + File.separator + i + ".txt");
            report.println("<tr>");
            report.println("<td><a href=\"" + i + ".txt\">" + i + "</td>");
            report.println("<td>" + temp[1] + "</td>");
            report.println("<td>" + temp[0] + "</td>");
            report.println("<td>" + temp[0] * 100 / temp[1] + "</td>");
            report.println("</tr>");
        }
        report.println("</table></body></html>");
        report.close();
    }

    //return the ratio of positive KASH protein/total protein;
    private long[] countKASH(String FileName) throws IOException {
        FastaReader proteinReader = new FastaReader(FileName);
        KASHTailHunter kashTailHunter = new KASHTailHunter(tmdFrameLength,
                tmdHydrophobicThreshold,
                maxKASHTailLength,
                minKASHTailLength,
                regex);
        Protein protein;

        long positives = 0, total = 0;
        while ((protein = proteinReader.readNext()) != null) {
            if (kashTailHunter.searchKASH(protein.sequence).length() > 0) {
                positives++;
            }
            total++;
        }
        proteinReader.closeFile();
        return new long[]{positives, total};
    }

    @Override
    public void process(List<String> logs) {
        for (String log : logs) {
            logList.append(log);
            logList.append("\n");
            logTextPrinter.println(log);
        }
    }

    @Override
    public void done() {
        //output log
        if (isCancelled()) {
            if (proc != null) {
                proc.destroy();
            }
            publish("***Final BLASTing nr cancelled***");
        } else {
            publish(">>>report generated, please see report.html in the nr_blastp_homology_group folder ");
            publish("***Final BLASTing nr finished***");
        }

    }
}
