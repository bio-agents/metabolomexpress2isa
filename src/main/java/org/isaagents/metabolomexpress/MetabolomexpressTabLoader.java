package org.isaagents.metabolomexpress;


import com.sun.agents.javac.util.Pair;
import org.apache.commons.lang.StringUtils;
import org.isaagents.io.FileType;
import org.isaagents.io.Loader;
import org.isaagents.magetoisatab.io.fileprocessing.CleanupUtils;
import org.isaagents.magetoisatab.io.fileprocessing.RemoveColumnUtil;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class MetabolomexpressTabLoader {

    public static final Character TAB_DELIM = '\t';

    public String[] columnNames;

    public Map<Integer, String[]> rowValues;

    public Map<Integer, String[]> samples;
    private String[] cleanedHeader;

    public MetabolomexpressTabLoader() {
        samples = new HashMap<Integer, String[]>();
    }

    public void loadMetabolomexpressTab(String url, String metabolomexpressFilename) throws IOException {

        try {

            File file = new File(url);

            System.out.println("URL: " + url);
            //this is to get the basename of the metabolomexpress file
            String[] tokens = metabolomexpressFilename.split("\\.(?=[^\\.]+$)");

            //trying to create a directory to host all files resulting from the conversion using the basename of the file as input to name the directory
            boolean success = (new File(DownloadUtils.CONVERTED_DIRECTORY + File.separator + tokens[0])).mkdirs();


            if (success) {
                System.out.println("Directory: " + tokens[0] + " created");
            }

            HashMap<String, String> studyAttributes = new HashMap<String, String>();

            studyAttributes.put("studyName", null);
            studyAttributes.put("studyDescription", null);
            studyAttributes.put("studyPMID", null);
            studyAttributes.put("studyPubAuthors", null);
            studyAttributes.put("studyPubTitle", null);
            studyAttributes.put("journal", null);

            String personFirstNames = "";
            String personLastNames = "";
            String studyName = "";
            String studyDescription = "";
            String studyPMID = "";
            String studyPubAuthors = "";
            String studyPubTitle = "";
            String journal = "";
            String studyDate = "";

            String PersonNames = "";
            String PersonMails = "";

            List<Map<String, List<String>>> protocolsByTypes = new ArrayList<Map<String, List<String>>>();

            List<String[]> samplePrepProtocol = new ArrayList<String[]>();
            List<String[]> extractionPrepProtocol = new ArrayList<String[]>();
            List<String[]> analyticalPrepProtocol = new ArrayList<String[]>();


            List<String[]> studyGenotype = new ArrayList<String[]>();
            List<String[]> studyEnvironment = new ArrayList<String[]>();
            List<String[]> studyTreatment = new ArrayList<String[]>();
            List<String[]> studyExplant = new ArrayList<String[]>();
            List<String[]> studyBiosource = new ArrayList<String[]>();
            List<String[]> studyExtract = new ArrayList<String[]>();
            List<String[]> studyAnalytical = new ArrayList<String[]>();
            List<String[]> studyRun = new ArrayList<String[]>();

            if (file.exists()) {

                //reading in the file if it exists
                Loader fileReader = new Loader();
                List<String[]> sheetData = fileReader.loadSheet(url, FileType.TAB);

                int iter = 0;

                int biosampleprep = 0;
                int biosampleprepEnd = 0;

                int extractionprep = 0;
                int extractionprepEnd = 0;

                int analyticalprep = 0;
                int analyticalprepEnd = 0;

                int explantBegin = 0;
                int explantEnd = 0;

                int biosourceBegin = 0;
                int extractBegin = 0;

                int analyticalBegin = 0;
                int analyticalEnd = 0;

                int runBegin = 0;

                int genotypeBegin = 0;
                int genotypeEnd = 0;

                int envBegin = 0;
                int envEnd = 0;

                int treatmentBegin = 0;
                int treatmentEnd = 0;


                String gcmsProtocolDesc = null;
                String metaboliteIdProtocolDesc = null;
                for (String[] row : sheetData) {

                    if (row[0].contains("Experiment Name:")) {
                        row[0] = row[0].toString().replace("Experiment Name:", "Study Name");
                        studyName = row[1];
                    }
//                    if (row.toString().startsWith("Biological Experimentalist Name:")) {
//                         studyName=(row.toString().replace("Biological Experimentalist Name:","Study Person Name");
//                    }
                    if (row[0].startsWith("Experimental Hypothesis:")) {
                        studyDescription = row[1];
                    }
                    if (row[0].startsWith("Brief Description of Experiment:")) {
                        studyDescription = studyDescription + row[1];
                    }
                    if (row[0].startsWith("Literature Reference:") && !(row[1].contains("?"))) {
                        List<String> temp = getPublicationTitleAndAuthors(row[1]);
                        // if (row[1].contains("(")) {

                        System.out.println("TEMP" + temp.get(0) + "|" + temp.get(1));
                        studyPubTitle = temp.get(1);
                        studyPubAuthors = temp.get(0);

                        // }


                    }
                    if (row[0].startsWith("Journal:") && !(row[1].contains("?"))) {
                        journal = row[1];
                    }

                    if (row[0].startsWith("Publication Date")&& !(row[1].contains("?"))) {
                        studyDate = row[1];
                    }
                    if (row[0].startsWith("PubMed")&& !(row[1].contains("?"))) {
                        studyPMID = row[1];
                    }

                    studyAttributes.put("studyName", studyName);
                    studyAttributes.put("studyDescription", studyDescription);
                    studyAttributes.put("studyPMID", studyPMID);
                    studyAttributes.put("studyPubAuthors", studyPubAuthors);
                    studyAttributes.put("studyPubTitle", studyPubTitle);

                    studyAttributes.put("studyDate", studyDate);
                    studyAttributes.put("journal", journal);


                    if (row[0].startsWith("Biological Experimentalist Name") && !(row[1].contains("?"))) {
                        Pair<List<String>, List<String>> names = checkFirstLastNames(row[1]);
                        personFirstNames = names.fst.get(0) + "\t";
                        personLastNames = names.snd.get(0) + "\t";

                        System.out.println("First Name: " + personFirstNames + "Last Name: " + personLastNames);
                    }
                    if (row[0].startsWith("Biological Experimentalist Email") && !(row[1].contains("?")) ) {
                        PersonMails = row[1] + "\t";
                    }
                    if (row[0].startsWith("Metabolome Analyst Name") && !(row[1].contains("?"))) {
                        Pair<List<String>, List<String>> names = checkFirstLastNames(row[1]);
                        // PersonNames = PersonNames+ row[1] + "\t";
                        personFirstNames = personFirstNames + names.fst.get(0);
                        personLastNames = personLastNames + names.snd.get(0);

                    }
                    if (row[0].startsWith("Metabolome Analyst Email") && !(row[1].contains("?"))) {
                        PersonMails = PersonMails + row[1];
                    }


                    //pulling protocols, protocol descriptions and protocol type
                    String protocolType = null;
                    String protocolDesc = null;

                    gcmsProtocolDesc = "";
                    metaboliteIdProtocolDesc = "";

                    if (row[0].startsWith("BioSample Processing Protocol ID")) {
                        biosampleprep = iter;
                        protocolType = "sample preprocessing";
                        protocolDesc = row[1];
                        System.out.println("sampleprocesStart: " + biosampleprep);
                    }

                    if (row[0].startsWith("Genotype ID")) {
                        genotypeBegin = iter;
                    }
                    if (row[0].startsWith("[ENVIRONMENTS]")) {
                        envEnd = iter;
                    }
                    if (row[0].startsWith("Environment ID")) {
                        envBegin = iter;
                    }
                    if (row[0].startsWith("[TREATMENTS]")) {
                        envEnd = iter;
                    }

                    if (row[0].startsWith("Treatment ID")) {
                        treatmentBegin = iter;
                    }
                    if (row[0].startsWith("[COLLECTIONS]")) {
                        treatmentEnd = iter;
                    }

                    if (row[0].contains("BIOSAMPLE EXTRACTION PROTOCOLS")) {
                        biosampleprepEnd = iter;
                        System.out.println("sampleprocesEnd: " + biosampleprepEnd);
                    }

                    if (row[0].startsWith("Extraction Protocol ID")) {
                        extractionprep = iter;
                        protocolType = "material extraction";
                        protocolDesc = row[1];
                        System.out.println("extractprocesStart: " + extractionprep);

                    }

                    if (row[0].contains("[ANALYTICAL SAMPLE PREPARATION PROTOCOLS]")) {
                        analyticalprep = iter;
                        extractionprepEnd = iter;
                        System.out.println("extractionprepEnd: " + extractionprepEnd);
                    }


                    if (row[0].startsWith("Analytical Sample Preparation Protocol ID")) {
                        protocolType = "analytical sample preparation";
                        protocolDesc = row[1];

                    }
                    if (row[0].startsWith("File Formats")) {
                        protocolType = "data preprocessing";
                        protocolDesc = row[1];
                        System.out.println("File formats: " + row.toString());

                    }

                    if (row[0].startsWith("Raw Data Pre-Processing")) {
                        protocolDesc = protocolDesc + row[1];

                    }
                    if (row[0].startsWith("Data Transformation")) {
                        protocolType = "data transformation";
                        protocolDesc = row[1];

                    }
                    if (row[0].startsWith("Statistics")) {
                        protocolType = "data transformation";
                        protocolDesc = row[1];

                    }


                    //locating where blocks begin and if finding Header rows and formatting them as ISA compatible header
                    if (row[0].startsWith("Explant ID")) {
                        //System.out.println(iter+"|"+Arrays.toString(row));
                        explantBegin = iter;
                        cleanHeader(row);
                        studyExplant.add(row);
                    }

                    if (row[0].startsWith("BioSample ID")) {
                        System.out.println(iter+"|"+Arrays.toString(row));
                        analyticalprepEnd = iter - 2;
                        biosourceBegin = iter;
                        cleanHeader(row);
                        studyBiosource.add(row);
                        System.out.println(iter+"|"+Arrays.toString(row));
                    }

                    if (row[0].startsWith("Extract ID")) {
                        //System.out.println(iter+"|"+Arrays.toString(row));
                        extractBegin = iter;
                        cleanHeader(row);
                        studyExtract.add(row);
                        System.out.println(iter + "**|****" + Arrays.toString(row));
                    }

                    if (row[0].startsWith("Analytical Sample ID")) {
                        //System.out.println(iter+"|"+Arrays.toString(row));
                        analyticalBegin = iter;
                        cleanHeader(row);
                        studyAnalytical.add(row);
                        System.out.println(iter + "|" + Arrays.toString(row));
                    }

                    if (row[0].startsWith("Run ID")) {
                        // System.out.println(iter+"|"+Arrays.toString(row));
                        runBegin = iter;
                        cleanHeader(row);
                        studyRun.add(row);
                        // System.out.println(iter+"|"+Arrays.toString(row));
                    }

                    //locating where record blocks end
                    if (row[0].startsWith("**********CHEMICAL ANALYSIS")) {
                        //System.out.println(iter+"|"+Arrays.toString(row));
                        explantEnd = iter;
                    }


                    HashMap<String, String> gcmsParamValues = new HashMap<String, String>();

                    if (row[0].startsWith("Auto-injector")) {
                        gcmsParamValues.put("auto-injector", row[1]);
                        gcmsProtocolDesc = row[1];

                    }
                    if (row[0].startsWith("Chromatography Instrument:")) {
                        gcmsParamValues.put("chromatography instrument", row[1]);
                        gcmsProtocolDesc = gcmsProtocolDesc + row[1];
                    }
                    if (row[0].startsWith("Instrument Control Software:")) {
                        gcmsParamValues.put("instrument control software", row[1]);
                        gcmsProtocolDesc = gcmsProtocolDesc + row[1];
                    }
                    if (row[0].startsWith("Separation Column:")) {
                        gcmsParamValues.put("chromatography column", row[1]);
                        gcmsProtocolDesc = gcmsProtocolDesc + row[1];
                    }
                    if (row[0].startsWith("Mass Spectrometry Instrument:")) {
                        gcmsParamValues.put("mass spectrometry instrument", row[1]);
                        gcmsProtocolDesc = gcmsProtocolDesc + row[1];
                    }
                    if (row[0].startsWith("Sample Introduction (MS)")) {
                        gcmsParamValues.put("mass spectrometry sample introduction", row[1]);
                        gcmsProtocolDesc = gcmsProtocolDesc + row[1];
                    }
                    if (row[0].startsWith("Ionization:")) {
                        gcmsParamValues.put("ion source", row[1]);
                        gcmsProtocolDesc = gcmsProtocolDesc + row[1];
                    }
                    if (row[0].startsWith("Data Acquisition:")) {
                        gcmsParamValues.put("detector", row[1]);
                        gcmsProtocolDesc = gcmsProtocolDesc + row[1];
                    }


                    if (row[0].startsWith("**********INSTRUMENTAL ANALYSIS")) {
                        analyticalEnd = iter;
                    }


                    int dataProcessing = 0;
                    if (row.toString().startsWith("**********DATA PROC")) {
                        dataProcessing = iter;
                    }
                    int metaboliteIdentBegin = 0;
                    if (row.toString().startsWith("**********METABOLITE IDENTIFICATION METAD")) {
                        metaboliteIdentBegin = iter;
                    }


                    if (row[0].startsWith("First Identification Parameter:")) {
                        metaboliteIdProtocolDesc = row[1];

                    }
                    if (row[0].startsWith("Second Identification Parameter:")) {
                        metaboliteIdProtocolDesc = metaboliteIdProtocolDesc + row[1];
                    }
                    if (row[0].startsWith("Description of Library Matching Algorithm Logic:")) {
                        metaboliteIdProtocolDesc = metaboliteIdProtocolDesc + row[1];
                    }
                    if (row[0].startsWith("Metabolite Library Description:")) {
                        metaboliteIdProtocolDesc = metaboliteIdProtocolDesc + row[1];
                    }
                    if (row[0].startsWith("**********INSTRUMENTAL ANALYSIS")) {
                        analyticalEnd = iter;
                    }


                    int metaboliteIdentEnd = 0;
                    if (row.toString().startsWith("***********QUALITY CONTROL")) {
                        metaboliteIdentBegin = iter;
                    }


                    //TODO create a method taking a block as input and return a protocolsbyType object
                    if (iter > biosampleprep && biosampleprep > 0 && biosampleprepEnd == 0) {
                        Map<String, List<String>> protocol = new HashMap<String, List<String>>();
                        List<String> protocolattributes = new ArrayList<String>();
                        if (row.toString() != null) {
                            samplePrepProtocol.add(row);
                            String sampleprepProtDesc = "";
                            protocolattributes.add(row[0]);
                            for (int j = 1; j < row.length; j++) {
                                if (row[j] != null) {
                                    sampleprepProtDesc = sampleprepProtDesc + row[j] + " ";
                                }
                            }
                            protocolattributes.add(sampleprepProtDesc);
                            protocol.put("sample preparation", protocolattributes);
                            protocolsByTypes.add(protocol);

                        }
                    }
                    //TODO create a method taking a block as input and return a protocolsbyType object
                    if (iter > extractionprep && extractionprep > 0 && extractionprepEnd == 0) {
                        Map<String, List<String>> protocol = new HashMap<String, List<String>>();
                        List<String> protocolattributes = new ArrayList<String>();
                        if (row.toString() != null) {
                            extractionPrepProtocol.add(row);
                            String extractionprepProtDesc = "";
                            protocolattributes.add(row[0]);
                            for (int j = 1; j < row.length; j++) {
                                if (row[j] != null) {
                                    extractionprepProtDesc = extractionprepProtDesc + row[j] + " ";
                                }
                            }
                            protocolattributes.add(extractionprepProtDesc);
                            protocol.put("extraction", protocolattributes);
                            protocolsByTypes.add(protocol);

                        }
                    }

                    //TODO create a method taking a block as input and return a protocolsbyType object
                    if (iter > analyticalprep && analyticalprep > 0 && analyticalprepEnd == 0) {
                        Map<String, List<String>> protocol = new HashMap<String, List<String>>();
                        List<String> protocolattributes = new ArrayList<String>();
                        if (!row.toString().equals(null) || !row.toString().contains("Description")) {

                            analyticalPrepProtocol.add(row);
                            String analyticalprepProtDesc = "";
                            protocolattributes.add(row[0]);
                            for (int j = 1; j < row.length; j++) {
                                if ((row[j] != null)) {

                                    analyticalprepProtDesc = analyticalprepProtDesc + row[j] + " ";
                                }
                            }
                            protocolattributes.add(analyticalprepProtDesc);
                            protocol.put("analytical sample preparation", protocolattributes);
                            protocolsByTypes.add(protocol);
                        }
                    }


                    //adding data records  to relevant arraylists

                    if (iter > genotypeBegin && genotypeBegin > 0 && genotypeEnd == 0) {
                        if (row.toString() != null) {
                            studyGenotype.add(row);
                            System.out.println(Arrays.toString(row));
                        }
                    }
                    if (iter > envBegin && envBegin > 0 && envEnd == 0) {
                        if (row.toString() != null) {
                            studyEnvironment.add(row);
                            System.out.println(Arrays.toString(row));
                        }
                    }
                    if (iter > treatmentBegin && treatmentBegin > 0 && treatmentEnd == 0) {
                        if (row.toString() != null) {
                            studyTreatment.add(row);
                            System.out.println(Arrays.toString(row));
                        }
                    }


                    //explant records
                    if (iter > explantBegin && explantBegin > 0 && explantEnd == 0) {
                        if (StringUtils.isBlank(row.toString())) {
                            System.out.println("Blank!" + Arrays.toString(row));
                        }
                        else if (row.toString() != null && StringUtils.isNotBlank(row.toString())) {
                            studyExplant.add(row);
                            System.out.println("That "+Arrays.toString(row));
                        }
                    }
                    //biosource records
                    if (iter > biosourceBegin && biosourceBegin > 0 && extractBegin == 0) {
                        if (StringUtils.isBlank(row.toString())) {
                            System.out.println("Blank!" + Arrays.toString(row));
                        }
                        else if (row.toString() != null && StringUtils.isNotBlank(row.toString())) {

                            studyBiosource.add(row);
                            System.out.println("this " + Arrays.toString(row));
                        }
                    }
                    //extract records
                    if (iter > extractBegin && extractBegin > 0 && analyticalBegin == 0) {
                        if (row.toString() != null) {
                            studyExtract.add(row);

                            System.out.println("extract: " + Arrays.toString(row));
                        }
                    }
                    //analytical sample records
                    if (iter > analyticalBegin && analyticalBegin > 0 && analyticalEnd == 0) {
                        if (row.toString() != null) {
                            studyAnalytical.add(row);
                            System.out.println("analytical: " + Arrays.toString(row));
                        }
                    }
                    //analytical run records
                    if (iter > runBegin && runBegin > 0 && metaboliteIdentBegin == 0) {
                        if (row.toString() != null) {
                            studyRun.add(row);
                            System.out.println("run: " + Arrays.toString(row));
                        }
                    }

                    iter++;
                }

                List<String[]> merged4Study = stitchStudyRecords(studyBiosource, studyExplant);

                for (int columnIndex = 0; columnIndex < merged4Study.get(0).length; columnIndex++) {
                    if (merged4Study.get(0)[columnIndex] == null) {
                        merged4Study.get(0)[columnIndex] = "";
                    }
                }


                System.out.println(Arrays.toString(merged4Study.get(0)));
                
                CleanupUtils removeColumns = new RemoveColumnUtil("",null);
                merged4Study = removeColumns.processSpreadsheet(merged4Study);

                List<String[]> merged4assays = stitchAnalyticalRecords(studyExtract, studyAnalytical);

                List<String[]> merged4assaysAll = stitchAssayRecords(merged4assays, studyRun);

                merged4assaysAll = removeColumns.processSpreadsheet(merged4assaysAll);


                int counter = 0;
                String protocolnames = "";
                String protocoltypes = "";
                String protocoldescriptions = "";
                String protocolparameters = "";

                for (Map<String, List<String>> method : protocolsByTypes) {

                    for (String key : method.keySet()) {

                        if ((!method.values().iterator().next().get(1).contains("Description")) && (!method.values().iterator().next().get(1).replaceAll("\\s", "").equals(""))) {
                            protocolnames = protocolnames + method.values().iterator().next().get(0) + "\t";
                            protocoltypes = protocoltypes + key + "\t";
                            protocoldescriptions = protocoldescriptions + method.values().iterator().next().get(1) + "\t";
                            protocolparameters = protocolparameters + "" + "\t";
                            //  System.out.println(counter+ " :" + key + ": name: "+method.values().iterator().next().get(0)+"description" + method.values().iterator().next().get(1));
                            counter++;
                        }
                    }
                }


                protocolnames = protocolnames + "mass spectrometry" + "\t" + "metabolite identification" + "\n";

                protocoltypes = protocoltypes + "mass spectrometry" + "\t" + "data transformation" + "\n";

                protocoldescriptions = protocoldescriptions + gcmsProtocolDesc + "\t" + metaboliteIdProtocolDesc + "\n";

                protocolparameters = protocolparameters + "auto-injector;chromatography instrument;chromatography instrument control software;chromatography column;chromatography parameters;mass spectrometry instrument;sample injection;ion source;detector" + "";

                String outputdir = DownloadUtils.CONVERTED_DIRECTORY;

                printFiles(metabolomexpressFilename, studyAttributes, merged4Study, merged4assaysAll, protocolnames, protocoltypes, protocoldescriptions, personFirstNames, personLastNames, PersonMails);
            } else {
                System.out.println("ERROR: file not found!");
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    private static List<String[]> stitchStudyRecords(List<String[]> studyBiosource, List<String[]> studyExplants) {

        List<String[]> records = new ArrayList<String[]>();
        if   (studyExplants.get(0).length == studyBiosource.get(0).length) {
        for (String[] explantrecord : studyExplants) {
            for (String[] biosourcerecord : studyBiosource) {

                String[] newrecord = new String[studyExplants.get(0).length + studyBiosource.get(0).length];
              //  "| " + biosourcerecord[1] +   // + explantrecord[0]
                System.out.println(studyBiosource.get(0).length +  "|" +studyExplants.get(0).length  + "|"  + "|" +   explantrecord[0] );   //+ biosourcerecord[1]

                               if (biosourcerecord[0].equals(explantrecord[0])) {         //&& !biosourcerecord[1].equals("")
                   // System.out.println(studyBiosource.get(0).length + "| " + biosourcerecord[0] + "| " + biosourcerecord[1] + "|" +studyExplants.get(0).length  + "|" + explantrecord[0]);

                    int recordlength = 0;

                    for (int j = 0; j < explantrecord.length; j++) {
                        newrecord[j] = explantrecord[j];
                        recordlength++;
                    }

                    for (int k = 0; k < biosourcerecord.length; k++) {
                        if (k == 0) {
                            System.out.println("err:" + biosourcerecord[k]);
                        } else if (k == 1) {
                            System.out.println("sheeesh: " + biosourcerecord[k]);
                        } else {
                            newrecord[recordlength] = biosourcerecord[k];
                            recordlength++;
                        }
                        //System.out.println("current state: " + Arrays.toString(newrecord));
                    }
                   // System.out.println("position: " + recordlength + ": " + biosourcerecord[0]);
                    newrecord[newrecord.length - 1] = biosourcerecord[0];
                    //System.out.println("current state: " + Arrays.toString(newrecord));

                    records.add(newrecord);
                }

            }
        }
        }
        return records;
    }


    private static List<String[]> stitchAssayRecords(List<String[]> studyAnalytical, List<String[]> studyRun) {

        List<String[]> newrecords = new ArrayList<String[]>();

        int recordcounter=0;
        for (String[] runrecord : studyRun) {
            for (String[] analyticalrecord : studyAnalytical) {
                String[] newrecord = new String[studyAnalytical.get(0).length + studyRun.get(0).length+2];
                if (runrecord[0].equals(analyticalrecord[0]) && !runrecord[0].equals("")) {

                    int recordlength = 0;

                    for (int j = 1; j < analyticalrecord.length; j++) {
                        newrecord[j] = analyticalrecord[j];
                        recordlength++;
                    }

                    for (int k = 0; k < runrecord.length; k++) {
                        if (k == 0) {
                            System.out.println("err:" + runrecord[k]);
                        }
                        else if (k == 1) {
                            if (recordcounter>0){
                            newrecord[recordlength] ="data acquisition"+"\t" +runrecord[k];
                            recordlength++;
                            }
                            else {
                                newrecord[recordlength] ="Protocol REF"+"\t"+ runrecord[k];
                                recordlength++;
                            }

                        }
                        else {
                            newrecord[recordlength] = runrecord[k];
                            recordlength++;
                        }

                    }
                    newrecord[newrecord.length - 2] = runrecord[0];
                    if (recordcounter>0){
                    newrecord[newrecord.length-1] = analyticalrecord[0]+".CDF";
                    newrecords.add(newrecord);
                    }
                    else {
                        newrecord[newrecord.length-1] = "Raw Spectral Data File";
                        newrecords.add(newrecord);
                    }

                    //  System.out.println("BBB:L " + runrecord.toString());
//                    System.arraycopy(analyticalrecord, 1, newrecord, 0, analyticalrecord.length - 1);
//                    System.arraycopy(analyticalrecord, 0, newrecord, analyticalrecord.length - 1, 1);  //  analyticalrecord.length
//                    System.arraycopy(runrecord, 0, newrecord, analyticalrecord.length, 1);
//                    System.arraycopy(runrecord, 2, newrecord, analyticalrecord.length + 1, runrecord.length - 2);
                   // newrecords.add(newrecord);
                }

            }
            recordcounter++;
        }

        return newrecords;

    }


    private static List<String[]> stitchAnalyticalRecords(List<String[]> studyExtract , List<String[]> studyAnalytical) {

        List<String[]> newrecords = new ArrayList<String[]>();

        for (String[] extractrecord : studyExtract) {
            for (String[] analyticalrecord : studyAnalytical) {
                String[] newrecord = new String[studyAnalytical.get(0).length + studyExtract.get(0).length];
                if (extractrecord[0].equals(analyticalrecord[0]) && !extractrecord[0].equals("")) {


                    int recordlength = 0;

                    for (int j = 0; j < extractrecord.length; j++) {
                        newrecord[j] = extractrecord[j];
                        recordlength++;
                    }

                    for (int k = 0; k < analyticalrecord.length; k++) {
                        if (k == 0) {
                            System.out.println("err:" + analyticalrecord[k]);
                        } else if (k == 1) {
                            System.out.println("sheeesh: " + analyticalrecord[k]);
                        } else {
                            newrecord[recordlength] = analyticalrecord[k];
                            recordlength++;
                        }

                    }
                    newrecord[newrecord.length - 1] = analyticalrecord[0];

                    newrecords.add(newrecord);

                }

            }
        }
        return newrecords;
    }


    private static List<String> getPublicationTitleAndAuthors(String element) {
        List<String> publicationDetails = new ArrayList<String>();

        Pattern regex = Pattern.compile("(.*?)(\\([0-9]{4,}\\))(.*.^?)");
        Matcher matcher = regex.matcher(element);
        if (matcher.find()) {
            System.out.println(matcher.group(1) + "|" + matcher.group(2) + "|" + matcher.group(3));
            publicationDetails.add(matcher.group(1).toString());
            publicationDetails.add(matcher.group(3).toString());
        }

        return publicationDetails;
    }

    private static Pair<List<String>, List<String>> checkFirstLastNames(String element) {

        List<String> FirstNames = new ArrayList<String>();
        List<String> LastNames = new ArrayList<String>();
        Pair<List<String>, List<String>> Names = new Pair(FirstNames, LastNames);

        Pattern regex = Pattern.compile("^([A-z\\'\\.-]+)*\\s([A-Z]{1,}.)\\s([A-z\\'\\.-]*)$");
        Matcher matcher = regex.matcher(element);
        if (matcher.find()) {
            System.out.println(matcher.group(1) + "-" + matcher.group(2) + "-" + matcher.group(3));
            FirstNames.add(matcher.group(1).toString());
            LastNames.add(matcher.group(3).toString());
        }
        return Names;
    }


    //a method to convert metabolome express headers into ISA compatible headers
    private void cleanHeader(String[] row) {
        for (int i = 0; i < row.length; i++) {
            if (row[i].contains("BioSample ID")) {
                row[i] = row[i].replace("BioSample", "Sample").replace("ID", "Name");
            }
            if (row[i].contains("Explant ID")) {
                row[i] = row[i].replace("Explant", "Source").replace("ID", "Name").replace("(s)", "");
            }
            if (row[i].contains("Genotype ID")) {
                row[i] = row[i].replace("Genotype ID", "Characteristics[genotype]").replace("(s)", "");
            }
            if (row[i].contains("Environment ID")) {
                row[i] = row[i].replace("Environment ID", "Characteristics[environment]").replace("(s)", "");
            }
            if (row[i].contains("Treatment ID")) {
                row[i] = row[i].replace("Treatment ID", "Characteristics[treatment]").replace("(s)", "");
            }
            if (row[i].contains("Collection ID")) {
                row[i] = row[i].replace("Collection ID", "Characteristics[collection]").replace("(s)", "");
            }
            if (row[i].contains("Individual ID")) {
                row[i] = row[i].replace("Individual ID", "Characteristics[replicate]").replace("(s)", "");
            }
            if (row[i].contains("Extract ID")) {
                row[i] = row[i].replace("Extract ID", "Sample Name").replace("(s)", "");
            }
            if (row[i].contains("Analytical Sample ID")) {
                row[i] = row[i].replace("Analytical Sample ID", "Sample Name").replace("(s)", "");
            }
            if (row[i].contains("Run ID (Filename without Extension)")) {
                row[i] = row[i].replace("Run ID (Filename without Extension)", "MS Assay Name");
            }
            if (row[i].contains("Batch Sequence ID")) {
                row[i] = row[i].replace("Batch Sequence ID", "Parameter Value[batch sequence identifier]");
            }

            if (row[i].contains("Protocol ID")) {
                row[i] = "Protocol REF";
            }
            if (row[i].endsWith("Mass/Volume")) {
                row[i] = "Parameter Value[Mass or Volume]";
            }
            if (row[i].equalsIgnoreCase("extraction volume (ml)")) {
                row[i] = "Parameter Value[extraction volume in ml]";
            }
            if (row[i].contains("Units")) {
                row[i] = "Unit";
            }
            if (row[i].contains("Preparation Date")) {
                row[i] = row[i].replace("Preparation ", "").replace(" (YYYY-MM-DD)", "");
            }
            if (row[i].contains("Run Date Time")) {
                row[i] = row[i].replace("Run Date Time", "Parameter Value[acquisition timestamp]").replace(" (YYYY-MM-DD HH:MM)", "");
            }
        }
    }


    //public void printFiles(String metabolomeXpressFilename, List<String[]> MaterialSheetDataSubset, List<String[]> NewArcAssaySheetDataSubset) {
    public void printFiles(String metabolomeXpressFilename,
                           HashMap<String, String> studyAttributes,
                           List<String[]> studydata,
                           List<String[]> assaydata,
                           String protocolnames,
                           String protocoltypes,
                           String protocoldescriptions,
                           String personFirstNames,
                           String personLastNames,
                           String PersonMails
                           //ArrayList<HashMap<String, ArrayList<String>>> ProtocolsbyTypes  ,
                           //HashMap<String,String> protocolParamValues

    ) {
        try {

            //this is to get the basename of the mirada file
            String[] tokens = metabolomeXpressFilename.split("\\.(?=[^\\.]+$)");

            //we print the investigation files
            PrintStream investigationPs = new PrintStream(new File(DownloadUtils.CONVERTED_DIRECTORY + File.separator + tokens[0] + "/i_" + tokens[0] + "_investigation.txt"));

            investigationPs.print("ONTOLOGY SOURCE REFERENCE\t\t\n" +
                    "Term Source Name\tENVO\tNCBI\n" +
                    "Term Source File\t\t\n" +
                    "Term Source Version\tv 1.26\tv 1.26\n" +
                    "Term Source Description\tEnvironmental Ontology\tNCBI Taxonomy\n" +
                    "INVESTIGATION\t\t\n" +
                    "Investigation Identifier\t\t\n" +
                    "Investigation Title\t\t\n" +
                    "Investigation Description\t\t\n" +
                    "Investigation Submission Date\t\t\n" +
                    "Investigation Public Release Date\t\t\n" +
                    "INVESTIGATION PUBLICATIONS\t\t\n" +
                    "Investigation PubMed ID\t\t\n" +
                    "Investigation Publication DOI\t\t\n" +
                    "Investigation Publication Author list\t\t\n" +
                    "Investigation Publication Title\t\t\n" +
                    "Investigation Publication Status\t\t\n" +
                    "Investigation Publication Status Term Accession Number\t\t\n" +
                    "Investigation Publication Status Term Source REF\t\t\n" +
                    "INVESTIGATION CONTACTS\t\t\n" +
                    "Investigation Person Last Name\t\t\n" +
                    "Investigation Person First Name\t\t\n" +
                    "Investigation Person Mid Initials\t\t\n" +
                    "Investigation Person Email\t\t\n" +
                    "Investigation Person Phone\t\t\n" +
                    "Investigation Person Fax\t\t\n" +
                    "Investigation Person Address\t\t\n" +
                    "Investigation Person Affiliation\t\t\n" +
                    "Investigation Person Roles\t\t\n" +
                    "Investigation Person Roles Term Accession Number\t\t\n" +
                    "Investigation Person Roles Term Source REF\t\t\n" +
                    "\t\t\n" +
                    "STUDY\t\t\n" +
                    "Study Identifier\tTEST\t\n" +
                    "Study Title\t" + studyAttributes.get("studyName") + "\n" +
                    "Study Submission Date\t\n" +
                    "Study Public Release Date\t" + studyAttributes.get("studyDate") + "\n" +
                    "Study Description\t" + studyAttributes.get("studyDescription") + "\t\n" +
                    "Study File Name\t" + "s_" + tokens[0] + "_study_sample.txt\n" +
                    "Comment [MetabolomeXpress ID]\n" +
                    "STUDY DESIGN DESCRIPTORS\t\t\n" +
                    "Study Design Type\t\n" +
                    "Study Design Type Term Accession Number\t\n" +
                    "Study Design Type Term Source REF\t\n" +
                    "STUDY PUBLICATIONS\t\t\n" +
                    "Study PubMed ID\t" + studyAttributes.get("studyPMID") + "\n" +
                    "Study Publication DOI\t\t\n" +
                    "Study Publication Author list\t" + studyAttributes.get("studyPubAuthors") + "\n" +
                    "Study Publication Title\t" + studyAttributes.get("studyPubTitle") + "\n" +
                    "Study Publication Status\t\t\n" +
                    "Study Publication Status Term Accession Number\t\t\n" +
                    "Study Publication Status Term Source REF\t\t\n" +
                    "STUDY FACTORS\t\t\n" +
                    "Study Factor Name\t\n" +
                    "Study Factor Type\t\n" +
                    "Study Factor Type Term Accession Number\t\n" +
                    "Study Factor Type Term Source REF\t\n" +
                    "STUDY ASSAYS\t\t\n" +
                    "Study Assay Measurement Type\tmetabolite profiling\t\n" +
                    "Study Assay Measurement Type Term Accession Number\t\t\n" +
                    "Study Assay Measurement Type Term Source REF\t\t\n" +
                    "Study Assay Technology Type\tmass spectrometry\t\n" +
                    "Study Assay Technology Type Term Accession Number\t\t\n" +
                    "Study Assay Technology Type Term Source REF\t\t\n" +
                    "Study Assay Technology Platform\n" +
                    "Study Assay File Name\t" + "a_" + tokens[0] + "_assay.txt\n" +
                    "STUDY PROTOCOLS\n" +
                    "Study Protocol Name\t" + protocolnames + "\n" +
                    "Study Protocol Type\t" + protocoltypes + "\n" +
                    "Study Protocol Type Term Accession Number\t\t\n" +
                    "Study Protocol Type Term Source REF\t\t\n" +
                    "Study Protocol Description\t" + protocoldescriptions + "\n" +
                    "Study Protocol URI\t\t\n" +
                    "Study Protocol Version\t\t\n" +
                    "Study Protocol Parameters Name\t\n" +
                    "Study Protocol Parameters Name Term Accession Number\t\n" +
                    "Study Protocol Parameters Name Term Source REF\t\n" +
                    "Study Protocol Components Name\n" +
                    "Study Protocol Components Type\n" +
                    "Study Protocol Components Type Term Accession Number\n" +
                    "Study Protocol Components Type Term Source REF\n" +
                    "STUDY CONTACTS\t\t\n" +
                    "Study Person Last Name\t" + personLastNames.toString() + "\n" +
                    "Study Person First Name\t" + personFirstNames + "\n" +
                    "Study Person Mid Initials\t\t\n" +
                    "Study Person Email\t" + PersonMails + "\n" +
                    "Study Person Phone\t\t\n" +
                    "Study Person Fax\t\t\n" +
                    "Study Person Address\t\t\n" +
                    "Study Person Affiliation\t\t\n" +
                    "Study Person Roles\texperimentalist\tdata analyst\n" +
                    "Study Person Roles Term Accession Number\t\n" +
                    "Study Person Roles Term Source REF\t\n");


            //We print the ISA-study-sample file as deduced from the VAMPS file used as input
            PrintStream studyPs = new PrintStream(new File(DownloadUtils.CONVERTED_DIRECTORY + File.separator + tokens[0] + "/s_" + tokens[0] + "_study_sample.txt"));

            for (String[] record : studydata) {
                String newline = "";
                for (int i = 0; i < record.length; i++) {
                    if (i < record.length - 1) {
                        if (record[i] != null) {
                            newline = newline + record[i] + TAB_DELIM;
                        }
                    } else {
                        newline = newline + record[i];
                    }
                }

                System.out.println("record: " + newline);
                studyPs.println(newline);

            }

            PrintStream assayPs = new PrintStream(new File(DownloadUtils.CONVERTED_DIRECTORY + File.separator + tokens[0] + "/a_" + tokens[0] + "_assay.txt"));
            for (String[] record : assaydata) {
                String newline = "";
                for (int i = 0; i < record.length; i++) {
                    if (i < record.length - 1) {
                        newline = newline + record[i] + TAB_DELIM;
                    } else {
                        newline = newline + record[i];
                    }
                }
                System.out.println("record: " + newline);
                assayPs.println(newline);
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


}
