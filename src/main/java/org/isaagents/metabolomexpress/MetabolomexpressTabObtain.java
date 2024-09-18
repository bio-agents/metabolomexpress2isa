package org.isaagents.metabolomexpress;

import java.io.File;
import java.io.IOException;
import java.lang.Exception;import java.lang.String;import java.lang.System;import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MetabolomexpressTabObtain {

    public String metabolomexpressFileUrl;

    public MetabolomexpressTabObtain() {
    }

    public void initialise() {
        DownloadUtils.createTmpDirectory();
        System.out.println("Enter The metabolomexpress file name: ");
        readInput();
    }


    /**
     * A Method that asks for an accession number and checks it is well formed.
     */
    public void readInput() {

        Scanner inputstr = new Scanner(System.in);

        try {
            doConversion(inputstr.nextLine(), "Data");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }



    public File doConversion(String metabolomexpressFilename, String saveDirectory) throws Exception {

        DownloadUtils.CONVERTED_DIRECTORY = saveDirectory;

        String  metabolomexpressFileUrl;

        metabolomexpressFileUrl = ("file:///Users/prs/Documents/ISA-Metabolights/" + metabolomexpressFilename);

        Pattern  metabolomexpressregex = Pattern.compile("_METADATA");
        Matcher  metabolomexpressmatcher =  metabolomexpressregex.matcher(metabolomexpressFilename);

        try {

            if ( metabolomexpressmatcher.find()) {

                System.out.println("MetabolomeXpress input at: " +  metabolomexpressFileUrl);

                DownloadUtils.createDirectory(DownloadUtils.TMP_DIRECTORY);

                String  metabolomexpressTabDownloadLocation = DownloadUtils.TMP_DIRECTORY  + metabolomexpressFilename  ;

                DownloadUtils.downloadFile( metabolomexpressFileUrl,  metabolomexpressTabDownloadLocation);

                System.out.println("where? " +  metabolomexpressTabDownloadLocation);

                MetabolomexpressTabLoader metabolomexpressloader = new MetabolomexpressTabLoader();

                metabolomexpressloader.loadMetabolomexpressTab( metabolomexpressTabDownloadLocation,  metabolomexpressFilename);


                //this is to get the basename of the metabolomeXpress file
                 String[] tokens =  metabolomexpressFilename.split("\\.(?=[^\\.]+$)");

                System.out.println("File is: "+ DownloadUtils.CONVERTED_DIRECTORY + File.separator + tokens[0]+ "." +tokens[1]);

                return new File(DownloadUtils.CONVERTED_DIRECTORY + File.separator + tokens[0]);

            } else {

                throw new Exception("Sorry, file " + metabolomexpressFilename + " not found !");

            }

        } catch (IOException ioe) {
            System.out.println("Caught an IO exception :-o");
            ioe.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }

        return null;
    }


    public static void main(String[] argv) {
        MetabolomexpressTabObtain metabolomexpressReadFunction = new MetabolomexpressTabObtain();
        metabolomexpressReadFunction.initialise();
    }

}