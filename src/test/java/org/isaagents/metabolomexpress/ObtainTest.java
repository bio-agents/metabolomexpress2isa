package org.isaagents.metabolomexpress;

import org.junit.Test;

import java.io.File;

import static org.junit.Assert.assertTrue;

/**
 * Created by IntelliJ IDEA.
 * User: prs
 * Date: 02/10/2012
 * Time: 13:59
 * To change this template use File | Settings | File Templates.
 */
public class ObtainTest {


    @Test
    public void testObtain() {
        MetabolomexpressTabObtain obtain  = new MetabolomexpressTabObtain();
        try {
            obtain.doConversion("metabolomeXpress-test.txt", "Data");
            assertTrue("Directory should not be empty.", new File("Data").listFiles().length > 0);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
