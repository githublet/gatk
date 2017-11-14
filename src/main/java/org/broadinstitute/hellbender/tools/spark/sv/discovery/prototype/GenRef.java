package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.IOException;

public class GenRef extends GATKSparkTool {
    @Override public boolean requiresReference() { return true; }
    @Override public void runTool( final JavaSparkContext ctx ) {
        final String insert1a_5817182 = "Ctggtgtgtggatacgggggattactggtgtgtggatacgggagattcctggtgtgtggatacaggggattactggtgtgtggatacaggggattactggtgtgtggacagaggggattactggtatgtggacagaggggattactggtgtgtggatacaggggattactggtgtgtggacagagggattactggtgtgtggatatgggggattactggtgtgtggatacgagggatta";
        final String insert1b_5817272 = "GgattactggtgtgtggacagaggggattactggtatgtggacagaggggattactggtgtgtggatacaggggattactggtgtgtggacagagggattactggtgtgtggatatgggggattactggtgtgtggatacgagggattactggtgtgtggatacgggggattaccggtgtgtggacacgggggattacgggtgcatggacagaggggattaccagtgtgtggatacaggggattgctggtgtgtggacagaggggattactggtgtgtggatacgggggattactggtgtgtggaTATGGA";
        final String insert13_5817257 = "Gtgtgtggatccaggggattactggtgtgtggacagaggggattactggtatgtggacagaggggattactggtgtgtggatacaggggattactggtgtgtggacagagggattactggtgtgtggatatgggggattactggtgtgtggatacgagggattactggtgtgtggatacgggggattaccggtgtgtggacacgggggattacgggtgcatggacagaggggattaccagtgtgtggatacaggggattgctggtgtgtggacagaggggattactggtgtgtggatacgggggattactg";
        final int start = 5816154;
        final SimpleInterval refInterval = new SimpleInterval("chr1",start,5818485);
        final String ref;
        try {
            ref = new String(getReference().getReferenceBases(refInterval).getBases());
        } catch ( final IOException ioe ) { throw new GATKException("Couldn't get ref", ioe); }
        final String chm1 =
                ref.substring(0, 5817183-start) + insert1a_5817182.substring(1) +
                ref.substring(5817183-start, 5817273-start) + insert1b_5817272.substring(1) +
                ref.substring(5817273-start);
        final String chm13 =
                ref.substring(0, 5817258-start) + insert13_5817257.substring(1) + ref.substring(5817258);
        dumpFASTA("chm1", chm1);
        dumpFASTA("chm13", chm13);
    }

    private void dumpFASTA( final String name, final String seq ) {
        System.out.println(">"+name);
        final int len = seq.length();
        final int lineLen = 80;
        for ( int idx = 0; idx < len; idx += lineLen )
            System.out.println(seq.substring(idx, Math.min(len,idx+lineLen)));
    }
}
