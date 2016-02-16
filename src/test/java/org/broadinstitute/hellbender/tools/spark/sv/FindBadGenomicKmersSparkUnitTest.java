package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.lang3.StringUtils;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.*;
import java.util.*;

/**
 * Unit tests for FindBadGenomicKmersSpark
 */
public class FindBadGenomicKmersSparkUnitTest extends BaseTest {

    private static final int KMER_SIZE = SVConstants.KMER_SIZE;
    private static final String REFERENCE_FILE_NAME = hg19MiniReference;

    @Test(groups = "spark")
    public void badKmersTest() throws IOException {
        final File refSequenceFile = File.createTempFile("genomic_sequence", "txt");
        refSequenceFile.deleteOnExit();

        final String polyA = StringUtils.repeat('A', KMER_SIZE);
        final String polyC = StringUtils.repeat('C', KMER_SIZE);
        final String polyT = StringUtils.repeat('T', KMER_SIZE);

        try ( PrintWriter writer = new PrintWriter(refSequenceFile) ) {
            long nTimes = FindBadGenomicKmersSpark.MAX_KMER_FREQ;
            // write polyA and polyC the max number of times possible while evading the trigger
            while ( nTimes-- > 0L ) {
                writer.println(polyA);
                writer.println(polyC);
            }
            // tip polyA over the edge by writing its reverse complement
            writer.println(polyT);
        }

        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final List<SVKmer> badKmers = FindBadGenomicKmersSpark.processRefFile(ctx, refSequenceFile);

        // should have just one bad kmer:  polyA
        Assert.assertEquals(badKmers.size(), 1);
        Assert.assertEquals(badKmers.get(0).toString(KMER_SIZE), polyA);

        if ( !refSequenceFile.delete() )
            throw new GATKException("Unable to delete file "+refSequenceFile);
    }

    @Test(groups = "spark")
    public void miniRefTest() throws IOException {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final ReferenceMultiSource ref = new ReferenceMultiSource((PipelineOptions)null,
                REFERENCE_FILE_NAME,
                ReferenceWindowFunctions.IDENTITY_FUNCTION);
        final SAMSequenceDictionary dict = ref.getReferenceSequenceDictionary(null);
        if ( dict == null ) throw new GATKException("No reference dictionary available.");

        final HashMap<SVKmer,Long> kmerMap = new HashMap<>();
        for ( final SAMSequenceRecord rec : dict.getSequences() ) {
            final SimpleInterval interval = new SimpleInterval(rec.getSequenceName(), 1, rec.getSequenceLength());
            final byte[] bases = ref.getReferenceBases(null, interval).getBases();
            final SVKmerizer kmerizer = new SVKmerizer(bases, KMER_SIZE);
            while ( kmerizer.hasNext() ) {
                final SVKmer kmer = kmerizer.next().canonical(KMER_SIZE);
                final Long currentCount = kmerMap.getOrDefault(kmer,0L);
                kmerMap.put(kmer,currentCount+1);
            }
        }
        final Iterator<Map.Entry<SVKmer,Long>> kmerIterator = kmerMap.entrySet().iterator();
        while ( kmerIterator.hasNext() ) {
            if ( kmerIterator.next().getValue() <= FindBadGenomicKmersSpark.MAX_KMER_FREQ )
                kmerIterator.remove();
        }

        final List<SVKmer> badKmers = FindBadGenomicKmersSpark.findBadGenomicKmers(ctx,ref,null,null);
        final HashSet<SVKmer> badKmerSet = new HashSet<>(badKmers.size());
        for ( final SVKmer kmer : badKmers ) {
            badKmerSet.add(kmer);
        }
        Assert.assertEquals(badKmers.size(), badKmerSet.size());
        Assert.assertEquals(badKmerSet, kmerMap.keySet());
    }
}