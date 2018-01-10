package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.python.StreamingPythonScriptExecutor;
import org.broadinstitute.hellbender.utils.runtime.AsynchronousStreamWriterService;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.*;
import java.util.*;
import java.util.concurrent.TimeUnit;

/**
 * Created by sam on 11/17/17.
 */
@CommandLineProgramProperties(
        summary = NeuralNetStreamingExecutor.USAGE_SUMMARY,
        oneLineSummary = NeuralNetStreamingExecutor.USAGE_ONE_LINE_SUMMARY,
        programGroup = VariantEvaluationProgramGroup.class
)

public class NeuralNetStreamingExecutor extends VariantWalker {
    private final static String NL = String.format("%n");
    static final String USAGE_ONE_LINE_SUMMARY = "Apply 1d Convolutional Neural Net to filter annotated variants";
    static final String USAGE_SUMMARY = "Annotate a VCF with scores from 1d Convolutional Neural Network (CNN)." +
            "The CNN will look at the reference sequence and variant annotations to determine a Log Odds Score for each variant.";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file")
    private File outputFile = null; // output file produced by Python code

    @Argument(fullName = "architecture", shortName = "a", doc = "Neural Net architecture and weights hd5 file", optional = false)
    private String architecture = null;

    @Argument(fullName = "keep-info", shortName = "ki", doc = "Keep info fields in the vcf-like score file.", optional = true)
    private boolean keepInfo = true;

    @Argument(fullName = "python-batch-size", shortName = "pbs", doc = "Size of batches for python to do inference.", optional = true)
    private int pythonBatchSize = 256;

    @Argument(fullName = "python-sync-frequency", shortName = "psf", doc = "Size of data to queue for Python streaming.", optional = true)
    private int pythonSyncFrequency = 512;

    // Create the Python executor. This doesn't actually start the Python process, but verifies that
    // the requestedPython executable exists and can be located.
    final StreamingPythonScriptExecutor pythonExecutor = new StreamingPythonScriptExecutor(true);

    private FileOutputStream fifoWriter;
    private VariantContextWriter vcfWriter;
    private AsynchronousStreamWriterService<String> asyncWriter = null;
    private List<String> batchList = new ArrayList<>(pythonBatchSize);

    private boolean noSamples = true;
    private int curBatchSize = 0;


    @Override
    public boolean requiresReference(){
        return true;
    }

    @Override
    public void onTraversalStart() {
        // setup the header fields
        VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> inputHeaders = inputHeader.getMetaDataInSortedOrder();
        final Set<VCFHeaderLine> hInfo = new HashSet<>(inputHeaders);
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.CNN_1D_KEY));
        VariantRecalibrationUtils.addVQSRStandardHeaderLines(hInfo);
        final TreeSet<String> samples = new TreeSet<>();
        samples.addAll(inputHeader.getGenotypeSamples());
        if(samples.size() > 0) noSamples = false;
        hInfo.addAll(getDefaultToolVCFHeaderLines());
        final VCFHeader vcfHeader = new VCFHeader(hInfo, samples);

        vcfWriter = createVCFWriter(outputFile);
        vcfWriter.writeHeader(vcfHeader);
        vcfWriter.close();

        // Start the Python process, and get a FIFO from the executor to use to send data to Python. The lifetime
        // of the FIFO is managed by the executor; the FIFO will be destroyed when the executor is destroyed.
        pythonExecutor.start(Collections.emptyList(), true);
        final File fifoFile = pythonExecutor.getFIFOForWrite();

        // Open the FIFO for writing. Opening a FIFO for read or write will block until there is reader/writer
        // on the other end, so before we open it, send an ASYNCHRONOUS command (one that doesn't wait for a
        // response) to the Python process to open the FIFO for reading. The Python process will then block until
        // we open the FIFO. We can then call getAccumulatedOutput.
        pythonExecutor.sendAsynchronousCommand(String.format("fifoFile = open('%s', 'r')" + NL, fifoFile.getAbsolutePath()));
        try {
            fifoWriter = new FileOutputStream(fifoFile);
        } catch ( IOException e ) {
            throw new GATKException("Failure opening FIFO for writing", e);
        }

        asyncWriter = pythonExecutor.getAsynchronousStreamWriterService(fifoWriter, AsynchronousStreamWriterService.stringSerializer);
        batchList = new ArrayList<>(pythonSyncFrequency);

        // Also, ask Python to open our output file, where it will write the contents of everything it reads
        // from the FIFO. <code sendSynchronousCommand/>
        pythonExecutor.sendSynchronousCommand(String.format("tempFile = open('%s', 'a')" + NL, outputFile.getAbsolutePath()));
        pythonExecutor.sendSynchronousCommand("from keras.models import load_model" + NL);
        pythonExecutor.sendSynchronousCommand("import vqsr_cnn" + NL);
        pythonExecutor.sendSynchronousCommand(String.format("model = load_model('%s')", architecture) + NL);
        logger.info("Loaded CNN architecture:"+architecture);
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        referenceContext.setWindow(63, 64);
        transferToPythonViaFifo(variant, referenceContext);
    }

    private void transferToPythonViaFifo(final VariantContext variant, final ReferenceContext referenceContext) {
        if (variant.isSNP() || variant.isIndel()) {
            try {
                String ref = new String(Arrays.copyOfRange(referenceContext.getBases(), 0, 128), "UTF-8");
                String varData = getVariantDataString(variant, referenceContext);
                String isSnp = variant.isSNP() ? "1" : "0";
                String genos = "\t.";
                if (noSamples) genos = "";

                final String outDat = String.format("%s|%s|%s|%s|%s|%s|\n",
                        ref,
                        getVariantInfoString(variant),
                        varData,
                        GATKVCFConstants.CNN_1D_KEY,
                        genos,
                        isSnp);

                if (curBatchSize == pythonSyncFrequency) {
                    // wait for the last batch to complete before we start a new one
                    asyncWriter.waitForPreviousBatchCompletion(1, TimeUnit.MINUTES);
                    final String pythonCommand = String.format(
                            "vqsr_cnn.score_and_write_batch(model, tempFile, fifoFile, %d, %d)", curBatchSize, pythonBatchSize) + NL;
                    pythonExecutor.sendAsynchronousCommand(pythonCommand);
                    asyncWriter.startAsynchronousBatchWrite(batchList);
                    curBatchSize = 0;
                    batchList = new ArrayList<>(pythonSyncFrequency);
                }

                // write summary data to the FIFO
                batchList.add(outDat);
                curBatchSize++;
            } catch (UnsupportedEncodingException e) {
                throw new GATKException("Trying to make string from reference, but unsupported encoding UTF-8.", e);
            }
        }
    }

    private String getVariantDataString(final VariantContext variant, final ReferenceContext referenceContext){
        String varInfo = keepInfo ? getVariantInfoString(variant) : "" ;
        String alts = variant.getAlternateAlleles().toString().replace(" ", "");
        alts = alts.substring(1, alts.length() - 1);

        String varData = String.format("%s\t%d\t.\t%s\t%s\t%.2f\t.\t%s",
                variant.getContig(),
                variant.getStart(),
                variant.getReference().getBaseString(),
                alts,
                variant.getPhredScaledQual(),
                varInfo
        );

        return varData;
    }

    private String getVariantInfoString(final VariantContext variant){
        String varInfo = "";
        for (final String attributeKey : variant.getAttributes().keySet()) {
            varInfo += attributeKey + "=" + variant.getAttribute(attributeKey).toString().replace(" ", "").replace("[", "").replace("]", "") + ";";
        }
        return varInfo;
    }

    @Override
    public Object onTraversalSuccess() {
        asyncWriter.waitForPreviousBatchCompletion(1, TimeUnit.MINUTES);
        if (curBatchSize > 0){
            final String pythonCommand = String.format(
                    "vqsr_cnn.score_and_write_batch(model, tempFile, fifoFile, %d, %d)", curBatchSize, pythonBatchSize) + NL;
            pythonExecutor.sendAsynchronousCommand(pythonCommand);
            asyncWriter.startAsynchronousBatchWrite(batchList);
            asyncWriter.waitForPreviousBatchCompletion(1, TimeUnit.MINUTES);
        }

        pythonExecutor.sendSynchronousCommand("tempFile.close()" + NL);
        pythonExecutor.sendSynchronousCommand("fifoFile.close()" + NL);
        pythonExecutor.terminate();

        return true;
    }

}
