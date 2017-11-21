package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.python.StreamingPythonScriptExecutor;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.io.*;
import java.util.*;

/**
 * Created by sam on 11/17/17.
 */
@CommandLineProgramProperties(
        summary = "Example/toy program that uses a Python script.",
        oneLineSummary = "Example/toy program that uses a Python script.",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)

public class NeuralNetStreamingExecutor extends VariantWalker {
    private final static String NL = String.format("%n");

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file")
    private File outputFile; // output file produced by Python code

    @Argument(fullName = "architecture", shortName = "a", doc = "Neural Net architecture and weights hd5 file", optional = false)
    private String architecture = null;


    // Create the Python executor. This doesn't actually start the Python process, but verifies that
    // the requestedPython executable exists and can be located.
    final StreamingPythonScriptExecutor pythonExecutor = new StreamingPythonScriptExecutor(true);

    private FileWriter fifoWriter;
    private VariantContextWriter vcfWriter;

    private PrintStream outputStream = null;
    private boolean noSamples = true;

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
            fifoWriter = new FileWriter(fifoFile);
        } catch ( IOException e ) {
            throw new GATKException("Failure opening FIFO for writing", e);
        }

        // Also, ask Python to open our output file, where it will write the contents of everything it reads
        // from the FIFO. <code sendSynchronousCommand/>
        pythonExecutor.sendSynchronousCommand(String.format("tempFile = open('%s', 'a')" + NL, outputFile.getAbsolutePath()));
        pythonExecutor.sendSynchronousCommand("from keras.models import load_model" + NL);
        pythonExecutor.sendSynchronousCommand("import vqsr_cnn" + NL);
        logger.info("Imported from vqsr cnn.");
        pythonExecutor.sendSynchronousCommand(String.format("model = load_model('%s')", architecture) + NL);
        logger.info("Loaded architecture:"+architecture);
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        referenceContext.setWindow(63, 64);
        transferToPythonViaFifo(variant, referenceContext);

        //transferVariantSummaryToPython(variant, referenceContext);
    }

    private void transferToPythonViaFifo(final VariantContext variant, final ReferenceContext referenceContext) {
        try {
            String varData = getVariantDataString(variant, referenceContext);
            String ref = new String(Arrays.copyOfRange(referenceContext.getBases(),0,128), "UTF-8");
            String genos = "\t.";
            if (noSamples) genos = "";
            // write summary data to the FIFO
            fifoWriter.write(String.format("%s|%s|%s|%s|%s\n",
                    ref,
                    variant.getAttributes().toString(),
                    varData,
                    GATKVCFConstants.CNN_1D_KEY,
                    genos));

            // NOTE: flush the stream since we're going to immediately issue a Python command to read it.
            // Failure to flush the stream could result in Python blocking on the read.
            fifoWriter.flush();
        }
        catch ( IOException e ) {
            throw new GATKException("Failure writing to FIFO", e);
        }

        // Send a synchronous command (wait for the response prompt) to Python to consume the line written to the FIFO
        if(variant.isSNP()) {
            pythonExecutor.sendSynchronousCommand("vqsr_cnn.write_snp_score(model, tempFile, fifoFile.readline())" + NL);
        } else if(variant.isIndel()){
            pythonExecutor.sendSynchronousCommand("vqsr_cnn.write_indel_score(model, tempFile, fifoFile.readline())" + NL);
        }
    }


    private void transferVariantSummaryToPython(final VariantContext variant, final ReferenceContext referenceContext) {
        try {
            String varData = getVariantDataString(variant, referenceContext);
            String ref = new String(Arrays.copyOfRange(referenceContext.getBases(),0,128), "UTF-8");
            String genos = "\t.";

            if(variant.isSNP()){
                pythonExecutor.sendSynchronousCommand(String.format("score = vqsr_cnn.snp_score_from_reference_annotations(model, '%s', '%s')\n", ref, variant.getAttributes().toString()));
            } else if(variant.isIndel()){
                pythonExecutor.sendSynchronousCommand(String.format("score = vqsr_cnn.indel_score_from_reference_annotations(model, '%s', '%s')\n", ref, variant.getAttributes().toString()));
            }
            pythonExecutor.sendSynchronousCommand(String.format("tempFile.write('%s'+'%s='+str(score)+'%s'+str('\\n'))\n", varData, GATKVCFConstants.CNN_1D_KEY, genos));

        } catch (UnsupportedEncodingException e) {
            e.printStackTrace();
        }


    }

    private String getVariantDataString(final VariantContext variant, final ReferenceContext referenceContext){
        String varInfo = "";
        for (final String attributeKey : variant.getAttributes().keySet()) {
            varInfo += attributeKey + "=" + variant.getAttribute(attributeKey).toString().replace(" ", "").replace("[", "").replace("]", "") + ";";
        }

        String alts = variant.getAlternateAlleles().toString().replace(" ", "");
        alts = alts.substring(1, alts.length() - 1);

        String varData = String.format("%s\t%d\t.\t%s\t%s\t%.3f\t.\t%s",
                variant.getContig(),
                variant.getStart(),
                variant.getReference().getBaseString(),
                alts,
                variant.getPhredScaledQual(),
                varInfo
        );

        return varData;
    }

    @Override
    public Object onTraversalSuccess() {
        pythonExecutor.sendSynchronousCommand("tempFile.close()\n");

        return true;
    }


}
