package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.AnnotationUtils;
import org.broadinstitute.hellbender.tools.walkers.contamination.PileupSummary;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by sam on 11/18/17.
 */
@CommandLineProgramProperties(
        summary = "Apply a score cutoff to filter variants based on a recalibration table",
        oneLineSummary = " Apply a score cutoff to filter variants based on a recalibration table",
        programGroup = VariantProgramGroup.class
)
@DocumentedFeature
public class AddScores extends MultiVariantWalker {


    /////////////////////////////
    // Inputs
    /////////////////////////////
    @Argument(fullName="score_file", shortName="scoreFile", doc="The input score file used to annotate the VCF", optional=false)
    private FeatureInput<VariantContext> scores;


    /////////////////////////////
    // Outputs
    /////////////////////////////

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output annotated VCF file in which each variant has score", optional=false)
    private String output;



    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private VariantContextWriter vcfWriter;
    final static private String scoreKey = "CNN_1D";


    //---------------------------------------------------------------------------------------------------------------
    //
    // onTraversalStart
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public void onTraversalStart() {

        // setup the header fields
        VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> inputHeaders = inputHeader.getMetaDataInSortedOrder();
        final Set<VCFHeaderLine> hInfo = new HashSet<>(inputHeaders);
        VariantRecalibrationUtils.addVQSRStandardHeaderLines(hInfo);
        final TreeSet<String> samples = new TreeSet<>();
        samples.addAll(inputHeader.getGenotypeSamples());
        hInfo.addAll(getDefaultToolVCFHeaderLines());
        final VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter = createVCFWriter(new File(output));
        vcfWriter.writeHeader(vcfHeader);
    }



    //---------------------------------------------------------------------------------------------------------------
    //
    // apply
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext ref, final FeatureContext featureContext) {
        final List<VariantContext> scoreList =  featureContext.getValues(scores, vc.getStart());
        VariantContext scoreDatum = getMatchingVC(vc, scoreList, null);
        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        final VariantContext outputVC = builder.make();

        if( scoreDatum == null ) {
            throw new UserException("Encountered input variant which isn't found in the input score file. First seen at: " + vc );
        }

        final String lodString = scoreDatum.getAttributeAsString(GATKVCFConstants.VQS_LOD_KEY, null);
        if( lodString == null ) {
            throw new UserException("Encountered a malformed record in the input score file. There is no score for the record at: " + vc );
        }
        final double score;
        try {
            score = Double.valueOf(lodString);
        } catch (NumberFormatException e) {
            throw new UserException("Encountered a malformed record in the input score file. The lod is unreadable for the record at: " + vc );
        }

        builder.attribute(scoreKey, score);


        vcfWriter.add( outputVC );

    }


    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }



    private VariantContext getMatchingVC(final VariantContext target, final List<VariantContext> recalVCs, final Allele allele) {
        for( final VariantContext recalVC : recalVCs ) {
            if ( target.getEnd() == recalVC.getEnd() ) {
                return recalVC;
            }
        }
        return null;
    }


}

