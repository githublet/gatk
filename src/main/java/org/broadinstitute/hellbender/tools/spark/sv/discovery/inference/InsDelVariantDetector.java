package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AnnotatedVariantProducer;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoveryInputData;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoveryUtils;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Stream;

@Deprecated
public final class InsDelVariantDetector implements VariantDetectorFromLocalAssemblyContigAlignments {

    @Override
    public void inferSvAndWriteVCF(final JavaRDD<AssemblyContigWithFineTunedAlignments> assemblyContigs,
                                   final SvDiscoveryInputData svDiscoveryInputData) {

        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputData.referenceSequenceDictionaryBroadcast;
        final String outputPath = svDiscoveryInputData.outputPath;
        final Logger toolLogger = svDiscoveryInputData.toolLogger;

        // TODO: 11/23/17 take insertion mappings from the input and add them to VC
        final JavaPairRDD<byte[], List<ChimericAlignment>> contigSeqAndChimeras =
                assemblyContigs
                        .map( AssemblyContigWithFineTunedAlignments::getSourceContig )
                        .mapToPair(tig ->
                                org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.InsDelVariantDetector
                                        .convertAlignmentIntervalToChimericAlignment(tig,
                                StructuralVariationDiscoveryArgumentCollection.
                                        DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD,
                                StructuralVariationDiscoveryArgumentCollection.
                                        DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH,
                                referenceSequenceDictionaryBroadcast.getValue()));

        final List<VariantContext> annotatedInsDels =
                org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.InsDelVariantDetector
                        .produceVariantsFromSimpleChimeras(contigSeqAndChimeras, svDiscoveryInputData);

        SVVCFWriter.writeVCF(annotatedInsDels, outputPath,
                referenceSequenceDictionaryBroadcast.getValue(), toolLogger);
    }

    /**
     * Very similar to {@link ChimericAlignment#parseOneContig(AlignedContig, SAMSequenceDictionary, boolean, int, int, boolean)}, except that
     * badly mapped (MQ < 60) 1st alignment is no longer skipped.
     */
    @Deprecated
    private static Tuple2<byte[], List<ChimericAlignment>> convertAlignmentIntervalToChimericAlignment (final AlignedContig contig,
                                                                                                        final int mapQualThresholdInclusive,
                                                                                                        final int minAlignmentBlockSize,
                                                                                                        final SAMSequenceDictionary referenceDictionary) {

        final List<AlignmentInterval> alignmentIntervals = contig.alignmentIntervals;
        final Iterator<AlignmentInterval> iterator = alignmentIntervals.iterator();
        AlignmentInterval current = iterator.next();
        final List<ChimericAlignment> results = new ArrayList<>(alignmentIntervals.size() - 1);
        final List<String> insertionMappings = new ArrayList<>();
        while ( iterator.hasNext() ) {
            final AlignmentInterval next = iterator.next();
            if (ChimericAlignment.nextAlignmentMayBeInsertion(current, next, mapQualThresholdInclusive, minAlignmentBlockSize, true)) {
                if (iterator.hasNext()) {
                    insertionMappings.add(next.toPackedString());
                    continue;
                } else {
                    break;
                }
            }

            final ChimericAlignment ca = new ChimericAlignment(current, next, insertionMappings, contig.contigName, referenceDictionary);
            final boolean validStateForThisPath = ca.isNeitherSimpleTranslocationNorIncompletePicture();
            if ( ! validStateForThisPath )
                throw new GATKException.ShouldNeverReachHereException("Mapped assembled contigs are sent down the wrong path: " +
                        "contig suggesting \"translocation\" is sent down the insert/deletion path.\n" + contig.toString());
            results.add(ca);
        }

        if( results.size() > 1)
            throw new GATKException("Expecting only no more than one simple chimera but got more!\n" + contig.toString() + "\n" + results.toString());
        return new Tuple2<>(contig.contigSequence,results);
    }

    @Deprecated
    public static List<VariantContext> produceVariantsFromSimpleChimeras(final JavaPairRDD<byte[], List<ChimericAlignment>> contigSeqAndChimeras,
                                                                         final SvDiscoveryInputData svDiscoveryInputData) {

        final Broadcast<ReferenceMultiSource> referenceBroadcast = svDiscoveryInputData.referenceBroadcast;
        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputData.referenceSequenceDictionaryBroadcast;
        final List<SVInterval> assembledIntervals = svDiscoveryInputData.assembledIntervals;
        final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast = svDiscoveryInputData.cnvCallsBroadcast;
        final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs = svDiscoveryInputData.discoverStageArgs;
        final String sampleId = svDiscoveryInputData.sampleId;
        final Logger toolLogger = svDiscoveryInputData.toolLogger;

        final JavaPairRDD<NovelAdjacencyReferenceLocations, Iterable<ChimericAlignment>> narlsAndSources =
                contigSeqAndChimeras
                        .flatMapToPair(tigSeqAndChimeras -> {
                            final byte[] contigSeq = tigSeqAndChimeras._1;
                            final List<ChimericAlignment> chimericAlignments = tigSeqAndChimeras._2;
                            final Stream<Tuple2<NovelAdjacencyReferenceLocations, ChimericAlignment>> novelAdjacencyAndSourceChimera =
                                    chimericAlignments.stream()
                                            .map(ca -> new Tuple2<>(
                                                    new NovelAdjacencyReferenceLocations(ca, contigSeq,
                                                            referenceSequenceDictionaryBroadcast.getValue()), ca));
                            return novelAdjacencyAndSourceChimera.iterator();
                        })
                        .groupByKey();                                                      // group the same novel adjacency produced by different contigs together

        narlsAndSources.cache();

        SvDiscoveryUtils.evaluateIntervalsAndNarls(assembledIntervals, narlsAndSources.map(Tuple2::_1).collect(),
                referenceSequenceDictionaryBroadcast.getValue(), discoverStageArgs, toolLogger);

        List<VariantContext> annotatedVariants =
                narlsAndSources
                        .mapToPair(noveltyAndEvidence -> new Tuple2<>(noveltyAndEvidence._1,
                                new Tuple2<>(SimpleNovelAdjacencyInterpreter.inferTypeFromNovelAdjacency(noveltyAndEvidence._1), noveltyAndEvidence._2)))       // type inference based on novel adjacency and evidence alignments
                        .map(noveltyTypeAndEvidence ->
                                {
                                    final NovelAdjacencyReferenceLocations novelAdjacency = noveltyTypeAndEvidence._1;
                                    final SimpleSVType inferredSimpleType = noveltyTypeAndEvidence._2._1;
                                    final Iterable<ChimericAlignment> evidence = noveltyTypeAndEvidence._2._2;
                                    return AnnotatedVariantProducer
                                            .produceAnnotatedVcFromInferredTypeAndRefLocations(
                                                    novelAdjacency.leftJustifiedLeftRefLoc,
                                                    novelAdjacency.leftJustifiedRightRefLoc.getStart(),
                                                    novelAdjacency.complication,
                                                    inferredSimpleType,
                                                    null, // TODO: 1/21/18 implement this for InsDel
                                                    evidence,
                                                    referenceBroadcast,
                                                    referenceSequenceDictionaryBroadcast,
                                                    cnvCallsBroadcast,
                                                    sampleId);
                                })
                        .collect();

        narlsAndSources.unpersist();

        return annotatedVariants;
    }
}
