package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.*;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import scala.Tuple2;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * This deals with the special case where a contig has exactly two alignments
 * and seemingly has the complete alt haplotype assembled.
 * See criteria in {@link AssemblyContigWithFineTunedAlignments#hasIncompletePictureFromMultipleAlignments()}.
 *
 * TODO: 1/19/18
 *      Exactly how the returned type in {@link SimpleNovelAdjacency} is treated (trusted, or updated, or re-interpreted),
 *      is to be developed.
 */
public final class SimpleNovelAdjacencyInterpreter {

    static final int MORE_RELAXED_ALIGNMENT_MIN_LENGTH = 30;
    static final int MORE_RELAXED_ALIGNMENT_MIN_MQ = 20;

    public JavaPairRDD<SimpleNovelAdjacency, List<SvType>> inferTypeFromSingleContigSimpleChimera(final JavaRDD<AssemblyContigWithFineTunedAlignments> assemblyContigs,
                                                                                                  final SvDiscoveryInputData svDiscoveryInputData) {

        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputData.referenceSequenceDictionaryBroadcast;
        final Broadcast<ReferenceMultiSource> referenceBroadcast = svDiscoveryInputData.referenceBroadcast;

        final JavaRDD<SimpleNovelAdjacency> simpleNovelAdjacencies =
                getSimpleNovelAdjacencyJavaRDD(assemblyContigs, svDiscoveryInputData);

        return simpleNovelAdjacencies
                        .mapToPair(simpleNovelAdjacency ->
                                new Tuple2<>(simpleNovelAdjacency,
                                        inferTypeFromNovelAdjacency(simpleNovelAdjacency,
                                                referenceBroadcast.getValue(), referenceSequenceDictionaryBroadcast.getValue())));
    }

    private JavaRDD<SimpleNovelAdjacency> getSimpleNovelAdjacencyJavaRDD(final JavaRDD<AssemblyContigWithFineTunedAlignments> assemblyContigs,
                                                                         final SvDiscoveryInputData svDiscoveryInputData) {
        final String sampleId = svDiscoveryInputData.sampleId;
        final Logger toolLogger = svDiscoveryInputData.toolLogger;
        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputData.referenceSequenceDictionaryBroadcast;
        final Broadcast<ReferenceMultiSource> referenceBroadcast = svDiscoveryInputData.referenceBroadcast;
        final DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs = svDiscoveryInputData.discoverStageArgs;
        final List<SVInterval> assembledIntervals = svDiscoveryInputData.assembledIntervals;
        final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast = svDiscoveryInputData.cnvCallsBroadcast;

        final JavaRDD<SimpleNovelAdjacency> simpleNovelAdjacencies =
                assemblyContigs
                        .filter(tig ->
                                ChimericAlignment.splitPairStrongEnoughEvidenceForCA(tig.getSourceContig().alignmentIntervals.get(0),
                                        tig.getSourceContig().alignmentIntervals.get(1),
                                        MORE_RELAXED_ALIGNMENT_MIN_MQ, MORE_RELAXED_ALIGNMENT_MIN_LENGTH))
                        .mapToPair(tig -> {
                            final ChimericAlignment simpleChimera = ChimericAlignment.extractSimpleChimera(tig,
                                    referenceSequenceDictionaryBroadcast.getValue());
                            final byte[] contigSequence = tig.getSourceContig().contigSequence;
                            final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations =
                                    new NovelAdjacencyReferenceLocations(simpleChimera, contigSequence,
                                            referenceSequenceDictionaryBroadcast.getValue());
                            final byte[] altHaplotypeSeq = null; // TODO: 1/21/18 force all subtypes of NovelAdjacencyReferenceLocations to extract alt haplotype
                            final SimpleNovelAdjacency.NovelAdjacencyAndInferredAltHaptype novelAdjacencyAndInferredAltHaptype =
                                    new SimpleNovelAdjacency.NovelAdjacencyAndInferredAltHaptype(novelAdjacencyReferenceLocations, altHaplotypeSeq);
                            return new Tuple2<>(novelAdjacencyAndInferredAltHaptype, simpleChimera);
                        })
                        .groupByKey()       // group the same novel adjacency produced by different contigs together
                        .map(noveltyAndEvidence ->
                                new SimpleNovelAdjacency(noveltyAndEvidence._1, Lists.newArrayList(noveltyAndEvidence._2)));

        SvDiscoveryUtils.evaluateIntervalsAndNarls(assembledIntervals,
                simpleNovelAdjacencies.map(SimpleNovelAdjacency::getNovelAdjacencyReferenceLocations).collect(),
                referenceSequenceDictionaryBroadcast.getValue(), discoverStageArgs, toolLogger);
        return simpleNovelAdjacencies;
    }

    static List<SvType> inferTypeFromNovelAdjacency(final SimpleNovelAdjacency simpleNovelAdjacency,
                                                    final ReferenceMultiSource reference, final SAMSequenceDictionary referenceDictionary) {

        // based on characteristic of simple chimera, infer type
        final List<ChimericAlignment> alignmentEvidence = simpleNovelAdjacency.getAlignmentEvidence();
        final boolean allIndicateSimpleTransLoc = alignmentEvidence.stream().allMatch(ChimericAlignment::isLikelySimpleTranslocation);
        final boolean allIndicateInvDup = alignmentEvidence.stream().allMatch(ChimericAlignment::isLikelyInvertedDuplication);
        final boolean allInvolvesStrandSwitch = alignmentEvidence.stream().map(ca -> ca.strandSwitch).noneMatch(ss -> ss.equals(StrandSwitch.NO_SWITCH));
        final boolean allIndicateInsDel = alignmentEvidence.stream().allMatch(ChimericAlignment::isNeitherSimpleTranslocationNorIncompletePicture) &&
                                            alignmentEvidence.stream().map(ca -> ca.strandSwitch).allMatch(ss -> ss.equals(StrandSwitch.NO_SWITCH));

        final List<SvType> inferredType;
        final NovelAdjacencyReferenceLocations novelAdjacency = simpleNovelAdjacency.getNovelAdjacencyReferenceLocations();
        if ( allIndicateSimpleTransLoc ) {
            final Tuple2<BreakEndVariantType, BreakEndVariantType> orderedMatesForTranslocSuspect =
                    BreakEndVariantType.TransLocBND.getOrderedMates(novelAdjacency,
                            reference, referenceDictionary);
            inferredType = Arrays.asList(orderedMatesForTranslocSuspect._1, orderedMatesForTranslocSuspect._2);
        } else if ( allIndicateInvDup ) {
            inferredType = Collections.singletonList( new SimpleSVType.DuplicationInverted(novelAdjacency) );
        } else if ( allInvolvesStrandSwitch ) {
            final Tuple2<BreakEndVariantType, BreakEndVariantType> orderedMatesForInversionSuspect =
                    BreakEndVariantType.InvSuspectBND.getOrderedMates(novelAdjacency, reference);
            inferredType = Arrays.asList(orderedMatesForInversionSuspect._1, orderedMatesForInversionSuspect._2);
        } else if ( allIndicateInsDel ){
            inferredType = Collections.singletonList( InsDelVariantDetector.inferTypeFromNovelAdjacency(novelAdjacency) );
        } else {
            throw new GATKException
                    .ShouldNeverReachHereException("novel adjacency has its supporting chimeric alignments showing inconsistent behavior\n" +
                    simpleNovelAdjacency.toString());
        }

        return inferredType;
    }
}
