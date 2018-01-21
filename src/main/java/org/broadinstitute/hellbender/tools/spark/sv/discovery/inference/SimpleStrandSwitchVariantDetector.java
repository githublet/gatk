package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AnnotatedVariantProducer;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.BreakEndVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoveryInputData;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.utils.RDDUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;
import scala.Tuple3;

import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

// TODO: 1/21/18 all methods, except inferInvDupRange(), in this class will do nothing more than what new centralized class SimpleNovelAdjacencyInterpreter, hence will be removed
public final class SimpleStrandSwitchVariantDetector implements VariantDetectorFromLocalAssemblyContigAlignments {


    static final int MORE_RELAXED_ALIGNMENT_MIN_LENGTH = 30;
    static final int MORE_RELAXED_ALIGNMENT_MIN_MQ = 20;

    @Override
    public void inferSvAndWriteVCF(final JavaRDD<AssemblyContigWithFineTunedAlignments> assemblyContigs,
                                   final SvDiscoveryInputData svDiscoveryInputData) {

        final Broadcast<ReferenceMultiSource> referenceBroadcast = svDiscoveryInputData.referenceBroadcast;
        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputData.referenceSequenceDictionaryBroadcast;
        final String sampleId = svDiscoveryInputData.sampleId;
        final String outputPath = svDiscoveryInputData.outputPath;
        final Logger toolLogger = svDiscoveryInputData.toolLogger;

        svDiscoveryInputData.toolLogger.info(assemblyContigs.count() +
                " chimeras indicating either 1) simple strand-switch breakpoints, or 2) inverted duplication.");

        // TODO: 11/23/17 take insertion mappings from the input and add them to VC
        // split between suspected inv dup VS strand-switch breakpoints
        // logic flow: split the input reads into two classes--those judged by IsLikelyInvertedDuplication are likely invdup and those aren't
        //             finally send the two split reads down different path, one for inv dup and one for BND records
        final Tuple2<JavaRDD<AlignedContig>, JavaRDD<AlignedContig>> invDupAndStrandSwitchBreakpoints =
                RDDUtils.split(assemblyContigs.map( AssemblyContigWithFineTunedAlignments::getSourceContig ),
                        contig -> ChimericAlignment.isLikelyInvertedDuplication(contig.alignmentIntervals.get(0),
                                contig.alignmentIntervals.get(1)), false);

        final JavaRDD<VariantContext> simpleStrandSwitchBkpts =
                dealWithSimpleStrandSwitchBkpts(invDupAndStrandSwitchBreakpoints._2, referenceBroadcast,
                        referenceSequenceDictionaryBroadcast, toolLogger, sampleId);

        SVVCFWriter.writeVCF(simpleStrandSwitchBkpts.collect(), outputPath.replace(".vcf", "_simpleSS.vcf"),
                referenceSequenceDictionaryBroadcast.getValue(), toolLogger);
        simpleStrandSwitchBkpts.unpersist();

        final JavaRDD<VariantContext> invDups =
                dealWithSuspectedInvDup(invDupAndStrandSwitchBreakpoints._1, referenceBroadcast,
                        referenceSequenceDictionaryBroadcast, toolLogger, sampleId);

        SVVCFWriter.writeVCF(invDups.collect(), outputPath.replace(".vcf", "_invDup.vcf"),
                referenceSequenceDictionaryBroadcast.getValue(), toolLogger);
    }

    // =================================================================================================================

    /**
     * @throws IllegalArgumentException if the assumption that the input aligned assembly contig has 2 alignments
     *                                  mapped to the same chr with strand switch is invalid
     */
    private static Tuple2<ChimericAlignment, byte[]> convertAlignmentIntervalsToChimericAlignment
    (final AlignedContig contigWith2AIMappedToSameChrAndStrandSwitch, final Broadcast<SAMSequenceDictionary> referenceDictionary) {
        Utils.validateArg(AssemblyContigAlignmentSignatureClassifier.indicatesIntraChrStrandSwitchBkpts(contigWith2AIMappedToSameChrAndStrandSwitch),
                "assumption that input aligned assembly contig has 2 alignments mapped to the same chr with strand switch is invalid.\n" +
                        contigWith2AIMappedToSameChrAndStrandSwitch.toString());

        final AlignmentInterval intervalOne = contigWith2AIMappedToSameChrAndStrandSwitch.alignmentIntervals.get(0),
                                intervalTwo = contigWith2AIMappedToSameChrAndStrandSwitch.alignmentIntervals.get(1);

        return new Tuple2<>(new ChimericAlignment(intervalOne, intervalTwo, EMPTY_INSERTION_MAPPINGS,
                contigWith2AIMappedToSameChrAndStrandSwitch.contigName, referenceDictionary.getValue()), contigWith2AIMappedToSameChrAndStrandSwitch.contigSequence);
    }

    // TODO: 1/21/18 to be replaced with corresponding method in new centralized class SimpleNovelAdjacencyInterpreter
    /**
     * Roughly similar to {@link ChimericAlignment#nextAlignmentMayBeInsertion(AlignmentInterval, AlignmentInterval, Integer, Integer, boolean)}:
     *  1) either alignment may have very low mapping quality (a more relaxed mapping quality threshold);
     *  2) either alignment may consume only a "short" part of the contig, or if assuming that the alignment consumes
     *     roughly the same amount of ref bases and read bases, has isAlignment that is too short
     */
    static boolean splitPairStrongEnoughEvidenceForCA(final AlignmentInterval intervalOne,
                                                      final AlignmentInterval intervalTwo,
                                                      final int mapQThresholdInclusive,
                                                      final int alignmentLengthThresholdInclusive) {

        if (intervalOne.mapQual < mapQThresholdInclusive || intervalTwo.mapQual < mapQThresholdInclusive)
            return false;

        final int overlap = AlignmentInterval.overlapOnContig(intervalOne, intervalTwo);

        final int x = intervalOne.endInAssembledContig - intervalOne.startInAssembledContig + 1,
                  y = intervalTwo.endInAssembledContig - intervalTwo.startInAssembledContig + 1;

        return Math.min(x - overlap, y - overlap) >= alignmentLengthThresholdInclusive;
    }

    // =================================================================================================================

    // TODO: 1/21/18 to be replaced with corresponding method in new centralized class SimpleNovelAdjacencyInterpreter
    // workflow manager for simple strand-switch alignment contigs
    private JavaRDD<VariantContext> dealWithSimpleStrandSwitchBkpts(final JavaRDD<AlignedContig> contigs,
                                                                    final Broadcast<ReferenceMultiSource> broadcastReference,
                                                                    final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary,
                                                                    final Logger toolLogger,
                                                                    final String sampleId) {

        final JavaPairRDD<ChimericAlignment, byte[]> simpleStrandSwitchBkpts =
                contigs
                        .filter(tig ->
                                splitPairStrongEnoughEvidenceForCA(tig.alignmentIntervals.get(0), tig.alignmentIntervals.get(1),
                                        MORE_RELAXED_ALIGNMENT_MIN_MQ,  MORE_RELAXED_ALIGNMENT_MIN_LENGTH))
                        .mapToPair(tig -> convertAlignmentIntervalsToChimericAlignment(tig, broadcastSequenceDictionary)).cache();

        toolLogger.info(simpleStrandSwitchBkpts.count() + " chimeras indicating simple strand-switch breakpoints.");

        return simpleStrandSwitchBkpts
                .mapToPair(pair -> new Tuple2<>(new NovelAdjacencyReferenceLocations(pair._1, pair._2, broadcastSequenceDictionary.getValue()), pair._1))
                .groupByKey()
                .mapToPair(noveltyAndEvidence -> inferBNDType(noveltyAndEvidence, broadcastReference.getValue()))
                .flatMap(noveltyTypeAndEvidence ->
                        AnnotatedVariantProducer
                                .produceAnnotatedBNDmatesVcFromNovelAdjacency(
                                        noveltyTypeAndEvidence._1,
                                        noveltyTypeAndEvidence._2._1,
                                        noveltyTypeAndEvidence._2._2,
                                        broadcastReference,
                                        broadcastSequenceDictionary,
                                        sampleId).iterator());
    }

    // TODO: 1/21/18 new centralized class SimpleNovelAdjacencyInterpreter just need to implement the correct logic to trigger calling the following logic
    private static Tuple2<NovelAdjacencyReferenceLocations, Tuple2<Tuple2<BreakEndVariantType, BreakEndVariantType>, Iterable<ChimericAlignment>>>
    inferBNDType(final Tuple2<NovelAdjacencyReferenceLocations, Iterable<ChimericAlignment>> noveltyAndEvidence,
                 final ReferenceMultiSource reference) {

        final NovelAdjacencyReferenceLocations novelAdjacency = noveltyAndEvidence._1;
        final Iterable<ChimericAlignment> chimericAlignments = noveltyAndEvidence._2;
        return new Tuple2<>(novelAdjacency,
                new Tuple2<>(BreakEndVariantType.InvSuspectBND.getOrderedMates(novelAdjacency, reference), chimericAlignments));
    }

    // =================================================================================================================

    // TODO: 1/21/18 to be replaced with corresponding method in new centralized class SimpleNovelAdjacencyInterpreter
    private JavaRDD<VariantContext> dealWithSuspectedInvDup(final JavaRDD<AlignedContig> contigs,
                                                            final Broadcast<ReferenceMultiSource> broadcastReference,
                                                            final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary,
                                                            final Logger toolLogger,
                                                            final String sampleId) {

        final JavaPairRDD<ChimericAlignment, byte[]> invDupSuspects =
                contigs
                        .filter(tig ->
                                splitPairStrongEnoughEvidenceForCA(tig.alignmentIntervals.get(0), tig.alignmentIntervals.get(1),
                                        MORE_RELAXED_ALIGNMENT_MIN_MQ,  MORE_RELAXED_ALIGNMENT_MIN_LENGTH))
                        .mapToPair(tig -> convertAlignmentIntervalsToChimericAlignment(tig, broadcastSequenceDictionary)).cache();

        toolLogger.info(invDupSuspects.count() + " chimera indicating inverted duplication");

        return invDupSuspects
                .mapToPair(pair -> new Tuple2<>(new NovelAdjacencyReferenceLocations(pair._1, pair._2, broadcastSequenceDictionary.getValue()), new Tuple2<>(pair._1, pair._2)))
                .groupByKey()
                .flatMapToPair(SimpleStrandSwitchVariantDetector::inferInvDupRange)
                .map(noveltyTypeAndAltSeqAndEvidence -> {
                            final NovelAdjacencyReferenceLocations novelAdjacency = noveltyTypeAndAltSeqAndEvidence._1._1();
                            final SimpleSVType.DuplicationInverted invDup = noveltyTypeAndAltSeqAndEvidence._1._2();
                            final byte[] altHaplotypeSeq = noveltyTypeAndAltSeqAndEvidence._1._3();

                            final List<ChimericAlignment> evidence = noveltyTypeAndAltSeqAndEvidence._2;
                            return AnnotatedVariantProducer
                                    .produceAnnotatedVcFromInferredTypeAndRefLocations(
                                            novelAdjacency.leftJustifiedLeftRefLoc,
                                            novelAdjacency.leftJustifiedRightRefLoc.getStart(),
                                            novelAdjacency.complication,
                                            invDup,
                                            altHaplotypeSeq,
                                            evidence,
                                            broadcastReference,
                                            broadcastSequenceDictionary,
                                            null,
                                            sampleId);

                });
    }

    // TODO: 1/21/18 the only method that needs to be moved into the new centralized class
    private static Iterator<Tuple2<Tuple3<NovelAdjacencyReferenceLocations, SimpleSVType.DuplicationInverted, byte[]>, List<ChimericAlignment>>>
    inferInvDupRange(final Tuple2<NovelAdjacencyReferenceLocations, Iterable<Tuple2<ChimericAlignment, byte[]>>> noveltyAndEvidence) {

        final NovelAdjacencyReferenceLocations novelAdjacency = noveltyAndEvidence._1;
        final SimpleSVType.DuplicationInverted duplicationInverted = new SimpleSVType.DuplicationInverted(novelAdjacency);

        // doing this because the same novel adjacency reference locations might be induced by different (probably only slightly) alt haplotypes, so a single group by NARL is not enough
        final Iterable<Tuple2<ChimericAlignment, byte[]>> chimeraAndContigSeq = noveltyAndEvidence._2;
        final Set<Map.Entry<byte[], List<ChimericAlignment>>> alignmentEvidenceGroupedByAltHaplotypeSequence =
                Utils.stream(chimeraAndContigSeq)
                        .collect(
                                Collectors.groupingBy(caAndSeq ->
                                                novelAdjacency.complication.extractAltHaplotypeForInvDup(caAndSeq._1, caAndSeq._2),
                                        Collectors.mapping(caAndSeq -> caAndSeq._1, Collectors.toList())
                                )
                        )
                        .entrySet();

        return alignmentEvidenceGroupedByAltHaplotypeSequence.stream()
                .map(entry -> new Tuple2<>(new Tuple3<>(novelAdjacency, duplicationInverted, entry.getKey()),
                        entry.getValue()))
                .iterator();
    }
}
