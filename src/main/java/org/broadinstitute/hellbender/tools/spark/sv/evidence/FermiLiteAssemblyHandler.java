package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexCache;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembler;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils.AutoDelete;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import scala.Tuple2;

import java.io.*;
import java.util.*;

/** This LocalAssemblyHandler aligns assembly contigs with BWA, along with some optional writing of intermediate results. */
public final class FermiLiteAssemblyHandler implements FindBreakpointEvidenceSpark.LocalAssemblyHandler {
    private static final long serialVersionUID = 1L;
    private final String alignerIndexFile;
    private final int maxFastqSize;
    private final String fastqDir;
    private final boolean writeGFAs;

    public FermiLiteAssemblyHandler( final String alignerIndexFile, final int maxFastqSize,
                                     final String fastqDir, final boolean writeGFAs ) {
        this.alignerIndexFile = alignerIndexFile;
        this.maxFastqSize = maxFastqSize;
        this.fastqDir = fastqDir;
        this.writeGFAs = writeGFAs;
    }

    @Override
    public AlignedAssemblyOrExcuse apply( final Tuple2<Integer, List<SVFastqUtils.FastqRead>> intervalAndReads ) {
        final int intervalID = intervalAndReads._1();
        final String assemblyName = AlignedAssemblyOrExcuse.formatAssemblyID(intervalID);
        final List<SVFastqUtils.FastqRead> readsList = intervalAndReads._2();

        // bail if the assembly will be too large
        final int fastqSize = readsList.stream().mapToInt(FastqRead -> FastqRead.getBases().length).sum();
        if ( fastqSize > maxFastqSize ) {
            return new AlignedAssemblyOrExcuse(intervalID, "no assembly -- too big (" + fastqSize + " bytes).");
        }

        // record the reads in the assembly as a FASTQ, if requested
        if ( fastqDir != null ) {
            final String fastqName = String.format("%s/%s.fastq", fastqDir, assemblyName);
            final ArrayList<SVFastqUtils.FastqRead> sortedReads = new ArrayList<>(readsList);
            sortedReads.sort(Comparator.comparing(SVFastqUtils.FastqRead::getHeader));
            SVFastqUtils.writeFastqFile(fastqName, sortedReads.iterator());
        }

        // assemble the reads
        final long timeStart = System.currentTimeMillis();
        final FermiLiteAssembly initialAssembly = new FermiLiteAssembler().createAssembly(readsList);
        final int secondsInAssembly = (int)((System.currentTimeMillis() - timeStart + 500)/1000);
        if ( initialAssembly.getNContigs() == 0 ) {
            return new AlignedAssemblyOrExcuse(intervalID, "no assembly -- no contigs produced by assembler.");
        }

        // patch up the assembly using read pairs to link contigs
        final FermiLiteAssembly assembly;
        if ( fastqDir == null ) {
                assembly = reviseAssembly(initialAssembly, assemblyName, readsList, null);
        } else {
            final String detailsFile = String.format("%s/%s.details", fastqDir, assemblyName);
            try ( final Writer writer = new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(detailsFile))) ) {
                assembly = reviseAssembly(initialAssembly, assemblyName, readsList, writer);
            }
            catch ( final IOException ioe ) {
                throw new GATKException("Can't write "+detailsFile, ioe);
            }
        }

        // record the assembly as a GFA, if requested
        if ( fastqDir != null && writeGFAs ) {
            final String gfaName =  String.format("%s/%s.gfa", fastqDir, assemblyName);
            try ( final OutputStream os = BucketUtils.createFile(gfaName) ) {
                assembly.writeGFA(os);
            }
            catch ( final IOException ioe ) {
                throw new GATKException("Can't write "+gfaName, ioe);
            }
        }

        // align the assembled contigs to the genomic reference
        final AlignedAssemblyOrExcuse result;
        try ( final BwaMemAligner aligner = new BwaMemAligner(BwaMemIndexCache.getInstance(alignerIndexFile)) ) {
            aligner.setIntraCtgOptions();
            final List<byte[]> sequences =
                    assembly.getContigs().stream()
                            .map(FermiLiteAssembly.Contig::getSequence)
                            .collect(SVUtils.arrayListCollector(assembly.getNContigs()));
            final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(sequences);
            result = new AlignedAssemblyOrExcuse(intervalID, assembly, secondsInAssembly, alignments);
        }

        return result;
    }

    @VisibleForTesting
    static FermiLiteAssembly reviseAssembly( final FermiLiteAssembly assembly,
                                             final String assemblyName,
                                             final List<SVFastqUtils.FastqRead> readsList,
                                             final Writer writer ) {
        final int nContigs = assembly.getNContigs();
        if ( nContigs == 0 ) return assembly;

        // alignments of reads in readsList onto assembly contigs
        final List<List<BwaMemAlignment>> alignments;

        final String tmpDir = IOUtils.tempDir("fbes", null).getPath();
        try ( final AutoDelete tmpDirAD = new AutoDelete(tmpDir) ) {
            final String fastaFile = String.format("%s/%s.fasta", tmpDirAD.getPath(), assemblyName);
            try ( final AutoDelete fastaFileAD = new AutoDelete(fastaFile) ) {

                // write a temporary FASTA file with all the assembled contigs
                try ( final BufferedOutputStream os =
                              new BufferedOutputStream(BucketUtils.createFile(fastaFileAD.getPath())) ) {
                    for ( int id = 0; id != nContigs; ++id ) {
                        final String seqName = assemblyName + "." + id;
                        final byte[] line1 = (">" + seqName + "\n").getBytes();
                        final byte[] seq = assembly.getContigs().get(id).getSequence();
                        os.write(line1);
                        os.write(seq);
                        os.write('\n');
                    }
                }
                catch ( final IOException ioe ) {
                    throw new GATKException("Unable to write fasta file of contigs from assembly " + assemblyName, ioe);
                }

                // create a BWA-mem index from the assembled contigs
                final String imageFile = String.format("%s/%s.img", tmpDirAD.getPath(), assemblyName);
                try ( final AutoDelete imageFileAD = new AutoDelete(imageFile) ) {
                    BwaMemIndex.createIndexImageFromFastaFile(fastaFileAD.getPath(), imageFileAD.getPath());

                    // align the reads that were assembled onto the assembled contigs
                    try ( final BwaMemIndex assemblyIndex = new BwaMemIndex(imageFileAD.getPath());
                          final BwaMemAligner aligner = new BwaMemAligner(assemblyIndex) ) {
                        aligner.setIntraCtgOptions();
                        alignments = aligner.alignSeqs(readsList, SVFastqUtils.FastqRead::getBases);
                    }
                }
            }
        }
        catch ( final IOException ioe ) {
            throw new GATKException("Unable to clean up temporary file.", ioe);
        }

        if ( writer != null ) {
            try {
                writer.write("Initial assembly:\n");
                writer.write("id\tlen\tnReads\tseq\n");
                final HashMap<FermiLiteAssembly.Contig, Integer> idMap = new HashMap<>(SVUtils.hashMapCapacity(nContigs));
                int id = 0;
                final List<FermiLiteAssembly.Contig> contigs = assembly.getContigs();
                for ( final FermiLiteAssembly.Contig contig : contigs ) {
                    idMap.put(contig, id++);
                }
                for ( final FermiLiteAssembly.Contig contig : contigs ) {
                    final int contigId = idMap.get(contig);
                    final String sequence = new String(contig.getSequence());
                    writer.write(contigId + "\t" + sequence.length() + "\t" + contig.getNSupportingReads() + "\t");
                    writer.write(sequence.substring(0,30));
                    writer.write("...");
                    writer.write(sequence.substring(sequence.length()-30));
                    writer.write('\n');
                    for ( final FermiLiteAssembly.Connection connection : contig.getConnections() ) {
                        final int targetId = idMap.get(connection.getTarget());
                        final int overlapLen = connection.getOverlapLen();
                        writer.write("\t" + (connection.isRC() ? "-" : "+") + contigId +
                                "--(" + overlapLen + ")-->" +
                                (connection.isTargetRC() ? "-" : "+") + targetId + "\n");
                    }
                    final int nReads = readsList.size();
                    for ( int idx = 0; idx < nReads; idx += 2 ) {
                        final List<BwaMemAlignment> alignList1 = alignments.get(idx);
                        final List<BwaMemAlignment> alignList2 = alignments.get(idx + 1);
                        for ( final BwaMemAlignment alignment1 : alignList1 ) {
                            if ( SAMFlag.READ_UNMAPPED.isSet(alignment1.getSamFlag()) ) continue;
                            for ( final BwaMemAlignment alignment2 : alignList2 ) {
                                if ( SAMFlag.READ_UNMAPPED.isSet(alignment2.getSamFlag()) ) continue;
                                if ( alignment1.getRefId() != contigId && alignment2.getRefId() != contigId ) continue;
                                final boolean isRC1 = SAMFlag.READ_REVERSE_STRAND.isSet(alignment1.getSamFlag());
                                // notice that the following flag is inverted -- we expect pairs to be on opposite strands
                                final boolean isRC2 = !SAMFlag.READ_REVERSE_STRAND.isSet(alignment2.getSamFlag());
                                if ( alignment1.getRefId() != alignment2.getRefId() || isRC1 != isRC2 ) {
                                    writer.write("\t" + readsList.get(idx).getName() + "\t" +
                                            (isRC1 ? "-" : "+") + alignment1.getRefId() + ":" + alignment1.getRefStart() +
                                            "\t" + alignment1.getCigar() + "\t" +
                                            (isRC2 ? "-" : "+") + alignment2.getRefId() + ":" + alignment2.getRefStart() +
                                            "\t" + alignment2.getCigar() + "\n");
                                }
                            }
                        }
                    }
                }
            }
            catch ( final IOException ioe ) {
                throw new GATKException("Can't write assembly details file for assembly "+assemblyName, ioe);
            }
        }

        // find read pairs that are NOT innies on the same contig (i.e., pairs that connect the graph of assembled contigs)
        // and count the number of pairs that corroborate each connection
        final Map<Link, Integer> linkCounts = new HashMap<>();
        final int nReads = readsList.size();
        for ( int idx = 0; idx < nReads; idx += 2 ) {
            final List<BwaMemAlignment> alignList1 = alignments.get(idx);
            final List<BwaMemAlignment> alignList2 = alignments.get(idx + 1);
            for ( final BwaMemAlignment alignment1 : alignList1 ) {
                if ( SAMFlag.READ_UNMAPPED.isSet(alignment1.getSamFlag()) ) continue;
                for ( final BwaMemAlignment alignment2 : alignList2 ) {
                    if ( SAMFlag.READ_UNMAPPED.isSet(alignment2.getSamFlag()) ) continue;
                    final boolean isRC1 = SAMFlag.READ_REVERSE_STRAND.isSet(alignment1.getSamFlag());
                    // notice that the following flag is inverted -- we expect pairs to be on opposite strands
                    final boolean isRC2 = !SAMFlag.READ_REVERSE_STRAND.isSet(alignment2.getSamFlag());
                    if ( alignment1.getRefId() != alignment2.getRefId() || isRC1 != isRC2 ) {
                        final Link link =
                                new Link(new ContigStrand(alignment1.getRefId(), isRC1),
                                         new ContigStrand(alignment2.getRefId(), isRC2));
                        linkCounts.merge(link, 1, Integer::sum);
                    }
                }
            }
        }

        // we'll use a magic number for this proof of concept
        final int MIN_LINK_COUNT = 8;
        // TODO: find a principled way of determining what this number should be
        final Iterator<Map.Entry<Link, Integer>> linkItr = linkCounts.entrySet().iterator();
        while ( linkItr.hasNext() ) {
            final Map.Entry<Link, Integer> entry = linkItr.next();
            if ( entry.getValue() < MIN_LINK_COUNT ) {
                linkItr.remove();
                continue;
            }

            // for this proof of concept, we'll connect only contigs that are corroborated by a
            // FermiLiteAssembly.Connection, and we'll assume that the assembler got the overlap correct.
            // we'll ignore underlaps (i.e., connections with negative overlaps) for now.
            final Link link = entry.getKey();
            final FermiLiteAssembly.Contig sourceContig = assembly.getContig(link.getSourceId());
            final FermiLiteAssembly.Contig targetContig = assembly.getContig(link.getTargetId());
            boolean foundConnection = false;
            for ( final FermiLiteAssembly.Connection connection : sourceContig.getConnections() ) {
                if ( connection.getTarget() == targetContig &&
                        connection.isRC() == link.isSourceRC() &&
                        connection.isTargetRC() == link.isTargetRC() &&
                        connection.getOverlapLen() >= 0 ) {
                    // ugly hack -- replace the count with the overlap length from the connection, which we'll need below
                    entry.setValue(connection.getOverlapLen());
                    foundConnection = true;
                    break;
                }

            }
            if ( !foundConnection ) linkItr.remove();
        }

        if ( linkCounts.isEmpty() ) {
            return assembly;
        }

        if ( writer != null ) {
            try {
                writer.write("\nJoining " + linkCounts.size() + " contig pairs:\n");
            }
            catch ( final IOException ioe ) {
                throw new GATKException("Can't write assembly details file for assembly "+assemblyName, ioe);
            }
        }
        final List<FermiLiteAssembly.Contig> contigList = new ArrayList<>(nContigs);
        final Set<Integer> representedContigs = new HashSet<>(SVUtils.hashMapCapacity(nContigs));
        for ( final Map.Entry<Link, Integer> entry : linkCounts.entrySet() ) {
            final Link link = entry.getKey();
            if ( writer != null ) {
                try {
                    writer.write("\t" + (link.isSourceRC() ? "-" : "+") + link.getSourceId() +
                            "--(" + entry.getValue() + ")-->" +
                            (link.isTargetRC() ? "-" : "+") + link.getTargetId() + "\n");
                }
                catch ( final IOException ioe ) {
                    throw new GATKException("Can't write assembly details file for assembly "+assemblyName, ioe);
                }
            }
            final int overlapLen = entry.getValue();
            final FermiLiteAssembly.Contig sourceContig =
                    link.isSourceRC() ?
                            rcContig(assembly.getContig(link.getSourceId())) :
                            assembly.getContig(link.getSourceId());
            final FermiLiteAssembly.Contig targetContig =
                    link.isTargetRC() ?
                            rcContig(assembly.getContig(link.getTargetId())) :
                            assembly.getContig(link.getTargetId());
            final int sourceLength = sourceContig.getSequence().length;
            final int targetLength = targetContig.getSequence().length - overlapLen;
            final int contigLength = sourceLength + targetLength;
            final byte[] sequence = new byte[contigLength];
            System.arraycopy(sourceContig.getSequence(), 0, sequence, 0, sourceLength);
            System.arraycopy(targetContig.getSequence(), overlapLen, sequence, sourceLength, targetLength);
            final byte[] perBaseCoverage = new byte[contigLength];
            System.arraycopy(sourceContig.getPerBaseCoverage(), 0, perBaseCoverage, 0, sourceLength);
            System.arraycopy(targetContig.getPerBaseCoverage(), overlapLen, perBaseCoverage, sourceLength, targetLength);
            final int nSupportingReads = sourceContig.getNSupportingReads() + targetContig.getNSupportingReads();
            contigList.add(new FermiLiteAssembly.Contig(sequence, perBaseCoverage, nSupportingReads));
            representedContigs.add(link.getSourceId());
            representedContigs.add(link.getTargetId());
        }

        for ( int id = 0; id != nContigs; ++id ) {
            if ( !representedContigs.contains(id) ) {
                final FermiLiteAssembly.Contig contig = assembly.getContig(id);
                final byte[] sequence = contig.getSequence();
                final byte[] perBaseCoverage = contig.getPerBaseCoverage();
                final int nSupportingReads = contig.getNSupportingReads();
                contigList.add(new FermiLiteAssembly.Contig(sequence, perBaseCoverage, nSupportingReads));
            }
        }

        return new FermiLiteAssembly(contigList);
    }

    private static FermiLiteAssembly.Contig rcContig( final FermiLiteAssembly.Contig contig ) {
        final byte[] sequence = BaseUtils.simpleReverseComplement(contig.getSequence());
        final byte[] perBaseCoverage = Arrays.copyOf(contig.getPerBaseCoverage(), contig.getPerBaseCoverage().length);
        SequenceUtil.reverse(perBaseCoverage, 0, perBaseCoverage.length);
        return new FermiLiteAssembly.Contig(sequence, perBaseCoverage, contig.getNSupportingReads());
    }

    private static final class ContigStrand implements Comparable<ContigStrand> {
        private final int contigId;
        private final boolean isRC;

        public ContigStrand( final int contigId, final boolean isRC ) {
            this.contigId = contigId;
            this.isRC = isRC;
        }

        public int getId() { return contigId; }
        public boolean isRC() { return isRC; }

        public ContigStrand rc() { return new ContigStrand(contigId, !isRC); }

        @Override public boolean equals( final Object obj ) {
            return obj instanceof ContigStrand && equals((ContigStrand) obj);
        }

        public boolean equals( final ContigStrand that ) {
            return contigId == that.contigId && isRC == that.isRC;
        }

        @Override public int hashCode() {
            int hashVal = 113;
            hashVal = 47 * (hashVal + contigId);
            return 47 * (hashVal + (isRC ? 3 : 5));
        }

        @Override public int compareTo( final ContigStrand that ) {
            int result = Integer.compare(contigId, that.contigId);
            if ( result == 0 ) result = Boolean.compare(isRC, that.isRC);
            return result;
        }
    }

    private static final class Link {
        private final ContigStrand contigStrandSource;
        private final ContigStrand contigStrandTarget;

        /** Canonicalizing constructor: will swizzle to guarantee that id1 <= id2 */
        public Link(final ContigStrand contigStrandSource, final ContigStrand contigStrandTarget) {
            if ( contigStrandSource.getId() < contigStrandTarget.getId() ) {
                this.contigStrandSource = contigStrandSource;
                this.contigStrandTarget = contigStrandTarget;
            } else {
                this.contigStrandSource = contigStrandTarget.rc();
                this.contigStrandTarget = contigStrandSource.rc();
            }
        }

        public int getSourceId() { return contigStrandSource.getId(); }
        public boolean isSourceRC() { return contigStrandSource.isRC(); }
        public ContigStrand getSource() { return contigStrandSource; }
        public int getTargetId() { return contigStrandTarget.getId(); }
        public boolean isTargetRC() { return contigStrandTarget.isRC(); }
        public ContigStrand getTarget() { return contigStrandTarget; }

        @Override public boolean equals( final Object obj ) {
            return obj instanceof Link && equals((Link)obj);
        }

        public boolean equals( final Link link ) {
            return this == link ||
                    (contigStrandSource.equals(link.contigStrandSource) && contigStrandTarget.equals(link.contigStrandTarget));
        }

        @Override public int hashCode() {
            int hashVal = 113;
            hashVal = 47 * (hashVal + contigStrandSource.hashCode());
            return 47 * (hashVal + contigStrandTarget.hashCode());
        }
    }
}
