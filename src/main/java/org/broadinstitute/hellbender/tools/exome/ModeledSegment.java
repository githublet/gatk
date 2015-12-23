package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.Arrays;

/**
 * Class for storing a segment that was generated by a model (e.g., GATK CNV), though we do not necessarily know which model ahead of time.
 *
 * @author lichtens &lt;lichtens@broadinstitute.org&gt;
 */
public class ModeledSegment implements Locatable {

    public static final String NO_CALL = "";
    protected SimpleInterval simpleInterval;

    /**
     * Segment mean in log copy ratio space as determined by whatever model generated this segment.
     */
    private double segmentMean;

    private String call;

    private long originalProbeCount;

    public ModeledSegment(final SimpleInterval interval, final String call, final long originalProbeCount, final double segmentMeanInLogCR) {
        this.simpleInterval = Utils.nonNull(interval, "The input interval cannot be null");
        this.call = Utils.nonNull(call, String.format("The input call cannot be null.  Use empty string, instead (\"%s\")", NO_CALL));
        this.segmentMean = ParamUtils.isFinite(segmentMeanInLogCR, "Segment Mean must be finite.");
        this.originalProbeCount = ParamUtils.isPositiveOrZero(originalProbeCount, "Number of original probes must be positive or zero.");
    }

    public ModeledSegment(final SimpleInterval interval, final long originalProbeCount, final double segmentMean) {
        this(interval, NO_CALL, originalProbeCount, segmentMean);
    }

    public double getSegmentMean() {
        return segmentMean;
    }

    public void setSegmentMean(final double segmentMean) {
        this.segmentMean = segmentMean;
    }

    public double getSegmentMeanInCRSpace() {
        return Math.pow(2, segmentMean);
    }

    public void setSegmentMeanInCRSpace(final double segmentMeanInCRSpace) {
        this.segmentMean = Math.log(segmentMeanInCRSpace)/Math.log(2);
    }

    @Override
    public String getContig() {return simpleInterval.getContig(); }

    @Override
    public int getStart() {return simpleInterval.getStart(); }

    @Override
    public int getEnd() {return simpleInterval.getEnd(); }

    public SimpleInterval getSimpleInterval() {
        return simpleInterval;
    }

    public void setSimpleInterval(final SimpleInterval simpleInterval) {
        this.simpleInterval = Utils.nonNull(simpleInterval, "The input interval cannot be null");
    }

    /**
     * Returns the call.  Returns null for uncalled segments.
     *
     * @return never {@code null}
     */
    public String getCall() {
        return this.call;
    }

    /**
     * Sets the call.
     */
    public void setCall(final String call) {
        this.call = Utils.nonNull(call, String.format("The input call cannot be null.  Use empty string, (\"%s\")", NO_CALL));
    }

    public long getOriginalProbeCount() {
        return originalProbeCount;
    }

    public void setOriginalProbeCount(final long originalProbeCount) {
        this.originalProbeCount = ParamUtils.isPositiveOrZero(originalProbeCount, "Number of Probes must be positive or zero.");
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof ModeledSegment)) {
            return false;
        }

        final ModeledSegment modeledSegment = (ModeledSegment) o;
        return simpleInterval.equals(modeledSegment.simpleInterval) && call.equals(modeledSegment.call)
                && Math.abs(segmentMean - modeledSegment.segmentMean) < 2 * Math.ulp(segmentMean)
                && originalProbeCount == modeledSegment.originalProbeCount;
    }

    @Override
    public int hashCode() {
        final int[] hashes = new int[]{simpleInterval.hashCode(), call.hashCode(), Double.hashCode(segmentMean),
                Long.hashCode(originalProbeCount)};
        return Arrays.hashCode(hashes);
    }
}
