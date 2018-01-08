package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

/**
 * Integration tests for {@link NeuralNetStreamingExecutor}.
 * Created by sam on 1/8/18.
 */
public class NeuralNetStreamingExecutorIntegrationTest extends CommandLineProgramTest {
    private static String architectureHD5 = packageMainResourcesDir + "tools/walkers/vqsr/cnn_1d_annotations.hd5";
    private static final String space = " ";
    /**
     * Run the tool on sample VCF.
     */
    @Test(groups = {"python"})
    public void testInference() {
        final String inputVCF = largeFileTestDir + "VQSR/recalibrated_chr20_start.vcf";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .addArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, "%s")
                .addArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .addArgument("architecture", architectureHD5)
                .addArgument(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false");

        final IntegrationTestSpec spec = new IntegrationTestSpec(argsBuilder.toString(),
                Arrays.asList(largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf"));
        try {
            spec.executeTest("testInference", this);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
