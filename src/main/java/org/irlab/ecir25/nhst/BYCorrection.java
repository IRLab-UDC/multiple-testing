package org.irlab.ecir25.nhst;

import org.irlab.ecir25.util.ParallelArrays;

import java.util.Arrays;
import java.util.stream.IntStream;

public class BYCorrection implements CorrectionProcedure {
  @Override
  public double[] correct(double[] unadjustedPvalues, double alpha) {
    int total_comparisons = unadjustedPvalues.length;

    double[] unadjusted = Arrays.copyOf(unadjustedPvalues, unadjustedPvalues.length);
    double[] positions = IntStream.range(0, unadjusted.length).mapToDouble(val -> (double) val).toArray();
    ParallelArrays.sort(unadjusted, positions, true);

    double harmonic = 0;
    for (int i = 1; i <= total_comparisons; i++) {
      harmonic += (double) 1 / i;
    }

    harmonic = harmonic * total_comparisons;

    double[] newPvalues = new double[total_comparisons];
    int k = -1;
    for (int i = total_comparisons - 1; i >= 0; i--) {
      double threshold = alpha * ((double) (i + 1) / harmonic);
      if (unadjusted[i] <= threshold) {
        k = i;
        break;
      }
    }
    if (k == -1) {
      for (int i = 0; i < total_comparisons; i++) newPvalues[(int) positions[i]] = 1;
    } else {
      for (int i = 0; i <= k; i++) newPvalues[(int) positions[i]] = 0;
      for (int i = k + 1; i < total_comparisons; i++) newPvalues[(int) positions[i]] = 1;
    }

    return newPvalues;
  }
}
