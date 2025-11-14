package org.irlab.ecir25.nhst;

import org.irlab.ecir25.util.ParallelArrays;

import java.util.Arrays;
import java.util.stream.IntStream;

public class BHCorrection implements CorrectionProcedure {
  @Override
  public double[] correct(double[] unadjustedPvalues, double alpha) {
    int total_comparisons = unadjustedPvalues.length;

    double[] unadjusted = Arrays.copyOf(unadjustedPvalues, unadjustedPvalues.length);
    double[] positions = IntStream.range(0, unadjusted.length).mapToDouble(val -> (double) val).toArray();
    ParallelArrays.sort(unadjusted, positions, true);

    double[] newPvalues = new double[total_comparisons];
    // Step-up: find largest k s.t. p_(k) <= (k/m)*alpha and reject all i <= k.
    int k = -1;
    for (int i = total_comparisons - 1; i >= 0; i--) {
      double threshold = alpha * ((double) (i + 1) / total_comparisons);
      if (unadjusted[i] <= threshold) {
        k = i;
        break;
      }
    }
    if (k == -1) {
      for (int i = 0; i < total_comparisons; i++) newPvalues[(int) positions[i]] = 1; // accept all
    } else {
      for (int i = 0; i <= k; i++) newPvalues[(int) positions[i]] = 0; // reject 1..k
      for (int i = k + 1; i < total_comparisons; i++) newPvalues[(int) positions[i]] = 1; // accept rest
    }
    return newPvalues;
  }
}
