#include <cmath>
#include "cblas.h"
#include "MixtapeException.hpp"
#include "DiagonalGaussianHMMSufficientStats.hpp"
#include "DiagonalGaussianHMMEstep.hpp"

namespace Mixtape {

void gaussian_loglikelihood_diag(const float* __restrict__ sequence,
                                 const double* __restrict__ means,
                                 const double* __restrict__ variances,
                                 const double* __restrict__ means_over_variances,
                                 const double* __restrict__ means2_over_variances,
                                 const double* __restrict__ log_variances,
                                 const int n_observations,
                                 const int n_states, const int n_features,
                                 float* __restrict__ loglikelihoods)
{
    int t, i, j;
    double temp;
    double x;
    static const double log_M_2_PI = 1.8378770664093453; // np.log(2*np.pi)

    for (t = 0; t < n_observations; t++) {
        for (j = 0; j < n_states; j++) {
            temp = 0.0f;
            for (i = 0; i < n_features; i++) {
                x = sequence[t*n_features + i];
                temp += means2_over_variances[j*n_features + i]              \
                        - 2.0 * x * means_over_variances[j*n_features + i]   \
                        + (x * x) / variances[j*n_features + i]              \
                        + log_variances[j*n_features + i];
            }
            loglikelihoods[t*n_states + j] = -0.5 * (n_features * log_M_2_PI + temp);
        }
    }
}

//-------------------------------------------------------------------------------------------------------------//

DiagonalGaussianHMMEstep::DiagonalGaussianHMMEstep(const int numStates,
                                                   const DoubleArray2D& transmat,
                                                   const DoubleArray1D& startProb,
                                                   const int numFeatures,
                                                   const DoubleArray2D& means,
                                                   const DoubleArray2D& variances)
    : HMMEstep(numStates, transmat, startProb)
    , numFeatures_(numFeatures)
    , means_(means)
    , variances_(variances)
    , logVariances_(boost::extents[numStates][numFeatures])
    , meansOverVariances_(boost::extents[numStates][numFeatures])
    , means2OverVariances_(boost::extents[numStates][numFeatures])
{
    if (numStates_ != means.shape()[0] || numFeatures_ != means.shape()[1])
        throw MixtapeException("means has wrong shape");
    if (numStates_ != variances.shape()[0] || numFeatures_ != variances.shape()[1])
        throw MixtapeException("variances has wrong shape");

    typedef DoubleArray2D::index index;
    for (index i = 0; i < numStates_; i++) {
        for (index j = 0; j < numFeatures_; j++) {
            logVariances_[i][j] = log(variances[i][j]);
            meansOverVariances_[i][j] = means[i][j] / variances[i][j];
            means2OverVariances_[i][j] = (means[i][j] * means[i][j]) / variances[i][j];
        }
    }
}

FloatArray2D DiagonalGaussianHMMEstep::emissionLogLikelihood(const FloatArray2D& X) {
    FloatArray2D frameLogProb(boost::extents[X.shape()[0]][numStates_]);
    gaussian_loglikelihood_diag(&X[0][0], &means_[0][0], &variances_[0][0], &meansOverVariances_[0][0],
                                &means2OverVariances_[0][0], &logVariances_[0][0], X.shape()[0], numStates_,
                                numFeatures_, &frameLogProb[0][0]);
    return frameLogProb;
}

DiagonalGaussianHMMSufficientStats* DiagonalGaussianHMMEstep::initializeSufficientStats() {
    DiagonalGaussianHMMSufficientStats* stats = new DiagonalGaussianHMMSufficientStats(numStates_);
    return stats;
}

void DiagonalGaussianHMMEstep::accumulateSufficientStats(
    DiagonalGaussianHMMSufficientStats& stats, const FloatArray2D& seq, const FloatArray1D& frameLogProb,
    const FloatArray2D& posteriors, const DoubleArray2D& fwdLattice, const DoubleArray2D& bwdLattice)
{
    HMMEstep::accumulateSufficientStats(stats, seq, frameLogProb, posteriors, fwdLattice, bwdLattice);
    
    
    const float alpha = 1.0;
    const float beta = 1.0;
    const int length = seq.shape()[0];
    assert(seq.shape()[1] == numFeatures_);

    FloatArray2D obs(boost::extents[numStates_][numFeatures_]);
    FloatArray2D obs2(boost::extents[numStates_][numFeatures_]);
    FloatArray2D seq2(boost::extents[length][numFeatures_]);

    typedef FloatArray2D::index index;
    for (index i = 0; i < length; i++)
        for (index j = 0; j < length; j++)
            seq2[i][j] = seq[i][j];
    
    sgemm_("N", "T", &numFeatures_, &numStates_, &length, &alpha, &seq[0][0],
           &numFeatures_, &posteriors[0][0], &numStates_, &beta, &obs[0][0],
           &numFeatures_);
    sgemm_("N", "T", &numFeatures_, &numStates_, &length, &alpha, &seq2[0][0],
           &numFeatures_, &posteriors[0][0], &numStates_, &beta, &obs2[0][0],
           &numFeatures_);
    //for (k = 0; k < n_states; k++)
    //for (j = 0; j < sequence_lengths[i]; j++)
    //seq_post[k] += posteriors[j*n_states + k];
}




} // namespace
