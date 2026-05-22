    def _approx_bound(self, X, doc_topic_distr, sub_sampling):
        """Estimate the variational bound.

        Estimate the variational bound over "all documents" using only the
        documents passed in as X. Since log-likelihood of each word cannot
        be computed directly, we use this bound to estimate it.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Document word matrix.

        doc_topic_distr : ndarray of shape (n_samples, n_components)
            Document topic distribution. In the literature, this is called
            gamma.

        sub_sampling : bool, default=False
            Compensate for subsampling of documents.
            It is used in calculate bound in online learning.

        Returns
        -------
        score : float

        """

        def _loglikelihood(prior, distr, dirichlet_distr, size):
            # calculate log-likelihood
            score = np.sum((prior - distr) * dirichlet_distr)
            score += np.sum(gammaln(distr) - gammaln(prior))
            score += np.sum(gammaln(prior * size) - gammaln(np.sum(distr, 1)))
            return score

        is_sparse_x = sp.issparse(X)
        n_samples, n_components = doc_topic_distr.shape
        n_features = self.components_.shape[1]
        score = 0

        dirichlet_doc_topic = _dirichlet_expectation_2d(doc_topic_distr)
        dirichlet_component_ = _dirichlet_expectation_2d(self.components_)
        doc_topic_prior = self.doc_topic_prior_
        topic_word_prior = self.topic_word_prior_

        if is_sparse_x:
            X_data = X.data
            X_indices = X.indices
            X_indptr = X.indptr

        # E[log p(docs | theta, beta)]
        for idx_d in range(0, n_samples):
            if is_sparse_x:
                ids = X_indices[X_indptr[idx_d] : X_indptr[idx_d + 1]]
                cnts = X_data[X_indptr[idx_d] : X_indptr[idx_d + 1]]
            else:
                ids = np.nonzero(X[idx_d, :])[0]
                cnts = X[idx_d, ids]
            temp = (
                dirichlet_doc_topic[idx_d, :, np.newaxis] + dirichlet_component_[:, ids]
            )
            norm_phi = logsumexp(temp, axis=0)
            score += np.dot(cnts, norm_phi)

        # compute E[log p(theta | alpha) - log q(theta | gamma)]
        score += _loglikelihood(
            doc_topic_prior, doc_topic_distr, dirichlet_doc_topic, self.n_components
        )

        # Compensate for the subsampling of the population of documents
        if sub_sampling:
            doc_ratio = float(self.total_samples) / n_samples
            score *= doc_ratio

        # E[log p(beta | eta) - log q (beta | lambda)]
        score += _loglikelihood(
            topic_word_prior, self.components_, dirichlet_component_, n_features
        )

        return score