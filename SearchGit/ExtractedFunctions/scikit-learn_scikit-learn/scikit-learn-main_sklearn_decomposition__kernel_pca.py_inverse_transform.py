    def inverse_transform(self, X):
        """Transform X back to original space.

        ``inverse_transform`` approximates the inverse transformation using
        a learned pre-image. The pre-image is learned by kernel ridge
        regression of the original data on their low-dimensional representation
        vectors.

        .. note:
            :meth:`~sklearn.decomposition.fit` internally uses a centered
            kernel. As the centered kernel no longer contains the information
            of the mean of kernel features, such information is not taken into
            account in reconstruction.

        .. note::
            When users want to compute inverse transformation for 'linear'
            kernel, it is recommended that they use
            :class:`~sklearn.decomposition.PCA` instead. Unlike
            :class:`~sklearn.decomposition.PCA`,
            :class:`~sklearn.decomposition.KernelPCA`'s ``inverse_transform``
            does not reconstruct the mean of data when 'linear' kernel is used
            due to the use of centered kernel.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_components)
            Training vector, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        Returns
        -------
        X_original : ndarray of shape (n_samples, n_features)
            Original data, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        References
        ----------
        `Bakır, Gökhan H., Jason Weston, and Bernhard Schölkopf.
        "Learning to find pre-images."
        Advances in neural information processing systems 16 (2004): 449-456.
        <https://papers.nips.cc/paper/2003/file/ac1ad983e08ad3304a97e147f522747e-Paper.pdf>`_
        """
        if not self.fit_inverse_transform:
            raise NotFittedError(
                "The fit_inverse_transform parameter was not"
                " set to True when instantiating and hence "
                "the inverse transform is not available."
            )

        K = self._get_kernel(X, self.X_transformed_fit_)
        return np.dot(K, self.dual_coef_)