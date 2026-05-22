    def setUp(self):
        super().setUp()
        positive_var = torch.randn(20, dtype=torch.double).exp()
        positive_var2 = torch.randn(20, dtype=torch.double).exp()
        random_var = torch.randn(20, dtype=torch.double)
        simplex_tensor = softmax(torch.randn(20, dtype=torch.double), dim=-1)
        cov_tensor = torch.randn(20, 20, dtype=torch.double)
        cov_tensor = cov_tensor @ cov_tensor.mT
        self.distribution_pairs = [
            (Bernoulli(simplex_tensor), scipy.stats.bernoulli(simplex_tensor)),
            (
                Beta(positive_var, positive_var2),
                scipy.stats.beta(positive_var, positive_var2),
            ),
            (
                Binomial(10, simplex_tensor),
                scipy.stats.binom(
                    10 * np.ones(simplex_tensor.shape), simplex_tensor.numpy()
                ),
            ),
            (
                Cauchy(random_var, positive_var),
                scipy.stats.cauchy(loc=random_var, scale=positive_var),
            ),
            (Dirichlet(positive_var), scipy.stats.dirichlet(positive_var)),
            (
                Exponential(positive_var),
                scipy.stats.expon(scale=positive_var.reciprocal()),
            ),
            (
                FisherSnedecor(
                    positive_var, 4 + positive_var2
                ),  # var for df2<=4 is undefined
                scipy.stats.f(positive_var, 4 + positive_var2),
            ),
            (
                Gamma(positive_var, positive_var2),
                scipy.stats.gamma(positive_var, scale=positive_var2.reciprocal()),
            ),
            (Geometric(simplex_tensor), scipy.stats.geom(simplex_tensor, loc=-1)),
            (
                Gumbel(random_var, positive_var2),
                scipy.stats.gumbel_r(random_var, positive_var2),
            ),
            (
                GeneralizedPareto(
                    loc=random_var, scale=positive_var, concentration=random_var / 10
                ),
                scipy.stats.genpareto(
                    c=random_var / 10, loc=random_var, scale=positive_var
                ),
            ),
            (HalfCauchy(positive_var), scipy.stats.halfcauchy(scale=positive_var)),
            (HalfNormal(positive_var2), scipy.stats.halfnorm(scale=positive_var2)),
            (
                InverseGamma(positive_var, positive_var2),
                scipy.stats.invgamma(positive_var, scale=positive_var2),
            ),
            (
                Laplace(random_var, positive_var2),
                scipy.stats.laplace(random_var, positive_var2),
            ),
            (
                # Tests fail 1e-5 threshold if scale > 3
                LogNormal(random_var, positive_var.clamp(max=3)),
                scipy.stats.lognorm(
                    s=positive_var.clamp(max=3), scale=random_var.exp()
                ),
            ),
            (
                LowRankMultivariateNormal(
                    random_var, torch.zeros(20, 1, dtype=torch.double), positive_var2
                ),
                scipy.stats.multivariate_normal(random_var, torch.diag(positive_var2)),
            ),
            (
                Multinomial(10, simplex_tensor),
                scipy.stats.multinomial(10, simplex_tensor),
            ),
            (
                MultivariateNormal(random_var, torch.diag(positive_var2)),
                scipy.stats.multivariate_normal(random_var, torch.diag(positive_var2)),
            ),
            (
                MultivariateNormal(random_var, cov_tensor),
                scipy.stats.multivariate_normal(random_var, cov_tensor),
            ),
            (
                Normal(random_var, positive_var2),
                scipy.stats.norm(random_var, positive_var2),
            ),
            (
                OneHotCategorical(simplex_tensor),
                scipy.stats.multinomial(1, simplex_tensor),
            ),
            (
                Pareto(positive_var, 2 + positive_var2),
                scipy.stats.pareto(2 + positive_var2, scale=positive_var),
            ),
            (Poisson(positive_var), scipy.stats.poisson(positive_var)),
            (
                StudentT(2 + positive_var, random_var, positive_var2),
                scipy.stats.t(2 + positive_var, random_var, positive_var2),
            ),
            (
                Uniform(random_var, random_var + positive_var),
                scipy.stats.uniform(random_var, positive_var),
            ),
            (
                VonMises(random_var, positive_var),
                scipy.stats.vonmises(positive_var, loc=random_var),
            ),
            (
                Weibull(
                    positive_var[0], positive_var2[0]
                ),  # scipy var for Weibull only supports scalars
                scipy.stats.weibull_min(c=positive_var2[0], scale=positive_var[0]),
            ),
            (
                # scipy var for Wishart only supports scalars
                # SciPy allowed ndim -1 < df < ndim for Wishar distribution after version 1.7.0
                Wishart(
                    (
                        20
                        if version.parse(scipy.__version__) < version.parse("1.7.0")
                        else 19
                    )
                    + positive_var[0],
                    cov_tensor,
                ),
                scipy.stats.wishart(
                    (
                        20
                        if version.parse(scipy.__version__) < version.parse("1.7.0")
                        else 19
                    )
                    + positive_var[0].item(),
                    cov_tensor,
                ),
            ),
        ]