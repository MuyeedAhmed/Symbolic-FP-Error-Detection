    def _perturb_tensor(self, value, constraint):
        if isinstance(constraint, constraints._IntegerGreaterThan):
            return value + 1
        if isinstance(constraint, constraints._LessThan):
            return value - torch.rand(value.shape)
        if isinstance(
            constraint, (constraints._GreaterThan, constraints._GreaterThanEq)
        ):
            return value + torch.rand(value.shape)
        if isinstance(
            constraint,
            (constraints._PositiveDefinite, constraints._PositiveSemidefinite),
        ):
            return value + torch.eye(value.shape[-1])
        if value.dtype in [torch.float, torch.double]:
            transform = transform_to(constraint)
            delta = value.new(value.shape).normal_()
            return transform(transform.inv(value) + delta)
        if value.dtype == torch.long:
            result = value.clone()
            result[value == 0] = 1
            result[value == 1] = 0
            return result
        raise NotImplementedError