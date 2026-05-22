    def __matmul__(self, other: Variable) -> Variable:
        result = Variable(self.value @ other.value)

        with GradientTracker() as tracker:
            # if tracker is enabled, computation graph will be updated
            if tracker.enabled:
                tracker.append(OpType.MATMUL, params=[self, other], output=result)
        return result