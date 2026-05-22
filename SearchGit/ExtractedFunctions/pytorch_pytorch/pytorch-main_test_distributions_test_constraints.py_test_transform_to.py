def test_transform_to(constraint_fn, args, is_cuda):
    constraint = build_constraint(constraint_fn, args, is_cuda=is_cuda)
    t = transform_to(constraint)
    if constraint_fn is constraints.corr_cholesky:
        # (D * (D-1)) / 2 (where D = 4) = 6 (size of last dim)
        x = torch.randn(6, 6, dtype=torch.double)
    else:
        x = torch.randn(5, 5, dtype=torch.double)
    if is_cuda:
        x = x.cuda()
    y = t(x)
    if not constraint.check(y).all():
        raise AssertionError(f"Failed to transform_to({constraint})")
    x2 = t.inv(y)
    y2 = t(x2)
    if not torch.allclose(y, y2):
        raise AssertionError(f"Error in transform_to({constraint}) pseudoinverse")