  def test_experts_choose_masked_router_dtype_shape(self):
    tf_keras.mixed_precision.set_global_policy('mixed_bfloat16')
    num_groups = 2
    tokens_per_group = 3
    hidden_dim = tokens_per_group
    num_experts = tokens_per_group
    expert_capacity = 2
    x = np.zeros([num_groups, tokens_per_group, hidden_dim])
    x[0, 0, 0] += 1
    x[0, :2, :2] += 1
    x[1, 1:, 1:] += 1
    x[1, -1, -1] += 1

    router = moe.ExpertsChooseMaskedRouter(
        num_experts=num_experts,
        jitter_noise=0.1,
        use_bias=True,
        kernel_initializer=tf_keras.initializers.get('identity'),
        bias_initializer=tf_keras.initializers.get('ones'))
    router_mask = router(x, expert_capacity=expert_capacity, training=False)

    self.assertDTypeEqual(router_mask.dispatch_mask, tf.bfloat16)
    self.assertDTypeEqual(router_mask.combine_array, tf.bfloat16)

    expect_shape = [num_groups, tokens_per_group, num_experts, expert_capacity]
    self.assertEqual(expect_shape, router_mask.dispatch_mask.shape)
    self.assertEqual(expect_shape, router_mask.combine_array.shape)

    # top_k call may not be sorted, so can't compare the output directly
    # Check that the output contains only 0s and 1s
    out_dm = router_mask.dispatch_mask.numpy()
    self.assertSetEqual({0, 1}, set(out_dm.flatten().astype(np.int32)))
    # Check that the right tokens for selected
    out_dm_indices = np.dot(
        out_dm.transpose((0, 2, 3, 1)), np.arange(tokens_per_group))
    # Shape [num_groups, num_experts, expert_capacity]
    self.assertSetEqual({0, 1}, set(out_dm_indices[0, 0, :].astype(np.int32)))
    self.assertSetEqual({1, 2}, set(out_dm_indices[0, 1, :].astype(np.int32)))
    self.assertSetEqual({1, 2}, set(out_dm_indices[0, 2, :].astype(np.int32)))
    self.assertSetEqual({0, 1}, set(out_dm_indices[1, 0, :].astype(np.int32)))
    self.assertSetEqual({0, 1}, set(out_dm_indices[1, 1, :].astype(np.int32)))
    self.assertSetEqual({1, 2}, set(out_dm_indices[1, 2, :].astype(np.int32)))

    out_ca = router_mask.combine_array.numpy()
    out_ca = np.dot(out_ca, np.ones((expert_capacity,)))

    expected_combine_array = np.array([[[0.66, 0.0, 0.0], [0.42, 0.42, 0.16],
                                        [0.0, 0.33, 0.33]],
                                       [[0.33, 0.33, 0.0], [0.16, 0.42, 0.42],
                                        [0.0, 0.0, 0.66]]])
    self.assertAllClose(expected_combine_array, out_ca, atol=1e-2)