    def _split_psi_buckets(self, nodes: list[_PsiTreeNode]) -> list[list[_PsiTreeNode]]:
        """Split large Psi inputs so exact pair ranking is bounded per bucket."""
        if len(nodes) <= self._psi_bucket_size:
            return [nodes]

        node_embeddings = self._normalize_embeddings(np.asarray([node.embedding for node in nodes], dtype=np.float64))
        groups = [np.arange(len(nodes), dtype=int)]
        buckets = []

        while groups:
            group = np.asarray(groups.pop(), dtype=int)
            if len(group) <= self._psi_bucket_size:
                buckets.append(group.tolist())
                continue

            fanout = min(max(2, int(np.ceil(len(group) / self._psi_bucket_size))), len(group), 32)
            group_embeddings = node_embeddings[group]
            center_idx = np.linspace(0, len(group_embeddings) - 1, num=fanout, dtype=int)
            centers = group_embeddings[center_idx].copy()

            for _ in range(5):
                labels = np.argmax(group_embeddings @ centers.T, axis=1)
                for center_id in range(fanout):
                    mask = labels == center_id
                    if not np.any(mask):
                        continue
                    center = group_embeddings[mask].mean(axis=0)
                    norm = np.linalg.norm(center)
                    centers[center_id] = center / norm if norm > 0 else center

            labels = np.argmax(group_embeddings @ centers.T, axis=1)
            split_groups = [group[labels == center_id].tolist() for center_id in range(fanout)]
            split_groups = [bucket for bucket in split_groups if bucket]
            if len(split_groups) <= 1:
                split_groups = [
                    group[start:start + self._psi_bucket_size].tolist()
                    for start in range(0, len(group), self._psi_bucket_size)
                ]
            groups.extend(split_groups)

        buckets = [bucket for bucket in buckets if bucket]
        buckets.sort(key=lambda bucket: (len(bucket), bucket[0]))
        return [[nodes[idx] for idx in bucket] for bucket in buckets]