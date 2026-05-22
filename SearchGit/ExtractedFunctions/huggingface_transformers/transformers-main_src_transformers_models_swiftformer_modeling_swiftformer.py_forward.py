    def forward(self, x):
        query = self.to_query(x)
        key = self.to_key(x)

        query = torch.nn.functional.normalize(query, dim=-1)
        key = torch.nn.functional.normalize(key, dim=-1)

        query_weight = query @ self.w_g
        scaled_query_weight = query_weight * self.scale_factor
        scaled_query_weight = scaled_query_weight.softmax(dim=-1)

        global_queries = torch.sum(scaled_query_weight * query, dim=1)
        global_queries = global_queries.unsqueeze(1).repeat(1, key.shape[1], 1)

        out = self.proj(global_queries * key) + query
        out = self.final(out)

        return out