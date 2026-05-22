    def forward(self, p, cv, pv):
        pN = p.shape[1]
        vN = cv.shape[1]
        p_shortcut = p

        p1 = (
            self.p(p)
            .reshape([-1, pN, self.num_heads, self.dim // self.num_heads])
            .transpose([0, 2, 1, 3])
        )
        cv1 = (
            self.cv(cv)
            .reshape([-1, vN, self.num_heads, self.dim // self.num_heads])
            .transpose([0, 2, 1, 3])
        )
        pv1 = (
            self.pv(pv)
            .reshape([-1, vN, self.num_heads, self.dim // self.num_heads])
            .transpose([0, 2, 1, 3])
        )

        edge = F.softmax(p1.matmul(pv1.transpose((0, 1, 3, 2))), -1)  # B h N N
        p_c = (edge @ cv1).transpose((0, 2, 1, 3)).reshape((-1, pN, self.dim))

        x1 = self.norm1(p_shortcut + self.drop_path1(self.p_proj(p_c)))

        x = self.norm2(x1 + self.drop_path1(self.mlp(x1)))
        return x