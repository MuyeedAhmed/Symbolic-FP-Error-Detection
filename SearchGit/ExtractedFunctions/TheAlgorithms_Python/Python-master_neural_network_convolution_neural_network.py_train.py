    def train(
        self, patterns, datas_train, datas_teach, n_repeat, error_accuracy, draw_e=bool
    ):
        # model training
        print("----------------------Start Training-------------------------")
        print((" - - Shape: Train_Data  ", np.shape(datas_train)))
        print((" - - Shape: Teach_Data  ", np.shape(datas_teach)))
        rp = 0
        all_mse = []
        mse = 10000
        while rp < n_repeat and mse >= error_accuracy:
            error_count = 0
            print(f"-------------Learning Time {rp}--------------")
            for p in range(len(datas_train)):
                # print('------------Learning Image: %d--------------'%p)
                data_train = np.asmatrix(datas_train[p])
                data_teach = np.asarray(datas_teach[p])
                data_focus1, data_conved1 = self.convolute(
                    data_train,
                    self.conv1,
                    self.w_conv1,
                    self.thre_conv1,
                    conv_step=self.step_conv1,
                )
                data_pooled1 = self.pooling(data_conved1, self.size_pooling1)
                shape_featuremap1 = np.shape(data_conved1)
                """
                print('  -----original shape   ', np.shape(data_train))
                print('  ---- after convolution  ',np.shape(data_conv1))
                print('  -----after pooling  ',np.shape(data_pooled1))
               """
                data_bp_input = self._expand(data_pooled1)
                bp_out1 = data_bp_input

                bp_net_j = np.dot(bp_out1, self.vji.T) - self.thre_bp2
                bp_out2 = self.sig(bp_net_j)
                bp_net_k = np.dot(bp_out2, self.wkj.T) - self.thre_bp3
                bp_out3 = self.sig(bp_net_k)

                # --------------Model Leaning ------------------------
                # calculate error and gradient---------------
                pd_k_all = np.multiply(
                    (data_teach - bp_out3), np.multiply(bp_out3, (1 - bp_out3))
                )
                pd_j_all = np.multiply(
                    np.dot(pd_k_all, self.wkj), np.multiply(bp_out2, (1 - bp_out2))
                )
                pd_i_all = np.dot(pd_j_all, self.vji)

                pd_conv1_pooled = pd_i_all / (self.size_pooling1 * self.size_pooling1)
                pd_conv1_pooled = pd_conv1_pooled.T.getA().tolist()
                pd_conv1_all = self._calculate_gradient_from_pool(
                    data_conved1,
                    pd_conv1_pooled,
                    shape_featuremap1[0],
                    shape_featuremap1[1],
                    self.size_pooling1,
                )
                # weight and threshold learning process---------
                # convolution layer
                for k_conv in range(self.conv1[1]):
                    pd_conv_list = self._expand_mat(pd_conv1_all[k_conv])
                    delta_w = self.rate_weight * np.dot(pd_conv_list, data_focus1)

                    self.w_conv1[k_conv] = self.w_conv1[k_conv] + delta_w.reshape(
                        (self.conv1[0], self.conv1[0])
                    )

                    self.thre_conv1[k_conv] = (
                        self.thre_conv1[k_conv]
                        - np.sum(pd_conv1_all[k_conv]) * self.rate_thre
                    )
                # all connected layer
                self.wkj = self.wkj + pd_k_all.T * bp_out2 * self.rate_weight
                self.vji = self.vji + pd_j_all.T * bp_out1 * self.rate_weight
                self.thre_bp3 = self.thre_bp3 - pd_k_all * self.rate_thre
                self.thre_bp2 = self.thre_bp2 - pd_j_all * self.rate_thre
                # calculate the sum error of all single image
                errors = np.sum(abs(data_teach - bp_out3))
                error_count += errors
                # print('   ----Teach      ',data_teach)
                # print('   ----BP_output  ',bp_out3)
            rp = rp + 1
            mse = error_count / patterns
            all_mse.append(mse)

        def draw_error():
            yplot = [error_accuracy for i in range(int(n_repeat * 1.2))]
            plt.plot(all_mse, "+-")
            plt.plot(yplot, "r--")
            plt.xlabel("Learning Times")
            plt.ylabel("All_mse")
            plt.grid(True, alpha=0.5)
            plt.show()

        print("------------------Training Complete---------------------")
        print((" - - Training epoch: ", rp, f"     - - Mse: {mse:.6f}"))
        if draw_e:
            draw_error()
        return mse