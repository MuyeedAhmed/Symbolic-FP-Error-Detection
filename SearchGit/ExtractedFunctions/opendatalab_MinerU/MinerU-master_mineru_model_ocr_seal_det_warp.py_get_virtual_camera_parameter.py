    def get_virtual_camera_parameter(self):
        vcam_thz = 0
        vcam_thx1 = 180
        vcam_thy = 180
        vcam_thx2 = 0

        vcam_x = 0
        vcam_y = 0
        vcam_z = 100

        radian = np.pi / 180

        angle_z = radian * vcam_thz
        angle_x1 = radian * vcam_thx1
        angle_y = radian * vcam_thy
        angle_x2 = radian * vcam_thx2

        optic_x = vcam_x
        optic_y = vcam_y
        optic_z = vcam_z

        fu = 100
        fv = 100

        matT = np.zeros((4, 4))
        matT[0, 0] = cos(angle_z) * cos(angle_y) - sin(angle_z) * sin(angle_x1) * sin(
            angle_y
        )
        matT[0, 1] = cos(angle_z) * sin(angle_y) * sin(angle_x2) - sin(angle_z) * (
            cos(angle_x1) * cos(angle_x2) - sin(angle_x1) * cos(angle_y) * sin(angle_x2)
        )
        matT[0, 2] = cos(angle_z) * sin(angle_y) * cos(angle_x2) + sin(angle_z) * (
            cos(angle_x1) * sin(angle_x2) + sin(angle_x1) * cos(angle_y) * cos(angle_x2)
        )
        matT[0, 3] = optic_x
        matT[1, 0] = sin(angle_z) * cos(angle_y) + cos(angle_z) * sin(angle_x1) * sin(
            angle_y
        )
        matT[1, 1] = sin(angle_z) * sin(angle_y) * sin(angle_x2) + cos(angle_z) * (
            cos(angle_x1) * cos(angle_x2) - sin(angle_x1) * cos(angle_y) * sin(angle_x2)
        )
        matT[1, 2] = sin(angle_z) * sin(angle_y) * cos(angle_x2) - cos(angle_z) * (
            cos(angle_x1) * sin(angle_x2) + sin(angle_x1) * cos(angle_y) * cos(angle_x2)
        )
        matT[1, 3] = optic_y
        matT[2, 0] = -cos(angle_x1) * sin(angle_y)
        matT[2, 1] = cos(angle_x1) * cos(angle_y) * sin(angle_x2) + sin(angle_x1) * cos(
            angle_x2
        )
        matT[2, 2] = cos(angle_x1) * cos(angle_y) * cos(angle_x2) - sin(angle_x1) * sin(
            angle_x2
        )
        matT[2, 3] = optic_z
        matT[3, 0] = 0
        matT[3, 1] = 0
        matT[3, 2] = 0
        matT[3, 3] = 1

        matS = np.zeros((4, 4))
        matS[2, 3] = 0.5
        matS[3, 2] = 0.5

        self.ifu = 1 / fu
        self.ifv = 1 / fv

        self.matT = matT
        self.matS = matS
        self.K = np.dot(matT.T, matS)
        self.K = np.dot(self.K, matT)