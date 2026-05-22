    def testBinaryOps(self):
        testme = AllTests()
        # Binary operations

        callLst[:] = []
        testme + 1
        self.assertCallStack([("__add__", (testme, 1))])

        callLst[:] = []
        1 + testme
        self.assertCallStack([("__radd__", (testme, 1))])

        callLst[:] = []
        testme - 1
        self.assertCallStack([("__sub__", (testme, 1))])

        callLst[:] = []
        1 - testme
        self.assertCallStack([("__rsub__", (testme, 1))])

        callLst[:] = []
        testme * 1
        self.assertCallStack([("__mul__", (testme, 1))])

        callLst[:] = []
        1 * testme
        self.assertCallStack([("__rmul__", (testme, 1))])

        callLst[:] = []
        testme @ 1
        self.assertCallStack([("__matmul__", (testme, 1))])

        callLst[:] = []
        1 @ testme
        self.assertCallStack([("__rmatmul__", (testme, 1))])

        callLst[:] = []
        testme / 1
        self.assertCallStack([("__truediv__", (testme, 1))])


        callLst[:] = []
        1 / testme
        self.assertCallStack([("__rtruediv__", (testme, 1))])

        callLst[:] = []
        testme // 1
        self.assertCallStack([("__floordiv__", (testme, 1))])


        callLst[:] = []
        1 // testme
        self.assertCallStack([("__rfloordiv__", (testme, 1))])

        callLst[:] = []
        testme % 1
        self.assertCallStack([("__mod__", (testme, 1))])

        callLst[:] = []
        1 % testme
        self.assertCallStack([("__rmod__", (testme, 1))])


        callLst[:] = []
        divmod(testme,1)
        self.assertCallStack([("__divmod__", (testme, 1))])

        callLst[:] = []
        divmod(1, testme)
        self.assertCallStack([("__rdivmod__", (testme, 1))])

        callLst[:] = []
        testme ** 1
        self.assertCallStack([("__pow__", (testme, 1))])

        callLst[:] = []
        1 ** testme
        self.assertCallStack([("__rpow__", (testme, 1))])

        callLst[:] = []
        testme >> 1
        self.assertCallStack([("__rshift__", (testme, 1))])

        callLst[:] = []
        1 >> testme
        self.assertCallStack([("__rrshift__", (testme, 1))])

        callLst[:] = []
        testme << 1
        self.assertCallStack([("__lshift__", (testme, 1))])

        callLst[:] = []
        1 << testme
        self.assertCallStack([("__rlshift__", (testme, 1))])

        callLst[:] = []
        testme & 1
        self.assertCallStack([("__and__", (testme, 1))])

        callLst[:] = []
        1 & testme
        self.assertCallStack([("__rand__", (testme, 1))])

        callLst[:] = []
        testme | 1
        self.assertCallStack([("__or__", (testme, 1))])

        callLst[:] = []
        1 | testme
        self.assertCallStack([("__ror__", (testme, 1))])

        callLst[:] = []
        testme ^ 1
        self.assertCallStack([("__xor__", (testme, 1))])

        callLst[:] = []
        1 ^ testme
        self.assertCallStack([("__rxor__", (testme, 1))])