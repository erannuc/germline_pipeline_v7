
class Bits:
    def __init__(self, size):
        self.data = bytearray(int(size/8) + 1)
        self.size = size
        self.powers2a = tuple([2**i for i in range(0, 8)])
        self.revpowers2a = tuple([255 - 2**i for i in range(0, 8)])

    def on(self, index):
        # turn on element index
        bytearrayindex = int(index/8)
        byteindex = index % 8
        self.data[bytearrayindex] |= self.powers2a[byteindex]

    def off(self, index):
        # turn off element index
        bytearrayindex = int(index/8)
        byteindex = index % 8
        self.data[bytearrayindex] &= self.revpowers2a[byteindex]

    def get(self, index):
        # return value of a requested index
        bytearrayindex = int(index/8)
        byteindex = index % 8
        if self.data[bytearrayindex] & self.powers2a[byteindex]:
            return 1
        else:
            return 0

    def any(self):
        # return a boolean if there is at least one on value
        return any(self.data)

    def print(self, **kargs):
        start = 0
        end = self.size
        for k in kargs:
            if k == 'start':
                start = kargs[k]
            if k == 'end':
                end = kargs[k]
        print([self.get(i) for i in range(start, end)])

    def clear(self):
        self.data.clear()

if __name__ == "__main__":
    mybitsar = Bits(10)

    print(mybitsar.any())

    mybitsar.on(1)
    print(mybitsar.get(1))
    print(mybitsar.any())

    mybitsar.off(1)
    print(mybitsar.any())
    exit(0)

    print(mybitsar.get(0))
    print(mybitsar.get(2))

    mybitsar.off(1)
    print(mybitsar.get(1))
    print(mybitsar.get(0))
    print(mybitsar.get(2))

    mybitsar.on(9)
    print(mybitsar.get(9))
    print(mybitsar.get(8))

    mybitsar.off(9)
    print(mybitsar.get(9))
    print(mybitsar.get(8))

    mybitsar.on(1)
    mybitsar.on(9)

    print(mybitsar.get(1))
    print(mybitsar.get(9))

    mybitsar.print(end=3)
