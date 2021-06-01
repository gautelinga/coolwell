import dolfin as df

class GenericBC(df.SubDomain):
    def __init__(self, Lx, b, h):
        self.Lx = Lx
        self.b = b
        self.h = h
        df.SubDomain.__init__(self)


class PeriodicBC(GenericBC):
    def inside(self, x, on_boundary):
        return bool(x[0] < -self.Lx/2 + df.DOLFIN_EPS_LARGE and
                    x[0] > -self.Lx/2 - df.DOLFIN_EPS_LARGE and
                    on_boundary)

    def map(self, x, y):
        y[0] = x[0] - self.Lx
        y[1] = x[1]


class TopWall(GenericBC):
    def inside(self, x, on_boundary):
        return on_boundary and x[1] > self.h - df.DOLFIN_EPS_LARGE

class BtmWall(GenericBC):
    def inside(self, x, on_boundary):
        return on_boundary and x[1] < -self.b + df.DOLFIN_EPS_LARGE

class Wall(GenericBC):
    def inside(self, x, on_boundary):
        return on_boundary

class NotWall(GenericBC):
    def inside(self, x, on_boundary):
        return on_boundary and (x[0] < -self.Lx/2 + df.DOLFIN_EPS_LARGE or
                                x[0] > self.Lx/2 - df.DOLFIN_EPS_LARGE)
