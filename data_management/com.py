import math


class vector:
    def __init__(self, x: (int, float) = 0, y: (int, float) = 0, z: (int, float) = 0):
        self.x = x
        self.y = y
        self.z = z

    def __add__(self, other):
        if not isinstance(other, vector):
            raise TypeError()
        return vector(self.x + other.x, self.y + other.y, self.z + other.z)

    def __iadd__(self, other):
        if not isinstance(other, vector):
            raise TypeError()
        self.x += other.x
        self.y += other.y
        self.z += other.z
        return self

    def __sub__(self, other):
        if not isinstance(other, vector):
            raise TypeError()
        return vector(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return vector(other*self.x, other*self.y, other*self.z)
        elif isinstance(other, vector):
            return self.x*other.x + self.y*other.y + self.z*other.z
        else:
            raise TypeError()

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            return vector(self.x/other, self.y/other, self.z/other)
        else:
            raise TypeError()

    def __rmul__(self, other):
        return self*other

    def __str__(self):
        return f"{self.x}i " + (("- " + str(-self.y)) if self.y < 0 else ("+ " + str(self.y))) + \
                               ((" - " + str(-self.z)) if self.z < 0 else (" + " + str(self.z)))


def com_quantity(masses: list, vectors: list):
    if len(masses) != len(vectors):
        raise ValueError
    vecsum = vector()
    for mass, vec in zip(masses, vectors):
        vecsum += mass*vec
    return vecsum/sum(masses)


def main():
    masses = [1174.5, 1050.76, 822.192] # kg
    # from entry 299:
    positions = [vector(0.448051, -0.483545, 0.539052),
                 vector(-0.959978, -0.0767199, -0.172722),
                 vector(0.586812, 0.788788, -0.549293)]
    velocities = [vector(-9.9106e-05, -0.000156659, 6.05595e-05),
                  vector(9.24332e-05, 8.71225e-05, -5.43853e-05),
                  vector(2.34428e-05, 0.000112444, -1.70046e-05)]
    com_pos = com_quantity(masses, positions)
    com_vel = com_quantity(masses, velocities)
    print(f"COM position = {com_pos}\nCOM velocity = {com_vel}")
    m1r1 = (positions[0])*(masses[0])
    m2r2 = (positions[1])*(masses[1])
    m3r3 = (positions[2])*(masses[2])
    totmass = masses[0] + masses[1] + masses[2]
    summ = m1r1 + m2r2 + m3r3
    compos = summ/(totmass)
    print(f"m1r1: {m1r1}")
    print(f"m2r2: {m2r2}")
    print(f"m3r3: {m3r3}")
    print(f"totmass: {totmass}")
    print(f"sum: {summ}")
    print(f"compos: {compos}")


if __name__ == "__main__":
    main()
