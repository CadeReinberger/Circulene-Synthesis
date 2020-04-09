import math as m

#utilities that come up in the course of the project.
class util:
    def linspace(lower, upper, num):
        if num <=1:
            return None
        dx = (upper - lower)/(num - 1)
        return [lower + i * dx for i in range(0, num)]

    def circular_linspace(n, is_degrees = False):
        if is_degrees:
            return list(map(util.rad2deg, util.circular_linspace(n)))
        return util.linspace(0, 2 * m.pi, n+1)[:-1]

    def deg2rad(deg):
        return deg * m.pi / 180

    def rad2deg(rad):
        return rad * 180 / m.pi

    #the number here is effectively a rounding constant universally for the TeX
    def round_off(num):
        return round(num, 6)

    #There's gotta be a better way to do this xd.
    def ang_add(a1, a2):
        return two_vec.from_angle(a1).rotate(a2).ang()

    k_epsilon = .000001

#run of the mill vector class. Due to a bug in my python setup, I don't have access to numpy.
class two_vec:
    def __init__(self, m_x, m_y):
        self.x = m_x
        self.y = m_y

    def id():
        return two_vec(0, 0)

    def rotate(self, angle_rad):
        c = m.cos(angle_rad)
        s = m.sin(angle_rad)
        x_prime = self.x * c - self.y * s
        y_prime = self.x * s + self.y * c
        return two_vec(x_prime, y_prime)

    #_lambda is the scale factor
    def scale(self, _lambda):
        return two_vec(self.x * _lambda, self.y * _lambda)

    def add(self, other):
        return two_vec(self.x + other.x, self.y + other.y)

    def subtract(self, other):
        return self.add(other.scale(-1))

    def dot(self, other):
        return self.x * other.x + self.y * other.y

    #ccw angle in radians with the positive x axis zero
    def from_angle(angle_rad):
        return two_vec(m.cos(angle_rad), m.sin(angle_rad))

    def mag(self):
        return m.sqrt(self.dot(self))\

    #ccw angle in radians from the positive x axis
    def ang(self):
        return m.atan2(self.y, self.x)

    #check first that these angles don't nearly add up to 0
    def ang_avg(a_one, a_two):
        form = lambda a : util.ang_add(a, 0)
        a_one = form(a_one)
        a_two = form(a_two)
        diff = util.ang_add(a_two, -a_one)
        return util.ang_add(a_one, -m.copysign(diff, a_two))

#Carbon because verything in Python here is just to generate skeletons. It's just quicker
#for my purposes to not add in hetero atom functionality and add in other atoms
#and functional groups manually ad-hoc from the structures. It's not tough to
#write, but it's even less to do by hand, especially since I adding in mechanism
#handling is definitely beyond the scope of this project, so I'm going to be
#thinking about and modifying these molecules manually anyway.
class carbon_atom: #for our purposes this is just a vector. However, this improves readability.
    def __init__(self, m_pos):
        self.pos = m_pos

class carbon_bond:
    single_str = '-'
    double_left_str = '=^'
    double_right_str = '=_'
    triple_str = '~'

    def __init__(self, m_bond_str, m_bond_len, m_bond_angle):
        self.bond_str = m_bond_str
        self.bond_len = m_bond_len
        self.bond_angle = m_bond_angle

    def from_carbons(first_carbon, second_carbon, bond_str):
        diff_vec = second_carbon.pos.subtract(first_carbon.pos)
        len = diff_vec.mag()
        ang = diff_vec.ang()
        return carbon_bond(bond_str, len, ang)

    def swap_str(self, new_str):
        return carbon_bond(new_str, self.bond_len, self.bond_angle)

    #like always. ang is ccw in radians from positive x axis
    def rotate(self, ang):
        return carbon_bond(self.bond_str, self.bond_len, util.ang_add(self.bond_angle, ang))

    #scales the bond length by _lambda
    def scale(self, _lambda):
        new_len = m.fabs(_lambda) * self.bond_len
        should_flip_ang = _lambda < 0
        new_ang = self.bond_angle if not should_flip_ang else util.ang_add(m.pi, self.bond_angle)
        return carbon_bond(self.bond_str, new_len, new_ang)

    def is_bond(self):
        return True

    def as_str(self):
        return self.bond_str + '[:' + str(util.round_off(util.rad2deg(self.bond_angle))) + ',' + str(util.round_off(self.bond_len)) + ']'

class carbon_branch:
    def __init__(self, m_start, m_bonds):
        self.start = m_start
        self.bonds = m_bonds #bonds can also include branches. I added incrementally and couldn't be bothered to refactor.

    #This is so that branches can be treated iteratively like bonds for string creation of bonds while branching
    def as_str(self):
        return '(' + self.to_text() + ')'

    def to_text(self):
        return ''.join(map(lambda bond : bond.as_str(), self.bonds))

    def is_bond(self):
        return False

    #geenrates a skeleton from the carbon atoms in the order given. Returns to start at the end if is_circ.
    def skel_from_cs(cs, is_circ = False):
        m_start = cs[0]
        num_bonds = len(cs) - 1 if not is_circ else len(cs)
        m_bonds = []
        for bond_num in range(0, num_bonds):
            m_bonds.append(carbon_bond.from_carbons(cs[bond_num], cs[(bond_num + 1) % len(cs)], carbon_bond.single_str))
        return carbon_branch(m_start, m_bonds)

    def add_dubs_left(self, dubs):
        skel = self
        for b in dubs:
            skel = skel.swap_bond(b, carbon_bond.double_left_str)
        return skel

    def insert_branch(self, pos, branch):
        before = self.bonds[:pos]
        during = [branch] #remember branches are simply stored like elements
        after = self.bonds[pos:]
        res_bonds = before + during + after
        return carbon_branch(self.start, res_bonds)

    #this will do it's best to insert this branch in the most natural direction
    def smart_insert_branch(self, pos, branch):
        before = self.bonds[:pos]
        after = self.bonds[pos:]
        assert len(before) > 0
        assert len(after) > 0
        pred = before[-1]
        succ = after[0]
        assert pred.is_bond()
        assert succ.is_bond()
        unangled_branch = branch.rotate(-branch.bonds[0].bond_angle)
        if m.fabs(util.ang_add(-pred.bond_angle, succ.bond_angle)) < util.k_epsilon:
            branch_ang = pred.bond_angle
        else:
            ang_one = util.ang_add(pred.bond_angle, m.pi)
            ang_two = succ.bond_angle
            avg_ang = util.ang_add(.5 * (ang_one + ang_two), 0)
            if m.fabs(util.ang_add(ang_one, -avg_ang)) < util.k_epsilon:
                branch_ang = util.ang_add(ang_one, .5 * m.pi)
            else:
                branch_ang = util.ang_add(avg_ang, 0)
        rotated_branch = unangled_branch.rotate(branch_ang)
        if m.fabs(pred.bond_len - succ.bond_len):
            rotated_branch = rotated_branch.scale(pred.bond_len / branch.bonds[0].bond_len)
        return self.insert_branch(pos, rotated_branch)

    def insert_branches(self, positions, branches, tries_smart = True):
        num_inserted = 0
        cur = self
        for ind in range(0, len(branches)):
            if tries_smart:
                try:
                    cur = cur.smart_insert_branch(positions[ind] + num_inserted, branches[ind])
                except:
                    cur = cur.insert_branch(positions[ind] + num_inserted, branches[ind])
            else:
                cur = cur.insert_branch(positions[ind] + num_inserted, branches[ind])
            num_inserted += 1
        return cur

    #bondnum is indexed from 0 for the first bond + 1 for every branch along the way. (But Not subbranches).
    def swap_bond(self, bond_num, new_bond):
        cur_bond = self.bonds[bond_num]
        new_bonds = self.bonds
        new_bonds[bond_num] = cur_bond.swap_str(new_bond)
        return carbon_branch(self.start, new_bonds)

    #ang in radians ccw from 0.
    def rotate(self, ang):
        return carbon_branch(self.start, list(map(lambda a : a.rotate(ang), self.bonds)))

    #multiplies all bond lengths by _lambda
    def scale(self, _lambda):
        return carbon_branch(self.start, list(map(lambda a : a.scale(_lambda), self.bonds)))

    #Note that despite the name, this is concatenation and therefore noncommutative
    def add(self, other):
        return carbon_branch(self.start, self.bonds + other.bonds)

class molecules:
    #unit ethane oriented horizontally left to right
    def ethane():
        return carbon_branch(carbon_atom(two_vec.id()), [carbon_bond(carbon_bond.single_str, 1.0, 0.0)])

    #mostly for testing, this is just a simple cylcoalkane
    def cycloalkane(c):
        carbons = list(map(lambda a : carbon_atom(two_vec.from_angle(a)), util.circular_linspace(c)))
        return carbon_branch.skel_from_cs(carbons, True).rotate(util.deg2rad(-90))

    def get_pentane():
        return molecules.cycloalkane(5).to_text()

    def get_hexane():
        return molecules.cycloalkane(6).to_text()

    def get_heptane():
        return molecules.cycloalkane(7).to_text()

    def benzene():
        return molecules.cycloalkane(6).add_dubs_left([0, 2, 4])

    def get_benzene():
        return molecules.benzene().to_text()


    #these aren't test or utils, these are the methods I'm actually using in my LaTeX source code for
    #my project for the synthesis of [7]Circulene
    class synth_skels:

        #use cases: 5,5'--dimethyl--2,2'--dinitrobiphenyl;
        def two_two_prime_five_five_prime_tetramethylbiphenyl():
            m_benzene = molecules.benzene().rotate(util.deg2rad(90))
            methyl = molecules.ethane()
            m_toluene = m_benzene.add(methyl)
            m_biphenyl = m_toluene.add(m_benzene.scale(-1))
            # true is to emphasize that this only works bc of smart insert
            two_two_prime_five_five_prime_tetramethylbiphenyl = m_biphenyl.insert_branches([1,2, 3, 4, 5, 8, 9, 10, 11, 12], [methyl] * 8, True)
            return two_two_prime_five_five_prime_tetramethylbiphenyl.to_text()

print(molecules.synth_skels.two_two_prime_five_five_prime_tetramethylbiphenyl())






