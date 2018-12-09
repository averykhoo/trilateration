import random
import time
from numpy import matrix
import numpy as np

from copy import deepcopy


def measure_distance(a, b, noise_variance=0.0):
    return np.sqrt(np.sum(np.power(a - b, 2))) + random.normalvariate(0, noise_variance)
    # return np.sqrt(np.sum(np.power(a - b, 2))) + random.random() * noise_variance * 2 - noise_variance


def F(A_, r_, C):
    return sum(np.abs(np.sqrt(np.sum(np.power(C - Ai, 2))) - ri) for Ai, ri in zip(A_, r_))


def delta(A_, r_, C, alpha=1):
    d_ = [np.sqrt(np.sum(np.power(C - Ai, 2))) for Ai in A_]

    # jTf := J(C).T * f(C)
    jTf_x = sum((C[0, 0] - Ai[0, 0]) * (di - ri) / di for Ai, di, ri in zip(A_, d_, r_))
    jTf_y = sum((C[1, 0] - Ai[1, 0]) * (di - ri) / di for Ai, di, ri in zip(A_, d_, r_))
    jTf_z = sum((C[2, 0] - Ai[2, 0]) * (di - ri) / di for Ai, di, ri in zip(A_, d_, r_))

    JTf = matrix([jTf_x, jTf_y, jTf_z]).T

    # JTJ := J(C).T * J(C)
    Jxx = sum((C[0, 0] - Ai[0, 0]) ** 2 / (di ** 2) for Ai, di in zip(A_, d_))
    Jyy = sum((C[1, 0] - Ai[1, 0]) ** 2 / (di ** 2) for Ai, di in zip(A_, d_))
    Jzz = sum((C[2, 0] - Ai[2, 0]) ** 2 / (di ** 2) for Ai, di in zip(A_, d_))

    Jxy = sum((C[0, 0] - Ai[0, 0]) * (C[1, 0] - Ai[1, 0]) / (di ** 2) for Ai, di in zip(A_, d_))
    Jxz = sum((C[0, 0] - Ai[0, 0]) * (C[2, 0] - Ai[2, 0]) / (di ** 2) for Ai, di in zip(A_, d_))
    Jyz = sum((C[1, 0] - Ai[1, 0]) * (C[2, 0] - Ai[2, 0]) / (di ** 2) for Ai, di in zip(A_, d_))

    JTJ = matrix([[Jxx, Jxy, Jxz], [Jxy, Jyy, Jyz], [Jxz, Jyz, Jzz]])

    return C - JTJ.I * JTf * alpha


def LSQ(A_, r_):
    """
    Linear least SQuares

    :param A_: station locations
    :param r_: distances from each station
    :return: most likely tag position
    """
    # di is distance between anchor i and reference anchor 0
    d_ = [np.sqrt(np.sum(np.power(Ai - A_[0], 2))) for Ai in A_[1:]]

    # Set up approximate Linear Least Squares model
    # bi == (x2 - x1)(x - x1) + (y2 - y1)(y - y1) + (z2 - z1)(z - z1)
    # bi ~=  0.5*(r1^2 - ri^2 + di^2)
    b_ = [0.5 * (r_[0] ** 2 - ri ** 2 + di ** 2) for ri, di in zip(r_[1:], d_)]

    # given Ax = b
    # assuming (A.T * A) is non-singular     (has an inverse)
    # assuming (A.T * A) is well-conditioned (similar input --> similar output):
    #                      A  * x ==                     b
    #                A.T * A  * x ==               A.T * b
    # (A.T * A).I * (A.T * A) * x == (A.T * A).I * A.T * b
    #                           x == (A.T * A).I * A.T * b
    A = np.concatenate([Ai - A_[0] for Ai in A_[1:]], axis=1).T
    b = matrix(b_).T
    c = matrix(A.T * A).I * A.T * b

    # because c is in the reference frame of A0, ie: c == [x - A0(x), y - A0(y), z - A0(z)]
    # we need to add back A0's coordinates to c to find tag coordinates in real world
    return c + A_[0]


def NLSQ(guess, f_error, f_update, max_iterations=100):
    """
    Non-linear Least SQuares

    :param A_: station locations
    :param r_: distances from each station
    :param guess: initial guess
    :param max_iterations: if NLSQ doesn't converge, stop computing
    :return: most likely tag position
    """
    # Aim of NLSQ is to minimise Error function F = sum(ei ** 2),
    # ei = sqrt((x - xi) ** 2+(y - yi) ** 2+(z - zi) ** 2)  -  r(i), i = 1 to n
    # LSQ_e = sum(np.abs(np.sqrt(np.sum(np.power(T - Ai, 2))) - ri) for Ai, ri in zip(A_, r_))
    best_error = curr_error = f_error(guess)
    # print curr_error, guess.T
    # best_guess = guess.copy()
    best_guess = deepcopy(guess)

    for _ in range(max_iterations):
        prev_err = curr_error
        # print guess
        # print df(guess)
        guess = f_update(guess)
        # guess = [a - b for a, b in zip(guess, df(guess))]
        curr_error = f_error(guess)
        # print curr_error,guess.T

        if curr_error < best_error:
            best_error = curr_error
            # best_guess = guess.copy()
            best_guess = deepcopy(guess)

        if abs(prev_err - curr_error) < 0.01:
            break

    # return best_guess
    return guess  # for testing against LSQ


if __name__ == '__main__':
    # locations of stations
    # must be 3-dimensional
    A_ = [
        # matrix([-2000, -2000, -2000]).T,
        matrix([1000, 1000, 1000]).T,
        # matrix([1000, 1000, -1000]).T,
        # matrix([1000, -1000, 1000]).T,
        # matrix([1000, -1000, -1000]).T,
        # matrix([-1000, 1000, 1000]).T,
        matrix([-1000, 1000, -1000]).T,
        matrix([-1000, -1000, 1000]).T,
        matrix([-1000, -1000, -1000]).T,
    ]

    # set up test
    t = time.time()
    nlsq_win = 0
    lsq_win = 0
    draw = 0
    errs = []

    # begin test
    for test_iter in range(10000):
        # generate random answer
        answer = [random.randrange(-1000000, 1000000, 1) / 1000.0 for _ in range(3)]

        # measurement with noise
        r_ = [measure_distance(matrix(answer).T, Ai) for Ai in A_]

        # get LSQ guess
        C = LSQ(A_, r_)


        def f_error(C):
            return F(A_, r_, C)


        def f_update(C):
            return delta(A_, r_, C, alpha=0.5)


        # get NLSQ guess
        C2 = NLSQ(C.copy(), f_error, f_update)

        # for getting average error later
        errs.append(measure_distance(matrix(answer).T, C2, 0))

        # compare errors of LSQ and NLSQ
        if F(A_, r_, C2) < F(A_, r_, C):
            nlsq_win += 1
        elif F(A_, r_, C2) == F(A_, r_, C):
            draw += 1
        else:
            lsq_win += 1

        # maybe print something
        if random.random() < .01:
            print('lsq ', C.T)
            print('nlsq', C2.T)
            print('ans ', matrix(answer))
            print('lerr', measure_distance(matrix(answer).T, C, 0))
            print('nerr', measure_distance(matrix(answer).T, C2, 0))
            print('ldst', F(A_, r_, C))
            print('ndst', F(A_, r_, C2))
            print('')

    # finished test, print results
    errs = sorted(errs)
    print('nlsq_win', nlsq_win)
    print('draw', draw)
    print('lsq_win ', lsq_win)
    print('time/tag', (time.time() - t) / 1000)
    print('avg offset', sum(errs) / 1000)
    print('med offset', errs[len(errs) // 2])

if 0 and __name__ == '__main__':
    # import autograd.numpy as np
    from autograd import grad

    A_ = [
        # [-2000, -2000, -2000],
        [1000, 1000, 1000],
        # [1000, 1000,   -1000],
        # [1000, -1000,  -1000],
        [-1000, 1000, -1000],
        [-1000, -1000, 1000],
        [-1000, -1000, -1000],
    ]

    # set up test
    e1 = 0
    e2 = 0

    # begin test
    for test_iter in range(1000):
        # generate random answer
        answer = [random.randrange(-1000000, 1000000, 1) / 1000.0 for _ in range(3)]

        # measurement with noise
        r_ = [measure_distance(matrix(answer).T, matrix(Ai).T) for Ai in A_]

        print('ans  ', answer)

        C = LSQ(map(lambda x: np.matrix(x).T, A_), r_)
        c = [i[0, 0] for i in list(C)]


        def f_error(C):
            return F(map(lambda x: np.matrix(x).T, A_), r_, C)


        def f_update(C):
            return delta(map(lambda x: np.matrix(x).T, A_), r_, C, alpha=0.5)


        print('lsq  ', C.T)
        C2 = NLSQ(C.copy(), f_error, f_update)
        print('alvin', C2.T)
        e1 += measure_distance(C2, matrix(answer).T, 0)


        def F2(A_, r_, C):
            # print [sum((c - a) ** 2 for c, a in zip(C, Ai)) ** 0.5 - ri for Ai, ri in zip(A_, r_)]
            return sum(abs(sum((c - a) ** 2 for c, a in zip(C, Ai)) ** 0.5 - ri) for Ai, ri in zip(A_, r_))


        def f2(C):
            return F2(A_, r_, C)


        g2 = grad(f2)


        def u2(c):
            dx, dy, dz = g2(c)
            x, y, z = c
            e = 0.5

            return [x - dx * e, y - dy * e, z - dz * e]


        # print 'lsq  ',c
        c2 = NLSQ(c, f2, u2)
        print('autog', c2)
        e2 += sum((c - a) ** 2 for c, a in zip(c2, answer)) ** .5

    print(e1 / 1000)
    print(e2 / 1000)
