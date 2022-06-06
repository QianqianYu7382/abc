import logging as log
import math
import itertools as it
import numpy as np
import scipy.special
import scipy.stats



def unit_vector(vec):
    """
    Returns unit vector
    """
    return vec / np.linalg.norm(vec)


def cos_sim(v1, v2):
    """
    Returns cosine of the angle between two vectors
    """

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.clip(np.tensordot(v1_u, v2_u, axes=(-1, -1)), -1.0, 1.0)


    ##return np.dot(v1, v2) / math.sqrt(np.dot(v1, v1) * np.dot(v2, v2))


def weat_association(W, A, B):
    """
    Returns association of the word w in W with the attribute for WEAT score.
    s(w, A, B)
    :param W: target words' vector representations
    :param A: attribute words' vector representations
    :param B: attribute words' vector representations
    :return: (len(W), ) shaped numpy ndarray. each rows represent association of the word w in W
    """
    return np.mean(cos_sim(W, A), axis=-1) - np.mean(cos_sim(W, B), axis=-1)


def weat_differential_association(X, Y, A, B):
    """
    Returns differential association of two sets of target words with the attribute for WEAT score.
    s(X, Y, A, B)
    :param X: target words' vector representations
    :param Y: target words' vector representations
    :param A: attribute words' vector representations
    :param B: attribute words' vector representations
    :return: differential association (float value)
    """
    return np.sum(weat_association(X, A, B)) - np.sum(weat_association(Y, A, B))


def weat_p_value(W, X, Y, A, B, parametric=True):
    """
    Returns one-sided p-value of the permutation test for WEAT score
    CAUTION: this function is not appropriately implemented, so it runs very slowly
    :param X: target words' vector representations
    :param Y: target words' vector representations
    :param A: attribute words' vector representations
    :param B: attribute words' vector representations
    :return: p-value (float value)
    """
    X = np.array(list(X), dtype=np.int)
    Y = np.array(list(Y), dtype=np.int)
    A = np.array(list(A), dtype=np.int)
    B = np.array(list(B), dtype=np.int)

    assert len(X) == len(Y)
    size = len(X)
    XY = np.concatenate((X, Y))

    if parametric:
        log.info('Using parametric test')
        s = weat_differential_association(X, Y, A, B)

        log.info('Drawing {} samples'.format(W))
        samples = []

        for _ in range(W):
            np.random.shuffle(XY)
            Xi = XY[:size]
            Yi = XY[size:]
            assert len(Xi) == len(Yi)
            si = weat_differential_association(Xi, Yi, A, B)
            samples.append(si)

        log.info('Inferring p-value based on normal distribution')
        (shapiro_test_stat, shapiro_p_val) = scipy.stats.shapiro(samples)
        log.info('Shapiro-Wilk normality test statistic: {:.2g}, p-value: {:.2g}'.format(
            shapiro_test_stat, shapiro_p_val))
        sample_mean = np.mean(samples)
        sample_std = np.std(samples, ddof=1)
        log.info('Sample mean: {:.2g}, sample standard deviation: {:.2g}'.format(
            sample_mean, sample_std))
        p_val = scipy.stats.norm.sf(s, loc=sample_mean, scale=sample_std)
        return p_val

    else:
        log.info('Using permutation test')
        s = weat_differential_association(X, Y, A, B)
        total_true = 0
        total_equal = 0
        total = 0
        ##print("hi")

        num_partitions = int(scipy.special.binom(2 * len(X), len(X)))
        if num_partitions > W:
            total_true += 1
            total += 1
            log.info('Drawing {} samples (and biasing by 1)'.format(W - total))
            for _ in range(W - 1):
                np.random.shuffle(XY)
                Xi = XY[:size]
                Yi = XY[size:]
                assert len(Xi) == len(Yi)
                si = weat_differential_association(Xi, Yi, A, B)
                if si > s:
                    total_true += 1
                elif si == s:  # use conservative test
                    total_true += 1
                    total_equal += 1
                total += 1
        """"
        else:
            print("xiaoyu")
            log.info('Using exact test ({} partitions)'.format(num_partitions))
            for Xi in it.combinations(XY, len(X)):
                Xi = np.array(Xi, dtype=np.int)
                assert 2 * len(Xi) == len(XY)
                si = weat_association(Xi, A, B)
                if si.any() > s.any():
                    total_true += 1
                elif si.any() == s.any():  # use conservative test
                    total_true += 1
                    total_equal += 1
                total += 1
        """
        if total_equal:
            log.warning('Equalities contributed {}/{} to p-value'.format(total_equal, total))
        print("total_true, total",total_true, total)
        return total_true / total


def weat_score(X, Y, A, B):
    """
    Returns WEAT score
    X, Y, A, B must be (len(words), dim) shaped numpy ndarray
    CAUTION: this function assumes that there's no intersection word between X and Y
    :param X: target words' vector representations
    :param Y: target words' vector representations
    :param A: attribute words' vector representations
    :param B: attribute words' vector representations
    :return: WEAT score
    """

    x_association = weat_association(X, A, B)
    y_association = weat_association(Y, A, B)


    tmp1 = np.mean(x_association, axis=-1) - np.mean(y_association, axis=-1)
    tmp2 = np.std(np.concatenate((x_association, y_association), axis=0))

    return tmp1 / tmp2


def wefat_p_value(W, A, B, parametric = True, num=1000):
    """
    Returns WEFAT p-value
    W, A, B must be (len(words), dim) shaped numpy ndarray
    CAUTION: not implemented yet
    :param W: target words' vector representations
    :param A: attribute words' vector representations
    :param B: attribute words' vector representations
    :return: WEFAT p-value
    """
    W = np.array(list(W), dtype=np.int)
    A = np.array(list(A), dtype=np.int)
    B = np.array(list(B), dtype=np.int)

    size = len(W)
    if parametric:
        log.info('Using parametric test')
        s = wefat_score(W, A, B)

        log.info('Drawing {} samples'.format(W))
        samples = []

        for _ in range(num):
            np.random.shuffle(W)
            Wi = W[:size]
            si = wefat_score(Wi, A, B)
            samples.append(si)

        log.info('Inferring p-value based on normal distribution')
        (shapiro_test_stat, shapiro_p_val) = scipy.stats.shapiro(samples)
        log.info('Shapiro-Wilk normality test statistic: {:.2g}, p-value: {:.2g}'.format(
            shapiro_test_stat, shapiro_p_val))
        sample_mean = np.mean(samples)
        sample_std = np.std(samples, ddof=1)
        log.info('Sample mean: {:.2g}, sample standard deviation: {:.2g}'.format(
            sample_mean, sample_std))
        p_val = scipy.stats.norm.sf(s, loc=sample_mean, scale=sample_std)
        return p_val

    else:
        log.info('Using permutation test')
        s = wefat_score(W, A, B)
        total_true = 0
        total_equal = 0
        total = 0

        #num_partitions = int(scipy.special.binom(2 * len(W)))
        #if num_partitions > num:
        total_true += 1
        total += 1
        log.info('Drawing {} samples (and biasing by 1)'.format(W - total))
        for _ in range(num - 1):
            np.random.shuffle(W)
            Wi = W[:size]
            si = wefat_score(Wi, A, B)
            if si > s:
                total_true += 1
            elif si == s:  # use conservative test
                total_true += 1
                total_equal += 1
            total += 1
        #print("total_true / total", total_true , total)
        return total_true / total


def wefat_score(W, A, B):
    """
    Returns WEFAT score
    W, A, B must be (len(words), dim) shaped numpy ndarray
    CAUTION: this function assumes that there's no intersection word between A and B
    :param W: target words' vector representations
    :param A: attribute words' vector representations
    :param B: attribute words' vector representations
    :return: WEFAT score
    """
    tmp1 = weat_association(W, A, B)

    tmp2 = np.std(np.concatenate((cos_sim(W, A), cos_sim(W, B)), axis=0))
    #print("tmp1/tmp2",tmp1,tmp2)

    return np.mean(tmp1 / tmp2)