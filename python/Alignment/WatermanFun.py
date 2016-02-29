import numpy
import copy


def chr_state(i):
    if i > 26:
        co = chr((i - 27) + ord('A'))
    else:
        co = chr((i - 1) + ord('a'))
    return co


def matchfun_naive(i1, i2):
    if i1 == i2:
        return 2
    else:
        return -1


def WatermanAligner(seq1, seq2, matchfun, alpha, beta):
    # this function assumes len(seq1)<=len(seq2)
    # gap penalty: k*alpha+beta, match reward: matchfun(char1, char2)
    # initialize
    n1 = len(seq1)
    n2 = len(seq2)
    matrix = []
    for i in range(0, n1 + 1):
        matrix.append([])
        for j in range(0, n2 + 1):
            matrix[i].append([])
        for j in range(0, n2 + 1):
            # matching score, no deletion.
            matrix[i][j].append(0)
            # matching score, deletion at the end of seq1.
            matrix[i][j].append(0)
            # matching score, deletion at the end of seq2.
            matrix[i][j].append(0)
            # matching score, deletion at both end.
            matrix[i][j].append(0)
            # matching status for each matrix
            matrix[i][j].extend([0, 0, 0, 0])
    # update the matrix
    for i in range(0, n1):
        for j in range(0, n2):
            #update matching case, matching_status=i if best score comes from i'th matrix
            m = [0, 0, 0, 0, 0]
            m[1] = matrix[i][j][0] + matchfun(seq1[i], seq2[j])
            m[2] = matrix[i][j][1] + matchfun(seq1[i], seq2[j])
            m[3] = matrix[i][j][2] + matchfun(seq1[i], seq2[j])
            m[4] = matrix[i][j][3] + matchfun(seq1[i], seq2[j])
            argm = numpy.argmax(m)
            matrix[i + 1][j + 1][0] = m[argm]
            matrix[i + 1][j + 1][4] = argm
            #update deletion at seq1, matching_status=i if best score comes from i'th matrix
            m = [0, 0, 0]
            m[1] = matrix[i+1][j][0] - alpha - beta
            m[2] = matrix[i+1][j][1] - alpha
            argm = numpy.argmax(m)
            matrix[i + 1][j + 1][1] = m[argm]
            matrix[i + 1][j + 1][5] = argm
            #update deletion at seq2, matching_status=i if best score comes from i'th matrix
            m = [0, 0, 0, 0]
            m[1] = matrix[i][j+1][0] - alpha - beta
            m[3] = matrix[i][j+1][2] - alpha
            argm = numpy.argmax(m)
            matrix[i + 1][j + 1][2] = m[argm]
            matrix[i + 1][j + 1][6] = argm
            #update deletion at both, matching_status=4/5 if best score comes from matrix 4 and with another deletion at seq1/2
            m = [0, 0, 0, 0, 0, 0]
            m[2] = matrix[i][j+1][1] - alpha - beta
            m[3] = matrix[i+1][j][2] - alpha - beta
            m[4] = matrix[i+1][j][3] - alpha
            m[5] = matrix[i][j+1][3] - alpha
            argm = numpy.argmax(m)
            matrix[i + 1][j + 1][3] = m[argm]
            matrix[i + 1][j + 1][7] = argm
    return matrix


def WatermanAligner_Even(seq1, seq2, matchfun, alpha):
    # this function assumes len(seq1)<=len(seq2)
    # gap penalty: k*alpha, match reward: matchfun(char1, char2)
    # initialize
    n1 = len(seq1)
    n2 = len(seq2)
    matrix = []
    for i in range(0, n1 + 1):
        matrix.append([])
        for j in range(0, n2 + 1):
            matrix[i].append([])
        for j in range(0, n2 + 1):
            # matching score
            matrix[i][j].append(0)
            # matching status, 1: match, 2: deletion in seq1, 3: deletion in seq2
            matrix[i][j].append(0)
    # update the matrix
    for i in range(0, n1):
        for j in range(0, n2):
            m = [0, 0, 0, 0]
            m[0] = 0
            m[1] = matrix[i][j][0] + matchfun(seq1[i], seq2[j])
            m[2] = matrix[i + 1][j][0] - alpha
            m[3] = matrix[i][j + 1][0] - alpha
            argm = numpy.argmax(m)
            matrix[i + 1][j + 1][0] = m[argm]
            matrix[i + 1][j + 1][1] = argm
    return matrix


def WatermanAligner_MD(seq1, seq2, matchfun, alpha, beta, gamma):
    # gap penalty: k*alpha+beta, matching gap: k*gamma, match reward: matchfun(char1, char2)
    return


def TraceBack(seq1, seq2, matrix, topn):
    return


def TraceBack_Even(seq1, seq2, matrix, topn):
    # this function assumes len(seq1)<=len(seq2)
    # coordinated by seq2, two matches are not allowed to have strong overlap.
    # initialize
    ol = 1.0 / 2
    n1 = len(seq1)
    n2 = len(seq2)
    # find matching regions
    score = [0] * n2
    index = [0] * n2
    for i in range(1, n2 + 1):
        for j in range(1, n1 + 1):
            if score[i - 1] < matrix[j][i][0]:
                score[i - 1] = matrix[j][i][0]
                index[i - 1] = j - 1
    v = numpy.argsort(score)
    v = list(v)
    v.reverse()
    scoreout=copy.deepcopy(score)
    t = 0
    matches = []
    for k in range(0, n2):
        jthis = v[k]
        ithis = index[v[k]]
        if score[jthis] == 0:
            continue
        # trace back
        coordthis = 0
        matchthis = [[score[jthis], ithis, jthis], '', '']
        while True:
            if matrix[ithis + 1][jthis + 1][0] == 0:
                break
            if matrix[ithis + 1][jthis + 1][1] == 1:
                matchthis[1] += chr_state(seq1[ithis])
                matchthis[2] += chr_state(seq2[jthis])
                ithis -= 1
                jthis -= 1
                coordthis += 1
            elif matrix[ithis + 1][jthis + 1][1] == 2:
                matchthis[1] += '-'
                matchthis[2] += chr_state(seq2[jthis])
                jthis -= jthis
                coordthis += 1
            elif matrix[ithis + 1][jthis + 1][1] == 3:
                matchthis[1] += chr_state(seq1[ithis])
                matchthis[2] += '-'
                ithis -= 1
        matchthis[1] = matchthis[1][::-1]
        matchthis[2] = matchthis[2][::-1]
        matches.append(matchthis)
        t += 1
        if t == topn:
            break
        # set scores near this match to zero
        coorddel = int(ol * coordthis)
        for i in range(max(0, v[k] - coorddel), min(n2, v[k] + coorddel + 1)):
            score[i] = 0
    return scoreout, matches


def TraceBack_MD(seq1, seq2, matrix, topn):
    return


def PrintMatrix(matrix, filename):
    f = open(''.join([filename, '.matrix']), 'w')
    n1 = len(matrix) - 1
    n2 = len(matrix[0]) - 1
    l = len(matrix[0][0]) / 2
    for i in range(1, n1 + 1):
        for j in range(1, n2 + 1):
            s = 0
            for k in range(0, l):
                if matrix[i][j][k] > s:
                    s = matrix[i][j][k]
            f.write(str(s))
            f.write('\t')
        f.write('\n')
    f.close()
    return


def VisualMap(seq2, matrix, filename, cutoff=float('inf')):
    # this function visualizes matching matrix to filename+'.html'
    f = open(''.join([filename, '.html']), 'w')
    n1 = len(matrix) - 1
    n2 = len(seq2)
    l = len(matrix[0][0]) / 2
    f.write('<html>\n<body>\n')
    for i in range(1, n2 + 1):
        s = 0
        for j in range(1, n1 + 1):
            for k in range(0, l):
                if matrix[j][i][k] > s:
                    s = matrix[j][i][k]
        c = chr_state(seq2[i - 1])
        f.write('<abbr title="')
        f.write(str(s))
        if s > cutoff:
            f.write('" style="color:#0101DF">')
        else:
            f.write('">')
        f.write(c)
        f.write('</abbr>\n')
    f.write('</body>\n</html>\n')
    f.close()
    return
