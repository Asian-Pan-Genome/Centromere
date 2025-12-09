def zeros(shape):
    retval = []
    for x in range(shape[0]):
        retval.append([])
        for y in range(shape[1]):
            retval[-1].append(0)
    return retval

match_award      = 10
mismatch_penalty = -5
gap_penalty      = -5 # both for opening and extanding


def match_score(alpha, beta):
    if alpha == beta:
        return match_award  # 奖励匹配
    else:
        return mismatch_penalty  # 惩罚错配

def finalize(align1, align2):
    align1 = align1[::-1]    # 反转序列1
    align2 = align2[::-1]    # 反转序列2
    
    identity = 0
    score = 0
    for i in range(len(align1)):
        if align1[i] == align2[i]:                
            identity += 1
            score += match_score(align1[i], align2[i])
        elif align1[i] != align2[i] and align1[i] != '-' and align2[i] != '-': 
            score += match_score(align1[i], align2[i])
        elif align1[i] == '-' or align2[i] == '-':          
            score += gap_penalty
    
    #identity = float(identity) / len(align1) * 100
    
    #print(f'Identity = {identity:.3f}%')
    #print(f'Score = {score}')
    #print("_".join(align1))
    #print('_'.join([' ' if a != b else a for a, b in zip(align1, align2)]))
    #print("_".join(align2))
    return identity

def edit_dist(list1, list2):
    assert len(list1) == len(list2)
    count = len(list1)

    temp1 = []
    temp2 = []
    pivot = 0
    while True:
        if list1[pivot] == '-' or list2[pivot] == '-':
            pivot += 1
        else:
            break
    for i in range(count):
        rot_i = (pivot + i) % count
        temp1.append(list1[rot_i])
        temp2.append(list2[rot_i])
    #print(temp1)
    #print(temp2)

    start = False
    rst = 0
    for i in range(count):
        t1 = temp1[i]
        t2 = temp2[i]
        if t1 != '-' and t2 != '-':
            if t1 != t2:
                rst += 1

            if start:
                start = False
                rst += 1
        else:
            if not start:
                start = True
    if start:
        rst += 1
    #print(rst)
    return rst

def needle(seq1, seq2):
    m, n = len(seq1), len(seq2)
    score = zeros((m+1, n+1))      # DP表格
    
    # 初始化边界
    for i in range(m + 1):
        score[i][0] = gap_penalty * i
    for j in range(n + 1):
        score[0][j] = gap_penalty * j
    
    # 填充DP表格
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1][j - 1] + match_score(seq1[i-1], seq2[j-1])
            delete = score[i - 1][j] + gap_penalty
            insert = score[i][j - 1] + gap_penalty
            score[i][j] = max(match, delete, insert)
    
    # 追溯并计算对齐
    align1, align2 = [], []
    i, j = m, n
    while i > 0 and j > 0:
        score_current = score[i][j]
        score_diagonal = score[i-1][j-1]
        score_up = score[i][j-1]
        score_left = score[i-1][j]
        
        if score_current == score_diagonal + match_score(seq1[i-1], seq2[j-1]):
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif score_current == score_left + gap_penalty:
            align1.append(seq1[i-1])
            align2.append('-')
            i -= 1
        elif score_current == score_up + gap_penalty:
            align1.append('-')
            align2.append(seq2[j-1])
            j -= 1
    
    while i > 0:
        align1.append(seq1[i-1])
        align2.append('-')
        i -= 1
    while j > 0:
        align1.append('-')
        align2.append(seq2[j-1])
        j -= 1
    
    #print("_".join(align1))
    #print("_".join(align2))
    identity = finalize(align1, align2)
    ed = edit_dist(align1, align2)
    return align1, align2, ed, identity

