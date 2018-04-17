import numpy as np
import copy
def two_compatible_or_not(first, last):
    if first >= last:
        return print('You have to guarantee j<k! Pls try again!')  ##可注释掉##
    True_False_test = []
    for x in Cj[last]:
        if x in Cj[first] and Start_time[last] >= End_time[first] and End_time[last] - Start_time[first] <= L:  # 有相同的可加工机器类别 & 满足L约束
            # print('Job %s and job %s are compatible in machine class %s ' % (first, last, x))  ##可注释掉##
            True_False_test.append(True)  # 可以相容的类别,True
    return any(True_False_test)  # 不相容的话返回False


def take_all_to_list(lst):  # make list中有list的list成为一个大list
    output_list = []
    for i in lst:
        if isinstance(i, list):
            output_list.extend(i)
        else:
            output_list.append(i)
    return output_list


def eva1_at_one_time(vector, index_list):  # 输入：（被修改的向量，修改位置索引向量）赋值被修改向量指定位置元素为1
    for ind in index_list:
        if isinstance(ind, list):
            for item in ind:
                i_tem = int(item)
                vector[i_tem] = 1
        else:
            i_ind = int(ind)
            vector[i_ind] = 1
    return vector


def to_int(lst):
    for i in range(len(lst)):
        lst[i] = int(lst[i])
    return lst

def to_flt(lst):
    for i in range(len(lst)):
        lst[i] = float(lst[i])
    return lst


def to_minus_1(a):
    for i in range(n):
        for j in range(n):
            a[i, j] = -1
    return a


# def find_ten_max(mat):  # 寻找矩阵中top10位置索引
#     mat_copy = copy.deepcopy(mat)
#     K = 0
#     ten_max_index = []
#     while 1:
#         for i in range(n):
#             for j in range(n):
#                 if mat_copy[i, j] == np.max(mat_copy):
#                     ten_max_index.append([i, j])
#                     mat_copy[i, j] = -100
#                     K += 1
#                     if K == 10:
#                         return ten_max_index
def find_ten_max(mat):  # 寻找矩阵中top10位置索引,且为正数
    mat_copy = copy.deepcopy(mat)
    K = 0
    ten_max_index = []
    while 1:
        for i in range(n):
            for j in range(n):
                if mat_copy[i, j] == np.max(mat_copy) and mat_copy[i,j] > 0:
                    ten_max_index.append([i, j])
                    mat_copy[i, j] = -100
                    K += 1
                    if K == 10:
                        return ten_max_index
                elif mat_copy[i, j] == np.max(mat_copy) and mat_copy[i,j] <= 0:  #最大值不是正数了，而且还并没有到10个
                    return ten_max_index



def backtrack_one_schedule(clas: int, head_tail_index: list, l_tensor: list) -> list:  # 输入一个类别如0，1，2；2个元素的list，如[2,11]
    j = head_tail_index[0]
    k = head_tail_index[1]

    L0 = int(l_tensor[clas][j, k])
    if L0 == j or L0 == -1:
        return head_tail_index
    # if L0==j:
    # return head_tail_index
    # elif L0==-1:
    # return head_tail_index
    else:
        mid_L = []
        while 1:
            mid_L.append(L0)
            L_tran = l_tensor[clas][j, L0]
            L_tran = int(L_tran)
            L0 = L_tran
            if L0 == -1:
                return head_tail_index
            elif L0 == j:
                mid_L.append(L0)
                mid_L.reverse()
                mid_L.append(k)
                return mid_L

def dedupe(items):         #去除list中重复元素，且保持顺序不变！
    seen = set()
    for item in items:
        if item not in seen:
            yield item
            seen.add(item)

def zero_one_vec_or_not(lliisstt):   #判断一个list是否全为0或1
    t_f_test = []
    for i in lliisstt:
        if i == 0 or i == 1:
            t_f_test.append(True)
        else :
            t_f_test.append(False)
    return all(t_f_test)

def del_none(lsst):
    for item in lsst:
        if None in item:
            item.remove(None)
    return lsst




########################################################################################################################
# 形成初始方案#
########################################################################################################################

n = 20  # jobs
c = 4  # classes
L = 100  # 连续工作时限
np.random.seed(99)
# w = to_int(np.random.choice(range(1,20),c,replace=False))
w = to_int( 20 * np.random.rand(c) ) # cost of different Classe(w0,w1,w2)
Start_time = 100 * np.random.rand(n)
Start_time = Start_time - Start_time.min()  # 最短开始时间设置为0
Start_time.sort()  # 从这里按开始时间给job排序，并重新定义index
Processing_time = 100 * np.random.rand(n)
End_time = Start_time + Processing_time

while 1:
    S_group_jobs = {}
    break_test = []
    for machine_class in range(c):
        rand_choice = list(np.random.choice(range(n), np.random.randint(1, n), replace=False))
        S_group_jobs[machine_class] = rand_choice
        break_test.extend(rand_choice)
    t_f_test = [False for x in list(range(n)) if x not in break_test]
    if t_f_test:  # [False]不是false也不是空集，会执行if语句
        pass
        # print(t_f_test)
        # print('pass')
    else :             # t_f_test成为[]是我们想要的结果，if不执行[],所以在else中执行
        break




for s in range(c):
    S_group_jobs[s].sort()

Cj = {}  # 可以处理第j个job的机器种类
for j in range(n):
    cj = []
    for cls in range(c):
        if j in S_group_jobs[cls]:
            cj.append(cls)
    Cj[j] = cj