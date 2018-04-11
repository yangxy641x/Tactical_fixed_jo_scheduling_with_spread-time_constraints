import numpy as np
import cplex
import time
from TFJ_Functions import *
tic=time.time()

S_num_Mi={}          # 不同class的机器各用几条
for cl in range(c):
    S_num_Mi[cl] = 0
#S_num_Mi = {0: 0, 1: 0, 2: 0}
S_group_index={}     # 记录初始方案——类别：job编号
for cl in range(c):
    S_group_index[cl]=[]
#S_group_index = {0: [], 1: [], 2: []}

for j in range(n):  # 第j个job
    group_of_j = Cj[j]  # 可以处理Jobj的class列表

    for cls in group_of_j:

        S_all_test = []

        if S_num_Mi[cls] == 0:
            S_num_Mi[cls] += 1
            S_group_index[cls].append(j)
            break

        else:
            for job in S_group_index[cls]:
                if isinstance(job, list):
                    # lenth=len(job)
                    # te=[]
                    jobb = job[-1]
                    test = Start_time[j] >= End_time[jobb] and End_time[j] - Start_time[jobb] <= L  # 多工作的话只检查之前的最近一个job
                    # te.append(test)
                    # for l in range(lenth):
                    # jobb=job[l]
                    # test=Start_time[j] >= End_time[jobb] and End_time[j]-Start_time[jobb] <= L
                    # te.append(test)
                    if test:
                        # if any(te):
                        S_all_test.append(True)
                    else:
                        S_all_test.append(False)
                else:
                    test = Start_time[j] >= End_time[job] and End_time[j] - Start_time[job] <= L
                    S_all_test.append(test)
            if any(S_all_test):
                same_line_job_index = S_all_test.index(True)
                same_line_job = S_group_index[cls][same_line_job_index]
                same_line = take_all_to_list([same_line_job] + [j])
                S_group_index[cls].append(same_line)
                del (S_group_index[cls][same_line_job_index])
                break
            else:
                S_num_Mi[cls] += 1
                S_group_index[cls].append(j)
                break

#########################################################################################################################
# 用CPLEX求初始方案#
#########################################################################################################################

A = {}  # 初始零矩阵(一般当n比较大时，S_num_Mi都会有正数的)
machine_num = []  # 初始不同class机器数量
for cls in range(c):
    Acls = np.zeros((n, S_num_Mi[cls]))
    A[cls] = Acls
    machine_num.append(S_num_Mi[cls])

for cls, mac_num in zip(range(c), machine_num):
    for col in range(mac_num):
        vec = A[cls][:, col]
        ind_lst = [S_group_index[cls][col]]
        new_vec = eva1_at_one_time(vec, ind_lst)
        A[cls][:, col] = new_vec

# 补齐剩余列(增加单位阵)                   #This is where I feel a bit of worried. 列相关？对应的n个系数w如何确定？#
appen_mat = np.eye(n)
A_matrix = np.hstack((A[0], A[1], A[2], appen_mat))

cpx = cplex.Cplex()

num_of_col = sum(machine_num) + n

my_obj_appen = [float(w[np.random.choice(Cj[x], 1)]) for x in range(n)]  # 因为每个job所属类别可能不止一个，所以这里采用随机选取可用类别的方法
initial_obj=[]
for objc in range(c):
    initial_obj += S_num_Mi[objc]*[float(w[objc])]
my_obj=initial_obj + my_obj_appen
#my_obj = [float(w[0])] * S_num_Mi[0] + [float(w[1])] * S_num_Mi[1] + [float(w[2])] * S_num_Mi[2] + my_obj_appen

my_rhs = [1] * n
my_sense = ['G'] * n
my_ub = [cplex.infinity] * num_of_col

rows = []
for i in range(n):
    rows.extend(num_of_col * [i])
# rows=np.array(rows)
# rows=np.ravel(rows)
cols = list(range(num_of_col)) * n
vals = np.ravel(A_matrix)

rows = list(rows)
cols = list(cols)
vals = list(vals)

rows = to_int(rows)
cols = to_int(cols)
vals = to_int(vals)

coe = zip(rows, cols, vals)
cpx.objective.set_sense(sense=cpx.objective.sense.minimize)
cpx.variables.add(obj=my_obj, ub=my_ub)
cpx.linear_constraints.add(rhs=my_rhs, senses=my_sense)
cpx.linear_constraints.set_coefficients(coe)
cpx.solve()
lmd = cpx.solution.get_dual_values()
optima_value = cpx.solution.get_values()
print('Opt ? ', cpx.solution.get_status_string())
print("Solution value  = ", cpx.solution.get_objective_value())
print('Primal variable = ', optima_value)
print('Dual variable = ', lmd)

#########################################################################################################################
                                                #每个工作紧后工作集合#
#########################################################################################################################
J_plus=[]               #每一个job后可以紧跟（immediate successor）的job
for j in range(n):
    J_plus_j=[]
    for k in range(j+1,n):
        if two_compatible_or_not(j,k):        #这里把两者(J_plus和compatible)混为等同
            J_plus_j.append(k)
    J_plus.append(J_plus_j)


#########################################################################################################################
# Pricing Algorithm__DP#
#########################################################################################################################
def pricing_algorithm(lmd, A_matrix):
    ######################每次循环用的都是一样的f_matrix l_tensor???除了lambda######################
    f_matrix = []  # 存放不同i的f(j,k)
    l_tensor = []  # 存放不同i的从j到k的倒数第二项（l）工作

    for cls in range(c):  # 本例中为0，1，2
        f_matrix_i = np.zeros((n, n))  # 第i类f(j,k)

        for x in S_group_jobs[cls]:  # 对所有的第cls类（机器能处理的）job：
            f_matrix_i[x, x] = lmd[x]

        initial_tensor_i = np.zeros((n, n))  # n*n  [j,k]的l
        tensor_i = to_minus_1(initial_tensor_i)  # 不存在的组合就初始化为-1

        for j in S_group_jobs[cls]:  #
            k_imm_succ = []  # j给定时，所有符合条件k的取值(**K必须是j的j_plus里面的吗？不一定？？**)
            for k in S_group_jobs[cls][S_group_jobs[cls].index(j) + 1:]:
                if k in J_plus[j]:
                    k_imm_succ.append(k)

            if len(k_imm_succ) != 0:  # 如果存在这样的k(此时j=j)
                for val_k in k_imm_succ:
                    l_value = []  # j k 给定时 所有符合条件的l
                    for l in S_group_jobs[cls][S_group_jobs[cls].index(j): S_group_jobs[cls].index(val_k)]:
                        # for l in range(j,val_k):
                        if val_k in J_plus[l]:
                            l_value.append(l)

                    if len(l_value) != 0:  # 如果存在这样的l(此时j=j k=val_k)
                        max_f_jl = -100

                        val_l_made_max = -1  # 记录使得fi(j,k)最大的那个l

                        for val_l in l_value:
                            f_jl = f_matrix_i[j, val_l] + lmd[val_k]
                            if f_jl >= max_f_jl:
                                max_f_jl = f_jl
                                val_l_made_max = val_l
                        f_matrix_i[j, val_k] = max_f_jl
                        tensor_i[j, val_k] = val_l_made_max
        f_matrix.append(f_matrix_i)
        l_tensor.append(tensor_i)

    ##################################################不安###############################################################
    fi_star = []  ##疑惑：加的这十列是不是同一类型的？？？##若是，前十有0怎么办？##若不是，还要在所有good_cls里找top10？？？
    for i in range(c):
        star = np.max(f_matrix[i])
        fi_star.append(star)

    good_cls = [x for x in range(c) if fi_star[x] > w[x]]  # 满足fi_star > wi 的类别，list
    #不满足就是达到最优条件，此时good_cls空集!!!

    Top_ten = {}  # 字典，形式为 good种类:[该种类矩阵中的top_10横纵坐标号]
    for cls in good_cls:
        top_ten_index_cls = find_ten_max(f_matrix[cls])
        Top_ten[cls] = top_ten_index_cls

    Top_ten_value = {}  # 所有Top_ten字典中对应的值,格式->>[类别，横坐标，纵坐标]：值

    for cls in good_cls:
        for i in range(10):
            ind = Top_ten[cls][i]
            val = f_matrix[cls][ind[0], ind[1]]
            Top_ten_value[cls, (ind[0], ind[1])] = val

    Top_ten_value_sorted = sorted(Top_ten_value.items(), key=lambda x: x[1], reverse=True)
    All_top_ten_value_sorted = Top_ten_value_sorted[:10]
    # All_top_ten_value_sorted
    # 所有good class中top10，格式为[((1, (5, 6)), 6.7227855863079178),
    # ((1, (6, 6)), 6.7227855863079178),
    # ((1, (8, 8)), 4.8807839924058367),
    # ((1, (2, 2)), 4.8807839924058367),
    # ((1, (2, 5)), 4.8807839924058367),
    # ((1, (2, 7)), 4.8807839924058367),
    # ((1, (3, 3)), 4.8807839924058367),
    # ((1, (5, 8)), 4.8807839924058367),
    # ((1, (5, 9)), 1.8420015939020811),
    # ((1, (7, 9)), 1.8420015939020811)]
    ##list(All_top_ten_value_sorted[0][0][1])   -->>[5,6]

    add_lst = []
    for i in range(len(All_top_ten_value_sorted)):
        clas = int(All_top_ten_value_sorted[i][0][0])
        adding_sch = backtrack_one_schedule(clas=clas, head_tail_index=list(All_top_ten_value_sorted[i][0][1]),
                                            l_tensor=l_tensor)
        add_lst.append(adding_sch)

    add_col = np.zeros((n, len(add_lst)))  # 加的10列column
    for i in range(len(add_lst)):
        eva1_at_one_time(vector=add_col[:, i], index_list=add_lst[i])
        # one_list=eva1_at_one_time(vector=add_col[:, i], index_list=add_lst[i])
        # add_col[:,i]=one_list

    AA_matrix = np.hstack((A_matrix, add_col))    #达到最优条件时，add_col.shape为(n,0),可以进行hstack运算且不改变A_matrix!!!
    A_matrix = AA_matrix
    # A_matrix.shape

    ##########################################################################################################
    # 这里按照同一种类取前十做
    # if len(good_cls)>1:
    # cmp=0
    # cls_star=None
    # for good in good_cls:
    # if fi_star[good]>cmp:
    # cmp=fi_star[good]
    # cls_star=good
    # elif len(good_cls)==0:
    # print('Achieved Optimal')

    # else :
    #  cls_star=good_cls[0]
    ##########################################################################################################
    return A_matrix, All_top_ten_value_sorted, fi_star


#########################################################################################################################
# Pricing Algorithm__Add cols#
#########################################################################################################################
def added_new_col_cplx(A_matrix, All_top_ten_value_sorted):
    cp=cplex.Cplex()
    num_of_cols = A_matrix.shape[1]

    obj_appen = []
    for i in range(10):
        clas = All_top_ten_value_sorted[i][0][0]
        obj_appen.append(w[clas])
    global my_obj
    my_obj = my_obj + obj_appen

    rhs = [1] * n
    sense = ['G'] * n
    ub = [cplex.infinity] * num_of_cols
    row = []
    for i in range(n):
        row.append(num_of_cols * [i])
    row = np.array(row)
    row = np.ravel(row)
    col = list(range(num_of_cols)) * n
    val = np.ravel(A_matrix)

    row = list(row)
    col = list(col)
    val = list(val)

    row = to_int(row)
    col = to_int(col)
    val = to_int(val)

    coes = zip(row, col, val)
    cp.objective.set_sense(sense=cp.objective.sense.minimize)
    cp.variables.add(obj=my_obj, ub=ub)
    cp.linear_constraints.add(rhs=rhs, senses=sense)
    cp.linear_constraints.set_coefficients(coes)
    cp.solve()
    optimal_value = cp.solution.get_values()
    lmd = cp.solution.get_dual_values()
    primal_x=optimal_value
    print("Solution value  = ", cp.solution.get_objective_value())
    print('Primal variable = ', optimal_value)
    print('Opt ? ', cp.solution.get_status_string())
    print('Num of cols :', num_of_cols)
    # print('Dual variable = ',lmd_reindex)
    return lmd, A_matrix,primal_x

while 1:
    A_matrix,All_top_ten_value_sorted,fi_star=pricing_algorithm(lmd,A_matrix)

    if all(fi_star<=w):
        nonzero_schedule_index=[]                #记录非零元素索引位置(A_matrix的列索引)
        for var_ind in range(len(primal_x)):
            if primal_x[var_ind]>0:
                nonzero_schedule_index.append(var_ind)

        schedule_machine_class=[]
        for pos_ind in nonzero_schedule_index:
            for clas in range(c):
                if my_obj[pos_ind]==w[clas]:
                    schedule_machine_class.append(clas)


        schedule=[]
        for nonzero_schedule in nonzero_schedule_index:
            schedule.append( A_matrix[: , nonzero_schedule] )               #取列

        real_schedule=[]                      #记录最优方案
        for s in schedule:
            schedule_s = []
            for i in range(n):
                if s[i] > 0:
                    schedule_s.append(i)
            real_schedule.append(schedule_s)

        print('\n')
        print('************************************************************************************')
        print('Congras! Relaxation IP achieved optimal...')
        print('************************************************************************************')
        print('\n')
        print('The currently optimal schedule portfolio is ：\n ', real_schedule)  #怎么记录类别……？
        print('The corresponding machine class are ： \n', schedule_machine_class)
        print('\n')
        toc=time.time()
        period=toc-tic
        print('%d jobs, %d machine classes, take time %s s'%(n,c,period))
        print('\n')
        break
    lmd,A_matrix,primal_x=added_new_col_cplx(A_matrix,All_top_ten_value_sorted)



