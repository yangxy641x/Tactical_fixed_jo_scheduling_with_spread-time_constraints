from TFJ import *
import copy
#########################################################################################################################
                                                    #B & B#
#########################################################################################################################
upper_bound=1e5
lower_bound,optimal_value,A_mat=price(lmd, A_matrix, J_plus, my_obj)  #实际上A_mat在这以后都没什么用...

if zero_one_vec_or_not(optimal_value):
    print('************************************************************************************')
    print('Nice! The sol. of relaxation ILP is EXACTLY the sol of ILP!')
    print('************************************************************************************')
    print('\n')

else :
    # update J_plus according to the positive components of Xl   到底是怎么更新J_plus的。。？？(这里是直接用最优分数方案的紧后一个)
    for update_job in range(n):
        J_plus[update_job].clear()
        for rel_sch in real_schedule:
            if update_job in rel_sch and rel_sch.index(update_job) != (len(rel_sch) - 1):
                immediate_job = rel_sch[rel_sch.index(update_job) + 1]
                # for immediate_job in rel_sch[rel_sch.index(update_job)+1 : ]:   #最优分数schedule的紧后工作
                #     if two_compatible_or_not(update_job, immediate_job):
                J_plus[update_job].append(immediate_job)
            elif update_job in rel_sch and rel_sch.index(update_job) == (len(rel_sch) - 1):
                J_plus[update_job].append(None)
        J_plus[update_job] = list(dedupe(J_plus[update_job]))

    Phi = []  # 存放Xs是分数的schedule
    Phi_column = []  # 存放所有分数s的0-1列向量
    for x in nonzero_schedule_index:
        if 0 < optimal_value[x] < 1:
            Phi_column.append(A_mat[:, x])
    for s in Phi_column:
        sch = []
        for i in range(n):
            if s[i] > 0:
                sch.append(i)
        Phi.append(sch)

    Jf = []  # job集，所有在‘分数schedule’中出现过的jobs(无序)
    Jf_raw = take_all_to_list(Phi)
    Jf = list(set(Jf_raw))

    Sj = defaultdict(list)  # 字典，格式为{分数schedule中的job:[出现过的分数schedule]}
    for job in Jf:
        for fraction_schedule in Phi:
            if job in fraction_schedule:
                Sj[job].append(fraction_schedule)

    # Is=[]           #所有分数schedule所属的机器种类
    # for fra_sch in Phi:
    #     for i in range(len(real_schedule)):
    #         if real_schedule[i] == fra_sch:
    #             Is.append(schedule_machine_class[i])

    # if Phi==[]:
    # print('This is REAL ILP optimal! Thanks for ur work!')
    for phi in Phi:  # Phi中schedule首尾添加None
        phi.insert(0, None)
        phi.append(None)

    separated_job = []  # 寻找separated job，即在分数schedule中紧后工作不止一个的job
    for job in Jf:
        job_immediate_successor = []
        for schedule in Sj[job]:
            imm_suc = schedule[(schedule.index(job)) + 1]
            job_immediate_successor.append(imm_suc)
        if len(set(job_immediate_successor)) > 1:
            separated_job.append(job)


    min_int_sol = None
    layer = 0
    while 1:
        layer += 1
        separated_job.sort()
        for sep_job in separated_job:  # 如[5,6,7,11]
            for j_plus_index in range(len(J_plus[sep_job]) + 1):  # 某一个separated job 的 descendant nodes
                if j_plus_index == len(J_plus[sep_job]):  # an immediate successor not from J^{+}
                    updata_J_plus = del_none(copy.deepcopy(J_plus))
                    updata_J_plus[sep_job].clear()
                    for imm in Compatible_set[sep_job]:
                        if imm not in J_plus[sep_job]:
                            updata_J_plus[sep_job].append(imm)
                    llmd, AA_matrix, my_obj, update_lmd, update_A_matrix, opt_val, obj_val = algorithm_gr(updata_J_plus)
                    update_obj_val, update_optimal_value = price(update_lmd, update_A_matrix, updata_J_plus, my_obj)[
                                                           : 2]
                    if len(take_all_to_list(real_schedule)) > n:
                        continue
                    elif zero_one_vec_or_not(update_optimal_value) and update_obj_val < upper_bound:  # 整数解且小于ub，作为上界
                        upper_bound = update_obj_val
                        min_int_sol = [(sep_job, J_plus[sep_job][j_plus_index]), upper_bound]
                    elif zero_one_vec_or_not(update_optimal_value) and update_obj_val >= upper_bound:  # 整数解且大于ub，舍掉
                        continue
                    elif (zero_one_vec_or_not(
                            update_optimal_value) == False) and update_obj_val >= upper_bound:  # 分数解，且大于ub，舍掉
                        continue
                    elif (zero_one_vec_or_not(
                            update_optimal_value) == False) and update_obj_val < lower_bound:  # 分数解，且小于lb，作为下界
                        lower_bound = update_obj_val
                else:
                    updata_J_plus = copy.deepcopy(J_plus)
                    updata_J_plus[sep_job].clear()
                    updata_J_plus[sep_job] = [J_plus[sep_job][j_plus_index]]
                    updata_J_plus = del_none(updata_J_plus)
                    llmd, AA_matrix, my_obj, update_lmd, update_A_matrix, opt_val, obj_val = algorithm_gr(updata_J_plus)
                    update_obj_val, update_optimal_value = price(update_lmd, update_A_matrix, updata_J_plus, my_obj)[
                                                           : 2]
                    if len(take_all_to_list(real_schedule)) > n:
                        continue
                    elif zero_one_vec_or_not(update_optimal_value) and update_obj_val < upper_bound:  # 整数解且小于ub，作为上界
                        upper_bound = update_obj_val
                        min_int_sol = [(sep_job, J_plus[sep_job][j_plus_index]), upper_bound]
                    elif zero_one_vec_or_not(update_optimal_value) and update_obj_val >= upper_bound:  # 整数解且大于ub，舍掉
                        continue
                    elif (zero_one_vec_or_not(
                            update_optimal_value) == False) and update_obj_val >= upper_bound:  # 分数解，且大于ub，舍掉
                        continue
                    elif (zero_one_vec_or_not(
                            update_optimal_value) == False) and update_obj_val < lower_bound:  # 分数解，且小于lb，作为下界
                        lower_bound = update_obj_val
        if min_int_sol is None:
            pass  #一般一层分支足以
        else :
            break


    sep_job=min_int_sol[0][0]
    updata_J_plus = copy.deepcopy(J_plus)
    updata_J_plus[sep_job].clear()
    updata_J_plus[sep_job] = [min_int_sol[0][1]]
    updata_J_plus = del_none(updata_J_plus)
    llmd, AA_matrix, my_obj, update_lmd, update_A_matrix, opt_val, obj_val = algorithm_gr(updata_J_plus)
    update_obj_val, update_optimal_value = price(update_lmd, update_A_matrix, updata_J_plus, my_obj)[: 2]
    print('************************************************************************************')
    print('Going through %s branching, the sol. of ILP is : \n %s'%(layer,real_schedule))
    print('************************************************************************************')
    print('The corresponding machine class are ： \n', schedule_machine_class)
    print('\t')
    print('Considering Branching, we take time :', time.time()-tic)
    print('\t')






















