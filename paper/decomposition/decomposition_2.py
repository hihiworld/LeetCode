# P0 = (1,2)
# AL = (1,3)
# AR = (2,6)
# AE = (4,5,7)
#
# side_task = [[0],[1],[0],[0,1],[0,1],[1],[0,1]]
# P = {1:set(),2:set(),3:[1],4:[2],5:[3,4],6:[4],7:[5,6]}
# S = {1:[3],2:[4],3:[5],4:[5,6],5:[7],6:[7],7:set()}
# Pa = {1:set(),2:set(),3:[1],4:[2],5:[1,2,3,4],6:[2,4],7:[1,2,3,4,5,6]}
# Sa = {1:[3,5,7],2:[4,5,6,7],3:[5,7],4:[5,6,7],5:[7],6:[7],7:set()}
# sides = [0,1]
# tasks = (1,2,3,4,5,6,7)
# products = (1,2,3)
# sequences = (1,2,3)
# stations = (1,2,3)
# time = [[1,3,1,3,2,3,1],
#         [2,1,3,2,1,2,3],
#         [2,1,1,2,3,1,2]]
P0 = (1,2,3)
AL = (1,4,6)
AR = (2,8,12)
AE = (3,5,7,9,10,11)

side_task = [[0],[1],[0,1],[0],[0,1],[0],[0,1],[1],[0,1],[0,1],[0,1],[1]]
P = {1:set(),2:set(),3:set(),4:[1],5:[2],6:[3],7:[4,5],8:[5],9:[5,6],10:[7,8],11:[9],12:[11]}
S = {1:[4],2:[5],3:[6],4:[7],5:[7,8,9],6:[9],7:[10],8:[10],9:[11],10:set(),11:[12],12:set()}
Pa = {1:set(),2:set(),3:set(),4:[1],5:[2],6:[3],7:[1,4,2,5],8:[2,5],9:[2,5,3,6],10:[1,2,4,5,7,8],11:[2,5,3,6,9],12:[2,5,3,6,9,11]}
Sa = {1:[4,7,10],2:[5,7,8,9,10,11,12],3:[6,9,11,12],4:[7,10],5:[7,8,9,10,11,12],6:[9,11,12],7:[10],8:[10],9:[11,12],10:set(),11:[12],12:set()}
tasks = (1,2,3,4,5,6,7,8,9,10,11,12)
products = (1,2,3)
stations = (1,2,3)
sides = (0,1)
sequences = (1,2,3)
time = [[2,3,2,3,1,1,3,3,3,2,2,1],
        [2,3,2,3,2,1,3,3,2,2,4,1],
        [2,6,2,3,1,1,3,3,2,2,2,1]]

from docplex.mp.model import Model
import json
model = Model()

def load_sol_json(filename):
    fp = json.load(open(filename, 'r'))['CPLEXSolution']['variables']
    P = {}
    for f in fp:
        if "X" in f['name']:
            P[f['name']] = int(1)
    return P

def C_oppo(i):
    if i in AR:
        return AL
    elif i in AL:
        return AR
    else:
        return set()

var_list_U = [(i,j,h,q) for i in tasks for j in products for h in tasks for q in products]

Y = model.binary_var_matrix(sequences, products, name='Y')
C = model.integer_var_cube(tasks, products,stations, name='CompletionTime')
U = model.integer_var_dict(var_list_U, name='U')
A = model.integer_var_matrix(products, stations, name='ArrivalTime')
D = model.integer_var_matrix(products, stations, name='DepartTime')
C_max = model.integer_var(lb=0, name='C_max')

M = 99999

filename = r'C:\Users\root\PycharmProjects\LeetCode\paper\decomposition\solution.json'
last_sol = load_sol_json(filename)
print(last_sol['X_1_1_1_0'])
list={}
for i in tasks:
    for j in products:
        for m in stations:
            for k in sides:
                if 'X_{}_{}_{}_{}'.format(i,j,m,k) in last_sol:
                    list[(i,j,m,k)] = int(1)

for i in tasks:
    for j in products:
        for m in stations:
            model.add_constraint(C[i, j, m] <= C_max)

for i in tasks:
    for j in products:
        model.add_constraint(model.sum(C[i,j,m] for m in stations) >= time[j-1][i-1])

# 没有安排在相同station上的C不存在
for i in  tasks:
    for j in products:
        for m in stations:
            b = list[(i,j,m,k) for k in side_task[i-1]]
            model.add_constraint(C[i,j,m] <= M * model.sum())
            print(model.add_constraint(C[i,j,m] <= M * model.sum_vars(last_sol['X_{}_{}_{}_{}'.format(i,j,m,k)] for k in side_task[i-1])))
# # 不同station上的有前后次序关系的task应满足
# for i in list(set(tasks).difference(set(P0))):
#     for r in P[i]:
#         for j in products:
#             for m in stations:
#                 for g in stations:
#                     model.add_constraint(C[i,j,m]-C[r,j,g] + M*(1-model.sum(last_sol['X_{}_{}_{}_{}'.format(i,j,m,k)] for k in side_task[i-1]))
#                                      +M*(1-model.sum(last_sol['X_{}_{}_{}_{}'.format(r, j, g, k)] for k in side_task[r-1]))
#                                      >= time[j-1][i-1])
#
# # 不同产品间task的关系，不能同时加工
# for j in products:
#     for q in products:
#         if j != q:
#             for i in tasks:
#                 for h in tasks:
#                     for m in stations:
#                         for k in list(set(side_task[i-1]).intersection(set(side_task[h-1]))):
#                             model.add_constraint(C[i,j,m]-C[h,q,m] +
#                                                  M * (1-last_sol['X_{}_{}_{}_{}'.format(i,j,m,k)])+
#                                                  M*(1-last_sol['X_{}_{}_{}_{}'.format(h,q,m,k)]) +
#                                                  M*U[i,j,h,q]
#                                          >= time[j-1][i-1])
# for j in products:
#     for q in products:
#         if j != q:
#             for i in tasks:
#                 for h in tasks:
#                     for m in stations:
#                         for k in list(set(side_task[i-1]).intersection(set(side_task[h-1]))):
#                             model.add_constraint(C[h,q,m]-C[i,j,m] + M * (1-last_sol['X_{}_{}_{}_{}'.format(i,j,m,k)])+
#                                                  M*(1-last_sol['X_{}_{}_{}_{}'.format(h,q,m,k)]) + M*(1-U[i,j,h,q])
#                                          >= time[j-1][h-1])
#
#
#
# # 同一产品内的task之间的约束，不能同时加工
# for i in list(set(tasks).difference(set(P0))):
#     for r in P[i]:
#         for j in products:
#             for m in stations:
#                 model.add_constraint(C[i,j,m]-C[r,j,m] + M*(1-model.sum(last_sol['X_{}_{}_{}_{}'.format(i,j,m,k)] for k in side_task[i-1]))
#                                      +M*(1-model.sum(last_sol['X_{}_{}_{}_{}'.format(r,j,m,k)] for k in side_task[r-1]))
#                                      >= time[j-1][i-1])
#
# for i in tasks:
#     r = list(set(tasks).difference(set(set(Pa[i]).union(set(Sa[i])).union(C_oppo(i)))))
#     for h in r:
#         if i < h:
#             for m in stations:
#                 for k in list(set(side_task[i-1]).intersection(set(side_task[h-1]))):
#                     for j in products:
#                         model.add_constraint(C[h,j,m]-C[i,j,m] + M*(1-last_sol['X_{}_{}_{}_{}'.format(i,j,m,k)])
#                                      +M*(1-last_sol['X_{}_{}_{}_{}'.format(h,j,m,k)])
#                                          +M*(1-U[i,j,h,j])
#                                      >= time[j-1][h-1])
#
# for i in tasks:
#     r = list(set(tasks).difference(set(set(Pa[i]).union(set(Sa[i])).union(C_oppo(i)))))
#     for h in r:
#         if i<h:
#             for m in stations:
#                 for k in list(set(side_task[i-1]).intersection(set(side_task[h-1]))):
#                     for j in products:
#                         model.add_constraint(C[i,j,m]-C[h,j,m] + M*(1-last_sol['X_{}_{}_{}_{}'.format(i,j,m,k)])
#                                      +M*(1-last_sol['X_{}_{}_{}_{}'.format(i,j,m,k)])
#                                         +M*U[i,j,h,j]
#                                      >= time[j-1][i-1])
#
# for s in sequences:
#     model.add_constraint(model.sum(Y[s,j] for j in products) == 1)
# for j in products:
#     model.add_constraint(model.sum(Y[s,j] for s in sequences) == 1)
#
# for s in sequences:
#     for j in products:
#         for q in products:
#             if s > 1 and q != j:
#                 for m in stations:
#                     model.add_constraint(D[j,m]-A[q,m] <= M*(2-Y[s-1,j]-Y[s,q]))
# for m in stations:
#     if m > 1:
#         for j in products:
#             model.add_constraint(A[j,m] >= D[j,m-1])
#
# for i in tasks:
#     for j in products:
#         for m in stations:
#             model.add_constraint(A[j,m]-C[i,j,m] <= time[j-1][i-1]*model.sum(last_sol['X_{}_{}_{}_{}'.format(i,j,m,k)] for k in side_task[i-1])
#                                  + M*(1-model.sum(last_sol['X_{}_{}_{}_{}'.format(i,j,m,k)] for k in side_task[i-1])))
#
# for j in products:
#     for m in stations:
#         for i in tasks:
#             model.add_constraint(D[j, m] >= C[i,j,m])
#
# for j in products:
#     for m in stations:
#         model.add_constraint(D[j, m] <= C_max)
#
# model.add_constraint(C_max>=14)
#
# model.minimize(C_max)
# sol = model.solve()
# print(sol)