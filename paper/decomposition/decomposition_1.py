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

model = Model()
def C_oppo(i):
    if i in AR:
        return AL
    elif i in AL:
        return AR
    else:
        return set()
var_list_X = [(i,p,m,k) for i in tasks for p in products for m in stations for k in sides]
var_list_U = [(i,j,h,q) for i in tasks for j in products for h in tasks for q in products]
X = model.binary_var_dict(var_list_X, name='X')
Z = model.binary_var_matrix(tasks,tasks, name='Z')
Y = model.binary_var_matrix(sequences, products, name='Y')
C = model.integer_var_cube(tasks, products,stations, name='CompletionTime')
U = model.integer_var_dict(var_list_U, name='U')
A = model.integer_var_matrix(products, stations, name='ArrivalTime')
D = model.integer_var_matrix(products, stations, name='DepartTime')
C_max = model.integer_var(lb=0, name='C_max')
W = model.integer_var(lb=0, name='W')
M = 99999

# ct1
for i in tasks:
    for j in products:
        model.add_constraint(model.sum(model.sum(
            X[i, j, m, k] for k in side_task[i-1]) for m in stations)  == 1, ctname='ct1_{}_{}'.format(i,j))

# 所有前置完成后才能进行下一个任务
for i in set(tasks).difference(P0):
    for r in P[i]:
        for j in products:
            model.add_constraint(model.sum(model.sum(g*X[r,j,g,k]for k in side_task[r-1]) for g in stations)
                                 <= model.sum(model.sum(m*X[i,j,m,k] for k in side_task[i-1]) for m in stations))


for m in stations:
    for k in sides:
        model.add_constraint(W >= model.sum(model.sum(time[j-1][i-1] * X[i,j,m,k] for j in products) for i in tasks))

obj = model.minimize(W)
sol = model.solve()
# print(model.print_information())
print(model.solve_details)
print(sol)