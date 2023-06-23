import numpy as np
import copy
import time
import random
from graphviz import Digraph
from IPython.display import Image

"""Define global variables"""
computed_table = dict()
unique_table = dict()
global_index_order = dict()
global_node_idx=0
add_find_time=0
add_hit_time=0
cont_find_time=0
cont_hit_time=0
epi=0.000001

root_of_unit = 8

rotate_angle = 2*np.pi/root_of_unit


def set_root_of_unit(k):
    global root_of_unit,rotate_angle
    root_of_unit = k
    rotate_angle = 2*np.pi/root_of_unit

class Index:
    """The index, here idx is used when there is a hyperedge"""
    def __init__(self,key,idx=0):
        self.key = key
        self.idx = idx
        
    def __eq__(self,other):
        if self.key == other.key and self.idx == other.idx:
            return True
        else:
            return False
        
    def __lt__(self,other):
        if global_index_order[self.key] < global_index_order[other.key]:
            return True
        elif self.key == other.key and self.idx<other.idx:
            return True
        
        return False
    
    def __str__(self):
        return str((self.key,self.idx))
    
    
class map_item:
    def __init__(self,x = None, rotate=0):
        self.x = x
        self.rotate = rotate%root_of_unit
        
    def __str__(self):
        s=''
        if self.x:
            s+='x'
        if self.rotate>0:
            s+=str(self.rotate)
        return s
    
    def self_copy(self):
        res = map_item(self.x,self.rotate)
        return res
    
class the_maps:    
    def __init__(self,data=dict()):
        self.data = data
        
    def self_copy(self):
        res=the_maps(dict())
        res.data={k:self.data[k].self_copy() for k in self.data}
        return res
        
    def add_item(self,item):
        self.data[item[0]] = map_item(item[1],item[2])
        
    def add(self,other):
        temp=other.self_copy()
        self.data.update(temp.data)
        
    def __mul__(self,other):
        """self*other;note that xt^kx=\omega^kt^{-k}"""
        if not self.data:
            return other.self_copy(),0
        if not other.data:
            return self.self_copy(),0
        
#         print(self.data,other.data)
        
        max_idx = max(max(self.data),max(other.data))
        
        the_phase = 0
        res_map = the_maps(dict())
        
        for k in range(max_idx,-1,-1):
            if not k in other.data:
                if k in self.data:
                    res_map.data[k] = self.data[k].self_copy()
            else:
                x_0 = other.data[k].x
                r_0 = other.data[k].rotate
                x_1 = None
                r_1 = 0
                
                if k in self.data:
                    x_1 = self.data[k].x
                    r_1 = self.data[k].rotate
                    
                if not x_0:
                    if not (x_1==None and (r_1+r_0)%root_of_unit==0):
                        res_map.data[k] = map_item(x_1,r_1+r_0)
                else:
                    the_phase += r_1
                    if x_1:
                        if not (r_0-r_1)%root_of_unit==0:
                            res_map.data[k] = map_item(None,r_0-r_1)
                    else:
                        res_map.data[k] = map_item('x',r_0-r_1)
                        
        return res_map,the_phase
        
    def __truediv__(self,other):
        """self/other"""
        
        if not other.data:
            return self.self_copy(),0
        else:
            max_idx = max(other.data)
        if self.data:
            max_idx = max(max(self.data),max_idx)
    
        the_phase = 0
        res_map = the_maps(dict())
    
        for k in range(max_idx,-1,-1):
            if not k in other.data:
                if k in self.data:
                    res_map.data[k] = self.data[k].self_copy()
            else:
                x_0 = other.data[k].x
                r_0 = other.data[k].rotate
                x_1 = None
                r_1 = 0
                
                if k in self.data:
                    x_1 = self.data[k].x
                    r_1 = self.data[k].rotate
                    
                if x_0==x_1:
                    if not (r_1-r_0)%root_of_unit==0:
                        res_map.data[k] = map_item(None,r_1-r_0)
                else:
                    res_map.data[k] = map_item('x',r_1+r_0)
                    the_phase += -r_0

        return res_map,the_phase
            
    def __str__(self):
        s = [(k,str(self.data[k])) for k in self.data]
        s.sort()
        s={k[0]:k[1] for k in s}
        return str(s)
        
        
    def self_hash(self):
        return str(self)
    
    def __eq__(self,other):
        return self.self_hash()==other.self_hash()
        
    
class Node:
    """To define the node of TDD"""
    def __init__(self,key,num=2):
        self.idx = 0
        self.key = key
        self.succ_num=num
        self.out_weight=[1]*num
        self.out_maps=[]
        self.successor=[None]*num
        self.meas_prob=[]
        self.isidentity = False #In the next step, we will add the special check anf operation for identity


class TDD:
    def __init__(self,node):
        """TDD"""
        self.weight = 1
        
        self.map = the_maps(dict())
        
        self.index_set=[]
        
        self.key_2_index=dict()
        self.index_2_key=dict()
        self.key_width=dict() #Only used when change TDD to numpy
        
        if isinstance(node,Node):
            self.node=node
        else:
            self.node=Node(node)
            
            
    def node_number(self):
        node_set=set()
        node_set=get_node_set(self.node,node_set)
        return len(node_set)
    
    def self_copy(self):
        temp = TDD(self.node)
        temp.weight = self.weight
        temp.map = self.map.self_copy()
        temp.index_set = copy.copy(self.index_set)
        temp.key_2_index=copy.copy(self.key_2_index)
        temp.index_2_key=copy.copy(self.index_2_key)
        return temp
    
    def show(self,real_label=True):
        edge=[]              
        dot=Digraph(name='reduced_tree')
        dot=layout(self.node,self.key_2_index,dot,edge,real_label)
        dot.node('-0','',shape='none')
        label1=str(complex(round(self.weight.real,2),round(self.weight.imag,2)))
        if self.map.data:
            label1+=str(self.map)
        dot.edge('-0',str(self.node.idx),color="blue",label=label1)
        dot.format = 'png'
        return Image(dot.render('output'))
    
    def to_array(self,var=[]):
        split_pos=0
        key_repeat_num=dict()
        var_idx=dict()       
        if var:
            for idx in var:
                if not idx.key in var_idx:
                    var_idx[idx.key]=1
                else:
                    var_idx[idx.key]+=1
        elif self.index_set:
            for idx in self.index_set:
                if not idx.key in var_idx:
                    var_idx[idx.key]=1
                else:
                    var_idx[idx.key]+=1
        if var:
            split_pos=len(var_idx)-1
        elif self.key_2_index:
            split_pos=max(self.key_2_index)
        else:
            split_pos=self.node.key
        orig_order=[]
        for k in range(split_pos+1):
            if k in self.key_2_index:
                if self.key_2_index[k] in var_idx:
                    key_repeat_num[k] = var_idx[self.key_2_index[k]]
            else:
                key_repeat_num[k]=1
            if k in self.key_2_index:
                for k1 in range(key_repeat_num[k]):
                    orig_order.append(self.key_2_index[k])
                     

        res = tdd_2_np(self,split_pos,key_repeat_num)
        
        return res

    def measure(self,split_pos=None):
        res=[]
        get_measure_prob(self)
        if split_pos==None:
            if self.key_2_index:
                split_pos=max(self.key_2_index)
            else:
                split_pos=self.node.key

        if split_pos==-1:
            return ''
        else:
            if split_pos!=self.node.key:
                l=random.randint(0,1)
                temp_res=self.measure(split_pos-1)
                res=str(l)+temp_res
                return res
            l=random.uniform(0,sum(self.node.meas_prob))
            if l<self.node.meas_prob[0]:
                temp_tdd=Slicing(self,self.node.key,0)
                temp_res=temp_tdd.measure(split_pos-1)
                res='0'+temp_res
            else:
                temp_tdd=Slicing(self,self.node.key,1)
                temp_res=temp_tdd.measure(split_pos-1)
                res='1'+temp_res
#         print(res)
        return res
    def get_amplitude(self,b):
        """b is the term for calculating the amplitude"""
        if len(b)==0:
            return self.weight
        
        if len(b)!=self.node.key+1:
            b.pop(0)
            return self.get_amplitude(b)
        else:
            temp_tdd=Slicing(self,self.node.key,b[0])
            b.pop(0)
            res=temp_tdd.get_amplitude(b)
            return res*self.weight
            
    def sampling(self,k):
        res=[]
        for k1 in range(k):
            temp_res=self.measure()
            res.append(temp_res )
#         print(res)
        return res
        
        
    def __eq__(self,other):
        if self.node==other.node and get_int_key(self.weight)==get_int_key(other.weight) and self.map==other.map:# and self.key_2_index==other.key_2_index
            return True
        else:
            return False
        
def layout(node,key_2_idx,dot=Digraph(),succ=[],real_label=True):
    col=['red','blue','black','green']
    if real_label and node.key in key_2_idx:
        if node.key==-1:
            dot.node(str(node.idx), str(1), fontname="helvetica",shape="circle",color="red")
        else:
            dot.node(str(node.idx), key_2_idx[node.key], fontname="helvetica",shape="circle",color="red")
    else:
        dot.node(str(node.idx), str(node.key), fontname="helvetica",shape="circle",color="red")
    for k in range(node.succ_num):
        if node.successor[k]:
            label1=str(complex(round(node.out_weight[k].real,2),round(node.out_weight[k].imag,2)))
            if node.out_maps[k].data:
                label1+=str(node.out_maps[k])
            if not node.successor[k] in succ:
                dot=layout(node.successor[k],key_2_idx,dot,succ,real_label)
                dot.edge(str(node.idx),str(node.successor[k].idx),color=col[k%4],label=label1)
                succ.append(node.successor[k])
            else:
                dot.edge(str(node.idx),str(node.successor[k].idx),color=col[k%4],label=label1)
    return dot        

        
def Ini_TDD(index_order=[]):
    """To initialize the unique_table,computed_table and set up a global index order"""
    global computed_table
    global unique_table
    global global_node_idx
    global add_find_time,add_hit_time,cont_find_time,cont_hit_time
    global root_of_unit,rotate_angle
    global_node_idx=0
    unique_table = dict()
    computed_table = dict()
    add_find_time=0
    add_hit_time=0
    cont_find_time=0
    cont_hit_time=0
    set_index_order(index_order)
    root_of_unit = 8
    rotate_angle = 2*np.pi/root_of_unit
    return get_identity_tdd()

def Clear_TDD():
    """To initialize the unique_table,computed_table and set up a global index order"""
    global computed_table
    global unique_table
    global global_node_idx
    global add_find_time,add_hit_time,cont_find_time,cont_hit_time
    global_node_idx=0
    unique_table.clear()
    computed_table.clear()
    add_find_time=0
    add_hit_time=0
    cont_find_time=0
    cont_hit_time=0
    global_node_idx=0


def get_identity_tdd():
    node = Find_Or_Add_Unique_table(-1)
    tdd = TDD(node)
    tdd.index_2_key={-1:-1}
    tdd.key_2_index={-1:-1}
    return tdd

def get_unique_table():
    return unique_table

def get_unique_table_num():
    return len(unique_table)

def set_index_order(var_order):
    global global_index_order
    global_index_order=dict()
    if isinstance(var_order,list):
        for k in range(len(var_order)):
            global_index_order[var_order[k]]=k
    if isinstance(var_order,dict):
        global_index_order = copy.copy(var_order)
    global_index_order[-1] = float('inf')
    
def get_index_order():
    global global_index_order
    return copy.copy(global_index_order)
    
def get_int_key(weight):
    """To transform a complex number to a tuple with int values"""
    global epi
    return (int(round(weight.real/epi)) ,int(round(weight.imag/epi)))

def get_node_set(node,node_set=set()):
    """Only been used when counting the node number of a TDD"""
    if not node in node_set:
        node_set.add(node)
        for k in range(node.succ_num):
            if node.successor[k]:
                node_set = get_node_set(node.successor[k],node_set)
    return node_set

def Find_Or_Add_Unique_table(x,weigs=[],succ_nodes=[],the_map2=[]):
    """To return a node if it already exist, creates a new node otherwise"""
    global global_node_idx,unique_table
    
#     print('c',x,weigs,the_map)
    
    if x==-1:
        if unique_table.__contains__(x):
            return unique_table[x]
        else:
            res=Node(x)
            res.idx=0
            unique_table[x]=res
        return res
    temp_key=[x]
    for k in range(len(weigs)):
        temp_key.append(get_int_key(weigs[k]))
        temp_key.append(succ_nodes[k])
    
    temp_key += [the_map2.self_hash()]
    
    temp_key=tuple(temp_key)
    
    if temp_key in unique_table:
        return unique_table[temp_key]
    else:
        res=Node(x,len(succ_nodes))
        global_node_idx+=1
        res.idx=global_node_idx
        res.out_weight=weigs
        res.successor=succ_nodes
        res.out_maps = [the_maps(dict()),the_map2]
        unique_table[temp_key]=res
        
    return res


def normalize(x,the_successors):
    """The normalize and reduce procedure"""
#     print('a',x,the_successors[0].weight,the_successors[0].map,the_successors[1].weight,the_successors[1].map)
    all_zero = True
    
#     for k in range(len(the_successors)):
#         if get_int_key(the_successors[k].weight)==(0,0):
#             the_successors[k].map=the_maps(dict())
    
    for k in range(len(the_successors)):
        if the_successors[k].node.key!=-1:
            all_zero = False
            break
        if the_successors[k].weight!=0:
            all_zero = False
            break        
            
    if all_zero:
        return the_successors[0]
    
    the_map = the_maps(dict())
    
    if abs(the_successors[0].weight) < abs(the_successors[1].weight)-epi/2:
        the_map.data[x] = map_item('x')
        the_successors = [the_successors[1],the_successors[0]]
    
    the_map.add(the_successors[0].map)
    if get_int_key(the_successors[1].weight)==(0,0):
        weigs=[1,0]  
        succ_nodes=[succ.node for succ in the_successors]
        node=Find_Or_Add_Unique_table(x,weigs,succ_nodes,the_maps(dict()))
        res=TDD(node)
        res.weight=the_successors[0].weight
        res.map = the_map
        return res
        
    the_map2,the_phase = the_successors[1].map/the_successors[0].map      
                
    weig_max = the_successors[0].weight
    w2 = the_successors[1].weight/the_successors[0].weight
    
    w2_rotate = int(np.angle(w2)//rotate_angle)
    
    if abs(np.angle(w2)%rotate_angle-rotate_angle)<epi:
        w2_rotate+=1
        w2 = np.abs(w2)
    elif w2_rotate != 0:
        w2 = np.abs(w2)*np.exp(1j*(np.angle(w2)%rotate_angle))
    
    the_phase = (the_phase + w2_rotate)%root_of_unit
    
    if the_phase != 0:
        if x in the_map.data:
            the_map.data[x] = map_item(the_map.data[x].x,the_phase)
        else:
            the_map.data[x] = map_item(None,the_phase)

        
    weigs=[1,w2]  
    succ_nodes=[succ.node for succ in the_successors]
#     print('b',x,weigs,the_map2,weig_max,the_map)
    node=Find_Or_Add_Unique_table(x,weigs,succ_nodes,the_map2)
    res=TDD(node)
    res.weight=weig_max
    res.map = the_map
    
#     print('c',x,weigs,the_map2,weig_max,the_map)
    return res

def get_count():
    global add_find_time,add_hit_time,cont_find_time,cont_hit_time
    print("add:",add_hit_time,'/',add_find_time,'/',add_hit_time/add_find_time)
    print("cont:",cont_hit_time,"/",cont_find_time,"/",cont_hit_time/cont_find_time)

def find_computed_table(item):
    """To return the results that already exist"""
    global computed_table,add_find_time,add_hit_time,cont_find_time,cont_hit_time
    if item[0]=='s':
        temp_key=item[1].index_2_key[item[2]]
        the_key=('s',get_int_key(item[1].weight),item[1].node,item[1].map.self_hash(),temp_key,item[3])
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = TDD(res[1])
            tdd.weight = res[0]
            tdd.map=res[2]
            return tdd
    elif item[0] == '+':
        the_key=('+',get_int_key(item[1].weight),item[1].node,item[1].map.self_hash(),get_int_key(item[2].weight),item[2].node,item[2].map.self_hash())
        add_find_time+=1
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = TDD(res[1])
            tdd.weight = res[0]
            tdd.map=res[2]
            add_hit_time+=1
            return tdd
        the_key=('+',get_int_key(item[2].weight),item[2].node,item[2].map.self_hash(),get_int_key(item[1].weight),item[1].node,item[1].map.self_hash())
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = TDD(res[1])
            tdd.weight = res[0]
            tdd.map=res[2]
            add_hit_time+=1
            return tdd
    else:
        the_key=('*',get_int_key(item[1].weight),item[1].node,item[1].map.self_hash(),get_int_key(item[2].weight),item[2].node,item[2].map.self_hash(),item[3][0],item[3][1],item[4])
        cont_find_time+=1
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = TDD(res[1])
            tdd.weight = res[0]
            tdd.map=res[2]
            cont_hit_time+=1            
            return tdd
        the_key=('*',get_int_key(item[2].weight),item[2].node,item[2].map.self_hash(),get_int_key(item[1].weight),item[1].node,item[1].map.self_hash(),item[3][1],item[3][0],item[4])
        if computed_table.__contains__(the_key):
            res = computed_table[the_key]
            tdd = TDD(res[1])
            tdd.weight = res[0]
            tdd.map=res[2]
            cont_hit_time+=1            
            return tdd
    return None

def insert_2_computed_table(item,res):
    """To insert an item to the computed table"""
    global computed_table,cont_time,find_time,hit_time
    if item[0]=='s':
        temp_key=item[1].index_2_key[item[2]]
        the_key = ('s',get_int_key(item[1].weight),item[1].node,item[1].map.self_hash(),temp_key,item[3])
    elif item[0] == '+':
        the_key = ('+',get_int_key(item[1].weight),item[1].node,item[1].map.self_hash(),get_int_key(item[2].weight),item[2].node,item[2].map.self_hash())
    else:
        the_key = ('*',get_int_key(item[1].weight),item[1].node,item[1].map.self_hash(),get_int_key(item[2].weight),item[2].node,item[2].map.self_hash(),item[3][0],item[3][1],item[4])
    computed_table[the_key] = (res.weight,res.node,res.map)
    
def get_index_2_key(var):
    var_sort=copy.copy(var)
    var_sort.sort()
    var_sort.reverse()
    idx_2_key={-1:-1}
    key_2_idx={-1:-1}
    n=0
    for idx in var_sort:
        if not idx.key in idx_2_key:
            idx_2_key[idx.key]=n
            key_2_idx[n]=idx.key
            n+=1
    return idx_2_key,key_2_idx
    
def get_tdd(U,var=[]):
    
#     if len(var)==0 or not isinstance(var[0],Index):
#         return np_2_tdd(U,var)
    
    idx_2_key, key_2_idx = get_index_2_key(var)
    order=[]
    for idx in var:
        order.append(idx_2_key[idx.key])
        
    tdd = np_2_tdd(U,order)
    tdd.index_2_key=idx_2_key
    tdd.key_2_index=key_2_idx
    tdd.index_set=var
    
#     if not order:
#         order=list(range(U_dim))
        
#     for k in range(max(order)+1):
#         split_pos=order.index(k)
#         tdd.key_width[k]=U.shape[split_pos]
        
    return tdd    
    
def get_tdd2(U,var,idx_2_key=None):
    #index is the index_set as the axis order of the matrix
    U_dim=U.ndim
    if sum(U.shape)==U_dim:
        node=Find_Or_Add_Unique_table(-1)
        res=TDD(node)
        for k in range(U_dim):
            U=U[0]
        res.weight=U
        return res
     
    if not idx_2_key:
        idx_2_key = get_index_2_key(var)
        
    min_index=min(var)
    x=min_index.key
    min_pos=var.index(min_index)
#     print(min_pos)
    new_var=copy.copy(var)
    new_var[min_pos]=Index(-1)
    split_U=np.split(U,2,min_pos)
    new_var_key=[idx.key for idx in new_var]
    
    while x in new_var_key:
        min_pos=new_var_key.index(x)
        split_U[0]=np.split(split_U[0],2,min_pos)[0]
        split_U[1]=np.split(split_U[1],2,min_pos)[1]
        new_var[min_pos]=Index(-1)
        new_var_key[min_pos]=-1
        
    low=get_tdd(split_U[0],new_var,idx_2_key)
    high=get_tdd(split_U[1],new_var,idx_2_key)
    tdd = normalize(idx_2_key[x], [low,high])
    for k in range(len(idx_2_key)):
        tdd.key_width[k]=2
    tdd.index_2_key=idx_2_key
    key_2_idx=dict()
    for k in idx_2_key:
        key_2_idx[idx_2_key[k]]=k
    tdd.key_2_index=key_2_idx
    return tdd

def np_2_tdd(U,order=[],key_width=True):
    #index is the index_set as the axis order of the matrix
    U_dim=U.ndim
    U_shape=U.shape
    if sum(U_shape)==U_dim:
        node=Find_Or_Add_Unique_table(-1)
        res=TDD(node)
        for k in range(U_dim):
            U=U[0]
        res.weight=U
        return res
    
    if not order:
        order=list(range(U_dim))
    
    if key_width:
        the_width=dict()
        for k in range(max(order)+1):
            split_pos=order.index(k)
            the_width[k]=U.shape[split_pos]
            
    x=max(order)
    split_pos=order.index(x)
    order[split_pos]=-1
    split_U=np.split(U,U_shape[split_pos],split_pos)
    
    while x in order:
        split_pos=order.index(x)
        for k in range(len(split_U)):
            split_U[k]=np.split(split_U[k],U_shape[split_pos],split_pos)[k]
        order[split_pos]=-1
    
    the_successors=[]
    for k in range(U_shape[split_pos]):
        res=np_2_tdd(split_U[k],copy.copy(order),False)
        the_successors.append(res)
    tdd = normalize(x,the_successors)
    
    if key_width:
        tdd.key_width=the_width

    return tdd
    
    
def np_2_tdd2(U,split_pos=None):
    #index is the index_set as the axis order of the matrix
    U_dim=U.ndim
    U_shape=U.shape
    if sum(U_shape)==U_dim:
        node=Find_Or_Add_Unique_table(-1)
        res=TDD(node)
        for k in range(U_dim):
            U=U[0]
        res.weight=U
        return res
    if split_pos==None:
        split_pos=U_dim-1
        
    split_U=np.split(U,U_shape[split_pos],split_pos)
    the_successors=[]
    for k in range(U_shape[split_pos]):
        res=np_2_tdd(split_U[k],split_pos-1)
        the_successors.append(res)
    tdd = normalize(split_pos,the_successors)
    for k in range(len(U_shape)):
        tdd.key_width[k]=U_shape[k]
    return tdd
    
def tdd_2_np(tdd,split_pos=None,key_repeat_num=dict()):
#     print(split_pos,key_repeat_num)
    if split_pos==None:
        split_pos=tdd.node.key
            
    if split_pos==-1:
        return tdd.weight
    else:
        the_succs=[]
        for k in range(tdd.key_width[split_pos]):
            succ=Slicing2(tdd,split_pos,k)
            succ.key_width=tdd.key_width
            temp_res=tdd_2_np(succ,split_pos-1,key_repeat_num)
            the_succs.append(temp_res)
        if not split_pos in key_repeat_num:
            r = 1
        else:
            r = key_repeat_num[split_pos]
            
        if r==1:
            res=np.stack(tuple(the_succs), axis=the_succs[0].ndim)
        else:
            new_shape=list(the_succs[0].shape)
            for k in range(r):
                new_shape.append(tdd.key_width[split_pos])
            res=np.zeros(new_shape)
            for k1 in range(tdd.key_width[split_pos]):
                f='res['
#                 print(the_succs[0].ndim,r-1)
                for k2 in range(the_succs[0].ndim):
                    f+=':,'
                for k3 in range(r-1):
                    f+=str(k1)+','
                f=f[:-1]+']'
                eval(f)[k1]=the_succs[k1]
        return res
    
    
def get_measure_prob(tdd):
    if tdd.node.meas_prob:
        return tdd
    if tdd.node.key==-1:
        tdd.node.meas_prob=[0.5,0.5]
        return tdd
    if not tdd.node.succ_num==2:
        print("Only can be used for binary quantum state")
        return tdd
    get_measure_prob(Slicing(tdd,tdd.node.key,0))
    get_measure_prob(Slicing(tdd,tdd.node.key,1))
    tdd.node.meas_prob=[0]*2
    tdd.node.meas_prob[0]=abs(tdd.node.out_weight[0])**2*sum(tdd.node.successor[0].meas_prob)*2**(tdd.node.key-tdd.node.successor[0].key-1)
    tdd.node.meas_prob[1]=abs(tdd.node.out_weight[1])**2*sum(tdd.node.successor[1].meas_prob)*2**(tdd.node.key-tdd.node.successor[1].key-1)
    return tdd

    
def cont(tdd1,tdd2):

    var_cont=[var for var in tdd1.index_set if var in tdd2.index_set]
    var_out1=[var for var in tdd1.index_set if not var in var_cont]
    var_out2=[var for var in tdd2.index_set if not var in var_cont]

    var_out=var_out1+var_out2
    var_out.sort()
    var_out_idx=[var.key for var in var_out]
    var_cont_idx=[var.key for var in var_cont]
    var_cont_idx=[var for var in var_cont_idx if not var in var_out_idx]
    
    idx_2_key={-1:-1}
    key_2_idx={-1:-1}
    
    n=0
    for k in range(len(var_out_idx)-1,-1,-1):
        if not var_out_idx[k] in idx_2_key:
            idx_2_key[var_out_idx[k]]=n
            key_2_idx[n]=var_out_idx[k]
            n+=1
        
    key_2_new_key=[[],[]]
    cont_order=[[],[]]
    
#     the_map_temp = dict()
    
#     tdd1_map_temp = copy.copy(tdd1.map)
#     tdd2_map_temp = copy.copy(tdd2.map)
    
    for k in range(len(tdd1.key_2_index)-1):
        v=tdd1.key_2_index[k]
        if v in idx_2_key:
            key_2_new_key[0].append(idx_2_key[v])
            
#             if k in tdd1.map:
#                 the_map_temp[idx_2_key[v]]=tdd1.map[k]
#                 tdd1.map.pop(k)
        else:
            key_2_new_key[0].append('c')

            
        cont_order[0].append(global_index_order[v])
        
    cont_order[0].append(float('inf'))
    
    
    
    for k in range(len(tdd2.key_2_index)-1):     
        v=tdd2.key_2_index[k]
        if v in idx_2_key:
            key_2_new_key[1].append(idx_2_key[v])
            
#             if k in tdd2.map:
#                 the_map_temp[idx_2_key[v]]=tdd2.map[k]
#                 tdd2.map.pop(k)            
        else:
            key_2_new_key[1].append('c')
        cont_order[1].append(global_index_order[v])
    cont_order[1].append(float('inf'))

    tdd=contract(tdd1,tdd2,key_2_new_key,cont_order,len(set(var_cont_idx)))
    tdd.index_set=var_out
    tdd.index_2_key=idx_2_key
    tdd.key_2_index=key_2_idx
    
    
    key_width=dict()
    for k1 in range(len(key_2_new_key[0])):
        if not key_2_new_key[0][k1]=='c' and not key_2_new_key[0][k1] ==-1:
            key_width[key_2_new_key[0][k1]]=tdd1.key_width[k1]
    for k2 in range(len(key_2_new_key[1])):
        if not key_2_new_key[1][k2]=='c' and not key_2_new_key[1][k2] ==-1:
            key_width[key_2_new_key[1][k2]]=tdd2.key_width[k2]             
   
    tdd.key_width=key_width
#     print(tdd1.key_width,tdd2.key_width,tdd.key_width)
    return tdd
    
def cont2(tdd1,tdd2,cont_var):
    """cont_var is in the form [[0],[3]]"""
    key_2_new_key=[[],[]]
    cont_order=[[],[]]
    cont_num=len(cont_var[0])
    num1=0
    num2=0
    cont_var[0].append(0)
    cont_var[1].append(0)
    for k in range(cont_num):
        for k1 in range(cont_var[0][k-1],cont_var[0][k]):
            key_2_new_key[0].append(num1)
            num1+=1
            cont_order[0].append(num2)
            num2+=1
        
        for k2 in range(cont_var[1][k-1],cont_var[1][k]):
            key_2_new_key[1].append(num1)
            num1+=1
            cont_order[1].append(num2)
            num2+=1
        cont_order[0].append(num2)
        cont_order[1].append(num2)
        num2+=1
        key_2_new_key[0].append('c')
        key_2_new_key[1].append('c')
        
    for k1 in range(cont_var[0][cont_num-1],tdd1.node.key):
        key_2_new_key[0].append(num1)
        num1+=1
        cont_order[0].append(num2)
        num2+=1
    for k2 in range(cont_var[1][cont_num-1],tdd2.node.key):
        key_2_new_key[1].append(num1)
        num1+=1
        cont_order[1].append(num2)
        num2+=1
    the_max=max(max(cont_order[0],cont_order[1]))
    cont_order[0]=[the_max-k for k in cont_order[0]]
    cont_order[1]=[the_max-k for k in cont_order[1]]
    cont_order[0].append(float('inf'))
    cont_order[1].append(float('inf'))
#     print(key_2_new_key,cont_order)
    tdd=contract(tdd1,tdd2,key_2_new_key,cont_order,cont_num)
    key_width=dict()
    for k1 in range(len(key_2_new_key[0])):
        if not key_2_new_key[0][k1]=='c' and not key_2_new_key[0][k1] ==-1:
            key_width[key_2_new_key[0][k1]]=tdd1.key_width[k1]
    for k2 in range(len(key_2_new_key[1])):
        if not key_2_new_key[1][k2]=='c' and not key_2_new_key[1][k2] ==-1:
            key_width[key_2_new_key[1][k2]]=tdd2.key_width[k2]             
   
    tdd.key_width=key_width
    return tdd

def contract(tdd1,tdd2,key_2_new_key,cont_order,cont_num):
    """The contraction of two TDDs, var_cont is in the form [[4,1],[3,2]]"""
    
    k1=tdd1.node.key
    k2=tdd2.node.key
    w1=tdd1.weight
    w2=tdd2.weight
    
    if k1==-1 and k2==-1:
        if w1==0:
            tdd=TDD(tdd1.node)
            tdd.weight=0
            return tdd
        if w2==0:
            tdd=TDD(tdd1.node)
            tdd.weight=0
            return tdd
        tdd=TDD(tdd1.node)
        tdd.weight=w1*w2
        if cont_num>0:
            tdd.weight*=2**cont_num
        return tdd

    if k1==-1:
        if w1==0:
            tdd=TDD(tdd1.node)
            tdd.weight=0
            return tdd
        if cont_num ==0 and key_2_new_key[1][k2]==k2:
            tdd=TDD(tdd2.node)
            tdd.weight=w1*w2
            tdd.map=tdd2.map.self_copy()
            return tdd
            
    if k2==-1:
        if w2==0:
            tdd=TDD(tdd2.node)
            tdd.weight=0
            return tdd        
        if cont_num ==0 and key_2_new_key[0][k1]==k1:
            tdd=TDD(tdd1.node)
            tdd.weight=w1*w2
            tdd.map=tdd1.map.self_copy()
            return tdd
    
    tdd1.weight=1
    tdd2.weight=1
    
    
    temp_key_2_new_key=[]
    temp_key_2_new_key.append(tuple([k for k in key_2_new_key[0][:(k1+1)]]))
    temp_key_2_new_key.append(tuple([k for k in key_2_new_key[1][:(k2+1)]]))    
    
    
    map1 = tdd1.map.self_copy()
    map2 = tdd2.map.self_copy()
    
#     print(994,k1,k2,w1,w2,map1,map2)
    
    remain_map = the_maps(dict())
    
    for k in range(len(temp_key_2_new_key[0])-1,-1,-1):
        if k in tdd1.map.data:
            if temp_key_2_new_key[0][k]!='c':
                remain_map.data[temp_key_2_new_key[0][k]]=tdd1.map.data.pop(k)
                
    for k in range(len(temp_key_2_new_key[1])-1,-1,-1):
        if k in tdd2.map.data:
            if temp_key_2_new_key[1][k]!='c':
                remain_map.data[temp_key_2_new_key[1][k]]=tdd2.map.data.pop(k)
            
#     print(1008,map1,map2,remain_map,tdd1.map,tdd2.map)
    
    
    tdd=find_computed_table(['*',tdd1,tdd2,temp_key_2_new_key,cont_num])
    if tdd:
        tdd.weight=tdd.weight*w1*w2
        tdd1.weight=w1
        tdd2.weight=w2
        tdd1.map = map1
        tdd2.map = map2
        if not get_int_key(tdd.weight)==(0,0):
            tdd.map,the_phase = remain_map*tdd.map
            tdd.weight*=np.exp(1j*the_phase*rotate_angle)
        else:
            tdd.map=the_maps(dict())
        return tdd
                
    if cont_order[0][k1]<cont_order[1][k2]:
        the_key=key_2_new_key[0][k1]
        if the_key!='c':
            the_successors=[]
            for k in range(tdd1.node.succ_num):
                res=contract(Slicing(tdd1,k1,k),tdd2,key_2_new_key,cont_order,cont_num)
                the_successors.append(res)
            tdd=normalize(the_key,the_successors)
        else:
            tdd=TDD(Find_Or_Add_Unique_table(-1))
            tdd.weight=0
            for k in range(tdd1.node.succ_num):
                res=contract(Slicing(tdd1,k1,k),tdd2,key_2_new_key,cont_order,cont_num-1)           
                tdd=add(tdd,res)
    elif cont_order[0][k1]==cont_order[1][k2]:
        the_key=key_2_new_key[0][k1]
        if the_key!='c':
            the_successors=[]
            for k in range(tdd1.node.succ_num):
                res=contract(Slicing(tdd1,k1,k),Slicing(tdd2,k2,k),key_2_new_key,cont_order,cont_num)
                the_successors.append(res)
            tdd=normalize(the_key,the_successors)
        else:
            tdd=TDD(Find_Or_Add_Unique_table(-1))
            tdd.weight=0
            for k in range(tdd1.node.succ_num):
                res=contract(Slicing(tdd1,k1,k),Slicing(tdd2,k2,k),key_2_new_key,cont_order,cont_num-1)           
                tdd=add(tdd,res)
    else:
        the_key=key_2_new_key[1][k2]
        if the_key!='c':
            the_successors=[]
            for k in range(tdd2.node.succ_num):
                res=contract(tdd1,Slicing(tdd2,k2,k),key_2_new_key,cont_order,cont_num)
                the_successors.append(res)
            tdd=normalize(the_key,the_successors)
        else:
            tdd=TDD(Find_Or_Add_Unique_table(-1))
            tdd.weight=0
            for k in range(tdd2.node.succ_num):
                res=contract(tdd1,Slicing(tdd2,k2,k),key_2_new_key,cont_order,cont_num-1)           
                tdd=add(tdd,res)
                
    if get_int_key(tdd.weight)==(0,0):
        tdd.map=the_maps(dict())
    insert_2_computed_table(['*',tdd1,tdd2,temp_key_2_new_key,cont_num],tdd)
    tdd.weight=tdd.weight*w1*w2
    tdd1.weight=w1
    tdd2.weight=w2
    tdd1.map = map1
    tdd2.map = map2

    tdd.map,the_phase = remain_map*tdd.map
    tdd.weight*=np.exp(1j*the_phase*rotate_angle)
    return tdd
    
def Slicing(tdd,x,c):
    """Slice a TDD with respect to x=c"""

    k=tdd.node.key
    
    if k==-1:
        return tdd.self_copy()
    
    if k<x:
        return tdd.self_copy()
    
    if k==x:
        if not k in tdd.map.data:
            res=TDD(tdd.node.successor[c])
            res.weight=tdd.node.out_weight[c]
            if not get_int_key(res.weight)==(0,0):
                res.map,the_phase=tdd.map*tdd.node.out_maps[c]
                res.weight*=np.exp(1j*rotate_angle*the_phase)
#             print('1131',tdd.map,tdd.node.out_maps[c],res.map)
            return res
        if not tdd.map.data[k].x:
            res=TDD(tdd.node.successor[c])
            res.weight=tdd.node.out_weight[c]
            if not get_int_key(res.weight)==(0,0):
                if c==1:
                    res.weight*=np.exp(1j*rotate_angle*tdd.map.data[k].rotate)
                temp_map,the_phase = tdd.map*tdd.node.out_maps[c]
                if k in temp_map.data:
                    temp_map.data.pop(k)
                res.map=temp_map
                res.weight*=np.exp(1j*rotate_angle*the_phase)
#                 print('1144',tdd.map,tdd.node.out_maps[c],res.map)
        else:
            res=TDD(tdd.node.successor[1-c])
            res.weight=tdd.node.out_weight[1-c]
            if not get_int_key(res.weight)==(0,0):
                if tdd.map.data[k].x:
                    if c==0:
                        res.weight*=np.exp(1j*rotate_angle*tdd.map.data[k].rotate)
                temp_map,the_phase = tdd.map*tdd.node.out_maps[1-c]
                if k in temp_map.data:
                    temp_map.data.pop(k)
                res.map=temp_map
                res.weight*=np.exp(1j*rotate_angle*the_phase)
#                 print('1156',tdd.map,tdd.node.out_maps[1-c],res.map)
        return res
    else:
        print("Not supported yet!!!")
        
        
def Slicing2(tdd,x,c):
    """Slice a TDD with respect to x=c"""

    k=tdd.node.key
    
    if k==-1:
        return tdd.self_copy()
    
    if k<x:
        return tdd.self_copy()
        
    if k==x:
        if not k in tdd.map.data:
            res=TDD(tdd.node.successor[c])
            res.weight=tdd.node.out_weight[c]*tdd.weight
            if not get_int_key(res.weight)==(0,0):
                res.map,the_phase=tdd.map*tdd.node.out_maps[c]
                res.weight*=np.exp(1j*rotate_angle*the_phase)
            return res
        if not tdd.map.data[k].x:
            res=TDD(tdd.node.successor[c])
            res.weight=tdd.node.out_weight[c]*tdd.weight
            if not get_int_key(res.weight)==(0,0):
                if c==1:
                    res.weight*=np.exp(1j*rotate_angle*tdd.map.data[k].rotate)
                temp_map,the_phase = tdd.map*tdd.node.out_maps[c]
                if x in temp_map.data:
                    temp_map.data.pop(x)
                res.map=temp_map
                res.weight*=np.exp(1j*rotate_angle*the_phase)                
        else:
            res=TDD(tdd.node.successor[1-c])
            res.weight=tdd.node.out_weight[1-c]*tdd.weight
            if not get_int_key(res.weight)==(0,0):
                if tdd.map.data[k].x:
                    if c==0:
                        res.weight*=np.exp(1j*rotate_angle*tdd.map.data[k].rotate)
                temp_map,the_phase = tdd.map*tdd.node.out_maps[1-c]
                if x in temp_map.data:
                    temp_map.data.pop(x)
                res.map=temp_map
                res.weight*=np.exp(1j*rotate_angle*the_phase)
        return res
    else:
        print("Not supported yet!!!")        
        
        

def add(tdd1,tdd2):
    """The apply function of two TDDs. Mostly, it is used to do addition here."""
    global global_index_order    
    
    k1=tdd1.node.key
    k2=tdd2.node.key
    
    if tdd1.weight==0:
        return tdd2.self_copy()
    
    if tdd2.weight==0:
        return tdd1.self_copy()
    
    if tdd1.node==tdd2.node and tdd1.map==tdd2.map:

        weig=tdd1.weight+tdd2.weight
        if get_int_key(weig)==(0,0):
            term=Find_Or_Add_Unique_table(-1)
            res=TDD(term)
            res.weight=0
            return res
        else:
            res=TDD(tdd1.node)
            res.weight=weig
            res.map=tdd1.map.self_copy()
            return res
        
    if find_computed_table(['+',tdd1,tdd2]):
        return find_computed_table(['+',tdd1,tdd2])
    the_successors=[]
    if k1>k2:
        x=k1
        for k in range(tdd1.node.succ_num):
            res=add(Slicing2(tdd1,x,k),tdd2)
            the_successors.append(res)
    elif k1==k2:
        x=k1
        for k in range(tdd1.node.succ_num):
            res=add(Slicing2(tdd1,x,k),Slicing2(tdd2,x,k))
            the_successors.append(res)        
    else:
        x=k2
        for k in range(tdd2.node.succ_num):
            res=add(tdd1,Slicing2(tdd2,x,k))
            the_successors.append(res)
            
    res = normalize(x,the_successors)
    if get_int_key(res.weight)==(0,0):
        res.map=the_maps(dict())
    insert_2_computed_table(['+',tdd1,tdd2],res)
    return res