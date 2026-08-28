import numpy as np
from TDD.TN import Index,Tensor,TensorNetwork
from TDD.TDD import Ini_TDD, set_root_of_unit
from qiskit.quantum_info.operators import Operator
import time

def is_diagonal(U):
    i, j = np.nonzero(U)
    return np.all(i == j)

def add_hyper_index(var_list,hyper_index):
    for var in var_list:
        if not var in hyper_index:
            hyper_index[var]=0
            
def reshape(U):
    if U.shape==(1,1):
        return U
    
    if U.shape[0]==U.shape[1]:
        split_U=np.split(U,2,1)
    else:
        split_U=np.split(U,2,0)
    split_U[0]=reshape(split_U[0])
    split_U[1]=reshape(split_U[1]) 
    return np.array([split_U])[0]            
            
def get_real_qubit_num(cir):
    """Return the number of qubits of a circuit.

    Uses the declared ``cir.num_qubits`` so that qubits no gate touches
    (identity wires) are not silently dropped.
    """
    return cir.num_qubits

def cir_2_tn(cir):
    """return the dict that link every quantum gate to the corresponding index"""
#     print(1)
#     t=time.time()
    
    
    hyper_index=dict()
    qubits_index = dict()
    start_tensors= dict()
    end_tensors = dict()
    
    qubits_num=get_real_qubit_num(cir)

    for k in range(qubits_num):
        qubits_index[k]=0
        
    tn=TensorNetwork([],tn_type='cir',qubits_num=qubits_num)
                
    gates=cir.data
    for k in range(len(gates)):
        g=gates[k]
        nam=g[0].name
        q = [q._index for q in g[1]]
        q.reverse()
        var=[]

        ts=Tensor([],[],nam,q)
        
        if nam=='reset':
            continue

        U=Operator(g[0]).data
        
        for k in q:
            var_in='x'+ str(k)+'_'+str(qubits_index[k])
            var_out='x'+ str(k)+'_'+str(qubits_index[k]+1)
            add_hyper_index([var_in,var_out],hyper_index)
            var+=[Index(var_in,hyper_index[var_in]),Index(var_out,hyper_index[var_out])]
            if qubits_index[k]==0 and hyper_index[var_in]==0:
                start_tensors[k]=ts
            end_tensors[k]=ts                
            qubits_index[k]+=1
        if len(q)>1:
            U=reshape(U)
            
        if len(q)==1:
            U=U.T
        ts.data=U
        ts.index_set=var
        tn.tensors.append(ts)
        
#         for k in ts.index_set:
#             print(k)
#         print(ts.data)         

    for k in range(qubits_num):
        if k in start_tensors:
            last1=Index('x'+str(k)+'_'+str(0),0)
            new1=Index('x'+str(k),0)            
            start_tensors[k].index_set[start_tensors[k].index_set.index(last1)]=new1
        if k in end_tensors:
            last2=Index('x'+str(k)+'_'+str(qubits_index[k]),hyper_index['x'+str(k)+'_'+str(qubits_index[k])])
            new2=Index('y'+str(k),0)            
            end_tensors[k].index_set[end_tensors[k].index_set.index(last2)]=new2
               
    for k in range(qubits_num):
        U=np.eye(2)
        if qubits_index[k]==0 and not 'x'+str(k)+'_'+str(0) in hyper_index:
            var_in='x'+str(k)
            var=[Index('x'+str(k),0),Index('y'+str(k),0)]
            ts=Tensor(U,var,'nu_q',[k])
            tn.tensors.append(ts)            
    
    all_indexs=[]
    for k in range(qubits_num):
        all_indexs.append('x'+str(k))
        for k1 in range(qubits_index[k]+1):
            all_indexs.append('x'+str(k)+'_'+str(k1))
        all_indexs.append('y'+str(k))
#     print(4)
#     print(time.time()-t)
    return tn,all_indexs

def _to_bit_list(bits, qubits_num, what):
    """Normalize a computational-basis bitstring/list to ``list[int]``.

    Accepts either a ``str`` of ``'0'``/``'1'`` (first char = qubit 0) or a
    ``list``/``tuple`` of ``0``/``1`` ints. Raises on length mismatch or
    non-binary entries instead of silently corrupting the tensor network.
    """
    if isinstance(bits, str):
        if len(bits) != qubits_num:
            raise ValueError(
                f"{what} length {len(bits)} does not match qubit count {qubits_num}"
            )
        out = []
        for c in bits:
            if c == "0":
                out.append(0)
            elif c == "1":
                out.append(1)
            else:
                raise ValueError(f"{what} contains non-binary character {c!r}")
        return out
    if isinstance(bits, (list, tuple)):
        if len(bits) != qubits_num:
            raise ValueError(
                f"{what} length {len(bits)} does not match qubit count {qubits_num}"
            )
        out = [int(b) for b in bits]
        if any(b not in (0, 1) for b in out):
            raise ValueError(f"{what} must contain only 0/1, got {bits!r}")
        return out
    raise TypeError(f"{what} must be str or list/tuple of 0/1, got {type(bits).__name__}")


def add_inputs(tn,input_s,qubits_num):
    U0=np.array([1,0])
    U1=np.array([0,1])
    bits = _to_bit_list(input_s, qubits_num, "input state")
    for k in range(qubits_num-1,-1,-1):
        ts=Tensor(U0 if bits[k] == 0 else U1,[Index('x'+str(k))],'in',[k])
        tn.tensors.insert(0,ts)

def add_outputs(tn,output_s,qubits_num):
    U0=np.array([1,0])
    U1=np.array([0,1])
    bits = _to_bit_list(output_s, qubits_num, "output state")
    for k in range(qubits_num):
        ts=Tensor(U0 if bits[k] == 0 else U1,[Index('y'+str(k))],'out',[k])
        tn.tensors.append(ts)

def add_trace_line(tn,qubits_num):
    U=np.eye(2)
    for k in range(qubits_num-1,-1,-1):
        var_in='x'+str(k)
        var=[Index('x'+str(k),0),Index('y'+str(k),0)]
        ts=Tensor(U,var,'tr',[k])
        tn.tensors.insert(0,ts)


def simulate(cir, initial_state=None, root_of_unit=2 ** 8):
    """Simulate a Qiskit circuit and return the statevector.

    Returns a ``(2**n,)`` complex :class:`numpy.ndarray` using **little-endian**
    basis ordering, matching the C++ ``test_state_output`` convention and
    Qiskit's ``Statevector.data``: ``state[i]`` has bit ``k`` = qubit ``k``
    (bit 0 is the least significant / qubit 0).

    Parameters
    ----------
    cir:
        Qiskit circuit to simulate.
    initial_state:
        Computational-basis product state, either a ``str`` of ``'0'``/``'1'``
        (first char = qubit 0) or a list/tuple of ``0``/``1`` ints. Defaults to
        the all-zero state.
    root_of_unit:
        Phase discretization (``rotate_angle = 2*pi/root_of_unit``). Default
        ``256`` is exact for Clifford+T gate phases.
    """
    n = cir.num_qubits
    if initial_state is None:
        initial_state = "0" * n

    tn, all_indexs = cir_2_tn(cir)
    add_inputs(tn, initial_state, n)
    Ini_TDD(index_order=all_indexs)
    set_root_of_unit(root_of_unit)
    tdd = tn.cont()

    state = np.array(
        [tdd.get_amplitude([(i >> k) & 1 for k in range(n)]) for i in range(2 ** n)],
        dtype=complex,
    )
    return state


def gen_cir(name=None,qubit_num = 1,gate_num = 1):
    from qiskit import QuantumCircuit
    import random
    cir=QuantumCircuit(qubit_num)
    
    if name=='Random_Clifford':
        gate_set = ['x','y','z','h','s','cx']
        
        for k in range(gate_num):
            g = gate_set[random.randint(0,len(gate_set)-1)]
            q = random.randint(0,qubit_num-1)
            if g=='cx':
                q2 = random.randint(0,qubit_num-1)
                while q2==q:
                    q2 = random.randint(0,qubit_num-1)
                eval('cir.'+g+str(tuple([q,q2])))
            else:
                eval('cir.'+g+str(tuple([q])))
                
        return cir
    
    if name=='Random_Clifford_T':
        gate_set = ['x','y','z','h','s','cx','t']
        
        for k in range(gate_num):
            g = gate_set[random.randint(0,len(gate_set)-1)]
            q = random.randint(0,qubit_num-1)
            if g=='cx':
                q2 = random.randint(0,qubit_num-1)
                while q2==q:
                    q2 = random.randint(0,qubit_num-1)
                eval('cir.'+g+str(tuple([q,q2])))
            else:
                eval('cir.'+g+str(tuple([q])))
                
        return cir    
    